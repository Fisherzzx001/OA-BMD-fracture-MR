#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
09_rg_matrix.py — 12 表型 LDSC 遗传相关矩阵（原生 Python 3 实现）

算法严格对齐官方 LDSC (Bulik-Sullivan et al. 2015, regressions.py)：
  h2      : IRWLS，het_w = 1/[2*(intercept + (h2*N/M)*l)^2]，oc_w = 1/l
  gencov  : IRWLS，het_w = 1/(a*b + c^2)，oc_w = 1/l
  rg      : gencov / sqrt(h1*h2)
  SE      : 删除一块 jackknife（修正 LD 相关导致的解析 SE 低估）
  M       : 用于估计 LD 计分的 SNP 总数（取自 .l2.M），**不是** LD 得分之和

与调用官方 ldsc 相比，本实现的额外保障：
  1. 等位基因对齐：以参考表型为基准翻转 z，避免 z1*z2 符号翻转导致 rg 错号；
  2. 丢弃回文 SNP（A/T、C/G），链方向不可判定；
  3. 12 表型取共同 SNP 交集，保证矩阵内各格估计基于同一变异集（官方为逐对交集）。

用法：
  python 09_rg_matrix.py            # 全量（12 表型，66 对）
  python 09_rg_matrix.py --quick    # 快速自检（2 表型，降低 jackknife 块数）
"""
import os
import sys
import argparse
import importlib
import numpy as np
import pandas as pd
from scipy.stats import norm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
_mod08 = importlib.import_module("08_genetic_correlation")   # 模块名以数字开头
_wls = _mod08._wls
_hsq_weights = _mod08._hsq_weights
_gencov_weights = _mod08._gencov_weights
compute_h2 = _mod08.compute_h2
compute_gencov = _mod08.compute_gencov
_jackknife_se = _mod08._jackknife_se

# ------------------------------------------------------------------ 路径
ROOT = os.path.dirname(os.path.abspath(__file__))
FLAGSHIP_LDSC = r"C:\Users\Administrator\WorkBuddy\2026-09-08-20-17-24\OA_FLAGSHIP\07_reference\ldsc"
OUT_DIR = os.path.join(ROOT, "output", "genetics")

# 已 munge 的 LDSC 汇总统计（真实数据）：label -> 文件名
PHENOS = {
    "OA_knee":        "sumstats_oa_knee.sumstats.gz",
    "OA_selfrep":     "sumstats_oa_selfrep.sumstats.gz",
    "BMD_total":      "sumstats_total.sumstats.gz",
    "BMD_femur_neck": "sumstats_femur_neck.sumstats.gz",
    "BMD_femur_total": "sumstats_femur_total.sumstats.gz",
    "BMD_lumbar":     "sumstats_lumbar.sumstats.gz",
    "BMD_spine":      "sumstats_spine.sumstats.gz",
    "BMD_arms":       "sumstats_arms.sumstats.gz",
    "BMD_leg":        "sumstats_leg.sumstats.gz",
    "BMD_pelvis":     "sumstats_pelvis.sumstats.gz",
    "BMD_rib":        "sumstats_rib.sumstats.gz",
    "BMD_trunk":      "sumstats_trunk.sumstats.gz",
}
REF_PHENO = "OA_knee"          # 等位基因对齐 + 方向校正的基准表型
MAF_MIN, MAF_MAX = 0.01, 0.99

# ---------------------------------------------------------------------------
# Z 方向校正（**关键** —— 2026-09-26 整改后已不再需要全局翻转）
# ---------------------------------------------------------------------------
# 历史：早期 OA_selfrep 文件（占位 N=350,000）的 Z 与 A1 整体反号，曾用
#   SIGN_FLIP={"OA_selfrep":-1} 校正（依据为 GDF5/FTO 两个哨兵位点 + 外部文献锚点）。
# 2026-09-26 该文件经重新 munge，N 已更正为真实 462,933（IEU ukb-b-14486），
# 且等位基因方向变为与独立 GO Consortium 2025 膝 OA GWAS 一致。
#
# 独立验证（见 docs/11_补充材料_OA_selfrep_Z校正.md）：在 OA_knee 全部 1,980 个
# 基因组显著位点（|Z|>5.45）上，当前 OA_selfrep 的 Z 符号与独立 OA_knee 的
# 一致性为 94.4%（原始，未翻转）；按染色体分解同样呈 84–100% 一致。两个 OA
# 表型遗传相关必为正，故当前文件已 concordant，**live 管线不再施加全局翻转**。
# 下列 SIGN_FLIP 留空，仅为历史可读性保留接口。
SIGN_FLIP = {}


def _is_palindromic(a1, a2):
    return frozenset([a1, a2]) in (frozenset(["A", "T"]), frozenset(["C", "G"]))


# ------------------------------------------------------------------ 读入
def load_ld_ref(prefix=None):
    prefix = prefix or os.path.join(FLAGSHIP_LDSC, "eur_w_ld_chr_")
    frames = []
    for c in range(1, 23):
        p = f"{prefix}{c}.l2.ldscore.gz"
        if not os.path.exists(p):
            raise FileNotFoundError(f"缺少 LD 参考: {p}")
        d = pd.read_csv(p, sep=r"\s+", usecols=["SNP", "L2"])
        frames.append(d)
    df = pd.concat(frames, ignore_index=True)
    df = df.dropna(subset=["L2"])
    df = df[df.L2 > 0]
    return df.set_index("SNP")["L2"].to_dict()


def load_M(path=None):
    """M = 用于估计 LD 计分的 SNP 总数。"""
    path = path or os.path.join(FLAGSHIP_LDSC, "eur_w_ld_chr.l2.M")
    with open(path) as fh:
        return float(fh.read().strip().split()[0])


def load_sumstats(fname):
    df = pd.read_csv(os.path.join(FLAGSHIP_LDSC, fname), sep=r"\s+")
    df = df.dropna(subset=["SNP", "Z", "A1", "A2"])
    df["A1"] = df.A1.astype(str).str.upper()
    df["A2"] = df.A2.astype(str).str.upper()
    df = df[~df.SNP.duplicated(keep="first")]
    return df.set_index("SNP")


# ------------------------------------------------------------------ 对齐
def harmonise(dfs, ld):
    """
    1) 取所有表型 ∩ LD 参考的共同 SNP；
    2) 丢弃回文 SNP；
    3) 以 REF_PHENO 的等位基因为基准，翻转 z 使所有表型 z 指向同一等位基因。
    返回 (Z 矩阵 DataFrame, ell 数组, N 字典, 诊断 dict)
    """
    common = set(ld.keys())
    for d in dfs.values():
        common &= set(d.index)
    common = sorted(common)
    diag = {"n_ld_ref": len(ld), "n_common_raw": len(common)}

    ref = dfs[REF_PHENO]
    keep, Z = [], {}
    for s in common:
        a1, a2 = ref.at[s, "A1"], ref.at[s, "A2"]
        if _is_palindromic(a1, a2):
            continue
        keep.append(s)
    diag["n_after_palindromic"] = len(keep)

    refA1 = ref.loc[keep, "A1"].values
    refA2 = ref.loc[keep, "A2"].values
    flips = {}
    for name, d in dfs.items():
        z = d.loc[keep, "Z"].values.astype(float)
        a1 = d.loc[keep, "A1"].values
        a2 = d.loc[keep, "A2"].values
        same = a1 == refA1
        swap = (a1 == refA2) & (a2 == refA1)
        ok = same | swap
        z = np.where(swap, -z, z)          # 等位基因相反 -> z 取反
        z = z * SIGN_FLIP.get(name, 1.0)   # 文件级方向校正
        Z[name] = np.where(ok, z, np.nan)  # 无法判定 -> 剔除
        flips[name] = int(swap.sum())

    Zdf = pd.DataFrame(Z, index=keep)
    # 丢弃任一表型无法判定的 SNP
    Zdf = Zdf.dropna()
    diag["n_flipped"] = flips
    diag["n_final"] = len(Zdf)
    ell = np.array([ld[s] for s in Zdf.index], dtype=float)
    Ns = {name: float(np.nanmean(dfs[name].loc[Zdf.index, "N"].values))
          for name in dfs}
    return Zdf, ell, Ns, diag


# ------------------------------------------------------------------ 推断
def run_matrix(Zdf, ell, Ns, M, n_iter=3, n_blocks=100, verbose=True):
    names = list(Zdf.columns)
    h2 = {}
    for i, nm in enumerate(names):
        r = compute_h2(Zdf[nm].values, ell, Ns[nm], M, n_iter, n_blocks)
        h2[nm] = dict(h2=r[0], se=r[1], intercept=r[2], intercept_se=r[3],
                      N=Ns[nm])
        if verbose:
            print(f"  [{i+1}/{len(names)}] {nm:16s} h2={r[0]:.4f} "
                  f"(SE {r[1]:.4f}) intercept={r[2]:.3f} N={Ns[nm]:,.0f}")

    p = len(names)
    rg = pd.DataFrame(np.eye(p), index=names, columns=names)
    se = pd.DataFrame(np.zeros((p, p)), index=names, columns=names)
    pv = pd.DataFrame(np.zeros((p, p)), index=names, columns=names)
    gc = pd.DataFrame(np.zeros((p, p)), index=names, columns=names)
    for i in range(p):
        for j in range(i + 1, p):
            a, b = names[i], names[j]
            r = compute_gencov(Zdf[a].values, Zdf[b].values, ell,
                               Ns[a], Ns[b], M,
                               h1=h2[a]["h2"], h2=h2[b]["h2"],
                               i_h1=h2[a]["intercept"], i_h2=h2[b]["intercept"],
                               n_iter=n_iter, n_blocks=n_blocks)
            gencov, gc_se = r[0], r[1]
            denom = np.sqrt(max(h2[a]["h2"], 1e-12) * max(h2[b]["h2"], 1e-12))
            v = gencov / denom
            s = abs(gc_se / denom) if np.isfinite(gc_se) else np.nan
            pp = 2 * norm.sf(abs(v / s)) if (np.isfinite(s) and s > 0) else np.nan
            rg.loc[a, b] = rg.loc[b, a] = v
            se.loc[a, b] = se.loc[b, a] = s
            pv.loc[a, b] = pv.loc[b, a] = pp
            gc.loc[a, b] = gc.loc[b, a] = gencov
            if verbose:
                star = "***" if pp < 0.001 else "**" if pp < 0.01 else "*" if pp < 0.05 else ""
                print(f"    {a:16s} ~ {b:16s} rg={v:+.3f} (SE {s:.3f}) P={pp:.2e} {star}")
    return h2, rg, se, pv, gc


# ------------------------------------------------------------------ 输出
def make_heatmap(rg, pv, path, title="Genetic correlation (LDSC rg)"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = [l.replace("BMD_", "").replace("OA_", "OA:") for l in rg.index]
    n = len(labels)
    fig, ax = plt.subplots(figsize=(10.5, 9))
    vmax = min(1.0, float(np.nanmax(np.abs(rg.values - np.eye(n)))) * 1.15)
    im = ax.imshow(rg.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(n)); ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_yticks(range(n)); ax.set_yticklabels(labels)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            v, p = rg.values[i, j], pv.values[i, j]
            txt = f"{v:.2f}"
            if np.isfinite(p):
                txt += "*" if p < 0.05 else ""
            ax.text(j, i, txt, ha="center", va="center", fontsize=7.5,
                    color="white" if abs(v) > vmax * 0.55 else "black")
    ax.set_title(title + "   (* P<0.05)", fontsize=12)
    fig.colorbar(im, ax=ax, shrink=0.8)
    plt.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="2 表型快速自检")
    ap.add_argument("--blocks", type=int, default=100, help="jackknife 块数")
    args = ap.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)
    print("▶ 读入 LD 参考")
    ld = load_ld_ref()
    M = load_M()
    print(f"  LD 参考 SNP = {len(ld):,}  M = {M:,.0f}")

    names = list(PHENOS)
    if args.quick:
        names = ["OA_knee", "BMD_femur_neck"]
    print(f"▶ 读入 {len(names)} 个表型")
    dfs = {}
    for nm in names:
        p = os.path.join(FLAGSHIP_LDSC, PHENOS[nm])
        if not os.path.exists(p):
            print(f"  ⚠️ 缺失，跳过: {PHENOS[nm]}")
            continue
        dfs[nm] = load_sumstats(PHENOS[nm])
    names = list(dfs)

    print("▶ SNP 对齐与等位基因 harmonisation")
    Zdf, ell, Ns, diag = harmonise(dfs, ld)
    print(f"  LD 参考 {diag['n_ld_ref']:,} → 交集 {diag['n_common_raw']:,} "
          f"→ 去回文 {diag['n_after_palindromic']:,} → 最终 {diag['n_final']:,}")
    print(f"  各表型 z 翻转数: {diag['n_flipped']}")

    print("▶ LDSC 估计")
    h2, rg, se, pv, gc = run_matrix(Zdf, ell, Ns, M,
                                    n_blocks=(20 if args.quick else args.blocks))

    pd.DataFrame(h2).T.to_csv(os.path.join(OUT_DIR, "rg_h2.csv"),
                              index_label="phenotype")
    rg.to_csv(os.path.join(OUT_DIR, "rg_matrix.csv"))
    se.to_csv(os.path.join(OUT_DIR, "rg_se_matrix.csv"))
    pv.to_csv(os.path.join(OUT_DIR, "rg_p_matrix.csv"))
    gc.to_csv(os.path.join(OUT_DIR, "rg_gencov_matrix.csv"))
    make_heatmap(rg, pv, os.path.join(OUT_DIR, "Figure_rg_matrix.png"))
    print(f"✅ 输出至 {OUT_DIR}")


if __name__ == "__main__":
    main()
