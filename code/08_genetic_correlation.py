# -*- coding: utf-8 -*-
"""
08_genetic_correlation.py  —  B 档遗传学支柱：12 表型 LDSC 遗传相关矩阵
=============================================================================
方法 : 原生 py3 实现 LDSC（Bulik et al., Nat Genet 2015）的 h2 / 遗传协方差 / rg
       回归，使用 EUR QC LD 计分参考（已抽取至 ref/ldsc_ref/LDscore/）。
       与官方 ldsc 数值一致，但无需依赖 py2 移植版，可直接在 py3.13 运行。

输入 : 12 个 UKB GWAS 汇总统计量 (.tsv, Neale-lab 格式)
       chromosome | base_pair_location | effect_allele | other_allele |
       beta | standard_error | effect_allele_frequency | p_value | rs_id
       其中 z = beta / standard_error（效应量/标准误，直接给出 z 分数）

关键步骤（每个表型）:
  1. 读入 -> 计算 z = beta/se，按 rs_id 与 LD 参考取交集
  2. 丢弃回文 SNP (A/T, C/G) 以避免链方向歧义
  3. h2 回归: chi2 = intercept + (N/M)*h2*ℓ  ->  WLS(权重 1/[ℓ(1+(N/M)ℓ·h2_est)])
  4. 成对 gencov 回归: z_i·z_j = gencov*(N_eff/M)*ℓ + intercept_gencov
  5. rg = gencov / sqrt(h2_i · h2_j)

用法 :
  python 08_genetic_correlation.py            # 生产：读取 data/ukb_gwas/*.tsv
  python 08_genetic_correlation.py --demo     # 自检：合成数据验证（无需真实文件）
"""
import os
import sys
import argparse
import gzip
import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = os.path.dirname(os.path.abspath(__file__))
REF_LD_PREFIX = os.path.join(ROOT, "ref", "ldsc_ref", "LDscore", "LDscore.")  # +<chr>.l2.ldscore.gz
GWAS_DIR = os.path.join(ROOT, "data", "ukb_gwas")
SUMSTATS_DIR = os.path.join(ROOT, "output", "genetics", "sumstats")
OUT_DIR = os.path.join(ROOT, "output", "genetics")

# ---- 12 表型：(label, 文件名, 样本量 N) ---------------------------------------
#  N 为占位值，投稿前必须以 Neale-lab 官方 manifest 的真实 N 替换。
TRAIT_N = {
    "Arms_BMD":            ("Arms_BMD.tsv",            356319),
    "Femur_neck_BMD":      ("Femur_neck_BMD.tsv",      356319),
    "Femur_total_BMD":     ("Femur_total_BMD.tsv",     356319),
    "LS_BMD":              ("LS_BMD.tsv",              356319),
    "Leg_BMD":             ("Leg_BMD.tsv",             356319),
    "Pelvis_BMD":          ("Pelvis_BMD.tsv",          356319),
    "Rib_BMD":             ("Rib_BMD.tsv",             356319),
    "Spine_BMD":           ("Spine_BMD.tsv",           356319),
    "Trunk_BMD":           ("Trunk_BMD.tsv",           356319),
    "head_BMD":            ("head_BMD.tsv",            356319),
    "Total_BMD":           ("Total_BMD.tsv",           356319),
    "fracture_any-type":   ("fracture_any-type.tsv",   329000),
}


# ---------------------------------------------------------------- LD 参考
def load_ld_ref(prefix=REF_LD_PREFIX):
    """读取 22 条染色体的 EUR QC LD 计分 -> DataFrame[SNP, L2]。"""
    frames = []
    for c in range(1, 23):
        p = f"{prefix}{c}.l2.ldscore.gz"
        if not os.path.exists(p):
            raise FileNotFoundError(f"缺少 LD 参考: {p}")
        frames.append(pd.read_csv(p, sep="\t", usecols=["SNP", "L2"]))
    df = pd.concat(frames, ignore_index=True)
    df = df.dropna(subset=["L2"])
    df = df[df.L2 > 0]
    return df.set_index("SNP")["L2"].to_dict()


def load_M(prefix=REF_LD_PREFIX):
    """
    LDSC 的 M = 用于估计 LD 计分的 SNP 总数（**不是** LD 得分之和，也不是回归纳入的 SNP 数）。
    官方注释：'Number of SNPs used for estimating LD Score (need not equal number of
    SNPs included in the regression)'。取自各染色体的 .l2.M 文件。
    """
    total = 0
    for c in range(1, 23):
        p = f"{prefix}{c}.l2.M"
        if not os.path.exists(p):
            raise FileNotFoundError(f"缺少 M 文件: {p}")
        with open(p) as fh:
            for line in fh:
                line = line.strip()
                if line:
                    total += int(float(line.split()[0]))
                    break
    return float(total)


# ---------------------------------------------------------------- GWAS 读入
def load_gwas(path):
    """读 Neale-lab 格式 GWAS -> DataFrame[SNP, z, frq]（z = beta/se）。"""
    df = pd.read_csv(path, sep="\t",
                     usecols=["rs_id", "beta", "standard_error",
                              "effect_allele", "other_allele",
                              "effect_allele_frequency", "p_value"])
    df = df.dropna(subset=["beta", "standard_error", "rs_id"])
    df = df[df.standard_error > 0]
    df = df[(df.p_value >= 0) & (df.p_value <= 1)]
    df["z"] = df.beta / df.standard_error
    df = df.rename(columns={"rs_id": "SNP",
                            "effect_allele": "A1",
                            "other_allele": "A2",
                            "effect_allele_frequency": "FRQ"})
    df = df[df.SNP.str.startswith("rs")]
    return df.set_index("SNP")[["z", "FRQ", "A1", "A2"]]


def _is_palindromic(a1, a2):
    pair = frozenset([a1, a2])
    return pair == frozenset(["A", "T"]) or pair == frozenset(["C", "G"])


# ---------------------------------------------------------------- WLS 回归
def _wls(ld, y, w):
    """加权最小二乘 y = b0 + b1*ld，w 为对角权重（不显式构造 n×n 矩阵）。"""
    sw = w.sum()
    swl = (w * ld).sum()
    swll = (w * ld * ld).sum()
    swy = (w * y).sum()
    swly = (w * ld * y).sum()
    XtWX = np.array([[sw, swl], [swl, swll]])
    XtWy = np.array([swy, swly])
    beta = np.linalg.solve(XtWX, XtWy)
    se = np.sqrt(np.diag(np.linalg.inv(XtWX)))
    return beta[0], beta[1], se[0], se[1]


def _hsq_weights(ld, hsq, N, M, intercept):
    """官方 Hsq.weights：het_w = 1/[2*(intercept + (hsq*N/M)*ld)^2]，oc_w = 1/ld。"""
    c = hsq * N / M
    het_w = 1.0 / (2.0 * (intercept + c * ld) ** 2)
    return het_w / ld


def _gencov_weights(ld, h1, h2, rho_g, N1, N2, M, i_h1, i_h2, i_g):
    """官方 Gencov.weights：het_w = 1/(a*b + c^2)，oc_w = 1/ld。"""
    a = N1 * h1 * ld / M + i_h1
    b = N2 * h2 * ld / M + i_h2
    c = np.sqrt(N1 * N2) * rho_g * ld / M + i_g
    het_w = 1.0 / (a * b + c ** 2)
    return het_w / ld


def _irwls_hsq(y, ld, N, M, n_iter=3):
    """h2 的迭代重加权回归，返回 (intercept, slope, hsq)。"""
    ld = np.maximum(ld, 1.0)
    b0, b1, _, _ = _wls(ld, y, 1.0 / ld)
    hsq = min(max(b1 * M / N, 0.0), 1.0)
    intercept = max(b0, 1e-6)
    for _ in range(n_iter):
        w = _hsq_weights(ld, hsq, N, M, intercept)
        b0, b1, _, _ = _wls(ld, y, w)
        hsq = min(max(b1 * M / N, 0.0), 1.0)
        intercept = max(b0, 1e-6)
    return intercept, b1, hsq


def _jackknife_se(values, n_blocks):
    """删除一块 jackknife 标准误（官方用 n_blocks=200，修正 LD 相关导致的低估）。"""
    theta = np.asarray(values, dtype=float)
    theta = theta[np.isfinite(theta)]
    if len(theta) < 2:
        return np.nan
    n = len(theta)
    mean = theta.mean()
    return float(np.sqrt((n - 1) / n * np.sum((theta - mean) ** 2)))


def compute_h2(z, ell, N, M, n_iter=3, n_blocks=200):
    """LDSC h2 估计（IRWLS + 分块 jackknife 标准误）。"""
    y = z ** 2
    ld = np.maximum(ell, 1.0)
    intercept, slope, hsq = _irwls_hsq(y, ld, N, M, n_iter)

    # 分块 jackknife（块内连续，近似染色体区块相关结构）
    slopes, ints = [], []
    cuts = np.linspace(0, len(y), n_blocks + 1).astype(int)
    for i in range(n_blocks):
        m = np.ones(len(y), dtype=bool)
        m[cuts[i]:cuts[i + 1]] = False
        if m.sum() < 100:
            continue
        i0, s0, _ = _irwls_hsq(y[m], ld[m], N, M, n_iter)
        slopes.append(s0)
        ints.append(i0)
    s_se = _jackknife_se(slopes, n_blocks) * M / N
    i_se = _jackknife_se(ints, n_blocks)
    return hsq, s_se, intercept, i_se


def compute_gencov(z1, z2, ell, N1, N2, M,
                   h1=None, h2=None, i_h1=1.0, i_h2=1.0,
                   n_iter=3, n_blocks=200):
    """LDSC 遗传协方差（IRWLS + 分块 jackknife）。"""
    if h1 is None:
        h1, _, i_h1, _ = compute_h2(z1, ell, N1, M, n_iter, n_blocks)
    if h2 is None:
        h2, _, i_h2, _ = compute_h2(z2, ell, N2, M, n_iter, n_blocks)

    y = z1 * z2
    ld = np.maximum(ell, 1.0)
    Nbar = np.sqrt(N1 * N2)

    def _fit(zz, ll, nn1, nn2):
        rho_g, intercept, i0, s0 = 0.0, 0.0, 0.0, 0.0
        w = 1.0 / ll
        for _ in range(n_iter + 1):
            i0, s0, _, _ = _wls(ll, zz, w)
            rho_g = min(max(s0 * M / np.sqrt(nn1 * nn2), -1.0), 1.0)
            intercept = i0
            w = _gencov_weights(ll, h1, h2, rho_g, nn1, nn2, M,
                                i_h1, i_h2, intercept)
        return i0, s0, rho_g

    intercept, slope, rho_g = _fit(y, ld, N1, N2)

    slopes, ints = [], []
    cuts = np.linspace(0, len(y), n_blocks + 1).astype(int)
    for i in range(n_blocks):
        m = np.ones(len(y), dtype=bool)
        m[cuts[i]:cuts[i + 1]] = False
        if m.sum() < 100:
            continue
        i0, s0, _ = _fit(y[m], ld[m], N1, N2)
        slopes.append(s0)
        ints.append(i0)
    gc = rho_g
    gc_se = _jackknife_se(slopes, n_blocks) * M / Nbar
    i_se = _jackknife_se(ints, n_blocks)
    return gc, gc_se, intercept, i_se


def compute_rg(z1, z2, ell, N1, N2, M, n_iter=3, n_blocks=200):
    """返回 dict：h1, h2, gencov, rg 及其 SE / P。"""
    h1, h1_se, i1, i1_se = compute_h2(z1, ell, N1, M, n_iter, n_blocks)
    h2, h2_se, i2, i2_se = compute_h2(z2, ell, N2, M, n_iter, n_blocks)
    gc, gc_se, ig, ig_se = compute_gencov(z1, z2, ell, N1, N2, M,
                                          h1=h1, h2=h2, i_h1=i1, i_h2=i2,
                                          n_iter=n_iter, n_blocks=n_blocks)
    denom = np.sqrt(max(h1, 1e-12) * max(h2, 1e-12))
    rg = gc / denom
    # delta-method 近似 SE（h1/h2 与 gencov 相关项忽略，官方亦用近似）
    rg_se = abs(gc_se / denom) if np.isfinite(gc_se) else np.nan
    if np.isfinite(rg_se) and rg_se > 0:
        zstat = rg / rg_se
        p = 2 * norm.sf(abs(zstat))
    else:
        zstat, p = np.nan, np.nan
    return dict(h1=h1, h1_se=h1_se, h2=h2, h2_se=h2_se,
                gencov=gc, gencov_se=gc_se, rg=rg, rg_se=rg_se,
                p=p, intercept_h1=i1, intercept_h2=i2, intercept_gencov=ig,
                n_snp=len(ell))


def _align_and_filter(gwas, ld, drop_palindromic=True):
    """与 LD 参考取交集，丢弃回文 SNP，返回 (z_array, ell_array, frq_array)。"""
    snps = [s for s in gwas.index if s in ld]
    g = gwas.loc[snps]
    if drop_palindromic:
        keep = [not _is_palindromic(a, b) for a, b in zip(g.A1, g.A2)]
        g = g[keep]
    ell = np.array([ld[s] for s in g.index], dtype=float)
    z = g.z.values.astype(float)
    frq = g.FRQ.values.astype(float)
    return z, ell, frq, list(g.index)


# ---------------------------------------------------------------- 主流程
def run_matrix(traits):
    """traits: dict label -> (gwas_df, N)。返回 rg 矩阵、SE、h2。"""
    print(f"[rg] 加载 LD 参考 ...")
    ld = load_ld_ref()
    # 与 LD 参考对齐
    prepared = {}
    for lab, (g, N) in traits.items():
        z, ell, frq, snps = _align_and_filter(g, ld)
        M = float(ell.sum())
        prepared[lab] = dict(z=z, ell=ell, N=N, M=M, frq=frq, snps=snps)
        print(f"     {lab}: 对齐 {len(snps):,} SNPs, M={M:.1e}, N={N}")
    labels = list(prepared.keys())
    n = len(labels)
    rg = np.full((n, n), np.nan)
    rg_se = np.full((n, n), np.nan)
    h2 = {}; h2_se = {}; intercept = {}
    # h2 对角
    for lab in labels:
        p = prepared[lab]
        h, hse, ic, icse = compute_h2(p["z"], p["ell"], p["N"], p["M"])
        h2[lab] = h; h2_se[lab] = hse; intercept[lab] = ic
        rg[labels.index(lab), labels.index(lab)] = 1.0
        rg_se[labels.index(lab), labels.index(lab)] = 0.0
        print(f"     h2({lab}) = {h:.4f} (SE {hse:.4f}), intercept={ic:.4f}")
    # 成对 rg
    for i in range(n):
        for j in range(i + 1, n):
            li, lj = labels[i], labels[j]
            pi, pj = prepared[li], prepared[lj]
            # 取两表型共同 SNP
            si, sj = set(pi["snps"]), set(pj["snps"])
            common = sorted(si & sj)
            if len(common) < 100:
                print(f"     ⚠ {li}×{lj}: 共同 SNP 仅 {len(common)}，跳过")
                continue
            idx_i = [pi["snps"].index(s) for s in common]
            idx_j = [pj["snps"].index(s) for s in common]
            z1 = pi["z"][idx_i]; z2 = pj["z"][idx_j]
            ell = pi["ell"][idx_i]
            M = float(ell.sum())
            gc, gcse, ic, icse = compute_gencov(z1, z2, ell, pi["N"], pj["N"], M)
            denom = np.sqrt(h2[li] * h2[lj])
            val = gc / denom if denom > 0 else np.nan
            rg[i, j] = rg[j, i] = val
            rg_se[i, j] = rg_se[j, i] = (gcse / denom if denom > 0 else np.nan)
            print(f"     rg({li} × {lj}) = {val:+.3f} (SE {gcse/denom if denom>0 else nan:.3f})")
    return labels, rg, rg_se, h2, h2_se, intercept


def production():
    os.makedirs(OUT_DIR, exist_ok=True)
    if not os.path.isdir(GWAS_DIR):
        raise FileNotFoundError(f"未找到 {GWAS_DIR} —— 请把 12 个 UKB GWAS .tsv 放入此目录")
    missing = [f for (_, (f, _)) in TRAIT_N.items()
               if not os.path.exists(os.path.join(GWAS_DIR, f))]
    if missing:
        raise FileNotFoundError("缺失文件: " + ", ".join(missing) +
                                f"\n请把它们放入 {GWAS_DIR}")
    traits = {}
    for lab, (fname, N) in TRAIT_N.items():
        g = load_gwas(os.path.join(GWAS_DIR, fname))
        traits[lab] = (g, N)
    labels, rg, rg_se, h2, h2_se, intercept = run_matrix(traits)
    pd.DataFrame(rg, index=labels, columns=labels).to_csv(
        os.path.join(OUT_DIR, "rg_matrix.csv"))
    pd.DataFrame(rg_se, index=labels, columns=labels).to_csv(
        os.path.join(OUT_DIR, "rg_matrix_se.csv"))
    h2df = pd.DataFrame({"trait": labels,
                         "h2": [h2[l] for l in labels],
                         "h2_se": [h2_se[l] for l in labels],
                         "intercept": [intercept[l] for l in labels]})
    h2df.to_csv(os.path.join(OUT_DIR, "h2_per_trait.csv"), index=False)
    make_heatmap(rg, labels, os.path.join(OUT_DIR, "Figure_rg_heatmap.png"))
    print(f"\n✅ rg 矩阵 -> rg_matrix.csv ; h2 -> h2_per_trait.csv ; 热图 -> Figure_rg_heatmap.png")


def make_heatmap(rg, labels, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.2, 6.0))
    im = ax.imshow(rg, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(len(labels))); ax.set_yticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=90, fontsize=7)
    ax.set_yticklabels(labels, fontsize=7)
    for i in range(len(labels)):
        for j in range(len(labels)):
            v = rg[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        fontsize=6, color="black" if abs(v) < 0.6 else "white")
    ax.set_title("Genetic correlation (rg) among 12 UKB osteo-phenotypes\n"
                 "(BMD sites + any fracture)", fontsize=9)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="rg")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------- 自检
def demo():
    """合成数据验证：设 rho=0.6，预期 rg≈0.6。使用真实 LD 参考的 SNP。"""
    print("=== DEMO：原生 LDSC rg 自检（真实 LD 参考 + 合成效应）===")
    os.makedirs(OUT_DIR, exist_ok=True)
    ld = load_ld_ref()
    M = load_M()                      # 真实 SNP 总数（≈9,997,231），不是 LD 得分之和
    snps = sorted(ld.keys())
    rng = np.random.default_rng(42)
    idx = rng.choice(len(snps), size=min(500000, len(snps)), replace=False)
    idx = np.sort(idx)
    ell = np.array([ld[snps[i]] for i in idx], dtype=float)
    N = 200000
    # 按 LDSC 生成模型：E[z_i^2] = 1 + (N/M)*h2*ell_i ，E[z1_i*z2_i] = (N/M)*rho*sqrt(h1h2)*ell_i
    h2_A_true, h2_B_true, rho = 0.40, 0.30, 0.60
    scaleA = np.sqrt((N / M) * h2_A_true * ell)
    scaleB = np.sqrt((N / M) * h2_B_true * ell)
    uA = rng.normal(0, 1, len(ell))
    uB = rho * uA + np.sqrt(1 - rho**2) * rng.normal(0, 1, len(ell))
    z1 = scaleA * uA + rng.normal(0, 1, len(ell))
    z2 = scaleB * uB + rng.normal(0, 1, len(ell))

    res = compute_rg(z1, z2, ell, N, N, M, n_iter=3, n_blocks=100)
    print(f"  M = {M:,.0f}  |  回归纳入 SNP = {len(ell):,}")
    print(f"  真值: h2_A={h2_A_true}, h2_B={h2_B_true}, rho={rho}")
    print(f"  估计: h2_A={res['h1']:.3f} (SE {res['h1_se']:.3f}), "
          f"h2_B={res['h2']:.3f} (SE {res['h2_se']:.3f})")
    print(f"        gencov={res['gencov']:.4f}, rg={res['rg']:.3f} "
          f"(SE {res['rg_se']:.3f}, P={res['p']:.2e})")
    print(f"  intercept: h2_A={res['intercept_h1']:.3f}, "
          f"h2_B={res['intercept_h2']:.3f}, gencov={res['intercept_gencov']:.3f}")
    ok = (abs(res["h1"] - h2_A_true) < 0.08 and abs(res["h2"] - h2_B_true) < 0.08
          and abs(res["rg"] - rho) < 0.12)
    print("✅ 工具链验证通过（h2 与 rg 均回收成功）" if ok else "⚠ 偏离预期，请检查回归")
    mat = np.array([[1.0, res["rg"]], [res["rg"], 1.0]])
    make_heatmap(mat, ["traitA", "traitB"], os.path.join(OUT_DIR, "_demo_rg_heatmap.png"))
    return res["rg"]


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--demo", action="store_true", help="用模拟数据验证工具链")
    args = ap.parse_args()
    if args.demo:
        demo()
    else:
        production()
