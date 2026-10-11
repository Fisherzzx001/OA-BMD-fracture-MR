#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
13_mr_analysis.py — OA ↔ BMD ↔ fracture 的两样本孟德尔随机化（Python，无 R）

遵循 local-mr-neale-gcst skill 的流程与坑位清单：
  · 用 chromosome+base_pair_location 建键（绝不用 rs_id）
  · 1-Mb 近似 clumping（每窗口取最显著 SNP）—— 正式发表需换 PLINK r²<0.001
  · 加权中位数权重 = 1/se_exposure²（不是 1/X²）
  · IVW + Egger(截距) + 加权中位数 + Cochran Q/I² + Radial 离群 + 留一法
  · 必做阳性对照：遗传性高 BMD → 低骨折

尺度说明（重要）：
  · 12 个已 munge 的 sumstats 只有 Z/N/MAF，按 SE = 1/sqrt(N·2p(1-p))、beta = Z·SE 重构。
    对二分类结局（OA），该 beta 在**线性概率尺度**而非 log-OR；因病例比例在全 SNP 恒定，
    差一个常数因子 → **不影响方向与显著性**，只影响绝对效应量。
  · fracture GWAS 有真实 beta/se（log-OR 尺度），直接使用。
"""
import os
import sys
import importlib
import numpy as np
import pandas as pd
from scipy.stats import norm, chi2

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
m09 = importlib.import_module("09_rg_matrix")
m10 = importlib.import_module("10_fracture_row")

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "output", "genetics")
P_THRESH = 5e-8
CLUMP_WINDOW = 1e6          # 1 Mb 近似 clumping


# ------------------------------------------------------------------ 载入
def load_munged(name):
    """已 munge 的 LDSC sumstats -> beta/se（重构）。"""
    df = pd.read_csv(os.path.join(m09.FLAGSHIP_LDSC, m09.PHENOS[name]), sep=r"\s+")
    df = df.dropna(subset=["SNP", "Z", "MAF", "N"])
    df = df[(df.MAF > 0) & (df.MAF < 1) & (df.N > 0)]
    se = 1.0 / np.sqrt(df.N.values * 2 * df.MAF.values * (1 - df.MAF.values))
    df = df.assign(beta=df.Z.values * se, se=se)
    return df[["CHR", "BP", "SNP", "A1", "A2", "beta", "se", "MAF", "P"]]


def load_fracture_df():
    """骨折 GWAS 有真实 beta/se（log-OR 尺度），直接用，不做 Z 重构。"""
    df = pd.read_csv(
        m10.FRACTURE, sep="\t",
        usecols=["chromosome", "base_pair_location", "effect_allele",
                 "other_allele", "beta", "standard_error",
                 "effect_allele_frequency", "p_value"],
        dtype={"chromosome": "int32", "base_pair_location": "int64"},
    )
    df = df.dropna(subset=["beta", "standard_error"])
    df = df[df.standard_error > 0]
    df = df.rename(columns={"chromosome": "CHR", "base_pair_location": "BP",
                            "effect_allele": "A1", "other_allele": "A2",
                            "standard_error": "se",
                            "effect_allele_frequency": "MAF", "p_value": "P"})
    df["A1"] = df.A1.astype(str).str.upper()
    df["A2"] = df.A2.astype(str).str.upper()
    return df[["CHR", "BP", "A1", "A2", "beta", "se", "MAF", "P"]]


def load_all():
    dfs = {nm: load_munged(nm) for nm in m09.PHENOS}
    dfs["fracture"] = load_fracture_df()
    return dfs


# ------------------------------------------------------------------ 工具变量
def select_instruments(df, p_thresh=P_THRESH, window=CLUMP_WINDOW):
    """p<p_thresh，按 1-Mb 窗口贪心 clumping（每窗口留最显著 SNP）。"""
    cand = df[df.P < p_thresh].sort_values("P")
    keep = []
    for chr_, grp in cand.groupby("CHR"):
        taken = []
        for _, r in grp.iterrows():
            if all(abs(r.BP - b) > window for b in taken):
                keep.append(r)
                taken.append(r.BP)
    return pd.DataFrame(keep).reset_index(drop=True)


def harmonise(inst, out):
    """按 chr:pos join，等位基因按名字对齐；EA==OA 则翻转，都不等则丢弃。"""
    o = out.set_index(["CHR", "BP"])
    rows = []
    for _, r in inst.iterrows():
        k = (r.CHR, r.BP)
        if k not in o.index:
            continue
        q = o.loc[k]
        if o.index.duplicated().any() and isinstance(q, pd.DataFrame):
            q = q.iloc[0]
        a1e, a2e, a1o, a2o = r.A1, r.A2, q.A1, q.A2
        if a1e == a1o and a2e == a2o:
            f = 1.0
        elif a1e == a2o and a2e == a1o:
            f = -1.0
        else:
            continue                       # 等位基因对不上 / 回文无法判定
        rows.append(dict(SNP=r.SNP, CHR=r.CHR, BP=r.BP,
                         beta_e=r.beta, se_e=r.se, p_e=r.P,
                         beta_o=f * q.beta, se_o=q.se, p_o=q.P,
                         eaf=r.MAF, flip=f))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ MR 估计
def _ivw(be, se_e, bo, so):
    w = 1.0 / (so ** 2)                     # 结局方差加权（标准 IVW）
    return float(np.sum(w * be * bo) / np.sum(w * be ** 2))


def _egger(be, bo, so):
    w = 1.0 / (so ** 2)
    x, y = be, bo
    mx = np.sum(w * x) / np.sum(w)
    my = np.sum(w * y) / np.sum(w)
    b = np.sum(w * (x - mx) * (y - my)) / np.sum(w * (x - mx) ** 2)
    a = my - b * mx
    # 截距 SE
    n = len(x)
    resid = y - (a + b * x)
    s2 = np.sum(w * resid ** 2) / max(n - 2, 1)
    se_b = np.sqrt(s2 / np.sum(w * (x - mx) ** 2))
    se_a = np.sqrt(s2 * (1 / np.sum(w) + mx ** 2 / np.sum(w * (x - mx) ** 2)))
    return float(b), float(a), float(se_b), float(se_a)


def _wmed(be, se_e, bo, so):
    """
    Burgess 2017 加权中位数：
      比值 β_j = b_out/b_exp；权重 w_j ∝ (b_exp/se_out)²  = 1/Var(β_j)；
      按 β_j 排序后取加权中位数，落在两个次序统计量之间时取均值。
    """
    ratio = bo / be
    w = (be / so) ** 2                      # = 1/Var(ratio)
    order = np.argsort(ratio)
    r, w = ratio[order], w[order]
    cw = np.cumsum(w)
    cw = cw / cw[-1]
    idx = int(np.searchsorted(cw, 0.5))
    idx = min(idx, len(r) - 1)
    if idx == 0:
        return float(r[0])
    if abs(cw[idx - 1] - 0.5) < 1e-12:
        return float(0.5 * (r[idx - 1] + r[idx]))
    return float(r[idx])


def _ivw_se(be, so):
    """IVW 一阶 delta 法 SE（固定效应）：1/sqrt(Σ b_e²/se_o²)。"""
    return float(1.0 / np.sqrt(np.sum(be ** 2 / so ** 2)))


def _ivw_re_se(be, bo, so):
    """随机效应 IVW SE：乘过离散因子 φ² = Q/df（高异质性时用）。"""
    b = _ivw(be, None, bo, so)
    resid = bo - b * be
    w = 1.0 / so ** 2
    Q = float(np.sum(w * resid ** 2))
    dfree = max(len(be) - 1, 1)
    phi2 = max(Q / dfree, 1.0)
    return float(np.sqrt(phi2 / np.sum(be ** 2 / so ** 2)))


def _q(be, bo, so):
    b = _ivw(be, None, bo, so)
    resid = bo - b * be
    w = 1.0 / (so ** 2)
    Q = float(np.sum(w * resid ** 2))
    dfree = max(len(be) - 1, 1)
    return Q, dfree, float(chi2.sf(Q, dfree))


def _radial(be, bo, so, z_thresh=3.0):
    """
    学生化残差剔除离群 SNP（MR-PRESSO 近似）。
    残差方差已知 = se_out²，故学生化残差 = resid/se_out（不要额外再除尺度）。
    """
    keep = np.ones(len(be), dtype=bool)
    for _ in range(10):
        if keep.sum() < 4:
            break
        b = _ivw(be[keep], None, bo[keep], so[keep])
        resid = bo - b * be
        z = np.where(so > 0, resid / so, 0.0)
        bad = keep & (np.abs(z) > z_thresh)
        if not bad.any():
            break
        keep &= ~bad
    return keep


def run_mr(be, se_e, bo, so, snps=None):
    be, se_e, bo, so = map(lambda a: np.asarray(a, float), (be, se_e, bo, so))
    res = dict(n_snp=len(be))
    if len(be) < 3:
        res["error"] = f"IV 数量不足 ({len(be)})"
        return res

    res["ivw"] = _ivw(be, se_e, bo, so)
    res["ivw_se"] = _ivw_se(be, so)
    res["ivw_z"] = res["ivw"] / res["ivw_se"] if res["ivw_se"] > 0 else np.nan
    res["ivw_p"] = float(2 * norm.sf(abs(res["ivw_z"]))) if np.isfinite(res["ivw_z"]) else np.nan

    if len(be) >= 3:
        sb, sa, se_b, se_a = _egger(be, bo, so)
        res["egger"] = sb; res["egger_se"] = se_b
        res["egger_p"] = float(2 * norm.sf(abs(sb / se_b))) if se_b > 0 else np.nan
        res["egger_intercept"] = sa; res["egger_intercept_se"] = se_a
        res["egger_intercept_p"] = float(2 * norm.sf(abs(sa / se_a))) if se_a > 0 else np.nan
    if len(be) >= 3:
        res["wmed"] = _wmed(be, se_e, bo, so)

    Q, dfree, Qp = _q(be, bo, so)
    res["Q"], res["Q_df"], res["Q_p"] = Q, dfree, Qp
    res["I2"] = float(max(Q - dfree, 0) / Q) if Q > 0 else 0.0
    res["ivw_re_se"] = _ivw_re_se(be, bo, so)
    res["ivw_re_p"] = float(2 * norm.sf(abs(res["ivw"] / res["ivw_re_se"]))) \
        if res["ivw_re_se"] > 0 else np.nan

    keep = _radial(be, bo, so)
    res["n_outlier"] = int((~keep).sum())
    if keep.sum() >= 3:
        res["ivw_radial"] = _ivw(be[keep], None, bo[keep], so[keep])
        res["ivw_radial_se"] = float(np.sqrt(np.sum((so[keep] / be[keep]) ** 2)
                                             / keep.sum()))

    # 留一法
    loo = []
    for i in range(len(be)):
        m = np.ones(len(be), bool); m[i] = False
        if m.sum() >= 2:
            loo.append(_ivw(be[m], None, bo[m], so[m]))
    if loo:
        res["loo_min"], res["loo_max"] = float(np.min(loo)), float(np.max(loo))
        res["loo_sign_flips"] = int(sum(1 for v in loo if np.sign(v) != np.sign(res["ivw"])))
    res["snps"] = list(snps) if snps is not None else None
    return res


def fmt(res):
    if "error" in res:
        return f"  ⚠️ {res['error']}"
    s = (f"  IVW(固定) β={res['ivw']:+.4f} (SE {res['ivw_se']:.4f}) "
         f"P={res['ivw_p']:.2e}  n_snp={res['n_snp']}")
    if "ivw_re_se" in res:
        s += (f"\n  IVW(随机) β={res['ivw']:+.4f} (SE {res['ivw_re_se']:.4f}) "
              f"P={res['ivw_re_p']:.2e}")
    if "egger" in res:
        s += (f"\n  Egger    β={res['egger']:+.4f} (SE {res['egger_se']:.4f}) "
              f"P={res.get('egger_p', float('nan')):.2e} | "
              f"截距={res['egger_intercept']:+.5f} (P={res.get('egger_intercept_p', float('nan')):.3f})")
    if "wmed" in res:
        s += f"\n  加权中位数 β={res['wmed']:+.4f}"
    s += (f"\n  Q={res['Q']:.1f} (df={res['Q_df']}, P={res['Q_p']:.3f}) "
          f"I²={res['I2']:.0%} | Radial 离群 {res['n_outlier']} 个")
    if "loo_min" in res:
        s += (f"\n  留一法   β ∈ [{res['loo_min']:+.4f}, {res['loo_max']:+.4f}] "
              f"符号翻转 {res['loo_sign_flips']}")
    return s


# ------------------------------------------------------------------ 主流程
PAIRS = [
    # (暴露, 结局, 说明)  —— 阳性对照放最前
    ("BMD_total", "fracture",      "阳性对照：高 BMD 应降低骨折"),
    ("BMD_femur_neck", "fracture", "阳性对照（股骨颈）"),
    ("BMD_total", "OA_knee",       "BMD → 膝 OA"),
    ("BMD_femur_neck", "OA_knee",  "股骨颈 BMD → 膝 OA"),
    ("BMD_lumbar", "OA_knee",      "腰椎 BMD → 膝 OA"),
    ("BMD_arms", "OA_knee",        "上肢 BMD → 膝 OA"),
    ("BMD_leg", "OA_knee",         "下肢 BMD → 膝 OA"),
    ("OA_knee", "fracture",        "膝 OA → 骨折"),
    ("OA_knee", "BMD_total",       "反向：膝 OA → BMD"),
    ("OA_knee", "BMD_femur_neck",  "反向：膝 OA → 股骨颈 BMD"),
]


def main():
    os.makedirs(OUT, exist_ok=True)
    print("▶ 载入 GWAS")
    dfs = load_all()
    print(f"  {len(dfs)} 个表型")

    print("▶ 选取工具变量（p<5e-8, 1-Mb 近似 clumping）")
    inst = {}
    for nm, d in dfs.items():
        iv = select_instruments(d)
        inst[nm] = iv
        print(f"  {nm:18s} GWS={int((d.P < P_THRESH).sum()):5d} → clumped={len(iv)}")

    print("\n▶ 两样本 MR")
    rows = []
    for exp, out, note in PAIRS:
        print(f"\n【{exp} → {out}】 {note}")
        h = harmonise(inst[exp], dfs[out])
        if len(h) < 3:
            print(f"  ⚠️ harmonise 后 IV 不足 ({len(h)})")
            continue
        res = run_mr(h.beta_e.values, h.se_e.values, h.beta_o.values, h.se_o.values,
                     snps=h.SNP.tolist())
        print(fmt(res))
        rows.append(dict(exposure=exp, outcome=out, note=note, **{
            k: v for k, v in res.items() if k != "snps"}))

    pd.DataFrame(rows).to_csv(os.path.join(OUT, "mr_results.csv"), index=False)
    print(f"\n✅ 输出至 {OUT}/mr_results.csv")


if __name__ == "__main__":
    main()
