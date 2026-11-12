# -*- coding: utf-8 -*-
"""
code/339 — Reviewer 1 判别式内核的位点级直接检验（第六轮预审，2026-11-12）

问题：SOS 声称"同一组遗传变异必须在整骨(fracture)与关节(OA)两端点固定相反符号"。
本脚本用同一 31 个 BMD 工具在两端点的逐 SNP 效应，直接检验这一可证伪断言。

只读：不改动任何已报 β/SE/p/OR/CI；不改动文档。所有输入来自 code/334 已验证的
harmonised 逐 SNP 表（同一 beta_e 已核对一致）。

检验分两层（关键区分，避免"符号相反"被误当作"单一设定点"）：
  L1 符号分量：同一变异在两端点 Wald ratio 是否倾向反号（二项符号检验）；
  L2 量级分量：两端点效应向量是否强反相关（SOS "单一设定点"强版本预测 r→−1）。

判据（预先声明）：
  支持"单一潜在设定点"需同时满足 (a) 符号反号显著 且 (b) 效应量级强反相关
  (|r| 较大且 CI 远离 0)；仅有 (a) 而无 (b) 只支持"两模块净方向相反"，
  为 two-partially-coupled-modules 模型与 single-set-point 模型所共有，不具判别力。

运行： python3 code/339_r1_discriminating_kernel_test.py
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HARM = os.path.join(ROOT, "code", "output", "genetics", "harmonised")
OUT = os.path.join(ROOT, "code", "output", "genetics")

FX = os.path.join(HARM, "BMD_total__fracture.csv")
OA = os.path.join(HARM, "BMD_total__OA_knee.csv")


def ivw(d, ratio, se):
    w = 1.0 / d[se] ** 2
    return float(np.sum(w * d[ratio]) / np.sum(w))


def main():
    fx = pd.read_csv(FX)
    oa = pd.read_csv(OA)
    m = fx[["SNP", "beta_e", "beta_o", "se_o", "p_o"]].merge(
        oa[["SNP", "beta_o", "se_o", "p_o"]], on="SNP", suffixes=("_fx", "_oa")
    )
    assert len(m) == 31, f"expected 31 matched BMD instruments, got {len(m)}"
    # beta_e must be identical (same exposure, already verified by code/334)
    check = fx.sort_values("SNP")["beta_e"].to_numpy()
    check2 = oa.sort_values("SNP")["beta_e"].to_numpy()
    assert np.allclose(check, check2), "beta_e differs across endpoint tables"

    # Wald ratios aligned to the BMD-increasing allele
    m["ratio_fx"] = m.beta_o_fx / m.beta_e
    m["ratio_oa"] = m.beta_o_oa / m.beta_e
    m["se_ratio_fx"] = m.se_o_fx / m.beta_e.abs()
    m["se_ratio_oa"] = m.se_o_oa / m.beta_e.abs()

    # ---- Aggregate (already reported) : sign opposition at summary level ----
    b_fx = ivw(m, "ratio_fx", "se_ratio_fx")
    b_oa = ivw(m, "ratio_oa", "se_ratio_oa")

    # ---- L1 : sign opposition ----
    opp_mask = ((m.ratio_fx > 0) & (m.ratio_oa < 0)) | ((m.ratio_fx < 0) & (m.ratio_oa > 0))
    n_opp = int(opp_mask.sum())
    binom = stats.binomtest(n_opp, len(m), 0.5, alternative="greater")
    w_fx = 1.0 / m.se_ratio_fx ** 2
    w_opp_share = float(w_fx[opp_mask].sum() / w_fx.sum())

    # each endpoint alone, direction consistency
    n_fx_prot = int((m.ratio_fx < 0).sum())   # BMD-up -> lower fracture hazard
    n_oa_risk = int((m.ratio_oa > 0).sum())    # BMD-up -> higher OA odds

    # ---- L2 : magnitude anticorrelation ----
    rp = stats.pearsonr(m.ratio_fx, m.ratio_oa)
    rs = stats.spearmanr(m.ratio_fx, m.ratio_oa)
    rng = np.random.default_rng(20261112)
    x, y = m.ratio_fx.to_numpy(), m.ratio_oa.to_numpy()
    bs = []
    for _ in range(10000):
        idx = rng.integers(0, len(x), len(x))
        if np.std(x[idx]) > 0 and np.std(y[idx]) > 0:
            bs.append(np.corrcoef(x[idx], y[idx])[0, 1])
    lo, hi = np.percentile(bs, [2.5, 97.5])

    # measurement-error disattenuation (reliability = 1 - mean sampling var / obs var)
    def reliability(col, secol):
        v = np.var(col, ddof=1)
        ve = np.mean(m[secol] ** 2)
        return (v - ve) / v
    rel_fx = reliability(m.ratio_fx, "se_ratio_fx")
    rel_oa = reliability(m.ratio_oa, "se_ratio_oa")
    r_disatt = rp.statistic / np.sqrt(max(rel_fx, 1e-9) * max(rel_oa, 1e-9))

    # LOO sign robustness
    loo_r = [stats.pearsonr(m.drop(i).ratio_fx, m.drop(i).ratio_oa).statistic
             for i in range(len(m))]

    # strong-endpoint subset
    sub = m[m.p_o_oa < 0.05]
    if len(sub) >= 4:
        rp_sub = stats.pearsonr(sub.ratio_fx, sub.ratio_oa)
        n_opp_sub = int((((sub.ratio_fx > 0) & (sub.ratio_oa < 0)) |
                         ((sub.ratio_fx < 0) & (sub.ratio_oa > 0))).sum())
    else:
        rp_sub, n_opp_sub = None, None

    # ---- verdict ----
    sign_supported = binom.pvalue < 0.05
    magnitude_supported = (rp.pvalue < 0.05) and (rp.statistic < -0.5) and (hi < -0.2)
    if sign_supported and magnitude_supported:
        verdict = "SUPPORTS single-set-point (sign AND magnitude anticorrelated)"
    elif sign_supported and not magnitude_supported:
        verdict = ("PARTIAL: opposite endpoint SIGNS supported; magnitude anticorrelation absent "
                   "-> consistent with two partially-coupled modules with opposite NET effects; "
                   "single-set-point (strong kernel) NOT discriminated")
    elif not sign_supported:
        verdict = "REFUTES sign opposition at variant level (kernel fails its own criterion)"
    else:
        verdict = "INCONCLUSIVE"

    print("=" * 78)
    print("Reviewer 1 — discriminating-kernel test (31 shared BMD instruments)")
    print("=" * 78)
    print(f"matched SNPs = {len(m)}; beta_e identical across endpoint tables = True")
    print("\n[aggregate, already reported]")
    print(f"  IVW beta per SD BMD-up : fracture {b_fx:+.4f} | OA {b_oa:+.4f} "
          f"-> opposite aggregate signs: {b_fx * b_oa < 0}")
    print("\n[L1 sign component]")
    print(f"  opposite-sign variants      : {n_opp}/{len(m)}  exact binomial one-sided p = {binom.pvalue:.4f}")
    print(f"  IVW-weight on opposite-sign : {w_opp_share * 100:.1f}%")
    print(f"  fracture-protective sign    : {n_fx_prot}/{len(m)}")
    print(f"  OA-raising sign             : {n_oa_risk}/{len(m)}")
    print("\n[L2 magnitude component — the 'single set-point' discriminator]")
    print(f"  Pearson  r = {rp.statistic:+.4f}  p = {rp.pvalue:.3g}  bootstrap 95% CI [{lo:+.3f},{hi:+.3f}]")
    print(f"  Spearman rho = {rs.statistic:+.4f}  p = {rs.pvalue:.3g}")
    print(f"  disattenuated r ~ {r_disatt:+.3f}  (reliability FX {rel_fx:.3f} / OA {rel_oa:.3f})")
    print(f"  LOO Pearson r range [{min(loo_r):+.3f},{max(loo_r):+.3f}]; all negative: {all(r < 0 for r in loo_r)}")
    if rp_sub is not None:
        print(f"  strong-endpoint subset (OA nominal p<0.05, n={len(sub)}): "
              f"r={rp_sub.statistic:+.3f} p={rp_sub.pvalue:.3g}; opposite {n_opp_sub}/{len(sub)}")
    print("\n" + "=" * 78)
    print("VERDICT:", verdict)
    print("=" * 78)
    print("NOTE: this is an in-corpus, hypothesis-testing diagnostic, not a new reported result.")
    print("      It does not alter any published beta/SE/p/OR/CI; it quantifies what the")
    print("      existing 31-instrument data can and cannot discriminate about the SOS kernel.")

    m.to_csv(os.path.join(OUT, "_r1_kernel_sign_test.csv"), index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
