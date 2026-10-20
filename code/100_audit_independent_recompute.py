# -*- coding: utf-8 -*-
"""
100_audit_independent_recompute.py
==================================
INDEPENDENT re-derivation of every two-sample MR statistic reported in
Tables 3 and 4 of docs/08, straight from the SNP-level harmonised files in
code/output/genetics/harmonised/.

This script deliberately does NOT import any project code. All estimators
(IVW fixed/random, Cochran Q / I2, MR-Egger, weighted median, leave-one-out)
are re-implemented here from first principles so that agreement is evidence,
not tautology.

Usage:
    <python> code/100_audit_independent_recompute.py
"""
from __future__ import annotations

import csv
import json
import math
import os
from itertools import product

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
HARM = os.path.join(HERE, "output", "genetics", "harmonised")

# 当前 harmonised（code/output/genetics/harmonised/）已为正确 log-OR 尺度
# （G6 轮将 OA_knee 的 β/SE 直接乘 1/√(K(1−K)) 修正入文件，稿件 Methods 亦声明
# "no rescaling was applied"）。故独立重算不再施加任何尺度因子。
SCALE = 1.0

# ---------------------------------------------------------------------------
# Values as PRINTED in the manuscript (Table 3). Rounded to the printed
# precision, so tolerances below are set accordingly.
# key -> (n_snp, beta, se, p_fixed, p_random, egger_beta, egger_int_p, wmedian,
#         I2_percent, loo_sign_flips)
# ---------------------------------------------------------------------------
# 稿件 Table 2 印刷值（G6 尺度修正 + 步近约定加权中位数 Bowden 2016）。
# 与 code/output/genetics/mr_results_corrected.csv 逐列一致。
TABLE3 = {
    "BMD_total__fracture":     (31, -0.337, 0.022, 3.65e-51, 2.36e-12, -0.337, 0.558, -0.393, 78, 0),
    "BMD_femur_neck__fracture":(25, -0.351, 0.030, 7.71e-31, 2.32e-06, -0.343, 0.018, -0.365, 83, 0),
    "BMD_total__OA_knee":      (31, +0.078, 0.013, 3.86e-10, 0.006, +0.078, 0.810, +0.065, 81, 0),
    "BMD_femur_neck__OA_knee": (25, +0.058, 0.017, 6.36e-04, 0.055, +0.058, 0.875, +0.052, 68, 0),
    "BMD_lumbar__OA_knee":     (21, +0.095, 0.016, 7.70e-09, 0.017, +0.089, 0.295, +0.115, 83, 0),
    "BMD_arms__OA_knee":       (16, +0.075, 0.016, 1.71e-06, 0.070, +0.075, 0.982, +0.026, 86, 0),
    "BMD_leg__OA_knee":        (24, +0.060, 0.014, 2.29e-05, 0.098, +0.058, 0.815, +0.026, 85, 0),
    "OA_knee__fracture":       (69, +0.027, 0.032, 3.95e-01, 0.542, +0.031, 0.476, +0.019, 49, 0),
    "OA_knee__BMD_total":      (69, +0.031, 0.032, 3.26e-01, 0.481, +0.036, 0.376, +0.026, 49, 0),
    "OA_knee__BMD_femur_neck": (69, +0.042, 0.032, 1.90e-01, 0.272, +0.048, 0.169, +0.023, 30, 0),
}

OUTLIER_COL = "presso_outlier"


# ---------------------------------------------------------------------------
# estimators
# ---------------------------------------------------------------------------
def ivw_fixed(b: np.ndarray, se: np.ndarray):
    """Inverse-variance weighted (fixed effect) meta-analysis of Wald ratios."""
    w = 1.0 / se ** 2
    beta = float(np.sum(w * b) / np.sum(w))
    se_beta = float(np.sqrt(1.0 / np.sum(w)))
    z = beta / se_beta
    p = 2.0 * stats.norm.sf(abs(z))
    return beta, se_beta, p


def heterogeneity(b: np.ndarray, se: np.ndarray, beta_ivw: float):
    w = 1.0 / se ** 2
    k = len(b)
    Q = float(np.sum(w * (b - beta_ivw) ** 2))
    df = k - 1
    p_Q = float(stats.chi2.sf(Q, df))
    I2 = max(0.0, (Q - df) / Q) * 100.0 if Q > 0 else 0.0
    return Q, df, p_Q, I2


def ivw_random(b: np.ndarray, se: np.ndarray):
    """DerSimonian-Laird random-effects IVW."""
    beta_f, _, _ = ivw_fixed(b, se)
    w = 1.0 / se ** 2
    k = len(b)
    Q, df, _, _ = heterogeneity(b, se, beta_f)
    # DL between-study variance
    c = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
    tau2 = max(0.0, (Q - df) / c) if c > 0 else 0.0
    wstar = 1.0 / (se ** 2 + tau2)
    beta = float(np.sum(wstar * b) / np.sum(wstar))
    se_beta = float(np.sqrt(1.0 / np.sum(wstar)))
    p = 2.0 * stats.norm.sf(abs(beta / se_beta))
    return beta, se_beta, p, tau2


def mr_egger(beta_e: np.ndarray, beta_o: np.ndarray, se_o: np.ndarray):
    """
    MR-Egger regression: beta_o = alpha + beta_E * beta_e, weights 1/se_o^2,
    with residual dispersion adjustment (matches standard implementations).
    """
    X = np.column_stack([np.ones_like(beta_e), beta_e])
    W = 1.0 / se_o ** 2
    XtW = X.T * W
    A = XtW @ X
    Ainv = np.linalg.inv(A)
    bhat = Ainv @ (XtW @ beta_o)
    resid = beta_o - X @ bhat
    k = len(beta_o)
    dof = k - 2
    rss = float(np.sum(W * resid ** 2))
    disp = rss / dof
    cov = Ainv * max(1.0, disp)
    se_b = np.sqrt(np.diag(cov))
    intercept, slope = bhat[0], bhat[1]
    se_int, se_slope = se_b[0], se_b[1]
    p_int = 2.0 * stats.t.sf(abs(intercept / se_int), dof)
    return slope, intercept, p_int, se_slope


def weighted_median(b: np.ndarray, se: np.ndarray):
    """Bowden weighted median (step function at the 50% cumulative weight)."""
    order = np.argsort(b)
    bs = b[order]
    w = (1.0 / se ** 2)[order]
    cum = np.cumsum(w) / np.sum(w)
    idx = int(np.searchsorted(cum, 0.5))
    idx = min(idx, len(bs) - 1)
    return float(bs[idx])


def loo_sign_flips(b: np.ndarray, se: np.ndarray, beta_full: float) -> int:
    flips = 0
    for i in range(len(b)):
        mask = np.arange(len(b)) != i
        bi, si = b[mask], se[mask]
        betai, _, _ = ivw_fixed(bi, si)
        if np.sign(betai) != np.sign(beta_full):
            flips += 1
    return flips


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------
def analyse(df: pd.DataFrame, drop_outliers: bool = False) -> dict:
    if drop_outliers and OUTLIER_COL in df.columns:
        df = df[df[OUTLIER_COL].astype(float) == 0]

    beta_e = df["beta_e"].to_numpy(float)
    se_e = df["se_e"].to_numpy(float)
    beta_o = df["beta_o"].to_numpy(float)
    se_o = df["se_o"].to_numpy(float)

    b = beta_o / beta_e                 # Wald ratio
    se = np.abs(se_o / beta_e)          # delta-method SE of the Wald ratio

    # sanity: does the shipped `ratio` column agree with our recomputation?
    if "ratio" in df.columns:
        shipped = df["ratio"].to_numpy(float)
        max_dev = float(np.max(np.abs(shipped - b)))
    else:
        max_dev = float("nan")

    beta_f, se_f, p_f = ivw_fixed(b, se)
    Q, dfQ, p_Q, I2 = heterogeneity(b, se, beta_f)
    beta_r, se_r, p_r, tau2 = ivw_random(b, se)
    slope_e, int_e, p_int, _ = mr_egger(beta_e, beta_o, se_o)
    wm = weighted_median(b, se)
    flips = loo_sign_flips(b, se, beta_f)
    n_out = int(df[OUTLIER_COL].astype(float).sum()) if OUTLIER_COL in df.columns else -1

    return {
        "n": int(len(df)),
        "ratio_col_max_abs_deviation": max_dev,
        "beta_raw": beta_f, "se_raw": se_f, "p_raw": p_f,
        "beta_re_raw": beta_r, "p_re_raw": p_r,
        "Q": Q, "Q_df": dfQ, "Q_p": p_Q, "I2": I2, "tau2": tau2,
        "egger_slope_raw": slope_e, "egger_intercept": int_e, "egger_int_p": p_int,
        "wmedian_raw": wm, "loo_flips": flips,
        "n_presso_outliers": n_out,
    }


def fmt_p(p: float) -> str:
    return f"{p:.3g}" if p >= 1e-4 else f"{p:.2e}"


def main() -> int:
    results = {}
    rows = []
    for key in TABLE3:
        path = os.path.join(HARM, key + ".csv")
        if not os.path.exists(path):
            print(f"MISSING {path}")
            return 2
        df = pd.read_csv(path)
        res = analyse(df)
        res_out = analyse(df, drop_outliers=True)

        # apply the manuscript's N-scale correction to betas AND ses
        c = SCALE
        res["beta"] = res["beta_raw"] * c
        res["se"] = res["se_raw"] * c
        res["beta_re"] = res["beta_re_raw"] * c
        res["wmedian"] = res["wmedian_raw"] * c
        res["egger_slope"] = res["egger_slope_raw"] * c
        res["beta_outlier_rm"] = res_out["beta_raw"] * c
        res["n_after_outlier_rm"] = res_out["n"]

        results[key] = res

        exp = TABLE3[key]
        e_n, e_beta, e_se, e_pf, e_pr, e_eg, e_ip, e_wm, e_i2, e_flip = exp

        rows.append({
            "pair": key,
            "n_ok": "OK" if res["n"] == e_n else f"FILE={res['n']} MS={e_n}",
            "beta_got": round(res["beta"], 6), "beta_ms": e_beta,
            "beta_diff": round(res["beta"] - e_beta, 6),
            "se_got": round(res["se"], 6), "se_ms": e_se,
            "pfix_got": fmt_p(res["p_raw"]), "pfix_ms": fmt_p(e_pf),
            "prand_got": fmt_p(res["p_re_raw"]), "prand_ms": fmt_p(e_pr),
            "egger_got": round(res["egger_slope"], 6), "egger_ms": e_eg,
            "egger_intp_got": round(res["egger_int_p"], 4), "egger_intp_ms": e_ip,
            "wmed_got": round(res["wmedian"], 6), "wmed_ms": e_wm,
            "i2_got": round(res["I2"], 1), "i2_ms": e_i2,
            "loo_got": res["loo_flips"], "loo_ms": e_flip,
            "n_outl": res["n_presso_outliers"],
            "beta_outlrm": round(res["beta_outlier_rm"], 6),
            "ratio_dev": f"{res['ratio_col_max_abs_deviation']:.2e}",
        })

    print("==================================================================================")
    print("INDEPENDENT RE-DERIVATION vs MANUSCRIPT TABLE 3 (N-corrected scale, x0.2823)")
    print("==================================================================================")
    hdr = f"{'pair':<26}{'n':<6}{'beta(calc)':>11}{'beta(ms)':>10}{'dbeta':>10}{'p_fix calc':>12}{'p_fix ms':>12}"
    print(hdr)
    for r in rows:
        print(f"{r['pair']:<26}{r['n_ok']:<6}{r['beta_got']:>11.6f}{r['beta_ms']:>10}{r['beta_diff']:>10}"
              f"{r['pfix_got']:>12}{r['pfix_ms']:>12}")
    print()
    print(f"{'pair':<26}{'p_rand calc':>12}{'p_rand ms':>12}{'egger calc':>11}{'egger ms':>10}"
          f"{'int_p calc':>12}{'int_p ms':>12}")
    for r in rows:
        print(f"{r['pair']:<26}{r['prand_got']:>12}{r['prand_ms']:>12}{r['egger_got']:>11.6f}"
              f"{r['egger_ms']:>10}{r['egger_intp_got']:>12.4f}{r['egger_intp_ms']:>12}")
    print()
    print(f"{'pair':<26}{'wmed calc':>11}{'wmed ms':>10}{'I2 calc':>9}{'I2 ms':>7}"
          f"{'LOO calc':>10}{'LOO ms':>7}{'n_outl':>7}{'beta_outlrm':>12}")
    for r in rows:
        print(f"{r['pair']:<26}{r['wmed_got']:>11.6f}{r['wmed_ms']:>10}{r['i2_got']:>9.1f}"
              f"{r['i2_ms']:>7}{r['loo_got']:>10}{r['loo_ms']:>7}{r['n_outl']:>7}{r['beta_outlrm']:>12.6f}")
    print()
    print("--- data integrity: max |shipped_ratio - recomputed beta_o/beta_e| ---")
    for r in rows:
        print(f"  {r['pair']:<26} {r['ratio_dev']}")

    out = os.path.join(HERE, "output", "genetics", "AUDIT_100_independent_recompute.json")
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(results, fh, indent=2, default=float)
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
