# -*- coding: utf-8 -*-
"""
101_audit_variants.py
=====================
Resolve the three items that did NOT reproduce exactly in
100_audit_independent_recompute.py:

  A) Is the reported beta the RAW file-derived IVW, or RAW x 0.2823 ?
     (decides whether the "x0.2823 N-correction" described in Methods was
      actually applied to the reported numbers)
  B) Which random-effects convention reproduces the reported p_random ?
     (decides the "only 2 of 5 survived random effects" claim)
  C) Which weighted-median convention reproduces the reported column ?

Also verifies the OR / CI columns, Table 4 power values, and back-checks
SE_e against the closed-form standardised-effect SE 1/sqrt(2 N f(1-f)),
which settles whether the exposure betas were standardised with the TRUE
N=31,873 or with the placeholder N=400,000.
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
HARM = os.path.join(HERE, "output", "genetics", "harmonised")
SCALE = 0.2823

TABLE3 = {
    "BMD_total__fracture":     (31, -0.337, 0.022, 3.7e-51, 2.4e-12, -0.337, 0.558, -0.364, 78, 0),
    "BMD_femur_neck__fracture":(25, -0.351, 0.030, 7.7e-31, 2.3e-06, -0.343, 0.018, -0.365, 83, 0),
    "BMD_total__OA_knee":      (31, +0.023, 0.004, 3.9e-10, 6.2e-03, +0.023, 0.810, +0.019, 81, 0),
    "BMD_lumbar__OA_knee":     (21, +0.028, 0.005, 7.7e-09, 1.7e-02, +0.026, 0.295, +0.037, 83, 0),
    "BMD_arms__OA_knee":       (16, +0.022, 0.005, 1.7e-06, 7.0e-02, +0.022, 0.982, +0.031, 86, 0),
    "BMD_leg__OA_knee":        (24, +0.018, 0.004, 2.3e-05, 9.8e-02, +0.017, 0.815, +0.025, 85, 0),
    "BMD_femur_neck__OA_knee": (25, +0.017, 0.005, 6.4e-04, 5.5e-02, +0.017, 0.875, +0.015, 68, 0),
    "OA_knee__fracture":       (69, +0.092, 0.108, 0.395,   0.542,   +0.106, 0.476, +0.118, 49, 0),
    "OA_knee__BMD_total":      (69, +0.106, 0.108, 0.326,   0.481,   +0.123, 0.376, +0.112, 49, 0),
    "OA_knee__BMD_femur_neck": (69, +0.142, 0.108, 0.190,   0.272,   +0.164, 0.169, +0.113, 30, 0),
}

# reported ORs / CIs fromthe results text
REPORTED_OR = {
    "BMD_total__fracture": (1.40, 1.34, 1.46),
    "BMD_femur_neck__fracture": (1.42, 1.34, 1.51),
}
GEFOS_OR = (1.55, 1.48, 1.63)

TABLE4_POWER = {
    "BMD_total__fracture": (15.05, 0.063, 1.00),
    "BMD_femur_neck__fracture": (11.55, 0.085, 1.00),
    "BMD_total__OA_knee": (6.26, 0.010, 1.00),
    "BMD_lumbar__OA_knee": (5.77, 0.014, 1.00),
    "BMD_arms__OA_knee": (4.78, 0.013, 1.00),
    "BMD_leg__OA_knee": (4.23, 0.012, 0.99),
    "BMD_femur_neck__OA_knee": (3.42, 0.014, 0.93),
    "OA_knee__fracture": (0.85, 0.304, 0.14),
    "OA_knee__BMD_total": (0.98, 0.303, 0.17),
    "OA_knee__BMD_femur_neck": (1.31, 0.303, 0.26),
}


def ivw(b, se):
    w = 1.0 / se ** 2
    beta = float(np.sum(w * b) / np.sum(w))
    se_beta = float(np.sqrt(1.0 / np.sum(w)))
    return beta, se_beta


def Q_parts(b, se, beta_f):
    w = 1.0 / se ** 2
    k = len(b)
    Q = float(np.sum(w * (b - beta_f) ** 2))
    return Q, k - 1, float(np.sum(w)), float(np.sum(w ** 2))


def re_variants(b, se):
    beta_f, se_f = ivw(b, se)
    Q, df, sw, sw2 = Q_parts(b, se, beta_f)
    c = sw - sw2 / sw
    tau2 = max(0.0, (Q - df) / c) if c > 0 else 0.0
    out = {}

    # (1) DerSimonian-Laird: additive tau^2
    ws = 1.0 / (se ** 2 + tau2)
    b1 = float(np.sum(ws * b) / np.sum(ws)); s1 = float(np.sqrt(1.0 / np.sum(ws)))
    out["DL_additive"] = (b1, s1, 2 * stats.norm.sf(abs(b1 / s1)))

    # (2) multiplicative heterogeneity: same point estimate, SE inflated
    phi = max(1.0, Q / df)
    b2, s2 = beta_f, se_f * np.sqrt(phi)
    out["multiplicative"] = (b2, s2, 2 * stats.norm.sf(abs(b2 / s2)))

    # (3) Paule-Mandel style tau^2 (iterative)
    lo, hi = 0.0, max(10.0, 10 * se_f ** 2)
    for _ in range(200):
        mid = (lo + hi) / 2.0
        wm_ = 1.0 / (se ** 2 + mid)
        beta_m = float(np.sum(wm_ * b) / np.sum(wm_))
        QM = float(np.sum(wm_ * (b - beta_m) ** 2))
        if QM > df:
            lo = mid
        else:
            hi = mid
    tau2_pm = (lo + hi) / 2.0
    wp = 1.0 / (se ** 2 + tau2_pm)
    b3 = float(np.sum(wp * b) / np.sum(wp)); s3 = float(np.sqrt(1.0 / np.sum(wp)))
    out["PauleMandel"] = (b3, s3, 2 * stats.norm.sf(abs(b3 / s3)))

    return out, tau2


def wmed_variants(b, se):
    out = {}
    order = np.argsort(b)
    bs, ses = b[order], se[order]

    # (a) step weighted median, weights 1/se^2
    w = 1.0 / ses ** 2
    cum = np.cumsum(w) / np.sum(w)
    i = min(int(np.searchsorted(cum, 0.5)), len(bs) - 1)
    out["step_w1se2"] = float(bs[i])

    # (b) kernel-weighted interpolated version
    #     (Bowden 2016 smooth version used by TwoSampleMR)
    delta = np.concatenate([[bs[0] - bs[-1]], np.diff(bs)])  # not used directly
    # bandwidth after Bowden: rho = (max-min)/... use simple kernel approach
    rng = bs.max() - bs.min()
    if rng > 0:
        bw = 0.253 * rng  # rough
        # kernel density of cumulative weight
        def F(x):
            return float(np.sum(w * stats.norm.cdf((x - bs) / bw)) / np.sum(w))
        grid = np.linspace(bs.min() - bw, bs.max() + bw, 4001)
        vals = np.array([F(x) for x in grid])
        out["kernel"] = float(grid[np.argmin(np.abs(vals - 0.5))])
    else:
        out["kernel"] = float("nan")

    # (c) simple unweighted median
    out["simple_median"] = float(np.median(b))

    # (d) weighted median with weights = beta_e^2 / se_o^2 (equivalent to 1/se_ratio^2)
    out["step_w1se2"] = out["step_w1se2"]  # keep
    return out


def main():
    print("=" * 100)
    print("A) WHICH SCALE ARE THE REPORTED BETAS ON?")
    print("=" * 100)
    print(f"{'pair':<26}{'IVW raw':>12}{'raw x .2823':>13}{'MS beta':>10}"
          f"{'|raw-MS|':>11}{'|scaled-MS|':>13}{'verdict':>14}")
    scale_verdict = {}
    for key, exp in TABLE3.items():
        df = pd.read_csv(os.path.join(HARM, key + ".csv"))
        b = df["beta_o"].to_numpy(float) / df["beta_e"].to_numpy(float)
        se = np.abs(df["se_o"].to_numpy(float) / df["beta_e"].to_numpy(float))
        beta_raw, se_raw = ivw(b, se)
        ms_beta = exp[1]
        d_raw = abs(beta_raw - ms_beta)
        d_sc = abs(beta_raw * SCALE - ms_beta)
        v = "RAW" if d_raw < d_sc else "SCALED"
        scale_verdict[key] = ("raw", beta_raw, se_raw) if v == "RAW" else ("scaled", beta_raw * SCALE, se_raw * SCALE)
        print(f"{key:<26}{beta_raw:>12.6f}{beta_raw*SCALE:>13.6f}{ms_beta:>10}"
              f"{d_raw:>11.6f}{d_sc:>13.6f}{v:>14}")

    print()
    print("=" * 100)
    print("B) RANDOM-EFFECTS p-VALUE: which convention reproduces the manuscript?")
    print("=" * 100)
    print(f"{'pair':<26}{'MS p_rand':>12}{'DL':>11}{'multiplic.':>13}{'PauleM':>11}"
          f"{'closest':>11}")
    for key, exp in TABLE3.items():
        df = pd.read_csv(os.path.join(HARM, key + ".csv"))
        b = df["beta_o"].to_numpy(float) / df["beta_e"].to_numpy(float)
        se = np.abs(df["se_o"].to_numpy(float) / df["beta_e"].to_numpy(float))
        variants, _ = re_variants(b, se)
        ms_pr = exp[4]
        dists = {k: abs(np.log10(v[2]) - np.log10(ms_pr)) for k, v in variants.items()}
        closest = min(dists, key=dists.get)
        print(f"{key:<26}{ms_pr:>12.3g}{variants['DL_additive'][2]:>11.3g}"
              f"{variants['multiplicative'][2]:>13.3g}{variants['PauleMandel'][2]:>11.3g}"
              f"{closest:>11}")

    print()
    print("=" * 100)
    print("C) WEIGHTED MEDIAN: which convention reproduces the manuscript?")
    print("=" * 100)
    exp_map = {k: v[7] for k, v in TABLE3.items()}
    for key in TABLE3:
        df = pd.read_csv(os.path.join(HARM, key + ".csv"))
        b = df["beta_o"].to_numpy(float) / df["beta_e"].to_numpy(float)
        se = np.abs(df["se_o"].to_numpy(float) / df["beta_e"].to_numpy(float))
        v = wmed_variants(b, se)
        print(f"{key:<26} MS={exp_map[key]:>8.3f}   " +
              "  ".join(f"{k}={val:>8.4f}" for k, val in v.items()))

    print()
    print("=" * 100)
    print("D) OR / CI reproduction and the 8.4% / 9.6% GEFOS comparison")
    print("=" * 100)
    for key, (or_ms, lo_ms, hi_ms) in REPORTED_OR.items():
        df = pd.read_csv(os.path.join(HARM, key + ".csv"))
        b = df["beta_o"].to_numpy(float) / df["beta_e"].to_numpy(float)
        se = np.abs(df["se_o"].to_numpy(float) / df["beta_e"].to_numpy(float))
        beta, seb = ivw(b, se)
        or_calc = float(np.exp(-beta))          # per 1 SD LOWER BMD
        lo = float(np.exp(-beta - 1.96 * seb))
        hi = float(np.exp(-beta + 1.96 * seb))
        gLO, gHI, gUP = GEFOS_OR
        gap_ratio = (gHI - or_calc) / gHI * 100
        gap_log = (np.log(gHI) - np.log(or_calc)) / np.log(gHI) * 100
        print(f"{key:<26} OR calc={or_calc:.4f} CI=[{lo:.3f},{hi:.3f}] | MS OR={or_ms} "
              f"CI=[{lo_ms},{hi_ms}]")
        print(f"{'':<26} below GEFOS: ratio-scale {gap_ratio:6.2f}%  |  log-scale {gap_log:6.2f}%")

    print()
    print("=" * 100)
    print("E) TABLE 4 power: abs z, MDE = 2.8*SE, achieved power = Phi(|z|-1.96)")
    print("=" * 100)
    print(f"{'pair':<26}{'z calc':>9}{'z ms':>7}{'MDE calc':>10}{'MDE ms':>8}"
          f"{'pow calc':>10}{'pow ms':>8}")
    for key, (z_ms, mde_ms, pow_ms) in TABLE4_POWER.items():
        df = pd.read_csv(os.path.join(HARM, key + ".csv"))
        b = df["beta_o"].to_numpy(float) / df["beta_e"].to_numpy(float)
        se = np.abs(df["se_o"].to_numpy(float) / df["beta_e"].to_numpy(float))
        beta, seb = ivw(b, se)
        z = beta / seb
        mde = 2.8 * seb
        power = float(stats.norm.cdf(abs(z) - 1.96))
        print(f"{key:<26}{z:>9.3f}{z_ms:>7}{mde:>10.4f}{mde_ms:>8}"
              f"{power:>10.3f}{pow_ms:>8}")

    print()
    print("=" * 100)
    print("F) WHICH N WERE THE EXPOSURE betas STANDARDISED WITH?")
    print("   closed form for standardised effect: SE = 1/sqrt(2 N f(1-f))")
    print("=" * 100)
    for key in ["BMD_total__OA_knee", "BMD_total__fracture", "OA_knee__BMD_total"]:
        df = pd.read_csv(os.path.join(HARM, key + ".csv")).head(6)
        # identify which column is the BMD one
        is_forward = key.startswith("BMD_")
        secol = "se_e" if is_forward else "se_o"
        print(f"-- {key} (checking column {secol})")
        for _, row in df.iterrows():
            se_obs = float(row[secol])
            f = float(row["eaf"])
            f = min(f, 1.0 - f)
            for N in (31873, 400000):
                se_theory = 1.0 / np.sqrt(2 * N * f * (1 - f))
                ratio = se_obs / se_theory
                flag = "  <== matches" if abs(ratio - 1) < 0.02 else ""
                print(f"     se_obs={se_obs:.6f}  N={N:<7} theory={se_theory:.6f}  ratio={ratio:.4f}{flag}")


if __name__ == "__main__":
    main()
