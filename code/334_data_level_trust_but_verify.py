#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""334 · Data-level trust-but-verify gate: re-derive reported statistics from
per-SNP instrument tables instead of comparing strings with a previous file.

Why this gate exists (2026-11-11 full-text final audit)
-------------------------------------------------------
code/149b checks that reported values in docs/08 are byte-identical to the
committed HEAD version. That is a *version-control* check: it proves
"numbers match the file", not "the numbers are correct". If the file itself
were wrong, 149b would still say SAFE. The audit brief explicitly requires
"数字与文件一致不等于数字对".

This gate recomputes the primary MR estimates from the harmonised per-SNP
instrument tables (code/output/genetics/harmonised/*.csv), which are one level
closer to the raw GWAS than the summary result CSVs, and compares them against
Table 2 in the manuscript. It also independently re-derives the four
closed-form quantities quoted in the abstract / Methods.

Independent means:
  * own IVW implementation (numpy, not the project's code/13 pipeline);
  * own Cochran Q / I^2 / random-effects SE / p-value;
  * own OR and CI from the recomputed beta/SE;
  * no import of any project result CSV except the per-SNP instrument table.

Hard constraint: this gate only READS. It never rewrites a reported number.
"""
from __future__ import annotations

import csv
import io
import math
import os
import re
import sys

import numpy as np
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HARM = os.path.join(ROOT, "code", "output", "genetics", "harmonised")
MANU = os.path.join(ROOT, "docs", "08_投稿稿_OA_BMD_fracture.md")

results: list[tuple[str, bool, str]] = []


def record(item: str, ok: bool, detail: str = "") -> None:
    results.append((item, bool(ok), detail))
    print(f"  [{'PASS' if ok else 'FAIL'}] {item}"
          + (f"  -> {detail}" if detail else ""))


def from_snps(pair: str) -> dict:
    """Own IVW re-derivation from per-SNP beta_e / beta_o / se_o."""
    path = os.path.join(HARM, pair + ".csv")
    if not os.path.isfile(path):
        raise FileNotFoundError(path)
    rows = list(csv.DictReader(io.open(path, encoding="utf-8")))
    be = np.array([float(r["beta_e"]) for r in rows])
    bo = np.array([float(r["beta_o"]) for r in rows])
    se_o = np.array([float(r["se_o"]) for r in rows])
    ratio = bo / be
    se = np.abs(se_o / be)
    w = 1.0 / se ** 2
    beta = float((w * ratio).sum() / w.sum())
    SE = float(1.0 / math.sqrt(w.sum()))
    Q = float((w * (ratio - beta) ** 2).sum())
    df = len(rows) - 1
    p = float(2 * stats.norm.sf(abs(beta / SE)))
    I2 = max(0.0, (Q - df) / Q) * 100 if Q > 0 else 0.0
    se_re = SE * math.sqrt(max(1.0, Q / df))
    p_re = float(2 * stats.norm.sf(abs(beta / se_re)))
    return dict(n=len(rows), beta=beta, SE=SE, p=p, Q=Q, df=df, I2=I2,
                se_re=se_re, p_re=p_re,
                OR=math.exp(abs(beta)),
                ci=(math.exp(beta - 1.96 * SE), math.exp(beta + 1.96 * SE)),
                # fracture rows are quoted per 1 SD LOWER BMD (protective)
                OR_low=math.exp(-beta),
                ci_low=(math.exp(-beta - 1.96 * SE), math.exp(-beta + 1.96 * SE)))


# Manuscript Table 2: (pair, n_snp, |beta|, SE, OR, CI_lo, CI_hi, I2)
TABLE2 = [
    ("BMD_total__fracture",       31, 0.337, 0.022, 1.401, 1.341, 1.464, 78),
    ("BMD_femur_neck__fracture",  25, 0.351, 0.030, 1.420, 1.338, 1.507, 83),
    ("BMD_total__OA_knee",        31, 0.078, 0.013, 1.081, 1.055, 1.108, 81),
    ("BMD_lumbar__OA_knee",       21, 0.095, 0.016, 1.099, 1.064, 1.135, 83),
    ("BMD_arms__OA_knee",         16, 0.075, 0.016, 1.078, 1.045, 1.112, 86),
    ("BMD_leg__OA_knee",          24, 0.060, 0.014, 1.062, 1.033, 1.092, 85),
    ("BMD_femur_neck__OA_knee",   25, 0.058, 0.017, 1.060, 1.025, 1.095, 68),
]

# Reverse pairs: quoted OR is present for OA knee -> fracture only.
REVERSE = ["OA_knee__BMD_total", "OA_knee__BMD_femur_neck", "OA_knee__fracture"]


def main() -> int:
    print("=" * 74)
    print("334 · DATA-LEVEL TRUST-BUT-VERIFY (own IVW from per-SNP instruments)")
    print("=" * 74)

    cache = {}
    for pair, *_ in TABLE2:
        try:
            cache[pair] = from_snps(pair)
        except Exception as e:
            record(f"per-SNP table readable: {pair}", False, str(e))
    for pair in REVERSE:
        try:
            cache[pair] = from_snps(pair)
        except Exception as e:
            record(f"per-SNP table readable: {pair}", False, str(e))

    print("\n(A) Table 2 primary estimates re-derived from per-SNP instruments")
    for pair, n, beta, se, orr, ci_lo, ci_hi, i2 in TABLE2:
        if pair not in cache:
            continue
        r = cache[pair]
        is_fx = pair.startswith("BMD") and "fracture" in pair
        orv = r["OR_low"] if is_fx else r["OR"]
        clo, chi = r["ci_low"] if is_fx else r["ci"]
        checks = {
            "n_snp": r["n"] == n,
            "|beta|": abs(abs(r["beta"]) - beta) < 0.0006,
            "SE": abs(r["SE"] - se) < 0.0011,   # manuscript shows 3 dp
            "OR": abs(orv - orr) < 0.002,
            "CI_lo": abs(clo - ci_lo) < 0.003,
            "CI_hi": abs(chi - ci_hi) < 0.003,
            "I2": abs(r["I2"] - i2) < 1.0,
        }
        bad = [k for k, v in checks.items() if not v]
        record(f"Table 2 {pair}", not bad,
               f"recomputed beta={r['beta']:+.6f} SE={r['SE']:.6f} "
               f"OR={orv:.4f} CI=({clo:.4f},{chi:.4f}) I2={r['I2']:.1f}% "
               f"p={r['p']:.3e}" + (f" | MISMATCH {bad}" if bad else ""))

    print("\n(B) Random-effects convention: 2/5 forward BMD→OA remain significant")
    fwd = [p for p, *_ in TABLE2 if "OA_knee" in p]
    nsig = sum(1 for p in fwd if p in cache and cache[p]["p_re"] < 0.05)
    record("random-effects IVW: exactly 2 of 5 forward estimates p<0.05",
           nsig == 2, f"observed {nsig}/5; "
           + ", ".join(f"{p.split('__')[0]}:{cache[p]['p_re']:.4f}" for p in fwd if p in cache))

    print("\n(C) Closed-form quantities quoted in Methods / Abstract")
    K = 0.0964
    f = 1 / math.sqrt(K * (1 - K))
    record("scale factor 1/sqrt(K(1-K)) = 3.389 (K=0.0964)",
           abs(f - 3.389) < 0.001, f"recomputed {f:.6f}")

    b, se = -0.337, 0.022
    OR = math.exp(-b)
    clo, chi = math.exp(-b - 1.96 * se), math.exp(-b + 1.96 * se)
    record("abstract fracture OR 1.40 (95% CI 1.34–1.46)",
           abs(OR - 1.40) < 0.005 and abs(clo - 1.34) < 0.005 and abs(chi - 1.46) < 0.005,
           f"recomputed {OR:.4f} ({clo:.4f}–{chi:.4f})")

    d, sd = 0.0188, 0.0200
    z = d / sd
    pz = 2 * stats.norm.sf(abs(z))
    record("F1 hip−knee adjudication z=+0.94, p=0.35",
           abs(z - 0.94) < 0.01 and abs(pz - 0.35) < 0.01,
           f"recomputed z={z:.4f}, p={pz:.4f}")

    bb, ss = 0.097, 0.016
    hor, hcl, hch = math.exp(bb), math.exp(bb - 1.96 * ss), math.exp(bb + 1.96 * ss)
    record("hip OA OR 1.102 (1.068–1.136)",
           abs(hor - 1.102) < 0.002 and abs(hcl - 1.068) < 0.002 and abs(hch - 1.136) < 0.002,
           f"recomputed {hor:.4f} ({hcl:.4f}–{hch:.4f})")

    print("\n(D) Manuscript strings still present after compliance edits")
    t = io.open(MANU, encoding="utf-8").read()
    for k in ["OR 1.060–1.099", "3.7×10⁻⁵¹", "3.389", "9/31", "31,873"]:
        record(f"docs/08 retained reported value {k!r}", k in t)
    record("L49 sample-size typo corrected", "N≈31,873 each" in t and "N≈1.873" not in t)

    n_fail = sum(1 for _, ok, _ in results if not ok)
    print("\n" + "=" * 74)
    print(f"TOTAL: {len(results) - n_fail}/{len(results)} PASS   FAIL: {n_fail}")
    print("VERDICT:", "DATA-LEVEL VERIFIED" if n_fail == 0 else "VERIFICATION FAILED")
    print("=" * 74)
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
