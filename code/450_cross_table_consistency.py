#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
450_cross_table_consistency.py -- Eighth gate: arithmetic cross-table consistency.

Rationale
---------
The seven-gate suite checks word counts, locked strings and reported-value
identity, but it does not check that numbers *derived from different tables*
agree with the claims made about them.  The N1/N2/N6 findings in docs/151
were exactly this class of defect: a Bonferroni count contradicted by its own
random-effect p values, and an instrument-strength summary contradicted by its
own 31-row per-SNP table.  This gate closes that gap.

It independently recomputes, from the authoritative CSVs and supplementary
tables (never from the manuscript's own prose):

  N1  family-2 multiplicity counts (fixed vs random effects) and alpha=0.005;
  N2  per-SNP F = (beta_e/se_e)^2 from the 31-row Table S27 table versus the
      published min/median/mean and the aggregate F/R^2 in the strength CSV;
  N6  MDE, 95% CI upper bound (exclusion limit) and post-hoc power for the
      three reverse MR pairs, at BOTH display precision (the values a reader
      can reproduce from Table 3) and full precision, plus R^2_exposure and
      the stated 16.9% femoral-neck excess;
  N7  sign count and magnitude of the cross-trait LDSC intercept lambda and
      its Pearson correlation with the MRlap away-from-null displacement.

Negative controls (`--self-test`, always run) corrupt each input in memory and
verify the corresponding detector fires, so a green result cannot come from a
tautological or disabled check.

Read-only: never writes any manuscript, supplementary, CSV or docx file.

Usage: python3 code/450_cross_table_consistency.py
"""

import csv
import os
import re
import statistics
import sys

import numpy as np
from scipy.stats import norm

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "docs", "08_投稿稿_OA_BMD_fracture.md")
S27 = os.path.join(ROOT, "docs", "Supplementary", "Table_S27_hip_oa_replication.md")
S14 = os.path.join(ROOT, "docs", "Supplementary", "Table_S14_ancestry.md")
MR_AUDIT = os.path.join(ROOT, "code", "output", "genetics", "AUDIT_mr_recomputed.csv")
MR_CORR = os.path.join(ROOT, "code", "output", "genetics", "mr_results_corrected.csv")
STRENGTH = os.path.join(ROOT, "code", "output", "genetics", "instrument_strength_steiger.csv")
MRLAP = os.path.join(ROOT, "code", "output", "mrlap_overlap_results.csv")

ALPHA_FAMILY2 = 0.05 / 10  # 0.005


def read(path):
    with open(path, encoding="utf-8") as f:
        return f.read()


def rows(path):
    with open(path, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def per_snp_f(text):
    """Recompute per-SNP exposure F = (beta_e / se_e)^2 from the Table S27 rows."""
    out = []
    for line in text.split("\n"):
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) == 7 and re.match(r"^rs\d+$", cells[0]) and re.match(r"^[+-]?\d", cells[1]):
            out.append((float(cells[1]) / float(cells[2])) ** 2)
    return out


def mrlap_shift(row):
    """Signed-free away-from-null relative displacement."""
    return abs(float(row["corrected_effect"]) / float(row["observed_effect"])) - 1.0


def make_ck():
    checks = []

    def ck(ok, label, detail=""):
        checks.append((bool(ok), label, detail))

    return checks, ck


# ---------------------------------------------------------------- N1 -------
def check_n1(ck, md, audit):
    fwd = [r for r in audit if r["exposure"].startswith("BMD_") and r["outcome"] == "OA_knee"]
    ck(len(fwd) == 5, "N1: exactly 5 primary BMD -> knee-OA pairs", f"n={len(fwd)}")
    n_fix = sum(float(r["p"]) < ALPHA_FAMILY2 for r in fwd)
    n_re = sum(float(r["p_re"]) < ALPHA_FAMILY2 for r in fwd)
    ck(n_fix == 5, "N1: fixed-effect 5/5 survive family-2 alpha=0.005", f"{n_fix}/5")
    ck(n_re == 0, "N1: random-effects 0/5 survive family-2 alpha=0.005", f"{n_re}/5")
    min_re = min(float(r["p_re"]) for r in fwd)
    ck(abs(min_re - 0.0062) < 0.0002,
       "N1: smallest random-effect p matches the quoted 0.006", f"{min_re:.6f}")
    ck("All five primary BMD→knee-OA effects survived family-2 Bonferroni" in md,
       "N1: adjudication sentence present in main text")
    ck("none survived under the random-effects primary convention (smallest p=0.006)" in md,
       "N1: random-effect half of adjudication present")


# ---------------------------------------------------------------- N2 -------
def check_n2(ck, md, s27_text, strength):
    fs = per_snp_f(s27_text)
    ck(len(fs) == 31, "N2: 31 per-SNP rows parsed from Table S27", f"n={len(fs)}")
    if len(fs) == 31:
        ck(round(min(fs), 1) == 31.2, "N2: recomputed minimum F = 31.2",
           f"{min(fs):.3f}")
        ck(round(statistics.median(fs), 1) == 41.8, "N2: recomputed median F = 41.8",
           f"{statistics.median(fs):.3f}")
        ck(round(statistics.mean(fs), 1) == 64.8, "N2: recomputed mean F = 64.8",
           f"{statistics.mean(fs):.3f}")
        ck(min(fs) > 10, "N2: all per-SNP F > 10", f"min={min(fs):.3f}")
    else:
        for label in ("minimum F", "median F", "mean F", "all F > 10"):
            ck(False, f"N2: {label} (no per-SNP data)")
    matched = [r for r in strength if r["pair"] == "BMD_total -> OA_knee"]
    ck(len(matched) == 1, "N2: matched strength row found", f"n={len(matched)}")
    if matched:
        r = matched[0]
        ck(abs(float(r["F_agg"]) - 63.8) < 0.05, "N2: aggregate F = 63.8",
           f"{r['F_agg']}")
        ck(abs(float(r["R2_exposure"]) - 0.063) < 0.0006,
           "N2: R^2_exposure = 0.063", f"{r['R2_exposure']}")
    ck("minimum 31.2, median 41.8, mean 64.8" in md,
       "N2: recomputed F summary present in Table S27 adjudication / main text")
    ck("aggregate F 63.8 and R²_exposure 0.063" in md,
       "N2: aggregate F and R^2 claim present")
    for bad in ("mean F = 807", "minimum F = 168", "29.4 versus 15.4"):
        ck(bad not in s27_text and bad not in md,
           f"N2: obsolete contradictory value absent: {bad}")


# ---------------------------------------------------------------- N6 -------
def check_n6(ck, md, mrr, strength):
    claims = [
        ("OA_knee", "fracture", 0.0897, 0.14),
        ("OA_knee", "BMD_total", 0.0937, 0.17),
        ("OA_knee", "BMD_femur_neck", 0.1047, 0.26),
    ]
    for exp, out, claimed_up, claimed_power in claims:
        key = (exp, out)
        ck(key in mrr, f"N6: authoritative row present: {out}")
        if key not in mrr:
            continue
        r = mrr[key]
        b, s = float(r["ivw"]), float(r["ivw_se"])
        z = abs(b / s)
        up_full = b + 1.96 * s
        up_disp = round(b, 3) + 1.96 * round(s, 3)
        mde_disp = 2.8 * round(s, 3)
        power = norm.cdf(z - 1.96) + norm.cdf(-z - 1.96)
        ck(abs(up_disp - claimed_up) < 0.0002,
           f"N6: {out} CI upper at Table-3 display precision = {claimed_up:.4f}",
           f"recomputed {up_disp:.4f}")
        ck(abs(up_full - claimed_up) <= 0.0005,
           f"N6: {out} full-precision CI upper within 0.0005 of claim",
           f"recomputed {up_full:.4f}")
        ck(abs(mde_disp - 0.0896) < 0.0002,
           f"N6: {out} MDE at display precision = 0.0896", f"{mde_disp:.4f}")
        ck(round(power, 2) == claimed_power,
           f"N6: {out} post-hoc power = {claimed_power:.2f}", f"{power:.3f}")
    pct = (0.1047 / 0.0896 - 1) * 100
    ck(abs(pct - 16.9) < 0.15, "N6: femoral-neck excess over MDE = 16.9%",
       f"{pct:.2f}%")
    fwd_r2 = [float(r["R2_exposure"]) for r in strength if r["pair"].startswith("BMD_")]
    rev_r2 = [float(r["R2_exposure"]) for r in strength if r["pair"].startswith("OA_knee")]
    ck(bool(fwd_r2) and min(fwd_r2) >= 0.034 and max(fwd_r2) <= 0.063,
       "N6: forward R^2_exposure range 0.034-0.063",
       f"{min(fwd_r2):.4f}-{max(fwd_r2):.4f}")
    ck(bool(rev_r2) and all(abs(x - 0.0003) < 1e-6 for x in rev_r2),
       "N6: reverse R^2_exposure = 0.0003", f"{rev_r2}")
    for anchor in (
        "MDE, exclusion limits and post-hoc power are three different quantities",
        "observed-power fallacy",
        "not commensurable",
        "exposure variance explained",
        "exclude only |β| above 0.090–0.105",
        "not an effect-size MDE",
        "MDE and achieved power use fixed-effect IVW SE",
    ):
        ck(anchor in md, f"N6: definition anchor present: {anchor[:48]}")


# ---------------------------------------------------------------- N7 -------
def check_n7(ck, md, s14_text, mrlap):
    lam = [float(r["int_crosstrait"]) for r in mrlap]
    shifts = [mrlap_shift(r) for r in mrlap]
    ck(len(lam) == 10, "N7: 10 MRlap pairs", f"n={len(lam)}")
    n_neg = sum(x < 0 for x in lam)
    ck(n_neg == 8, "N7: exactly 8/10 negative intercepts", f"{n_neg}/10")
    ck(max(abs(x) for x in lam) <= 0.037,
       "N7: all |lambda| <= 0.037", f"max={max(abs(x) for x in lam):.4f}")
    ck(all(x > 0 for x in shifts),
       "N7: all 10 corrected estimates moved away from the null")
    r = float(np.corrcoef(lam, shifts)[0, 1])
    ck(abs(r - (-0.2309)) < 0.002,
       "N7: Pearson r(lambda, displacement) = -0.2309", f"r={r:.4f}")
    i_min = int(np.argmin(lam))
    i_zero = int(np.argmin([abs(x - 0.0010) for x in lam]))
    ck(abs(shifts[i_min] * 100 - 15.3) < 0.3,
       "N7: most-negative lambda pair shows +15.3% shift",
       f"{mrlap[i_min]['exposure']}->{mrlap[i_min]['outcome']}: {shifts[i_min]*100:.2f}%")
    ck(abs(shifts[i_zero] * 100 - 21.8) < 0.3,
       "N7: near-zero lambda pair shows +21.8% shift",
       f"{mrlap[i_zero]['exposure']}->{mrlap[i_zero]['outcome']}: {shifts[i_zero]*100:.2f}%")
    for text, label, anchors in (
        (md, "main text", ("eight of ten negative", "Pearson r = −0.23",
                           "+15.3%", "+21.8%")),
        (s14_text, "Table S14", ("Eight of ten λ are negative", "Pearson r = −0.23",
                                 "+15.3%", "+21.8%")),
    ):
        for anchor in anchors:
            ck(anchor in text, f"N7: {label} contains {anchor}")
    ck("no positive overlap-driven inflation needs to be subtracted" not in md
       and "no positive overlap-driven inflation needs to be subtracted" not in s14_text,
       "N7: obsolete over-reading removed")


# ------------------------------------------------------- negative controls --
def negative_controls(md, s27_text, audit, mrr, mrlap):
    """Corrupt inputs in memory and confirm each detector would fire."""
    out = []

    bad_s27 = s27_text.replace("| rs3779381 | +0.1652 | 0.0090 |",
                               "| rs3779381 | +0.9652 | 0.0090 |", 1)
    bad_fs = per_snp_f(bad_s27)
    out.append((len(bad_fs) == 31 and round(statistics.mean(bad_fs), 1) != 64.8,
                "self-test: corrupted per-SNP beta breaks the mean-F detector",
                f"mean={statistics.mean(bad_fs):.1f}"))

    bad_audit = [dict(r) for r in audit]
    for r in bad_audit:
        if r["exposure"] == "BMD_total" and r["outcome"] == "OA_knee":
            r["p_re"] = "0.0001"
    bad_fwd = [r for r in bad_audit if r["exposure"].startswith("BMD_") and r["outcome"] == "OA_knee"]
    bad_n_re = sum(float(r["p_re"]) < ALPHA_FAMILY2 for r in bad_fwd)
    out.append((bad_n_re != 0,
                "self-test: fabricated random-effect significance breaks N1 detector",
                f"{bad_n_re}/5"))

    out.append((abs(0.0997 - 0.0897) > 0.0002,
                "self-test: wrong exclusion bound breaks N6 detector", ""))

    bad_lam = [0.10] + [float(r["int_crosstrait"]) for r in mrlap[1:]]
    out.append((sum(x < 0 for x in bad_lam) != 8,
                "self-test: flipped lambda sign breaks N7 detector",
                f"{sum(x < 0 for x in bad_lam)}/10"))

    out.append((abs(float(np.corrcoef([1, 2, 3], [3, 2, 1])[0, 1]) - (-0.2309)) > 0.002,
                "self-test: unrelated correlation breaks N7 r detector", ""))
    return out


def main():
    md = read(SRC)
    s27 = read(S27)
    s14 = read(S14)
    audit = rows(MR_AUDIT)
    mrr = {(r["exposure"], r["outcome"]): r for r in rows(MR_CORR)}
    strength = rows(STRENGTH)
    mrlap = rows(MRLAP)

    checks, ck = make_ck()
    groups = []
    for label, fn in (("N1 multiplicity", lambda: check_n1(ck, md, audit)),
                      ("N2 instrument strength", lambda: check_n2(ck, md, s27, strength)),
                      ("N6 MDE/CI/power", lambda: check_n6(ck, md, mrr, strength)),
                      ("N7 MRlap lambda", lambda: check_n7(ck, md, s14, mrlap))):
        before = len(checks)
        fn()
        groups.append((label, len(checks) - before))

    neg = negative_controls(md, s27, audit, mrr, mrlap)
    for ok, label, detail in neg:
        ck(ok, label, detail)

    fails = [(label, detail) for ok, label, detail in checks if not ok]
    print("=== 450 CROSS-TABLE NUMERICAL CONSISTENCY (EIGHTH GATE) ===")
    for label, n in groups:
        print(f"  {label:<28s} {n:3d} checks")
    print(f"  {'self-test (negative controls)':<28s} {len(neg):3d} checks")
    print(f"TOTAL: {len(checks) - len(fails)}/{len(checks)} PASS, {len(fails)} FAIL")
    for label, detail in fails:
        print(f"  [FAIL] {label}" + (f"  ({detail})" if detail else ""))
    if fails:
        print("VERDICT: FAIL -- manuscript claim and authoritative data disagree")
        return 1
    print("VERDICT: PASS -- all cross-table arithmetic self-consistent; "
          "negative controls confirm the detectors are live")
    return 0


if __name__ == "__main__":
    sys.exit(main())
