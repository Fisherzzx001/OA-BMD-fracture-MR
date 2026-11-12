#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""P1-2 — Report Burgess overlap bias on the TWO denominators side by side,
and quantify the MRlap-corrected positive control vs GEFOS on a comparable scale.

docs/148 P1-2: the manuscript body (L55) reports only the small 8.4%-of-fracture
denominator ("worst-case <=4.5%") while S14 already discloses the ~12% on the
70%-of-BMD basis.  This script makes the dual-denominator arithmetic explicit and
computes the MRlap relative displacement transferred onto the Table 2 IVW scale
(MRlap uses 48/34 instruments and its own exposure scale, so its raw corrected
beta is NOT the Table 2 corrected OR -- we transfer the SCALE-FREE ratio instead).

Self-checks:
  * overlap counts reproduce 22330/31873=70.1% and 22330/264973=8.4%
  * Burgess bias = overlap_fraction * |r_g| (Burgess 2016 approximation)
  * MRlap ratio corrected/observed > 1 (moves away from null) for both controls
  * transferred positive-control OR moves TOWARD GEFOS 1.55 (the honest finding)
Read-only; writes a JSON audit record.  Dry by default, --apply to emit.
"""
import json, sys, os, numpy as np, pandas as pd

APPLY = "--apply" in sys.argv
HERE = os.path.dirname(os.path.abspath(__file__))
MR = os.path.join(HERE, "output", "mrlap_overlap_results.csv")
OUT = os.path.join(HERE, "output", "genetics", "P1_2_overlap_dual_denominator.json")

# Verified UKB overlap counts (manuscript S14 / docs/148 P1-2)
N_BMD = 31873          # 100% UK Biobank
N_FRX = 264973         # GEFOS any-fracture
N_OVR = 22330          # 1,832 cases + 20,498 controls in UKB
FRAC_OF_BMD = N_OVR / N_BMD
FRAC_OF_FRX = N_OVR / N_FRX
# Symmetric ("geometric-mean") overlap: the Burgess overlap.prop is a symmetric
# between-sample quantity when the two GWAS differ greatly in size. sqrt(f1*f2).
OVERLAP_GEO = float(np.sqrt(FRAC_OF_BMD * FRAC_OF_FRX))
# per-pair r_g from LDSC (mrlap csv 'rg' col is cross-trait rg for the pair)
rg_bmd_frx_total = -0.486795020610176
rg_bmd_frx_femur = -0.508643754850951
rg_bmd_oa_total = 0.162084331881118

m = pd.read_csv(MR)
def row(ex, out):
    return m[(m.exposure == ex) & (m.outcome == out)].iloc[0]

checks = []
def chk(name, cond, detail=""):
    checks.append((name, bool(cond), detail))

chk("overlap 70.1% of BMD", abs(FRAC_OF_BMD - 0.701) < 0.002, f"{FRAC_OF_BMD:.4f}")
chk("overlap 8.4% of fracture", abs(FRAC_OF_FRX - 0.084) < 0.002, f"{FRAC_OF_FRX:.4f}")

# ---- Burgess bias on both denominators -----------------------------------
def burgess(frac, rg):
    return frac * abs(rg)

burg = {
    "BMD_total->kneeOA": {
        "overlap_frac_pct": 2.6, "r_g": rg_bmd_oa_total,
        "bias_on_2.6pct_pct": burgess(0.026, rg_bmd_oa_total) * 100,
    },
    "BMD_total->fracture": {
        "overlap_of_fracture_pct": FRAC_OF_FRX * 100,
        "overlap_of_BMD_pct": FRAC_OF_BMD * 100,
        "overlap_geometric_mean_pct": OVERLAP_GEO * 100,
        "r_g": rg_bmd_frx_total,
        "bias_on_8.4pct_pct": burgess(FRAC_OF_FRX, rg_bmd_frx_total) * 100,
        "bias_on_70pct_pct": burgess(FRAC_OF_BMD, rg_bmd_frx_total) * 100,
        "bias_on_geometric_mean_pct": burgess(OVERLAP_GEO, rg_bmd_frx_total) * 100,
    },
    "BMD_femur_neck->fracture": {
        "overlap_of_fracture_pct": FRAC_OF_FRX * 100,
        "overlap_of_BMD_pct": FRAC_OF_BMD * 100,
        "overlap_geometric_mean_pct": OVERLAP_GEO * 100,
        "r_g": rg_bmd_frx_femur,
        "bias_on_8.4pct_pct": burgess(FRAC_OF_FRX, rg_bmd_frx_femur) * 100,
        "bias_on_70pct_pct": burgess(FRAC_OF_BMD, rg_bmd_frx_femur) * 100,
        "bias_on_geometric_mean_pct": burgess(OVERLAP_GEO, rg_bmd_frx_femur) * 100,
    },
}
# The 8.4%-of-fracture basis uses the SAME 22330 people; the 70%-of-BMD basis
# is the correct "fraction of the exposure sample that is shared" for the
# positive control.  Burgess 2016 bias acts in the direction of r_g (here away
# from null: r_g<0, MR beta<0 => overlap biases the protective estimate down).
chk("total Burgess 8.4% ~ 4.1%", abs(burg["BMD_total->fracture"]["bias_on_8.4pct_pct"] - 4.09) < 0.15,
    f"{burg['BMD_total->fracture']['bias_on_8.4pct_pct']:.2f}%")
chk("total Burgess 70% ~ 34% (NOT 12%)", 33 < burg["BMD_total->fracture"]["bias_on_70pct_pct"] < 35,
    f"{burg['BMD_total->fracture']['bias_on_70pct_pct']:.2f}%")
chk("geometric-mean overlap ~ 24%", 23 < OVERLAP_GEO * 100 < 25,
    f"{OVERLAP_GEO*100:.2f}%")
chk("geometric-mean Burgess ~ 12% = the mislabeled '70%' figure",
    11 < burg["BMD_total->fracture"]["bias_on_geometric_mean_pct"] < 13,
    f"{burg['BMD_total->fracture']['bias_on_geometric_mean_pct']:.2f}%")

# ---- MRlap relative displacement transferred to Table 2 IVW ---------------
# Table 2 primary IVW betas (verified N=31,873 scale); per SD LOWER BMD the OR
# is exp(-beta) because beta is for higher BMD and is negative.
PRIMARY = {"BMD_total": -0.337, "BMD_femur_neck": -0.351}
GEFOS_OR, GEFOS_LOG = 1.55, float(np.log(1.55))
transfer = {}
for ex in PRIMARY:
    r = row(ex, "fracture")
    ratio = r.corrected_effect / r.observed_effect
    b0 = PRIMARY[ex]
    b1 = b0 * ratio
    or0, or1 = float(np.exp(-b0)), float(np.exp(-b1))
    transfer[ex] = {
        "mrlap_n_iv": int(r.m_IVs),
        "mrlap_observed": float(r.observed_effect),
        "mrlap_corrected": float(r.corrected_effect),
        "relative_displacement_pct": float((ratio - 1) * 100),
        "egger_intercept_p": float(r.egger_intercept_p),
        "table2_ivw_beta_perSD_higher": b0,
        "table2_or_perSD_lower": or0,
        "transferred_beta": float(b1),
        "transferred_or_perSD_lower": or1,
        "gap_to_GEFOS_OR_before_pct": (GEFOS_OR - or0) / GEFOS_OR * 100,
        "gap_to_GEFOS_OR_after_pct": (GEFOS_OR - or1) / GEFOS_OR * 100,
        "gap_to_GEFOS_log_before_pct": (GEFOS_LOG - (-b0)) / GEFOS_LOG * 100,
        "gap_to_GEFOS_log_after_pct": (GEFOS_LOG - (-b1)) / GEFOS_LOG * 100,
    }
    chk(f"{ex} MRlap moves away from null", ratio > 1.0, f"ratio={ratio:.4f}")
    chk(f"{ex} transfer moves toward GEFOS", transfer[ex]["gap_to_GEFOS_log_after_pct"] < transfer[ex]["gap_to_GEFOS_log_before_pct"],
        f"{transfer[ex]['gap_to_GEFOS_log_before_pct']:.1f}% -> {transfer[ex]['gap_to_GEFOS_log_after_pct']:.1f}%")

result = {
    "task": "docs/148 P1-2",
    "overlap_counts": {"N_BMD": N_BMD, "N_fracture": N_FRX, "N_overlap_UKB": N_OVR,
                       "frac_of_BMD_pct": FRAC_OF_BMD * 100, "frac_of_fracture_pct": FRAC_OF_FRX * 100},
    "burgess_bias": burg,
    "mrlap_transferred_positive_control": transfer,
    "caveats": [
        "Burgess bias = overlap fraction x |r_g| is an approximation (Burgess 2016); "
        "the 70%-of-BMD denominator is the fraction of the EXPOSURE sample shared, "
        "appropriate for asking whether the two-sample assumption holds for the "
        "BMD instrument -> fracture path; the 8.4%-of-fracture denominator is the "
        "fraction of the OUTCOME sample shared. Both describe the SAME 22,330 people.",
        "LABEL ERROR CORRECTED: manuscript S14 and docs/148 P1-2 previously wrote "
        "'up to ~12% on the 70%-of-BMD basis'. But 70% x |r_g|=0.4868 = 34.1%, NOT 12%. "
        "The ~12% figure actually corresponds to the symmetric geometric-mean overlap "
        "sqrt(0.701 x 0.084)=24.3% x |r_g| ~ 11.8-12.4% (code/121 used overlap~0.23, "
        "consistent with this). The three denominators (8.4% / geometric-mean 24.3% / "
        "70%) now each carry their correctly-derived bias (4.1% / ~12% / 34%). The "
        "literal Burgess formula is olsbias x overlap.prop x (1/F); the manuscript's "
        "overlap x r_g omits the 1/F attenuation and is a conservative upper bound, "
        "so even the 34% is a worst-case ceiling, not an expected bias.",
        "MRlap re-clumps to 48 (total) / 34 (femur-neck) instruments, not the primary "
        "31/25, and reports effects on its own exposure scale; its raw corrected beta "
        "is therefore NOT the Table 2 corrected OR. We transfer the scale-free ratio "
        "corrected/observed onto the Table 2 IVW beta, which assumes the RELATIVE "
        "overlap displacement is approximately instrument-set invariant.",
        "This transfer is a sensitivity analysis, not a substitute point estimate; "
        "the primary Table 2 OR (1.40/1.42) remains the reported value.",
    ],
    "checks": [{"name": n, "pass": p, "detail": d} for n, p, d in checks],
}

print("=" * 74)
print("P1-2  Overlap dual denominator + MRlap positive-control transfer")
print("=" * 74)
print(f"UKB overlap {N_OVR}: {FRAC_OF_BMD*100:.1f}% of BMD | {FRAC_OF_FRX*100:.1f}% of fracture")
for k, v in burg.items():
    print(f"\n{k}: r_g={v['r_g']:+.4f}")
    for kk in v:
        if kk.startswith("bias_") or kk.endswith("_pct"):
            print(f"   {kk}: {v[kk]:.2f}%")
print("\n[label check] Burgess bias on 70%-of-BMD = "
      f"{burg['BMD_total->fracture']['bias_on_70pct_pct']:.1f}% (NOT the 12% written in S14); "
      f"12% is the geometric-mean-overlap basis ({burg['BMD_total->fracture']['bias_on_geometric_mean_pct']:.1f}%)")
for ex, t in transfer.items():
    print(f"\n{ex}->fracture  MRlap n={t['mrlap_n_iv']} (primary 31/25), Egger-int p={t['egger_intercept_p']:.3f}")
    print(f"   MRlap {t['mrlap_observed']:+.4f}->{t['mrlap_corrected']:+.4f} ({t['relative_displacement_pct']:+.1f}%)")
    print(f"   Table2 OR {t['table2_or_perSD_lower']:.3f} -> transferred {t['transferred_or_perSD_lower']:.3f} (GEFOS 1.55)")
    print(f"   gap(log-odds) {t['gap_to_GEFOS_log_before_pct']:.1f}% -> {t['gap_to_GEFOS_log_after_pct']:.1f}%")
nfail = sum(1 for _, p, _ in checks if not p)
print(f"\nchecks: {len(checks)-nfail}/{len(checks)} passed")
for n, p, d in checks:
    print(f"  [{'OK' if p else 'FAIL'}] {n} {d}")
if nfail:
    print("FAILED -- not written"); sys.exit(1)
if APPLY:
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(result, open(OUT, "w", encoding="utf-8"), indent=2, ensure_ascii=False)
    print(f"written {OUT}")
else:
    print("DRY-RUN OK -- re-run with --apply")
