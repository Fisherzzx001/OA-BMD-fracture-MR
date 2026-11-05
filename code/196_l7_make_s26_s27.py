#!/usr/bin/env python3
"""Generate Supplementary Tables S26 (F1: prospective hip-OA MR) and S27 (F2: geometry arm)."""
import json, math, os, csv

D = "G:/OA_B_PROJECT_cache/l7"
PROJ = "C:/Users/Administrator/WorkBuddy/2026-09-08-22-17-40/OA_B_PROJECT"
SUP = os.path.join(PROJ, "docs", "Supplementary")
os.makedirs(SUP, exist_ok=True)

f1 = json.load(open(os.path.join(D, "l7_f1_final.json")))
f2 = json.load(open(os.path.join(D, "l7_f2_hsm.json")))

def fmt(x, n=4):
    return f"{x:+.{n}f}" if isinstance(x, float) else str(x)

# ---------------- S26 ----------------
h, k = f1["hip"], f1["knee"]
adj = f1["adjudication"]
lines26 = []
lines26.append("## Supplementary Table S26. Prospective hip-osteoarthritis replication (falsifying test F1): two-sample MR of total-body BMD on hip OA, and the conjunctive adjudication against the knee-OA estimate\n")
lines26.append("F1 was pre-registered as the prospective replication of the BMD->OA finding at the hip, with a conjunctive criterion: the hip effect must be positive AND of knee-like magnitude. The hip outcome is the Hatzikotoulas 2025 mixed-ancestry hip-OA GWAS (GCST90566798; 95,395 European cases, 1,011,691 European controls plus East Asian, African, South Asian, Hispanic and admixed groups; maximum per-variant 97,328 cases and 1,055,375 controls), and the knee outcome is the same study's knee-OA GWAS (GCST90566800; maximum per-variant 172,256 cases and 1,144,248 controls) — the file from which the manuscript's own knee-OA estimate was derived, so that the two rows below are produced by an identical pipeline on an identical instrument set (31 total-body BMD instruments, Qian 2025).\n")
lines26.append("| Outcome | n_snp | IVW beta (fixed) | SE | p | IVW p (random) | OR per SD (95% CI) | MR-Egger beta (SE) | Egger p | Egger intercept p | Weighted median (95% CI) | Cochran Q | I2 | LOO sign flips |")
lines26.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
for nm, r in (("hip OA (GCST90566798)", h), ("knee OA (GCST90566800)", k)):
    iv, eg, wm = r["ivw_fixed"], r["egger"], r["wmedian"]
    lines26.append(f"| {nm} | {r['n']} | {iv['beta']:+.3f} | {iv['se']:.3f} | {iv['p']:.1e} | "
                   f"{r['ivw_random']['p']:.1e} | {iv['OR']:.3f} ({iv['lo']:.3f}-{iv['hi']:.3f}) | "
                   f"{eg['slope']:+.3f} ({eg['se']:.3f}) | {eg['p']:.1e} | {eg['p_int']:.2f} | "
                   f"{wm['beta']:+.3f} ({wm['lo']:+.3f} to {wm['hi']:+.3f}) | {r['Q']:.1f} ({r['Q_df']} df) | "
                   f"{r['I2']:.0f}% | {r['loo_sign_flips']} |")
lines26.append("")
lines26.append(f"**Adjudication.** Difference in IVW log-odds (hip minus knee) = {adj['diff']:+.4f}; SE under the conservative assumption of independent outcomes = {adj['se_indep']:.4f}; z = {adj['z']:+.2f}, two-sided p = {adj['p_two']:.2f}. Ratio of hip to knee log-odds = {adj['ratio']:.2f}. The 95% intervals overlap from {max(h['ivw_fixed']['lo'], k['ivw_fixed']['lo']):.4f} to {min(h['ivw_fixed']['hi'], k['ivw_fixed']['hi']):.4f}. Instrument strength: mean F = {h['F_mean']:.0f}, minimum F = {h['F_min']:.0f}; SNP-level R2 for the exposure {h['R2_exposure']:.1f} versus {h['R2_outcome']:.1f} for the outcome, Steiger-type directionality {h['steiger']:+.3f}. Leaving out the {h['n_outliers']} variants with absolute studentised residual above 2 leaves the estimate positive and strengthened (beta = {h['after_outlier']['beta']:+.4f}, SE {h['after_outlier']['se']:.4f}, p = {h['after_outlier']['p']:.1e}); leave-one-out produced no sign reversal (range {h['loo_min']:+.4f} to {h['loo_max']:+.4f}).")
lines26.append("")
lines26.append("**Verification against the published knee-OA row.** Re-deriving the knee-OA estimate from the raw GWAS Catalog file reproduces the manuscript's Table 2 values: IVW beta +0.078 (this table +0.078), OR 1.081 (1.055-1.108), I2 81% (80.8%) and weighted median +0.065 (+0.027 to +0.177) (+0.066, +0.027 to +0.174). The MR-Egger intercept is non-significant in both (this table p = 0.63; manuscript p = 0.810).")
lines26.append("")
lines26.append("**Decision.** F1 passes its conjunctive criterion. The prospective same-pipeline hip estimate is positive (p = 6.6e-10) and statistically indistinguishable from the knee estimate (z = 0.94, p = 0.35); the hip estimate lies 24% above the knee estimate, within the magnitude band spanned by the only prior same-pipeline comparison (Hartley et al., hip 27% below knee). P6 is therefore restored to its pre-registered directional form.")
lines26.append("")
# per-SNP table
lines26.append("**Per-instrument values (hip OA outcome; the 31 total-body BMD instruments).** beta_e and se_e are the total-body BMD exposure effects per allele; hip_lnOR is the hip-OA log-odds on the BMD-increasing allele with its standard error.")
lines26.append("")
lines26.append("| SNP | beta_e (BMD) | se_e | hip ln(OR) | se | hip OR | N (cases+controls) |")
lines26.append("| --- | --- | --- | --- | --- | --- | --- |")
per = json.load(open(os.path.join(D, "l7_results.json")))["per_snp"]
for r in sorted(per, key=lambda x: x["rs"]):
    lines26.append(f"| {r['rs']} | {r['beta_e']:+.4f} | {r['se_e']:.4f} | {r['hip_lnOR']:+.4f} | "
                   f"{r['hip_se']:.4f} | {r['hip_OR']:.3f} | {int(r['hip_N']):,} |")
lines26.append("")
s26 = "\n".join(lines26) + "\n"
open(os.path.join(SUP, "Table_S26_f1_hip_replication.md"), "w", encoding="utf-8").write(s26)

# ---------------- S27 ----------------
lines27 = []
lines27.append("## Supplementary Table S27. Geometry arm of the mass-versus-geometry decomposition (falsifying test F2): two-sample MR of ten hip-shape modes on hip OA and knee OA, plus a conditional BMD-plus-hip-shape multivariable MR\n")
lines27.append("Areal BMD fuses bone mass, size and geometry, so the BMD->OA effect could in principle be transmitted by bone geometry rather than by mass or formation. Height indexes stature but not bone width or geometry, leaving this split to F2. The ten statistical hip-shape modes (HSM1-HSM10) of Faber et al. 2025 (PMID 39574169; deposited in the GWAS Catalog as GCST90482704-GCST90482713) are derived from DXA imaging and encode geometry only, with no density component; they are used here as pure-geometry instruments in the same pipeline and with the same hip and knee outcomes as Table S26. Instruments are the genome-wide-significant lead variants of each mode (P < 5e-9, 193 distinct positions); exposure effects are the authors' reported betas with direction, and standard errors are derived from the reported P values. Hip-shape mode 5 yielded only one instrument and is not analysable. A conditional BMD-plus-hip-shape multivariable MR extending this screen is reported in the second part of this table.")
lines27.append("")
lines27.append("### Part A. Single-variable two-sample MR of each hip-shape mode on hip OA and knee OA\n")
lines27.append("| Mode | Outcome | n_snp | IVW beta | SE | p (fixed) | OR (95% CI) | p (random) | Cochran Q | I2 | MR-Egger beta | Egger intercept p | Weighted median | Mean F |")
lines27.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
def mode_label(acc):
    return f"HSM{int(acc[-2:]) - 3}"

for m in f2["modes"]:
    for onm, key in (("hip OA", "hip"), ("knee OA", "knee")):
        r = m[key]
        if "ivw" not in r:
            lines27.append(f"| {mode_label(m['acc'])} | {onm} | {r['n']} | not analysable (insufficient instruments) | | | | | | | | | | |")
            continue
        i, e = r["ivw"], r["egger"]
        lines27.append(f"| {mode_label(m['acc'])} | {onm} | {r['n']} | {i['beta']:+.4f} | {i['se']:.4f} | {i['p']:.3g} | "
                       f"{i['OR']:.3f} ({i['lo']:.3f}-{i['hi']:.3f}) | {i['p_random']:.3g} | {r['Q']:.1f} ({r['Q_df']} df) | "
                       f"{r['I2']:.0f}% | {e['beta']:+.4f} | {e['p_int']:.3g} | {r['wmedian']['beta']:+.4f} | {r['F_mean']:.0f} |")
lines27.append("")
ag = f2["aggregate"]
# Bonferroni under the manuscript's random-effects convention
for nm in ("hip", "knee"):
    ps = [m[nm]["ivw"]["p_random"] for m in f2["modes"] if "ivw" in m[nm]]
    ag[nm]["any_significant_bonf_random"] = any(p < 0.005 for p in ps)
    ag[nm]["n_bonf_fixed"] = sum(1 for m in f2["modes"] if "ivw" in m[nm] and m[nm]["ivw"]["p"] < 0.005)
lines27.append(f"**Pooled across modes.** Hip OA: pooled beta {ag['hip']['beta']:+.4f} (SE {ag['hip']['se']:.4f}), OR {ag['hip']['OR']:.3f} ({ag['hip']['lo']:.3f}-{ag['hip']['hi']:.3f}); across-mode heterogeneity Q = {ag['hip']['Q']:.1f} ({ag['hip']['Q_df']} df), I2 = {ag['hip']['I2']:.0f}%. Knee OA: pooled beta {ag['knee']['beta']:+.4f} (SE {ag['knee']['se']:.4f}), OR {ag['knee']['OR']:.3f} ({ag['knee']['lo']:.3f}-{ag['knee']['hi']:.3f}); Q = {ag['knee']['Q']:.1f} ({ag['knee']['Q_df']} df), I2 = {ag['knee']['I2']:.0f}%. The largest single-mode estimate is {ag['hip']['max_abs_beta']:.3f} (hip) and {ag['knee']['max_abs_beta']:.3f} (knee) in absolute log-odds, against 0.097 and 0.078 for total-body BMD on the same outcomes; the smallest per-mode P value is {ag['hip']['min_p']:.2g} (hip) and {ag['knee']['min_p']:.2g} (knee), and under the random-effects convention used throughout the manuscript no mode reaches the Bonferroni threshold of 0.005 across the ten modes (hip {ag['hip']['any_significant_bonf_random']}, knee {ag['knee']['any_significant_bonf_random']}); {ag['hip']['n_bonf_fixed']} of ten modes (hip) and {ag['knee']['n_bonf_fixed']} of ten (knee) would qualify on fixed-effect P values alone, which is why the random-effects convention governs.")
lines27.append("")
lines27.append("**Interpretation of Part A.** The pooled geometry arm is an order of magnitude smaller than the BMD arm and carries no consistent direction (hip and knee pooled estimates have opposite signs), so the BMD->OA effect is not accounted for by the hip-geometry axes tested. Individual modes, however, are heterogeneous (I2 75-95%) and several show significant MR-Egger intercepts, so the per-mode estimates are pleiotropy-affected and are not read causally; the single-variable design is informative rather than a clean falsification, and it is reported as such. The conditional BMD-plus-hip-shape MVMR in Part B now supplies the formal test and reaches the same qualitative conclusion for knee OA (and is inconclusive for hip OA). The literal volumetric-instrument form of F2 (areal BMD replaced by a volume- or area-conditioned BMD instrument) remains untested because no public summary-statistics file for a volumetric BMD trait was found in the GWAS Catalog; the manuscript's height conditioning (9.5-33.2% attenuation) and BMI conditioning already cover the stature and adiposity axes.")

# ---------------- S27 Part B: conditional multivariable MR (BMD + hip-shape) ----------------
l8 = json.load(open(os.path.join(D, "L8_F2_conditional_mvmr.json")))
lines27.append("")
lines27.append("### Part B. Conditional multivariable MR (BMD + hip-shape): formal geometry-conditioning of falsifying test F2")
lines27.append("")
lines27.append("**Design.** A shared-instrument MVMR is infeasible because the 31 total-body BMD instruments and the HSM lead variants are disjoint (zero overlap; see Table S17), so we condition within the HSM-instrument frame: each HSM lead variant is used as the instrument, and GEFOS total-body BMD beta/SE (tidy_femur_total.tsv, 754,224 SNPs) are retrieved on the same variants. Of the 193 HSM lead variants, 59 carried a usable GEFOS BMD effect and entered the conditional MR; n_snp in Part B is therefore the GEFOS-overlapping subset of each mode's lead variants and is smaller than n_snp in Part A. IVW-MVMR follows Sanderson & Windmeijer (2015) and estimates theta_BMD|HSM (BMD's direct effect on OA holding hip-shape geometry constant) and theta_HSM|BMD (geometry's direct effect holding BMD constant), with only outcome SEs used for weighting. Egger-MVMR gives the horizontal-pleiotropy intercept; the Sanderson-Smith conditional F assesses instrument strength for the added BMD coefficient. Modes with <4 overlapping instruments are reported as insufficient (need >=4 for stable identification).")
lines27.append("")
lines27.append("| Mode | Outcome | n_snp | theta_BMD|HSM (SE) | p_BMD | theta_HSM|BMD (SE) | p_HSM | cond. F | Egger-int p |")
lines27.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
def _pf(p):
    if p == 0 or p < 1e-15:
        return "<1e-15"
    if p < 1e-3:
        return f"{p:.1e}"
    return f"{p:.3f}"
for mk in [f"GCST904827{i:02d}" for i in range(4, 14)]:
    hsm = int(mk[-2:]) - 3
    for onm, key in (("hip OA", "hip"), ("knee OA", "knee")):
        r = l8[key][mk]
        if "theta_BMD_given_HSM" not in r:
            lines27.append(f"| HSM{hsm} | {onm} | {r['n']} | insufficient (need >=4) | | | | | |")
            continue
        lines27.append(
            f"| HSM{hsm} | {onm} | {r['n']} | {r['theta_BMD_given_HSM']:+.3f} ({r['se_BMD']:.3f}) | "
            f"{_pf(r['p_BMD'])} | {r['theta_HSM_given_BMD']:+.3f} ({r['se_HSM']:.3f}) | {_pf(r['p_HSM'])} | "
            f"{r['cond_F_BMD_given_HSM']:.2f} | {_pf(r['p_intercept'])} |")
lines27.append("")
lines27.append("**Reading of Part B.** In the knee-OA models - where the marginal BMD->OA effect is robustly positive (Table 2) - the conditional BMD effect given geometry (theta_BMD|HSM) is positive in four of the five estimable modes (HSM2 +0.133, HSM4 +0.283, HSM8 +1.065, HSM10 +0.539) and negative but non-significant in HSM1 (-0.127, p=0.55); the two best-powered modes (HSM8, HSM10, n=5 and 4 with the largest conditional F) are both clearly positive and significant (p=0.0013 and 3.4e-06). Because conditioning on hip-shape geometry leaves the BMD->knee-OA effect intact and positive, bone geometry does not absorb (mediate) the BMD signal - this formally upgrades F2 from a single-variable null screen to a conditional test. In the hip-OA models the conditional estimates are heterogeneous and sign-flipping (HSM6: marginal BMD->OA +0.19 -> conditional -0.92; HSM1 and HSM6 negative, HSM2/HSM4/HSM10 positive), giving no consistent direction; this is consistent with the generally weaker hip-OA evidence and leaves the hip-OA geometry conclusion inconclusive rather than confirmed.")
lines27.append("")
lines27.append("**Why the magnitudes are not quoted as a single coefficient.** The HSM lead variants are strong instruments for geometry but weak for the added BMD coefficient: the conditional F is far below the conventional 10 (range 0.11-11.87; only HSM10 reaches ~11.9, all others <4 and most <1), per-mode n is 4-9, I2 is 48-98%, and Egger-MVMR intercepts are significant in several modes (e.g. HSM1/HSM2 hip, p<1e-7). The conditional BMD coefficients are therefore direction-indicative, not magnitude-precise; cross-mode pooling is not attempted because the per-mode estimates are direction-heterogeneous and the instrument sets are not independent samples. theta_HSM|BMD (geometry's direct OA effect given BMD) is significant in several modes but direction-heterogeneous (e.g. HSM2 hip +0.134, p=4.8e-4; HSM6 hip -0.754, p<1e-15; HSM8 knee +0.177, p=6.2e-4), corroborating the Part A finding that no hip-shape axis reproducibly drives OA.")
lines27.append("")
lines27.append("**Conclusion for F2.** The formal conditional MVMR corroborates the single-variable geometry arm: for knee OA the BMD effect survives geometry conditioning, so bone geometry is not the transmission path; for hip OA the conditional test is underpowered and inconclusive. The literal volumetric-instrument form of F2 (areal BMD replaced by a volume- or area-conditioned BMD instrument) remains untested because no public summary-statistics file for a volumetric BMD trait was found in the GWAS Catalog; the manuscript's height conditioning (9.5-33.2% attenuation) and BMI conditioning already cover the stature and adiposity axes.")
s27 = "\n".join(lines27) + "\n"
open(os.path.join(SUP, "Table_S27_f2_geometry_arm.md"), "w", encoding="utf-8").write(s27)

print("Wrote:")
print(" ", os.path.join(SUP, "Table_S26_f1_hip_replication.md"))
print(" ", os.path.join(SUP, "Table_S27_f2_geometry_arm.md"))
print("S26 words:", len(s26.split()), " S27 words:", len(s27.split()))
