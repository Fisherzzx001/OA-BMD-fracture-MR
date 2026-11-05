#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Second compression pass on docs/08 (prose only; numbers untouched).
Asserts each old string matches exactly once before replacing."""
import io

PATH = "docs/08_投稿稿_OA_BMD_fracture.md"

PAIRS = [
    # ---- Results ----
    ("The two definitions also agreed at the variant level: across 1,980 genome-wide-significant knee-OA loci, 94.4% had concordant signed Z with self-reported OA, confirming they relate to bone in the same genetic direction.",
     "The two definitions also agreed at the variant level (94.4% concordant signed Z across 1,980 genome-wide-significant knee-OA loci), confirming they relate to bone in the same genetic direction."),

    ("The estimate held under random effects (p=2.4×10⁻¹²), the MR-Egger intercept was null (p=0.558), weighted median agreed (−0.364), and leave-one-out produced no sign reversal.",
     "The estimate held under random effects (p=2.4×10⁻¹²); MR-Egger intercept null (p=0.558); weighted median agreed (−0.364); leave-one-out produced no sign reversal."),

    ("The conservative 1-Mb set gave β=−0.366 (OR 1.44); these agree with the published GEFOS estimate (OR 1.55, 1.48–1.63) [10].",
     "The conservative 1-Mb set gave β=−0.366 (OR 1.44), agreeing with GEFOS (OR 1.55, 1.48–1.63) [10]."),

    ("MR supported a positive effect of BMD on knee OA across all five exposures (Table 2): total-body β=+0.078 (p=3.9×10⁻¹⁰), lumbar-spine +0.095 (p=7.7×10⁻⁹), arms +0.075 (p=1.7×10⁻⁶), legs +0.060 (p=2.3×10⁻⁵), femoral-neck +0.058 (p=6.4×10⁻⁴) — OR 1.060–1.099 per SD higher BMD, against OR 1.40 (1.20–1.63) per SD in the largest one-sample estimate [22], a 3.6–5.8-fold log-odds difference (adjudicated §4).",
     "MR supported a positive BMD→knee-OA effect at all five sites (Table 2: β=+0.058 to +0.095, p=6.4×10⁻⁴ to 3.9×10⁻¹⁰; OR 1.060–1.099 per SD), against OR 1.40 (1.20–1.63) per SD in the largest one-sample estimate [22] (3.6–5.8-fold log-odds difference, adjudicated §4)."),

    ("Heterogeneity was high (I²=68–86%); under random-effects IVW only two of five exposures remained significant (total-body p=0.006, lumbar-spine p=0.017; a second random-effects convention additionally retains arms, p=0.033), so the count depends on the random-effects convention and we base the causal reading on estimator concordance rather than fixed-effect p-values. All five effects remained positively causal after MR-PRESSO outlier removal, with four of five strengthening. The positive control was comparably heterogeneous (I² 78% and 83%), indicating high I² is a property of the BMD instrument sets rather than an OA-specific IVW failure; a multi-estimator robustness overview of all ten pairs is in Supplementary Fig. S11.",
     "Heterogeneity was high (I²=68–86%); under random-effects IVW only two of five remained significant (total-body p=0.006, lumbar-spine p=0.017; a second convention additionally retains arms, p=0.033), so the count depends on the random-effects convention and the causal reading rests on estimator concordance. All five remained positively causal after MR-PRESSO outlier removal (four strengthened). The positive control was comparably heterogeneous (I² 78% and 83%), indicating high I² is a BMD-instrument-set property not an OA-specific IVW failure; a multi-estimator overview is in Fig. S11."),

    ("MR-Egger intercepts were all p>0.24 and weighted medians matched the IVW sign at all five sites, with bootstrap 95% CIs excluding the null at three of five; with 16–31 instruments the intercept test has limited power, so a non-significant intercept excludes detectable unbalanced pleiotropy rather than pleiotropy itself.",
     "MR-Egger intercepts were all p>0.24 and weighted medians matched the IVW sign at all five sites, with bootstrap 95% CIs excluding the null at three of five; with 16–31 instruments the intercept test has limited power, so a non-significant intercept excludes detectable unbalanced pleiotropy."),

    ("Leave-one-out produced no sign reversals and directions agreed with r_g (+0.12 to +0.25).",
     "Leave-one-out produced no sign reversals; directions agreed with r_g (+0.12 to +0.25)."),

    ("Removing the MR-PRESSO outliers flagged for these pairs moved two point estimates substantially (total-body BMD to +0.051, fracture to +0.059) without altering their non-significance, so the nulls reflect both low power and sensitivity to instrument selection.",
     "Removing the flagged MR-PRESSO outliers moved two estimates substantially (total-body BMD to +0.051, fracture to +0.059) without altering non-significance, so the nulls reflect both low power and instrument-selection sensitivity."),

    ("The Steiger directionality test supported BMD as the causally prior trait for every forward pair (z = 3.05–6.11, one-sided p = 4.8×10⁻⁴ to 5.0×10⁻¹⁰) and was non-significant for the three reverse pairs (|z|<0.1, p≈0.48–0.49).",
     "The Steiger test supported BMD as the causally prior trait for every forward pair (z = 3.05–6.11, p = 4.8×10⁻⁴ to 5.0×10⁻¹⁰) and was non-significant for the three reverse pairs (|z|<0.1, p≈0.48–0.49)."),

    ("The MR-PRESSO distortion test was significant for 9/10 pairs (expected from the widespread horizontal pleiotropy common in complex-trait MR [8]) and identified 0–3 outliers per pair; removing outliers did not reverse any direction, and strict-clumped estimates reproduced the conservative 1-Mb clumping.",
     "The MR-PRESSO distortion test was significant for 9/10 pairs (expected from widespread horizontal pleiotropy [8]) and identified 0–3 outliers per pair; removing outliers reversed no direction and strict-clumped estimates reproduced the conservative 1-Mb clumping."),

    ("The MRlap correction moved all ten estimates further from the null (median +20.5%; difference p<0.05 for 9/10), the direction expected when the towards-null weak-instrument and winner's-curse components dominate the away-from-null overlap component, and it reversed no direction or significance conclusion; the cross-trait intercept of the correction was |λ| ≤ 0.037 for all ten pairs.",
     "The MRlap correction moved all ten estimates further from the null (median +20.5%; difference p<0.05 for 9/10), the direction expected when towards-null weak-instrument and winner's-curse components dominate the away-from-null overlap component, reversing no direction or significance; the cross-trait intercept was |λ| ≤ 0.037."),

    ("At WNT16, SOX6 and GDF5, loci nominated as shared effectors, colocalisation with knee OA was not supported (PP.H4 ≤0.005 at each; e.g. WNT16 knee-OA z=0.94, p=0.35, PP.H4=0.002; GDF5 knee-OA z=−13.7 but no BMD genome-wide signal within ±500 kb, PP.H4 0.003–0.017).",
     "At WNT16, SOX6 and GDF5, nominated shared effectors, colocalisation with knee OA was not supported (PP.H4 ≤0.005 at each; e.g. WNT16 knee-OA z=0.94, p=0.35, PP.H4=0.002; GDF5 knee-OA z=−13.7 with no BMD signal within ±500 kb, PP.H4 0.003–0.017)."),

    ("The knee-OA null reflects power — under the MR estimate the expected knee-OA z at the strongest instrument is 0.85, and reaching genome-wide significance would require ~5×10⁷ participants (~42× the present sample).",
     "The knee-OA null reflects power: the expected knee-OA z at the strongest instrument is 0.85, and genome-wide significance would require ~5×10⁷ participants (~42× the sample)."),

    ("How much of the BMD→OA effect is skeletal size? Conditioning the BMD→OA effect on height attenuated it by 9.5–33.2% across the two tabulated univariable baselines (13.5–33.2% against the primary Table 2 baselines and 9.5–21.2% against the MVMR-set baselines; primary n_h=100/EUR; Supplementary Tables S9 and S10 and Supplementary Fig. S10), leaving 67–91% independent of height.",
     "How much of the BMD→OA effect is skeletal size? Height conditioning attenuated it by 9.5–33.2% across the two univariable baselines (13.5–33.2% vs Table 2 and 9.5–21.2% vs the MVMR set; Tables S9–S10, Fig. S10), leaving 67–91% independent of height."),

    ("The positive control confirmed that the established BMD→fracture effect survives the same conditioning: BMD→fracture moved only −1.7% to −3.0% across the same two baselines (−2.8% against the MVMR-set baseline), the control being insensitive to baseline choice.",
     "The positive control confirmed the BMD→fracture effect survives conditioning (moved −1.7% to −3.0%; −2.8% vs MVMR set), insensitive to baseline choice."),

    ("A parallel BMI-conditioned MVMR (Table S16) left θ_BMD positive and fixed-effect-significant at all five sites (total-body +0.069, p=2.6×10⁻⁸), attenuating −3.7% to 46.0% (the estimate rose at the femoral-neck site) — below the height attenuation at four of five sites; only the non-weight-bearing arms site exceeded it (46.0% vs 13.5%), losing significance under random effects.",
     "A parallel BMI-conditioned MVMR (Table S16) left θ_BMD positive and significant at all five sites (total-body +0.069, p=2.6×10⁻⁸), attenuating −3.7% to 46.0% — below height attenuation at four of five sites; only non-weight-bearing arms exceeded it (46.0% vs 13.5%), losing significance under random effects."),

    ("Because height indexes stature but not bone width or geometry [35], this residual is not identified as bone mass; the mass-versus-geometry split is resolved only by test F2.",
     "Because height indexes stature not bone geometry [35], this residual is not bone mass; the mass–geometry split is resolved only by F2."),

    # ---- Discussion ----
    ("Jointly analysing 13 skeletal phenotypes gave a consistent architecture: BMD sites form a strongly inter-correlated cluster, OA correlates positively with every BMD site, fracture negatively with every BMD site, and the two OA definitions correlate strongly. On pleiotropy-robust evidence, genetically higher BMD raises knee-OA risk; the reverse (OA→BMD) and OA→fracture estimates were non-significant and underpowered.",
     "Jointly analysing 13 skeletal phenotypes gave a consistent architecture: BMD sites inter-correlate strongly, OA correlates positively with every BMD site, fracture negatively with every BMD site, and the two OA definitions correlate strongly. On pleiotropy-robust evidence, genetically higher BMD raises knee-OA risk; the reverse (OA→BMD) and OA→fracture estimates were non-significant and underpowered."),

    ("These contradict the directional form of P6 (motivated by Zhou's protective finding; P6 and F1 in Table S11).",
     "These contradict the directional form of P6 (Table S11)."),

    ("We therefore restate P6 as a magnitude claim (a hip effect positive and at least as large as the knee's would fail P6), deviating from the pre-registered directional form prompted by contradictory hip evidence, scoring only prospective hip estimates.",
     "We therefore restate P6 as a magnitude claim (a hip effect positive and ≥ the knee's would fail P6), deviating from the pre-registered directional form given contradictory hip evidence and scoring only prospective hip estimates."),

    ("We adjudicate rather than cite selectively: the hip evidence disagrees with itself (protective [5], harmful [22,36]), and decisively for F1 the criterion is conjunctive — the hip effect must be positive AND of knee-like magnitude.",
     "We adjudicate rather than cite selectively: the hip evidence disagrees with itself (protective [5], harmful [22,36]), and for F1 the criterion is conjunctive — the hip effect must be positive AND knee-like in magnitude."),

    ("A systemic, polygenic bone-forming set-point readable as areal BMD would protect at the whole-bone endpoint and raise OA risk at the joint — the Systemic Osteogenic Set-point (SOS [34]), operationalised as areal BMD as a heuristic, not a fitted model.",
     "A systemic, polygenic bone-forming set-point readable as areal BMD would protect at the whole-bone endpoint and raise OA risk at the joint — the Systemic Osteogenic Set-point (SOS [34]), a heuristic not a fitted model."),

    ("If the causal quantity is volumetric density, an areal-BMD estimate is diluted by skeletal size, so volumetric-BMD (or area-conditioned) instruments should yield a BMD→OA effect at least as large as the areal one — the quantitative form of F2, which areal data cannot test.",
     "If the causal quantity is volumetric density, an areal-BMD estimate is diluted by size, so volumetric-BMD (or area-conditioned) instruments should yield a BMD→OA effect at least as large as the areal one — the quantitative form of F2, which areal data cannot test."),

    ("The reading generates six predictions (P1–P6) and ten falsifying tests (F1–F10), set out below and in Table S11.",
     "The reading generates six predictions (P1–P6) and ten falsifying tests (F1–F10) below and in Table S11."),

    ("**P2** systemic not local-mechanical — systemic clause supported (non-weight-bearing arms β = +0.075, p = 1.7×10⁻⁶), ordering clause under-powered (arms − femoral neck +0.017, one-sided p = 0.228).",
     "**P2** systemic not local-mechanical — systemic clause supported (non-weight-bearing arms β=+0.075, p=1.7×10⁻⁶); ordering under-powered (arms−femoral neck +0.017, p=0.228)."),

    ("**P3** diffuse polygenic — 0/31 loci colocalise (PP.H4 ≥ 0.8); direction-consistent, non-discriminating.",
     "**P3** diffuse polygenic — 0/31 loci colocalise (PP.H4 ≥ 0.8); direction-consistent."),

    ("**P4** size a separable partial component — height attenuation 9.5–33.2% and BMI attenuation −3.7% to 46.0% (fracture controls −1.7% to −3.0% and 0.9%); 67–91% residual unresolved into mass versus geometry (F2).",
     "**P4** size a separable partial component — height attenuation 9.5–33.2% and BMI −3.7% to 46.0% (fracture controls −1.7% to −3.0% and 0.9%); 67–91% residual unresolved into mass versus geometry (F2)."),

    ("**F4** eQTL colocalisation (cross-trait arm executed: fracture 9/31 versus knee-OA 0/31 loci at PP.H4≥0.8, paired Wilcoxon p=0.0043, §3; eQTLGen whole-blood arm executed, p=0.001156, but blood is neither bone nor cartilage, so the tissue-literal arm remains unexecuted);",
     "**F4** eQTL colocalisation (cross-trait arm executed: fracture 9/31 vs knee-OA 0/31 at PP.H4≥0.8, paired Wilcoxon p=0.0043, §3; eQTLGen whole-blood arm executed, p=0.001156, but blood is neither bone nor cartilage, so the tissue-literal arm is unexecuted);"),

    ("**F5** osteophyte-stratified outcome — F1–F3 and F5 prospective, requiring external data; F5 pre-registered here as the immediate next step and not performed.",
     "**F5** osteophyte-stratified outcome — F1–F3 and F5 prospective (external data); F5 is the immediate next step, not performed."),

    ("**F6** cross-site homogeneity (executed: Q = 6.92, 4 df, p = 0.14, τ² = 0) and **F7** joint-locus exclusion (executed: GDF5 not a BMD instrument at any site; excluding WNT16 or SOX6 changes no estimate's direction or significance) both passed.",
     "**F6** cross-site homogeneity (Q=6.92, 4 df, p=0.14, τ²=0) and **F7** joint-locus exclusion (GDF5 not a BMD instrument anywhere; excluding WNT16/SOX6 changes no estimate) both passed."),

    ("**F9** site-gradient ordering (executed: arms − femoral neck +0.017, one-sided p = 0.228).",
     "**F9** site-gradient ordering (arms−femoral neck +0.017, p=0.228)."),

    ("Our BMD→OA finding is consistent with Alhassan et al. [16], who identified BMD among exposures with MR support for a causal OA association; the ten-site design extends that signal with site specificity and a quantified size decomposition.",
     "Our BMD→OA finding is consistent with Alhassan et al. [16] and the ten-site design extends it with site specificity and a quantified size decomposition."),

    # ---- Methods ----
    ("**LDSC.** Genetic correlation (r_g) was estimated by LDSC [6] natively in Python from the Bulik-Sullivan implementation, on a common set of 690,150 HapMap3 SNPs with 1000 Genomes Phase 3 European LD scores, using delete-one-block jackknife (100 blocks) for standard errors. Palindromic SNPs were excluded; alleles were harmonised to one reference so every matrix cell uses the same variants. Both Bonferroni (α=0.05/78) and Benjamini–Hochberg FDR were applied across 78 phenotype pairs.",
     "**LDSC.** Genetic correlation (r_g) was estimated by LDSC [6] in Python (Bulik-Sullivan) on 690,150 HapMap3 SNPs with 1000 Genomes Phase 3 European LD scores, using delete-one-block jackknife (100 blocks) for SEs. Palindromic SNPs were excluded; alleles were harmonised to one reference so every matrix cell uses the same variants. Both Bonferroni (α=0.05/78) and Benjamini–Hochberg FDR were applied across 78 phenotype pairs."),

    ("Instruments were selected at p<5×10⁻⁸ and clumped at r²<0.001 over 10 Mb against 1000 Genomes Phase 3 EUR (strict PLINK-equivalent clumping benchmarked against the PLINK binary; 20–35 instruments retained per analysed BMD site, 116 for knee OA).",
     "Instruments were selected at p<5×10⁻⁸ and clumped at r²<0.001 over 10 Mb against 1000 Genomes Phase 3 EUR (strict PLINK-equivalent clumping benchmarked against PLINK; 20–35 retained per BMD site, 116 for knee OA)."),

    ("The BMD exposure set followed a prespecified rule: of the 11 DXA-derived BMD sites reported by Qian et al. 2025, the cranial head site — the only non-load-bearing region, whose instruments carry non-skeletal pleiotropy (head BMD is genetically correlated with intracranial aneurysm in the source GWAS [2]) — was excluded a priori, leaving ten prespecified canonical DXA regions analysed under the same rule: a site enters MR only if strict clumping completes.",
     "The BMD exposure set followed a prespecified rule: of the 11 DXA-derived BMD sites in Qian et al. 2025, the cranial head site — the only non-load-bearing region, whose instruments carry non-skeletal pleiotropy (head BMD is genetically correlated with intracranial aneurysm [2]) — was excluded a priori, leaving ten prespecified canonical DXA regions; a site enters MR only if strict clumping completes."),

    ("Random-effects IVW used the multiplicative convention SE_random = SE_fixed × √[max(1, Q/df)] (as in TwoSampleMR), not an additive DerSimonian–Laird τ²; where the convention changes an inference, both are reported (§3, Table 2 note).",
     "Random-effects IVW used the multiplicative convention SE_random = SE_fixed × √[max(1, Q/df)] (TwoSampleMR), not additive DerSimonian–Laird τ²; where the convention changes an inference, both are reported (§3, Table 2 note)."),

    ("We report MR-Egger, weighted median and MR-PRESSO alongside IVW and applied MRlap [15], which simultaneously corrects weak-instrument bias and winner's curse while accounting for participant overlap, to all ten pairs.",
     "We report MR-Egger, weighted median and MR-PRESSO alongside IVW, and applied MRlap [15] (corrects weak-instrument bias and winner's curse while accounting for overlap) to all ten pairs."),

    ("The joint correction moved every estimate further from the null (median +20.5%; observed-versus-corrected difference p<0.05 for 9/10 pairs), showing that the towards-null components dominate the away-from-null overlap component, and it reversed no direction or significance conclusion.",
     "The joint correction moved every estimate further from the null (median +20.5%; observed-versus-corrected difference p<0.05 for 9/10 pairs), showing towards-null components dominate the away-from-null overlap component, reversing no direction or significance."),

    ("The primary IVW estimates are therefore conservative with respect to the combined bias, but not overlap alone; we report MRlap alongside, rather than instead of, IVW.",
     "The primary IVW estimates are conservative for the combined bias, not overlap alone; we report MRlap alongside IVW."),

    ("We did not assume K: of the SNPs carrying both a deposited log-odds effect and the reconstructed Z, 39 matched, with |Z_deposited/Z_reconstructed| = 1.000000 (SD 5×10⁻¹⁴), confirming the two files describe one association, and implied K(1−K)=N_eff/N median 0.0871 (IQR 0.0835–0.0881), i.e. K=0.0964; we applied the factor 1/√[K(1−K)] = 3.389.",
     "We did not assume K: of SNPs carrying both a deposited log-odds effect and the reconstructed Z, 39 matched, with |Z_deposited/Z_reconstructed| = 1.000000 (SD 5×10⁻¹⁴), confirming one association, and implied K(1−K)=N_eff/N median 0.0871 (IQR 0.0835–0.0881), i.e. K=0.0964; we applied 1/√[K(1−K)] = 3.389."),

    ("Height was chosen as the primary size proxy: BMI is itself a major knee-OA risk factor and part of the adiposity pathway through which a BMD-raising process could raise knee-OA risk, whereas height indexes skeletal size without adiposity.",
     "Height was the primary size proxy: BMI is a major knee-OA risk factor and part of the adiposity pathway through which a BMD-raising process could raise knee-OA risk, whereas height indexes skeletal size without adiposity."),

    ("First, the exposures differ: Hartley instrumented heel-ultrasound eBMD (scaling with bone size) whereas we instrumented DXA areal BMD (area-normalised).",
     "First, the exposures differ: Hartley instrumented heel-ultrasound eBMD (bone size) whereas we instrumented DXA areal BMD (area-normalised)."),

    ("Second, the designs differ: Hartley's estimate is one-sample 2SLS (overlap 100%) whereas ours is bounded at ≤2.6% and MRlap moved every estimate a median 20.5% further from the null.",
     "Second, the designs differ: Hartley's is one-sample 2SLS (overlap 100%) whereas ours is bounded at ≤2.6% and MRlap moved every estimate a median 20.5% further from the null."),

    # ---- Highlights (H3 char fix) ----
    ("- BMD instruments colocalise with fracture (9/31), not knee OA (0/31): endpoint-specific.",
     "- BMD instruments colocalise with fracture (9/31) not knee OA (0/31): endpoint-specific."),
]

def main():
    with io.open(PATH, "r", encoding="utf-8") as fh:
        text = fh.read()
    for i, (old, new) in enumerate(PAIRS, 1):
        cnt = text.count(old)
        if cnt != 1:
            raise SystemExit(f"[FAIL] pair #{i} matched {cnt} times. old[:70]={old[:70]!r}")
        text = text.replace(old, new, 1)
    with io.open(PATH, "w", encoding="utf-8") as fh:
        fh.write(text)
    print(f"Applied {len(PAIRS)} replacements OK.")

if __name__ == "__main__":
    main()
