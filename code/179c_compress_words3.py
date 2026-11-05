#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Remaining compression pairs: F9 (corrected) + Methods cuts + H3 char fix.
(Alhassan already done in code/179 first batch.)"""
import io
PATH = "docs/08_投稿稿_OA_BMD_fracture.md"
PAIRS = [
    ("**F9** site-gradient ordering (executed: arms − femoral neck +0.017, one-sided p=0.228).",
     "**F9** site-gradient ordering (arms−femoral neck +0.017, p=0.228)."),

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
    print(f"Applied {len(PAIRS)} remaining replacements OK.")
if __name__ == "__main__":
    main()
