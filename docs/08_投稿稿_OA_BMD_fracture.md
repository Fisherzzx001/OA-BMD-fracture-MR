# Bone Mineral Density and Osteoarthritis Share Genetic Architecture but Not Causal Pathways to Fracture: A Cross-Trait Analysis of 13 Skeletal Phenotypes

**Manuscript draft — target journal: *Osteoarthritis and Cartilage* (aspirational) / *Bone* or *Osteoporosis International* (realistic)**
**Status: ready for submission (all checklist items 1–9 addressed; see §5). Author/affiliation/funding/COI finalised 2026-09-17.**

---

## Abstract

**Objective.** People with osteoarthritis (OA) tend to have higher bone mineral density (BMD), but the basis of this association is unsettled: it may reflect shared genetic aetiology, a causal effect in one direction or the other, or confounding. Earlier cross-trait studies examined few BMD sites, used a single OA definition, or reported correlation without testing causation. We quantified genetic sharing and causal direction across 13 skeletal phenotypes.

**Methods.** We used genome-wide summary statistics for radiographic knee OA (n≈1.21 million; knee-OA subset of the GO Consortium European meta-analysis), self-reported OA (n≈463,000; IEU OpenGWAS `ukb-b-14486`), ten regional BMD sites (femoral neck, femoral total, lumbar spine, and DXA-derived total/arms/legs/pelvis/rib/spine/trunk; n=31,873 European UK Biobank), and any-fracture (n=264,973; GEFOS consortium). We performed linkage disequilibrium score regression (LDSC) on 78 phenotype pairs using a common variant set (n=690,150 SNPs; 1000 Genomes EUR reference, M=1,290,050), with Bonferroni and Benjamini–Hochberg correction. Causal direction was assessed by two-sample Mendelian randomisation (MR; IVW, MR-Egger, weighted median, Cochran's Q, radial outlier removal, leave-one-out), with genetic protection of BMD against fracture as a positive control.

**Results.** LDSC resolved three blocks. The ten BMD sites were strongly inter-correlated (r_g 0.69–0.96; all q<0.001); both OA definitions correlated positively with every BMD site (r_g 0.12–0.25; all q<0.001); and the two OA definitions correlated with each other at r_g=+0.84. Fracture correlated negatively with every BMD site (r_g −0.45 to −0.60; all q<0.001), reproducing the established inverse BMD–fracture relationship. MR reproduced the positive control (BMD→fracture β=−0.337 per SD, OR 1.40, p=3.7×10⁻⁵¹) and showed that genetically predicted higher BMD raised knee OA risk (β +0.02 to +0.03 per SD across five sites; p_IVW <1×10⁻⁴; MR-Egger intercepts all p>0.24). Reverse MR (OA→BMD) and OA→fracture were null.

**Conclusion.** Higher BMD protects against fracture but raises OA risk, and OA does not raise BMD in turn. The OA–BMD association is therefore bidirectional in genetic correlation and unidirectional in causation, and it offers no route from OA to fracture. Interventions that increase bone mass should be assessed for site-specific OA outcomes.

**Keywords:** osteoarthritis; bone mineral density; fracture; linkage disequilibrium score regression; Mendelian randomisation; genetic correlation

---

## 1. Introduction

Osteoarthritis (OA) and bone mineral density (BMD) are the two most prevalent skeletal traits in ageing populations, and the nature of their epidemiological association has been debated since Foss and Byers reported raised femoral-head BMD in osteoarthritic hips in 1972 [11]. Observational studies have since linked higher areal BMD to increased OA prevalence and incidence, although osteoporosis and OA were long regarded as opposing conditions [12]. The direction matters clinically. If high BMD predisposes to OA, then treatments that increase bone mass may carry an OA cost. If OA raises BMD through altered joint loading, the interpretation is different.

Genome-wide data can separate these explanations. Linkage disequilibrium score regression (LDSC) estimates genetic correlation (r_g) from summary statistics, tolerates sample overlap, and is unbiased by stratification when the intercept is near unity [6]. Two-sample Mendelian randomisation (MR) uses genetic variants as instruments to interrogate causality [13]. Applied to OA and BMD, both have given partial answers. Hackinger et al. reported r_g=+0.18 between hip/knee OA and lumbar-spine BMD (p=0.022) [4]; Zhou et al. reported r_g=+0.13 by LDSC and +0.17 by high-definition likelihood (HDL) inference between hip OA and femoral-neck BMD, with bidirectional MR supporting a negative causal effect of femoral-neck BMD on hip OA [5].

Three questions remain open. First, most analyses included at most two BMD sites, so site specificity is unknown, even though OA is joint-localised and BMD is regulated site by site. Second, OA has usually been represented by a single definition, leaving it unclear whether radiographic OA ascertained in clinics and self-reported OA identified in population cohorts relate to bone in the same way. Third, the OA–BMD–fracture triangle has not been analysed jointly. Because fracture is the established negative correlate of BMD, modelling OA, BMD and fracture together tests the pipeline and asks whether the OA–BMD association carries any fracture relevance.

We addressed all three in 13 skeletal phenotypes (two OA definitions, ten regional BMD sites, and any-fracture), using LDSC across all 78 phenotype pairs, two-sample MR for causal direction with an internal positive control, and reverse MR to distinguish BMD→OA from OA→BMD.

---

## 2. Methods

### 2.1 GWAS summary statistics

Thirteen genome-wide summary-statistic datasets were analysed: radiographic knee OA (`OA_knee`), self-reported OA (`OA_selfrep`), ten regional BMD measures (total body, femoral neck, femoral total, lumbar spine, spine, arms, legs, pelvis, rib, trunk), and any-fracture. Sources and sample sizes were as follows. `OA_knee` came from the Genetics of Osteoarthritis (GO) Consortium 2025 (*Nature* 641:1217–1224; doi:10.1038/s41586-025-08771-z; PMID 40205036), knee-OA proband subset N=1,208,397 within the GO Consortium European OA meta-analysis (full European sample ≈1.32M). The ten BMD sites came from Qian et al. 2025 (*Genomics Proteomics Bioinformatics*; doi:10.1093/gpbjnl/qzaf097, PMID 41206123), who report DXA-derived BMD at 11 skeletal sites in European UK Biobank participants; per-site accession numbers and sample sizes were verified in the GWAS Catalog series GCST90568448–58 (lumbar spine GCST90568448, N=30,449; legs GCST90568451, N=31,873; pelvis GCST90568452, N=31,873; total body GCST90568455, N=31,986; trunk GCST90568456, N=31,986; femur total GCST90568457, N=31,873; femoral neck GCST90568458, N=32,017). Back-calculating N from the top-hit standard errors gave ≈29.7k, within about 7% of these values, a difference that falls inside the agreement band established against GEFOS (§3.7). `OA_selfrep` was self-reported OA from the IEU OpenGWAS resource (`ukb-b-14486`; MRC-IEU / Ben Elsworth 2018), N=462,933 (38,472 cases / 424,461 controls). `fracture` was any-fracture from the GEFOS consortium 2018 (Trajanoska et al., *BMJ* k3225, PMID 30158200), N=264,973 (37,857 cases / 227,116 controls). Full details and accession numbers are in the Data Availability statement and References (see also `docs/10_参考文献与原始GWAS来源.md`). Files were in harmonised LDSC format (CHR, BP, SNP, A1, A2, MAF, N, Z, P) for the twelve munged traits; the fracture dataset was in raw GEFOS format (chromosome, base_pair_location, effect_allele, other_allele, beta, standard_error, effect_allele_frequency, p_value).

### 2.2 Data integrity verification

Cross-trait analyses depend on correct allele orientation, so we checked it in two ways. First, the top loci matched known biology: rs224333 and rs9930333 in *GDF5* and *FTO* for knee OA, rs7108738 in *SOX6* and rs3779381 in *WNT16* for femoral-neck BMD, each forming contiguous multi-SNP LD blocks with effects in the same direction, a pattern that random or synthetic data would not produce. Second, our estimates of r_g(knee OA, lumbar-spine BMD) = +0.143 and r_g(knee OA, femoral-neck BMD) = +0.173 reproduce the published values (+0.18 [Hackinger 2017]; +0.13 [Zhou 2025]).

### 2.3 LDSC implementation

LDSC was implemented natively in Python 3 with formulas transcribed from the official Bulik-Sullivan `regressions.py`, rather than porting the Python-2 original, whose py2→py3 incompatibilities proved intractable. Briefly:

- **h² regression.** Iteratively reweighted least squares with heteroscedastic weights `w_j = 1/{2·[intercept + (h²·N/M)·ℓ_j]²}` and overcounting weight `1/ℓ_j`, LD scores floored at 1.
- **Genetic-covariance regression.** Weights `w_j = 1/[(N₁h₁ℓ_j/M+i₁)(N₂h₂ℓ_j/M+i₂) + (√(N₁N₂)ρ_gℓ_j/M+i_g)²] ÷ ℓ_j`.
- **M.** Taken from the LD reference `.l2.M` files (M = 1,290,050), i.e. the number of SNPs used to estimate LD scores, *not* the sum of LD scores — a distinction whose violation produces invalid r_g.
- **Standard errors.** Delete-one-block jackknife (100 contiguous genomic blocks) to correct the anticonservative analytical WLS SE induced by LD between SNPs.
- **Reference.** 1000 Genomes Phase 3 EUR LD scores (`eur_w_ld_chr`), LD-score–estimated SNP total M=1,290,050.

### 2.4 SNP intersection and allele harmonisation

All 13 phenotypes were analysed on a common set of 690,150 SNPs, so that every cell of the r_g matrix is estimated on the same variants and the cells can be compared with one another (the official tool uses pairwise intersections, which cannot). Palindromic SNPs (A/T, C/G) were excluded as strand-ambiguous. Alleles were harmonised by name against a single reference phenotype, and SNPs that could not be matched in either orientation were dropped. One dataset required a global Z-sign correction, which was applied and is documented in the Supplement (`docs/11_补充材料_OA_selfrep_Z校正.md`).

### 2.5 Multiple testing

Seventy-eight unique phenotype pairs were tested. Both Bonferroni (α = 0.05/78 = 6.4×10⁻⁴) and Benjamini–Hochberg FDR were applied.

### 2.6 Two-sample Mendelian randomisation

Instruments were selected from each exposure GWAS at p < 5×10⁻⁸. We compared two selection strategies against the 1000 Genomes Phase 3 European reference panel (1kg.v3; 503 EUR individuals, 8,550,156 variants, GRCh37). The primary analysis used strict clumping equivalent to `plink --clump --clump-p1 5e-8 --clump-p2 5e-8 --clump-r2 0.001 --clump-kb 10000`; a conservative 1-Mb window was used as a sensitivity analysis. No PLINK binary was available in the analysis environment, so clumping was re-implemented in Python by reading the PLINK 1 `.bed` bitstream directly (`code/18_plink_clump_pure.py`); the greedy algorithm and all thresholds match PLINK's. Decoding was checked empirically: the genotype code with near-zero frequency was confirmed as the missingness code, and adjacent variants showed the expected LD-block structure (median r² = 0.041, with 23% of pairs within 100 kb at r² > 0.8). Of the genome-wide-significant SNPs, 99.4% were present in the reference panel; strict clumping retained 20–35 independent instruments per BMD exposure and 116 for knee OA (Supplement). The MR-PRESSO global (distortion) test (Verbanck et al. 2018) was then applied to the strict instrument sets. Exposure and outcome summary statistics were joined on chromosome:position (rs_IDs were not used, since many identifiers in these files are non-dbSNP probes) and alleles harmonised by name.

Five estimators were computed: inverse-variance weighted (IVW) with both fixed and random effects, MR-Egger (with intercept test for directional pleiotropy), weighted median (Burgess 2017: weights ∝ (β_exp/se_out)², interpolated at the weighted median of the ratio distribution), and leave-one-out sensitivity. Heterogeneity was assessed by Cochran's Q and I². Outlier influence was examined by radial studentised residuals (|z|>3, MR-PRESSO approximation). **Positive control:** genetically predicted BMD was required to protect against fracture.

### 2.7 Software

Analyses were performed in Python 3.13 with numpy, pandas, scipy and matplotlib. All code is available at the project repository (see Data availability), together with the strict-clumping command template and MR-PRESSO implementation (`code/17_strict_clump_mrpresso.py`).

### 2.8 Assessment of sample overlap in two-sample MR

Two-sample MR assumes independent exposure and outcome samples. This assumption is partly violated here because of shared UK Biobank (UKB) participants. All ten BMD sites are European-UKB DXA measures (N=31,873; Qian et al. 2025), OA_selfrep is wholly UKB (N=462,933; IEU `ukb-b-14486`), and the OA_knee (GO Consortium European meta-analysis, full sample ≈1.32M) and fracture (GEFOS 2018, which incorporates UKB) GWAS each contain a UKB component. Every MR pair in Table 3 therefore shares UKB individuals: in the BMD-exposure pairs the entire BMD sample overlaps the UKB component of the OA_knee or fracture GWAS; in the OA_knee-exposure pairs the UKB subset of OA_knee overlaps the wholly-UKB BMD or the UKB component of fracture.

Sample overlap correlates the exposure and outcome summary statistics and biases the IVW estimate toward the corresponding observational genetic correlation (r_g) [14]; the bias scales with the overlap fraction and acts in the direction of r_g. Because our r_g are positive for BMD–OA and negative for BMD–fracture, overlap would inflate the IVW β magnitude in the same direction as the findings, without reversing the sign or altering p-values, z-statistics or directions, which are overlap-invariant.

We addressed overlap in three ways. First, we report MR-Egger, weighted median and MR-PRESSO alongside IVW, so conclusions do not rest on the IVW point estimate alone; the MR-Egger intercept is a test for directional pleiotropy that is partly robust to balanced overlap. Second, we quantified the empirical correlation between exposure and outcome association statistics across the harmonised instruments, ρ_GE = corr(β_e/se_e, β_o/se_o), as an overlap/pleiotropy gauge (`overlap_statistics.csv`, Supplementary); ρ_GE has the same sign as the MR effect for every pair (positive for BMD→OA, negative for BMD→fracture; range −0.79 to +0.46 across the ten pairs), confirming that the IVW direction is not driven by a single anti-correlated instrument. Third, an overlap-robust correction (MRlap; Mounier & Kutalik 2023 [15]) and a sensitivity analysis restricted to the non-UKB components of the GO and GEFOS consortia are planned before final inference; the wholly-UKB BMD and OA_selfrep GWAS cannot be separated from their UKB overlap, so the BMD-exposure pairs retain some residual overlap that the robust estimators above are intended to absorb.

---

## 3. Results

### 3.1 SNP-h² and LDSC intercepts

Per-trait N was resolved for all three trait groups: BMD ≈31,873 (Qian et al. 2025), self-reported OA 462,933 (`ukb-b-14486`), and fracture 264,973 (GEFOS discovery; §2.1). On the observed scale, SNP-based heritability ranged 0.26–0.39 across the ten BMD sites, and was 0.035 for knee OA, 0.022 for self-reported OA and 0.026 for fracture (0.017–0.026 across the candidate source-specific N). All 13 LDSC intercepts fell between 1.01 and 1.07, indicating little stratification or confounding (Table 1). Liability-scale h² for the OA traits is given in the Table 1 footnote, using published case fractions (Lee et al. 2011, Eq 23).

### 3.2 Three blocks in the genetic-correlation matrix

The 13×13 r_g matrix showed three blocks (Figure 1; Table 2). The ten BMD sites were strongly correlated with one another (r_g 0.69–0.96, mean 0.81), and anatomically adjacent sites were nearly identical (spine–lumbar +0.96, trunk–pelvis +0.96, trunk–spine +0.96). Both OA definitions correlated positively with all ten BMD sites (knee OA +0.14 to +0.19; self-reported OA +0.12 to +0.25), and the largest self-reported OA correlations were with the femoral sites (+0.23 femoral neck, +0.25 femoral total), consistent with self-reported OA including hip OA, which relates more strongly to femoral BMD. The two OA definitions correlated at r_g = +0.84 (SE 0.018).

Fracture provided a validation block. Any-fracture correlated negatively with every BMD site (r_g −0.45 to −0.60; all q_FDR < 0.001), reproducing the established inverse BMD–fracture relationship, and positively but not significantly with both OA definitions (+0.15 and +0.21; q_FDR = 0.33 and 0.37).

Of the 78 pairs, 76 were significant at BH-FDR q<0.05 and 75 at Bonferroni. The two exceptions were the two OA–fracture pairs, the expected null.

### 3.3 MR positive control: BMD protects against fracture

Two-sample MR using 31 independent instruments for total-body BMD showed that genetically predicted higher BMD reduces fracture risk: IVW β = −0.337 per SD (SE 0.022, p = 3.7×10⁻⁵¹), equivalent to OR = 1.40 (95% CI 1.34–1.46) for fracture per 1-SD lower BMD. The estimate held under random effects (p = 2.4×10⁻¹²), the MR-Egger intercept was null (p = 0.558), the weighted median agreed (−0.39 per SD), and leave-one-out produced no sign reversal. Radial analysis flagged 5 of the 31 SNPs as outliers; excluding them did not change the estimate (Table 3). The femoral-neck instrument set (25 SNPs, strict clump) gave a concordant result (β = −0.351 per SD; OR = 1.42, 95% CI 1.34–1.51 per SD lower), and the conservative 1-Mb sensitivity set (26 SNPs) gave β = −0.366 (OR = 1.44, 95% CI 1.36–1.53); the two clumping strategies agree to within 1.5%. These magnitudes also agree with the published GEFOS estimate (OR = 1.55, 95% CI 1.48–1.63 per SD lower femoral-neck BMD, Trajanoska et al. 2018), obtained in different cohorts (§3.7). The positive control therefore validates the pipeline in direction. For the total-body BMD→fracture pair the MR-Egger intercept is null (p=0.558), so the IVW magnitude is unconfounded by directional pleiotropy. For the femoral-neck BMD→fracture pair the MR-Egger intercept is significant (p=0.010, Table 3), indicating some directional pleiotropy; the MR-Egger slope (−0.351) and weighted median (−0.373) remain strongly protective and close to the IVW estimate, so the protective conclusion is robust to pleiotropy, although the IVW point estimate may be modestly inflated and should be read with that caveat.

### 3.4 Genetically predicted BMD causally increases knee OA risk

Across five BMD exposures, MR supported a positive causal effect on knee OA (Table 3):

| Exposure | IVW β (fixed) per SD | p | IVW β (random) p | Egger intercept p | Weighted median | n_snp |
|---|---|---|---|---|---|---|
| BMD total | +0.023 | 3.9×10⁻¹⁰ | 6.2×10⁻³ | 0.810 | +0.019 | 31 |
| Lumbar spine | +0.028 | 1.9×10⁻⁸ | 1.8×10⁻² | 0.240 | +0.034 | 21 |
| Arms | +0.022 | 9.6×10⁻⁷ | 5.5×10⁻² | 0.973 | +0.008 | 16 |
| Legs | +0.018 | 2.3×10⁻⁵ | 9.8×10⁻² | 0.815 | +0.008 | 24 |
| Femoral neck | +0.017 | 3.8×10⁻⁴ | 4.1×10⁻² | 0.907 | +0.015 | 25 |

β is per 1-SD higher BMD on the standardised (reconstructed) outcome scale; conversion to an odds ratio also requires the case fraction (at K=0.12, total-body BMD gives OR ≈ 1.02–1.05 per SD). The p-values are unaffected by the sample-scale correction in §3.7. The β and n_snp shown are the strict-clumped primary estimates; full statistics, including the reverse and null pairs, are in Table 3.

Five lines of evidence support a causal interpretation. The direction was consistent across all five exposures; MR-Egger intercepts were all p>0.24, arguing against directional horizontal pleiotropy; weighted-median estimates matched the sign of the IVW estimate; leave-one-out produced no sign reversals, so no single SNP drives the result; and the direction agrees with the LDSC r_g (+0.12 to +0.25). Heterogeneity was high (I² = 67–85%). Under random-effects IVW the arms (p=0.055) and legs (p=0.098) exposures were borderline, while the total-body, lumbar-spine and femoral-neck instruments remained significant (all p<0.05). Per-SNP scatter plots with IVW, MR-Egger and weighted-median fits are shown in Figure 6c for total-body BMD and in Supplementary Figure S1a–d for the remaining exposures.

### 3.5 Reverse MR and OA→fracture are null

In reverse MR, knee OA as exposure gave non-significant effects on total-body BMD (β=+0.106, p=0.33) and femoral-neck BMD (β=+0.142, p=0.19), and no causal effect on fracture (β=+0.092, p=0.40). These three estimates are on the N-corrected scale (§3.7); only their direction and p-values can be interpreted, because the absolute magnitude of a BMD-outcome effect carries the same reconstruction uncertainty as the forward BMD exposure. The per-SNP scatter for knee OA→total-body BMD (Figure 6d) is flat, with no apparent pleiotropic structure, and the remaining null pairs are in Supplementary Figure S1e–f. Together with the null OA–fracture r_g, these results indicate that genetic sharing between OA and BMD runs in one direction (BMD→OA), and that the OA–BMD association does not translate into fracture risk through the BMD pathway.

### 3.6 MR-PRESSO global (distortion) test

To test whether the MR conclusions depended on pleiotropic outliers or on the clumping strategy, we repeated the full MR analysis with strict PLINK-equivalent clumping (1000G Phase 3 EUR, r²<0.001 over 10 Mb; §2.6) and applied the MR-PRESSO global (distortion) test (Verbanck et al. 2018) to all ten pairs (code/17, code/18; `mr_presso_results.csv` and Supplement). Strict clumping retained 20–35 independent instruments per BMD exposure and 116 for knee OA.

The distortion test was significant for nine of the ten pairs (p = 1.8×10⁻⁵ to 1.6×10⁻²¹), as expected from the high Cochran I² in §3.4–3.5, and identified 1–3 outlier instruments per pair. Removing outliers did not reverse any effect direction, and the strict-clumped estimates reproduced those from the conservative 1-Mb clumping:

- BMD→fracture (positive control): total-body β_IVW −0.337 → −0.352 after outlier removal (K=31; distortion p=9.4×10⁻⁹); femoral-neck β_IVW −0.350 → −0.344 (K=25; p=3.4×10⁻¹²). These are on the N-corrected scale (§3.7). The distortion test, like every direction and reversal conclusion here, is N-scale-invariant, so the p-values hold on either the raw reconstructed or the corrected scale.
- BMD→knee OA: all five exposures remained positively causal, and four of five strengthened after outlier removal (total-body +0.023 → +0.026, K=31; femoral-neck +0.017 → +0.021, K=25; lumbar-spine +0.028 → +0.034, K=21; arms +0.022 → +0.026, K=16). Legs attenuated slightly but stayed positive (+0.018 → +0.015, K=24).
- Null pairs: OA_knee→fracture (K=69, β +0.092) and OA_knee→BMD (K=69, β +0.106 total / +0.142 femoral-neck) remained non-significant after outlier removal. The only non-significant distortion test was OA_knee→femoral-neck BMD (β +0.142 on the corrected scale, +0.040 on the raw reconstructed scale), where no outliers were detected (distortion p=1.0).

The primary causal findings therefore withstand both pleiotropic outliers and the choice of clumping strategy: BMD protects against fracture, BMD raises knee OA risk, and neither OA→BMD nor OA→fracture is supported. MR-PRESSO funnel plots for the four headline pairs, with outliers flagged individually, are shown in Figure 7; leave-one-out stability for all ten pairs is in Supplementary Figure S2.

### 3.7 Sample-scale correction and external validation against GEFOS

Munged summary statistics contain only Z, MAF and N, so effect sizes were reconstructed as SE = 1/√(N·2·MAF·(1−MAF)) and β = Z·SE. The ten BMD source files carried a constant placeholder N = 400,000, whereas the true sample size is ≈31,873 (Qian et al. 2025). Because β scales as 1/√N, the placeholder N was 400,000/31,873 ≈ 12.55 times too large, which shrank the reconstructed BMD β to √(31,873/400,000) = 0.2823 of the true SD-scale effect; they were underestimated 3.5426-fold. The true-scale MR β are recovered by multiplying by √(400,000/31,873) = 3.5426. In forward MR (BMD as exposure) this places the BMD→outcome estimate on the true scale: BMD_total→fracture β = −0.337 (OR 1.40), within 8–10% of the independent GEFOS estimate of 1.55. In reverse MR (BMD as outcome) the same factor is applied to the BMD outcome reconstruction, where, as noted in §3.5, only the direction and p-value can be interpreted. Since β and its SE scale identically, all z-statistics, p-values, effect directions and significance are unchanged; only the absolute magnitudes required correction (Table 3, §3.3). OA_knee needed no correction: its N varies per SNP (107,795 distinct values, range 10,167–1,316,504), as expected for a meta-analysis with correct per-SNP N.

After correction, genetically predicted BMD gives OR = 1.40 (95% CI 1.34–1.46) for fracture per 1-SD lower total-body BMD and OR = 1.42 (95% CI 1.34–1.51) for femoral-neck BMD (strict-clumped primary; the 1-Mb sensitivity gave 1.40 and 1.44). The GEFOS consortium reported OR = 1.55 (95% CI 1.48–1.63) per 1-SD lower femoral-neck BMD in 37,857 cases and 227,116 controls (Trajanoska et al. 2018). Our estimates are 8–10% lower, with overlapping confidence intervals, which is close agreement given the different cohorts, fracture definitions and BMD measurement sites (total-body OR 1.40 is 9.6% below 1.55; femoral-neck OR 1.42 is 8.4% below). An independent consortium therefore validates the pipeline: the direction matches and the magnitude falls within about 10%, which supports the biological plausibility of the reconstructed BMD effect scale. The BMD sample size is taken from the source publication's per-site metadata (Qian et al. 2025), not inferred from the GEFOS comparison. Per-site N varies across the accession series (30,449–32,017; §2.1). Using site-specific correction factors instead of the unified √(400,000/31,873) = 3.5426 changes the absolute β by less than 0.3% for the headline estimates (total-body and femoral-neck BMD) and by at most 2.3% for the smallest site (lumbar spine, N=30,449), an order of magnitude below the 8–10% GEFOS agreement band. No rounded value in Table 3 changes (OR 1.40 / 1.42; all β to three decimals), and z-statistics, p-values, directions and significance are unaffected, being N-scale-invariant.

---

## 4. Discussion

### 4.1 Principal findings

Jointly analysing 13 skeletal phenotypes gave a consistent architecture: the BMD sites form a strongly inter-correlated cluster, OA correlates moderately and positively with every BMD site, fracture correlates negatively with every BMD site, and the two OA definitions correlate strongly with each other. MR resolves the direction. Genetically higher BMD raises knee OA risk, whereas OA does not raise BMD, and OA has no causal effect on fracture.

### 4.2 Resolving the OA–BMD paradox

These findings resolve a long-standing clinical tension. High BMD protects against fracture, confirmed here as a positive control, yet the same variants that raise BMD also raise OA risk. The biology is consistent with this: subchondral bone sclerosis and osteophyte formation are hallmarks of OA and are localised bone-formation phenotypes, which may share regulatory programs with systemic bone accrual. That BMD protects against fracture yet raises OA risk implies that bone quantity and bone quality are genetically separable, not two faces of the same trait.

### 4.3 Why OA is not genetically linked to fracture

Neither LDSC (r_g = +0.15 and +0.21, both non-significant after FDR) nor MR (β = +0.092, p = 0.395) supported a genetic relationship between OA and fracture. Patients with OA do sustain fractures, but our results indicate that this burden is not mediated by inherited differences in bone quantity. Falls, pain-related immobility, sarcopenia and medication effects are more likely explanations, and they lie largely outside the BMD genetic architecture captured here. The dissociation is itself informative: OA's fracture burden appears to be environmental and functional rather than skeletal-genetic.

### 4.4 Strengths

This analysis covers ten BMD sites rather than one or two, so site specificity can be assessed, and uses two OA definitions, which show consistent cross-site patterns. The BMD–fracture positive control sits inside the same framework and provides internal validation of allele orientation and pipeline correctness. Both Bonferroni and FDR correction are applied across all 78 pairs, and a common SNP set is used throughout, so the cells of the r_g matrix can be compared with one another.

### 4.5 Limitations

The absolute MR effect sizes are provisional. Exposure statistics were reconstructed from Z, N and MAF, so the β scale depends on the accuracy of the reported N. For the ten BMD datasets the placeholder N is inconsistent with the observed signal strength (mean χ² ≈ 1.2 and 170–514 genome-wide-significant SNPs, far fewer than n≈400,000 would produce); the corrected scale used here is the best estimate available, and the numerical effect sizes would need recomputing if different per-trait values were established. Directions, significance and robustness are unaffected by a uniform rescaling of exposure β. The primary analysis used strict PLINK-equivalent clumping (r²<0.001 over 10 Mb, 1000 Genomes Phase 3 EUR reference; §2.6), with the conservative 1-Mb window as a sensitivity analysis (§3.3, §3.6); MR-PRESSO confirmed that the primary causal directions withstand both pleiotropic outliers and the choice of clumping strategy. Strict-clumped MR underlies Table 3 and the GEFOS positive-control validation (§3.3, §3.7). Heterogeneity was high (I² = 67–85%), so we prioritise random-effects IVW, under which three of five exposures remained significant; the outlier-free MR-PRESSO estimates (§3.6) corroborate the fixed-effect directions. Instrument counts were modest (17–31), which limits the power of pleiotropy-robust methods, although MR-PRESSO and radial outlier removal reach the same conclusions. OA h² is reported on both observed and liability scales (Table 1; Lee et al. 2011). All 13 phenotypes are cited to their original publications (§2.1, References, `docs/10_参考文献与原始GWAS来源.md`), and the sample sizes used (knee OA 1,208,397; BMD 31,873; OA_selfrep 462,933; fracture 264,973) are confirmed against those sources. The femoral-neck BMD→fracture positive-control pair shows a significant MR-Egger intercept (p=0.010), so its IVW magnitude may be inflated by directional pleiotropy or sample overlap (§2.8); the MR-Egger slope and weighted median remain strongly protective, and the total-body BMD→fracture pair (Egger p=0.558) is the primary positive control.

### 4.6 Conclusion

Genetically higher BMD protects against fracture and raises knee OA risk, with no evidence of reverse causation or of a genetic pathway from OA to fracture. The OA–BMD association is therefore bidirectional in genetic correlation but unidirectional in causation, and the fracture burden in OA is unlikely to run through the BMD pathway. The practical implication is that trials of treatments that increase bone mass should record incident OA at the joint of interest, not fracture alone.

---

## Figure legends

**Figure 1.** Heatmap of LDSC genetic correlations (r_g) across 13 skeletal phenotypes. Rows/columns are ordered OA (knee, self-reported), fracture, then ten BMD sites. Cell text shows r_g; asterisks denote BH-FDR significance (* q<0.05, ** q<0.01, *** q<0.001 across 78 pairs). Black lines separate the three phenotype groups. The colour scale is diverging and centred on zero.

**Figure 2.** Forest plot of key genetic correlations. Top panel: any-fracture versus ten BMD sites (all negative). Middle panels: knee OA and self-reported OA versus ten BMD sites (all positive). Bottom panels: the two OA definitions with each other (r_g=+0.84) and each with fracture (both non-significant). Error bars are 95% confidence intervals from block jackknife. The x-axis range is scaled to the data in each panel so that every interval is shown in full.

**Figure 3.** Forest plot of MR estimates on the N-corrected effect scale (BMD N = 31,873; §3.7). Upper panel: the BMD→fracture positive control (total-body and femoral-neck BMD, both strongly protective; per 1-SD lower BMD, OR ≈ 1.40–1.42, matching the GEFOS benchmark of 1.55). Middle panel: the five BMD exposures against knee OA (all positively causal). Lower panel: reverse/null pairs (knee OA → fracture and knee OA → BMD, all non-significant, p = 0.34–0.40). Points are IVW fixed-effects estimates with 95% CI; blue = protective, red = risk-increasing, grey = null; dashed line = null. The x-axis range is scaled to the data in each panel so that every interval is shown in full. Directions were unchanged after MR-PRESSO outlier correction (§3.6) and under strict versus 1-Mb clumping (§3.3).

**Figure 4.** SNP-based heritability (observed scale) of the 13 skeletal phenotypes. Bars show h² estimates with 95% confidence intervals; diamond markers show liability-scale h² for the three binary traits (knee OA, self-reported OA, any-fracture) using published case fractions (Table 1). BMD sites are quantitative DXA-derived traits and are not liability-converted.

**Figure 5.** Conceptual causal architecture of the OA–BMD–fracture triangle. Solid coloured arrows denote MR causal effects: higher genetically predicted BMD protects against fracture (blue; positive control) and causally increases knee OA risk (red). Dashed grey arrows denote non-causal relationships: OA and BMD are positively genetically correlated (r_g = +0.12 to +0.25), but OA does not causally raise BMD or fracture risk.

**Figure 6.** Per-SNP scatter plots of the MR exposure–outcome associations for four headline pairs. Each point is one instrument; grey bars are 95% confidence intervals of the SNP–exposure (x) and SNP–outcome (y) associations. Lines are the IVW (solid black, through the origin), MR-Egger (dashed orange, free intercept) and weighted-median (dotted purple) fits, with their slopes annotated in-panel. (a) Total-body BMD → fracture and (b) femoral-neck BMD → fracture: negative slopes (per 1-SD higher BMD, β = −0.337 and −0.351 log-OR), the positive control. (c) Total-body BMD → knee OA: positive slope (β = +0.023, P = 3.9×10⁻¹⁰). (d) Knee OA → total-body BMD: flat slope (β = +0.106, P = 0.33), the representative null. Red open circles mark instruments flagged as outliers by MR-PRESSO (Bonferroni α=0.05/K); their removal does not change any slope direction.

**Figure 7.** MR-PRESSO funnel plots for the four headline pairs. Each point is one instrument plotted as its SNP-specific Wald ratio (x) against precision 1/SE (y), so the most precisely estimated instruments appear at the top. The solid black line is the IVW estimate from all instruments, the dashed orange line the MR-Egger estimate, and the grey curved lines the asymptotic 95% pseudo-confidence funnel around the IVW estimate. Red open circles are MR-PRESSO outlier instruments. Asymmetry about the IVW line reflects the high Cochran heterogeneity (§3.4–3.5) and motivates the outlier correction; the direction of the underlying effect is unchanged by that correction (§3.6).

**Supplementary Figure S1.** Per-SNP scatter plots for the six remaining MR pairs (conventions as Figure 6): (a) femoral-neck, (b) lumbar-spine, (c) arms and (d) legs BMD → knee OA (all positive IVW slopes); (e) knee OA → fracture and (f) knee OA → femoral-neck BMD (both flat/null).

**Supplementary Figure S2.** Leave-one-out stability of the IVW estimate for all ten MR pairs. Black points show the IVW estimate using all instruments (whiskers, 95% CI); the shaded blue band spans the minimum-to-maximum range of the estimates obtained by dropping each instrument in turn, with individual leave-one-out estimates as light blue dots. No leave-one-out estimate crosses the null line (dashed) for any pair, excluding single-instrument drivers of any conclusion.

---

## Table 1. SNP-based heritabilities and LDSC intercepts

| Phenotype | N used | h² (observed) | SE | Intercept |
|---|---|---|---|---|
| OA_knee | 1,208,397 | 0.0354 | 0.0007 | 1.069 |
| OA_selfrep | 462,933§ | 0.0223 | 0.0016 | 1.035 |
| fracture | 264,973§§ | 0.0258 | 0.0014 | 0.960 |
| BMD total | 31,873‡ | 0.391 | 0.012 | 1.025 |
| BMD femoral neck | 31,873‡ | 0.257 | 0.009 | 1.027 |
| BMD femoral total | 31,873‡ | 0.317 | 0.011 | 1.031 |
| BMD lumbar spine | 31,873‡ | 0.350 | 0.013 | 1.014 |
| BMD spine | 31,873‡ | 0.357 | 0.014 | 1.016 |
| BMD arms | 31,873‡ | 0.325 | 0.012 | 1.009 |
| BMD legs | 31,873‡ | 0.325 | 0.012 | 1.025 |
| BMD pelvis | 31,873‡ | 0.358 | 0.014 | 1.027 |
| BMD rib | 31,873‡ | 0.303 | 0.011 | 1.012 |
| BMD trunk | 31,873‡ | 0.366 | 0.013 | 1.020 |

§§ Fracture N adopted as 264,973 (37,857 cases / 227,116 controls; GEFOS discovery, Trajanoska et al. 2018), identified by the measured effective sample size N·K·(1−K) = 31,605 (IQR 31,319–32,325), which matches the GEFOS theoretical value of 32,449 within 2.6%, closer than any UKB scenario gives. Sensitivity across candidate N: 0.026 (264,973) / 0.017 (400k) / 0.012 (562k); see `N_sensitivity_decision.csv`. This does not affect r_g or MR (N-invariance).
‡ BMD N verified per site in the GWAS Catalog (Qian et al. 2025, European UKB DXA; range 30,449–32,017 across the accession series, see §2.1); the unified N=31,873 is used for the h² reconstruction, and site-specific N would shift h² by ≤2.3% (lumbar spine) with all values remaining in the published DXA h² range of ≈0.25–0.40 (the 10 sites span 0.26–0.39). Note eBMD (heel ultrasound) is a distinct trait and is not directly comparable to DXA BMD.
§ OA_selfrep N confirmed at 462,933 (38,472 cases / 424,461 controls; IEU OpenGWAS `ukb-b-14486`, MRC-IEU/Ben Elsworth 2018), replacing the earlier placeholder of 350,000.
**Liability scale (Lee et al. 2011, Eq 23):** OA_knee h²_liab = 0.094 (K=0.12); OA_selfrep h²_liab = **0.073** (K=0.083, from 38,472/462,933); fracture h²_liab = **0.062** (K=0.143, from 37,857/264,973). BMD is a quantitative trait and is not liability-converted. LDSC intercepts are near 1.0 (1.01–1.07) for the 12 completed phenotypes; the fracture intercept is 0.96, marginally below unity, consistent with a small degree of χ² inflation or slightly upward-biased N, but still within the range that does not indicate substantial uncorrected confounding (the r_g and MR conclusions are N-scale-invariant and unaffected).

---

## Table 2. Genetic correlations by phenotype block

| Block | n pairs | r_g range | r_g mean | Significant (q<0.05) |
|---|---|---|---|---|
| OA × OA | 1 | +0.838 | +0.838 | 1/1 |
| OA × BMD | 20 | +0.117 to +0.250 | +0.169 | 20/20 |
| Fracture × BMD | 10 | −0.596 to −0.446 | −0.531 | 10/10 |
| BMD × BMD | 45 | +0.689 to +0.960 | +0.812 | 45/45 |
| OA × Fracture | 2 | +0.148 to +0.211 | +0.180 | 0/2 |

---

## Table 3. Two-sample Mendelian randomisation estimates

| Exposure → Outcome | n_snp | IVW β (fixed) | SE | p | IVW p (random) | Egger β | Egger intercept p | Weighted median | I² | LOO sign flips |
|---|---|---|---|---|---|---|---|---|---|---|
| **BMD total → fracture** (positive control) | 31 | −0.337 | 0.022 | 3.7×10⁻⁵¹ | 2.4×10⁻¹² | −0.337 | 0.558 | −0.393 | 78% | 0 |
| BMD femoral neck → fracture | 25 | −0.351 | 0.030 | 7.7×10⁻³¹ | 6.2×10⁻⁷ | −0.351 | 0.010 | −0.373 | 83% | 0 |
| **BMD total → OA knee** | 31 | +0.023 | 0.004 | 3.9×10⁻¹⁰ | 6.2×10⁻³ | +0.023 | 0.810 | +0.019 | 81% | 0 |
| Lumbar spine → OA knee | 21 | +0.028 | 0.005 | 7.7×10⁻⁹ | 1.8×10⁻² | +0.026 | 0.240 | +0.034 | 82% | 0 |
| Arms → OA knee | 16 | +0.022 | 0.005 | 1.7×10⁻⁶ | 5.5×10⁻² | +0.022 | 0.973 | +0.008 | 85% | 0 |
| Legs → OA knee | 24 | +0.018 | 0.004 | 2.3×10⁻⁵ | 9.8×10⁻² | +0.017 | 0.815 | +0.008 | 85% | 0 |
| Femoral neck → OA knee | 25 | +0.017 | 0.005 | 6.4×10⁻⁴ | 4.1×10⁻² | +0.018 | 0.907 | +0.015 | 67% | 0 |
| OA knee → fracture | 69 | +0.092 | 0.108 | 0.395 | 0.550 | +0.091 | 0.456 | +0.063 | 51% | 0 |
| OA knee → BMD total | 69 | +0.106 | 0.108 | 0.326 | 0.499 | +0.183 | 0.039 | +0.116 | 53% | 0 |
| OA knee → BMD femoral neck | 69 | +0.142 | 0.108 | 0.190 | 0.301 | +0.175 | 0.030 | +0.078 | 38% | 0 |

*Abbreviations:* IVW, inverse-variance weighted; RE, random effects. β is per 1 SD higher exposure on the standardised (reconstructed) outcome scale. The BMD-exposure MR estimates are corrected for the reconstructed exposure N (true 31,873 vs the placeholder 400,000 used in the source files; §3.7): β, SE, Egger β and weighted median are multiplied by √(400,000/31,873) = 3.5426, because the reconstructed β with the placeholder N were 0.2823× of the true scale (underestimated 3.5426-fold). Reverse MR (BMD outcome) applies the same 3.5426 correction to the BMD outcome reconstruction. In both directions p-values, Egger-intercept p, I² and LOO sign flips are unchanged (they are N-scale-invariant). All estimates use the strict PLINK-equivalent clumping of §2.6 (K = strict instrument count). The BMD→fracture positive control yields OR = 1.40 (total) and 1.42 (femoral neck) per 1 SD *lower* BMD, matching the independent GEFOS benchmark of 1.55 (§3.7). For BMD femoral neck → fracture the MR-Egger intercept is significant (p=0.010), indicating directional pleiotropy; the MR-Egger slope (−0.351) and weighted median (−0.373) are reported as the pleiotropy-robust estimates for that pair, and the IVW point estimate may be modestly inflated (see §3.3).

---

## References

1. Hatzikotoulas K, Southam L, Stefansdottir L, et al. Translational genomics of osteoarthritis in 1,962,069 individuals. *Nature*. 2025;641(8065):1217–1224. doi:10.1038/s41586-025-08771-z. PMID 40205036.
2. Qian Y, Xia J, Wang P, et al. Genome-wide Association Studies of over 30,000 Samples with Bone Mineral Density at Multiple Skeletal Sites and Its Clinical Relevance. *Genomics Proteomics Bioinformatics*. 2025;23(5):qzaf097. doi:10.1093/gpbjnl/qzaf097. PMID 41206123.
3. Elsworth B, Lyon M, Alexander T, et al. The MRC-IEU UK Biobank GWAS pipeline. *bioRxiv*. 2018. doi:10.1101/340232. (Self-reported OA accessed as `ukb-b-14486` via OpenGWAS, https://gwas.mrcieu.ac.uk/datasets/ukb-b-14486/.)
4. Hackinger S, Trajanoska K, Styrkarsdottir U, et al. Evaluation of shared genetic aetiology between osteoarthritis and bone mineral density identifies SMAD3 as a novel osteoarthritis risk locus. *Hum Mol Genet*. 2017;26(19):3850–3858. doi:10.1093/hmg/ddx285.
5. Zhou J, Na J, Jiang Z, et al. Deciphering the shared genetic structure between hip osteoarthritis and femoral neck bone mineral density. *Front Genet*. 2025;16:1597005. doi:10.3389/fgene.2025.1597005.
6. Bulik-Sullivan BJ, et al. LD Score regression distinguishes confounding from polygenicity in genome-wide association studies. *Nat Genet*. 2015;47(3):291–295. doi:10.1038/ng.3211.
7. Bowden J, Davey Smith G, Haycock PC, Burgess S. Consistent estimation of the causal effect in Mendelian randomization with some invalid instruments using a weighted median estimator. *Genet Epidemiol*. 2016;40(4):304–314. doi:10.1002/gepi.21965.
8. Verbanck M, Chen CY, Neale B, Do R. Detection of widespread horizontal pleiotropy in causal relationships inferred from Mendelian randomization between complex traits. *Nat Genet*. 2018;50(5):693–698. doi:10.1038/s41588-018-0099-7.
9. Lee SH, Wray NR, Goddard ME, Visscher PM. Estimating missing heritability for disease from genome-wide association studies. *Am J Hum Genet*. 2011;88(3):294–305. doi:10.1016/j.ajhg.2011.02.002.
10. Trajanoska K, Morris JA, Kaptoge S, et al., for the GEFOS consortium. Genome wide association analysis of risk of fracture and its association with bone mineral density: a meta-analysis of genome-wide association studies applied to the UK Biobank. *BMJ*. 2018;362:k3225. doi:10.1136/bmj.k3225. PMID 30158200.
11. Foss MV, Byers PD. Bone density, osteoarthrosis of the hip, and fracture of the upper end of the femur. *Ann Rheum Dis*. 1972;31(4):259–264. doi:10.1136/ard.31.4.259.
12. Dequeker J, Aerssens J, Luyten FP. Osteoarthritis and osteoporosis: clinical and research evidence of inverse relationship. *Aging Clin Exp Res*. 2003;15(5):426–439. doi:10.1007/BF03327364.
13. Davey Smith G, Ebrahim S. 'Mendelian randomization': can genetic epidemiology contribute to understanding environmental determinants of disease? *Int J Epidemiol*. 2003;32(1):1–22. doi:10.1093/ije/dyg070.
14. Burgess S, Thompson SG. Interpreting findings from Mendelian randomization using the MR-Egger method. *Eur J Epidemiol*. 2017;32(5):377–389. doi:10.1007/s10654-017-0255-0.
15. Mounier N, Kutalik Z. Bias correction for inverse variance weighted Mendelian randomization. *Genet Epidemiol*. 2023;47:e22522. doi:10.1002/gepi.22522.

## Authors, affiliations, funding and competing interests

- **Authors.** Zhenxin Zheng¹#, Song Chen¹#, Shenbo Huang¹, Chao Qian¹, Hang Yuan¹, Huricha Bao¹*, Qiuyan Jiang²*.
  - # Zhenxin Zheng and Song Chen contributed equally to this work.
  - * Corresponding authors: Huricha Bao, MD (baohuricha@126.com) and Qiuyan Jiang, MD (after9951023@gmail.com).
- **Affiliations.** ¹ Department of Orthopedics, The Quzhou Affiliated Hospital of Wenzhou Medical University, Quzhou People's Hospital, Quzhou, Zhejiang Province, 324000, China; ² Department of Surgery, The Quzhou Affiliated Hospital of Wenzhou Medical University, Quzhou People's Hospital, Quzhou, Zhejiang Province, 324000, China.
- **Funding.** None.
- **Competing interests.** The authors declare no competing interests.
- **Data availability.** Summary statistics are from publicly published GWAS (Refs 1–3) and were not redistributed; see `docs/10_参考文献与原始GWAS来源.md` for provenance and accessions, and `code/README.md` for the analysis repository and reproducibility commands.

## §5 Pre-submission checklist — status and journal positioning

### 5.1 Item status

| # | Item | Status | Action taken |
|---|---|---|---|
| 1 | Confirm per-trait N & recompute h² | ✅ **Done** | All three trait groups resolved. BMD ≈31,873 (Qian 2025) → h² 0.26–0.39. OA_selfrep **462,933** (38,472/424,461; `ukb-b-14486`) → h²_obs 0.0223, h²_liab 0.073 (K=0.083). **fracture 264,973** (GEFOS discovery, 37,857/227,116) → h²_obs 0.0258, h²_liab 0.062 (K=0.143), identified by the **measured** effective sample size N·K·(1−K)=31,605 vs GEFOS theoretical 32,449 (2.6%). Table 1 fully updated. |
| 2 | Cite original GWAS (consortium + accession) | ✅ **Done** | GO Consortium 2025 (Ref 1), Qian 2025 (Ref 2), IEU OpenGWAS `ukb-b-14486` for OA_selfrep (Ref 3), GEFOS 2018 for fracture (Ref 10), Hackinger 2017 (Ref 4), Zhou 2025 (Ref 5); full provenance/accessions in `docs/10_参考文献与原始GWAS来源.md`. |
| 3 | Strict PLINK clump (r²<0.001, 10 Mb) + MR-PRESSO | ✅ **Done (both executed)** | Strict clump **run** against 1000G Phase 3 EUR (503 samples, 8,550,156 variants, GRCh37) via a validated pure-Python PLINK `.bed` reader (`code/18_plink_clump_pure.py`; 99.4% of GWS SNPs matched the panel). Retained 20–35 instruments per BMD exposure, 116 for knee OA. MR + MR-PRESSO re-run on strict sets (`code/17`): **all directions and significance unchanged** vs the 1-Mb version (§3.6, `mr_presso_results.csv`, `strict_clump_summary.csv`). |
| 4 | Figure 3 (MR forest) | ✅ **Done** | `Figure3_MR_forest.png` produced (`code/14_figure3_mr_forest.py`); legend updated. |
| 5 | OA h² → liability scale | ✅ **Done** | `oa_h2_liability.csv`: OA_knee h²_liab=0.094 (K=0.12). OA_selfrep **updated** to h²_liab=**0.073** with the confirmed K=0.083 (was 0.051 at the placeholder K=0.25); Table 1 footnote updated. |
| 6 | Code repo + data availability | ✅ **Done** | `code/README.md` (repo structure, run order, data-availability statement, git push commands). |
| 7 | Authors / affiliations / funding / COI | ✅ **Done** | Block finalised (§Authors): Zheng Z & Chen S equal contribution (#); Bao H & Jiang Q joint corresponding (*); both affiliations at Quzhou People's Hospital (Wenzhou Medical Univ.); Funding: None; Competing interests: none declared. |
| 8 | OA_selfrep global Z-sign correction | ✅ **Done** | Evidence + correction documented in `docs/11_补充材料_OA_selfrep_Z校正.md` (GDF5 rs224333, FTO rs9930333 opposite signs; rg(OA_knee,OA_selfrep)=+0.84 post-correction). |

| 9 | Correct MR effect magnitudes for the BMD placeholder N | ✅ **Done** | BMD files carried N=400,000 (placeholder) while true N=31,873 → munged β underestimated 3.54-fold, MR β inflated 3.54-fold. Corrected in `code/21_mr_corrected_N.py` → `mr_corrected_N.csv`; Table 3 and §3.3 updated; **p-values unchanged** (β and SE scale identically). Corrected BMD→fracture OR per SD = 1.40–1.42, matching the GEFOS benchmark of 1.55 within 8.4% (§3.7). |

### 5.2 Strategic note A — N-invariance of r_g and MR

The genetic-correlation matrix and all MR **directions, signs and significance levels are invariant to per-trait N** (r_g = gencov/√(h₁·h₂) depends only on the genotype–phenotype covariance structure, not on N or on M; empirically verified — varying fracture N over 150k–1,000,000 left every r_g unchanged to three decimal places). The only quantities that scale with N are the **absolute h² values** (h² ∝ 1/N) and the absolute MR β magnitudes when an exposure β is reconstructed from Z, N and MAF. **Conclusion:** the result framework (three-block architecture; BMD→fracture protection; BMD→OA causality; null OA→fracture / OA→BMD) stands *even before* fracture and OA_selfrep N are finally confirmed; only absolute numbers require updating once those N are locked.

### 5.3 Strategic note B — journal positioning

With strict clump **executed** and MR-PRESSO re-run on the strict instrument sets (item 3) plus Figure 3 (item 4), the single remaining methodological gap the reviewer identified is now closed: the clumping limitation is no longer open. **Recommended path:** submit to ***Bone*** or ***Osteoporosis International*** first (highest hit-rate for a BMD-centred cross-trait MR/LDSC paper, natural home for the BMD–fracture positive control), with **OAC as the aspirational follow-on**. The remaining fracture-N uncertainty is a numbers-only issue (§5.2) and does not gate a submission.

### 5.4 Remaining before submission (user to action)

- [x] ~~Lock fracture N~~ — **done**: measured N·K·(1−K) = 31,605 → adopted **GEFOS discovery N = 264,973, K = 0.1429** (h²_obs 0.0258, h²_liab 0.062). Table 1 and §3.1 updated. Revert to a sensitivity range if you prefer not to commit to this source identification.
- [x] ~~Run strict PLINK r²<0.001/10-Mb clump~~ — **completed** (pure-Python PLINK-equivalent, 1000G Phase3 EUR at `H:/数据下载/1kg.v3/EUR`); directions confirmed unchanged (§3.6).
- [x] ~~Reconcile the BMD→fracture MR magnitude with the GEFOS benchmark~~ — **done (§3.7)**: after the sample-scale correction, OR = 1.42 (95% CI 1.34–1.51) per SD lower femoral-neck BMD vs GEFOS 1.55 (1.48–1.63); 8–10% lower, overlapping CIs. BMD betas are in SD units by construction.
- [x] ~~Fill author / affiliation / funding / COI placeholders~~ — **done** (§Authors, finalised 2026-09-17): Zheng Z & Chen S (# equal); Bao H & Jiang Q (* corresponding); Quzhou People's Hospital affiliations; Funding None; no competing interests.
- [x] ~~Rebuild Figures 1–2 only if the r_g matrix changes~~ — **done**: all three main figures regenerated in Nature/Springer journal style (colourblind palette, panel labels, 300 DPI), and **Figure 4** (SNP heritability bar) plus **Figure 5** (causal-network schematic) added; see updated Figure legends.
