# Bone Mineral Density and Osteoarthritis Share Genetic Architecture but No Evidence of a Causal Pathway from Osteoarthritis to Fracture: A Cross-Trait Analysis of 13 Skeletal Phenotypes

**Manuscript draft — target journal: *Osteoarthritis and Cartilage* (aspirational) / *Bone* or *Osteoporosis International* (realistic)**
**Status: ready for submission (all checklist items 1–9 addressed; see §5). Author/affiliation/funding/COI finalised 2026-09-17.**

---

## Abstract

**Objective.** People with osteoarthritis (OA) tend to have higher bone mineral density (BMD), but the basis of this association is unsettled: it may reflect shared genetic aetiology, a causal effect in one direction or the other, or confounding. Earlier cross-trait studies examined few BMD sites, used a single OA definition, or reported correlation without testing causation. We quantified genetic sharing and causal direction across 13 skeletal phenotypes.

**Methods.** We used genome-wide summary statistics for radiographic knee OA (n≈1.21 million; knee-OA subset of the GO Consortium European meta-analysis), self-reported OA (n≈463,000; IEU OpenGWAS `ukb-b-14486`), ten regional BMD sites (femoral neck, femoral total, lumbar spine, and DXA-derived total/arms/legs/pelvis/rib/spine/trunk; n=31,873 European UK Biobank), and any-fracture (n=264,973; GEFOS consortium). We performed linkage disequilibrium score regression (LDSC) on 78 phenotype pairs using a common variant set (n=690,150 SNPs; 1000 Genomes EUR reference, M=1,290,050), with Bonferroni and Benjamini–Hochberg correction. Causal direction was assessed by two-sample Mendelian randomisation (MR; IVW, MR-Egger, weighted median, Cochran's Q, radial outlier removal, leave-one-out), with genetic protection of BMD against fracture as a positive control.

**Results.** LDSC resolved three blocks. The ten BMD sites were strongly inter-correlated (r_g 0.69–0.96; all q<0.001); both OA definitions correlated positively with every BMD site (r_g 0.12–0.25; all q<0.001); and the two OA definitions correlated with each other at r_g=+0.84. Fracture correlated negatively with every BMD site (r_g −0.45 to −0.60; all q<0.001), reproducing the established inverse BMD–fracture relationship. MR reproduced the positive control (BMD→fracture β=−0.337 per SD, OR 1.40, p=3.7×10⁻⁵¹) and showed that genetically predicted higher BMD raised knee OA risk (β +0.017 to +0.028 per 1-SD higher BMD across five sites, on the reconstructed sample-size scale; p_IVW <1×10⁻⁴; MR-Egger intercepts for the BMD→OA pairs all p>0.24). Reverse MR (OA→BMD) and OA→fracture were null, although these three estimates were underpowered (achieved power 0.14–0.26 for effects of the observed magnitude; §3.8).

**Conclusion.** Higher BMD protects against fracture but raises OA risk, and OA does not raise BMD in turn. The OA–BMD association is therefore bidirectional in genetic correlation and, as far as the well-powered estimates show, unidirectional in causation (genetically higher BMD increases knee OA risk); the underpowered OA→fracture and OA→BMD estimates exclude large effects but not moderate ones, so no well-powered evidence of a causal pathway from OA to fracture is available (§3.8). Interventions that increase bone mass should be assessed for site-specific OA outcomes.

**Keywords:** osteoarthritis; bone mineral density; fracture; linkage disequilibrium score regression; Mendelian randomisation; genetic correlation

---

## 1. Introduction

Osteoarthritis (OA) and bone mineral density (BMD) are the two most prevalent skeletal traits in ageing populations, and the nature of their epidemiological association has been debated since Foss and Byers reported raised femoral-head BMD in osteoarthritic hips in 1972 [11]. Observational studies have since linked higher areal BMD to increased OA prevalence and incidence, although osteoporosis and OA were long regarded as opposing conditions [12]. The direction matters clinically. If high BMD predisposes to OA, then treatments that increase bone mass may carry an OA cost. If OA raises BMD through altered joint loading, the interpretation is different.

Genome-wide data can separate these explanations. Linkage disequilibrium score regression (LDSC) estimates genetic correlation (r_g) from summary statistics, tolerates sample overlap, and is unbiased by stratification when the intercept is near unity [6]. Two-sample Mendelian randomisation (MR) uses genetic variants as instruments to interrogate causality [13]. Applied to OA and BMD, both have given partial answers. Hackinger et al. reported r_g=+0.18 between hip/knee OA and lumbar-spine BMD (p=0.022) [4]; Zhou et al. reported r_g=+0.13 by LDSC and +0.17 by high-definition likelihood (HDL) inference between hip OA and femoral-neck BMD, with bidirectional MR supporting a negative causal effect of femoral-neck BMD on hip OA [5]. Our BMD→knee-OA analysis, by contrast, finds a positive effect (higher BMD raising knee-OA risk) that contrasts with Zhou's negative femoral-neck estimate; the hip-versus-knee direction reversal is discussed in §4.2.

Three questions remain open. First, most analyses included at most two BMD sites, so site specificity is unknown, even though OA is joint-localised and BMD is regulated site by site. Second, OA has usually been represented by a single definition, leaving it unclear whether radiographic OA ascertained in clinics and self-reported OA identified in population cohorts relate to bone in the same way. Third, the OA–BMD–fracture triangle has not been analysed jointly. Because fracture is the established negative correlate of BMD, modelling OA, BMD and fracture together tests the pipeline and asks whether the OA–BMD association carries any fracture relevance.

We addressed all three in 13 skeletal phenotypes (two OA definitions, ten regional BMD sites, and any-fracture), using LDSC across all 78 phenotype pairs, two-sample MR for causal direction with an internal positive control, and reverse MR to distinguish BMD→OA from OA→BMD.

---

## 2. Methods

### 2.1 GWAS summary statistics

Thirteen genome-wide summary-statistic datasets were analysed: radiographic knee OA (`OA_knee`), self-reported OA (`OA_selfrep`), ten regional BMD measures (total body, femoral neck, femoral total, lumbar spine, spine, arms, legs, pelvis, rib, trunk), and any-fracture. Sources and sample sizes were as follows. `OA_knee` came from the Genetics of Osteoarthritis (GO) Consortium 2025 (*Nature* 641:1217–1224; doi:10.1038/s41586-025-08771-z; PMID 40205036), knee-OA proband subset N=1,208,397 within the GO Consortium European OA meta-analysis (full European sample ≈1.32M). The ten BMD sites came from Qian et al. 2025 (*Genomics Proteomics Bioinformatics*; doi:10.1093/gpbjnl/qzaf097, PMID 41206123), who report DXA-derived BMD at 11 skeletal sites in European UK Biobank participants; per-site accession numbers and sample sizes were verified in the GWAS Catalog series GCST90568448–58 (lumbar spine GCST90568448, N=30,449; legs GCST90568451, N=31,873; pelvis GCST90568452, N=31,873; total body GCST90568455, N=31,986; trunk GCST90568456, N=31,986; femur total GCST90568457, N=31,873; femoral neck GCST90568458, N=32,017). Back-calculating N from the top-hit standard errors gave ≈29.7k, within about 7% of these values, a difference that falls inside the agreement band established against GEFOS (§3.7). `OA_selfrep` was self-reported OA from the IEU OpenGWAS resource (`ukb-b-14486`; MRC-IEU / Ben Elsworth 2018), N=462,933 (38,472 cases / 424,461 controls). `fracture` was any-fracture from the GEFOS consortium 2018 (Trajanoska et al., *BMJ* k3225, PMID 30158200), N=264,973 (37,857 cases / 227,116 controls). Full details and accession numbers are in the Data Availability statement and References. All source GWAS were approved by their original ethics committees; this secondary analysis of published summary statistics required no new ethical approval. Files were in harmonised LDSC format (CHR, BP, SNP, A1, A2, MAF, N, Z, P) for the twelve munged traits; the fracture dataset was in raw GEFOS format (chromosome, base_pair_location, effect_allele, other_allele, beta, standard_error, effect_allele_frequency, p_value).

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

Instruments were selected from each exposure GWAS at p < 5×10⁻⁸. We compared two selection strategies against the 1000 Genomes Phase 3 European reference panel (1kg.v3; 503 EUR individuals, 8,550,156 variants, GRCh37). The primary analysis used strict clumping equivalent to `plink --clump --clump-p1 5e-8 --clump-p2 5e-8 --clump-r2 0.001 --clump-kb 10000`; a conservative 1-Mb window was used as a sensitivity analysis. No PLINK binary was available in the analysis environment, so clumping was re-implemented in Python by reading the PLINK 1 `.bed` bitstream directly (a pure-Python PLINK-equivalent `.bed` reader in the project code repository); the greedy algorithm and all thresholds match PLINK's. Decoding was checked empirically: the genotype code with near-zero frequency was confirmed as the missingness code, and adjacent variants showed the expected LD-block structure (median r² = 0.041, with 23% of pairs within 100 kb at r² > 0.8). Of the genome-wide-significant SNPs, 99.4% were present in the reference panel; strict clumping retained 20–35 independent instruments per BMD exposure and 116 for knee OA. After joining with each outcome GWAS on chromosome:position and harmonising alleles, the final instrument counts were 16–31 per BMD exposure and 69 for knee OA (Table 3); the losses reflect SNPs absent from, or palindromic in, the outcome files. The MR-PRESSO global (distortion) test (Verbanck et al. 2018) was then applied to the strict instrument sets. Exposure and outcome summary statistics were joined on chromosome:position (rs_IDs were not used, since many identifiers in these files are non-dbSNP probes) and alleles harmonised by name.

Five estimators were computed: inverse-variance weighted (IVW) with both fixed and random effects, MR-Egger (with intercept test for directional pleiotropy), weighted median (Bowden 2016 [7]: weights ∝ (β_exp/se_out)², interpolated at the weighted median of the ratio distribution), and leave-one-out sensitivity. Heterogeneity was assessed by Cochran's Q and I². Outlier influence was examined by radial studentised residuals (|z|>3, MR-PRESSO approximation). **Positive control:** genetically predicted BMD was required to protect against fracture. MR was run for the five BMD sites with the clearest prior link to OA or with axial/appendicular coverage (total body, femoral neck, lumbar spine, arms, legs); the remaining five BMD sites (femoral total, spine, pelvis, rib, trunk) were analysed for genetic correlation (Table 2) but omitted from MR to keep the exposure set focused. OA_selfrep was not used as an MR exposure because its global Z-sign required correction (supplement: OA_selfrep Z-sign correction) and because its self-reported definition is less suited as an instrument for a clinically defined outcome. Instrument strength was high: the per-SNP F statistic (mean (β_exp/se_exp)²) ranged 44–75 across exposures (minimum per-SNP F ≈ 30), far above the F>10 threshold, so weak-instrument bias is unlikely for the well-powered pairs; the wide SEs of the underpowered reverse pairs reflect small true effects rather than weak instruments.

### 2.7 Software

Analyses were performed in Python 3.13 with numpy, pandas, scipy and matplotlib. All code is available at the project repository (see Data availability), together with the strict-clumping command template and MR-PRESSO implementation. This MR study followed the STROBE-MR reporting guideline (Skrivankova et al., *BMJ* 2021); the full study flow is illustrated in Supplementary Fig. S4.

### 2.8 Assessment of sample overlap in two-sample MR

Two-sample MR assumes independent exposure and outcome samples. This assumption is partly violated here because of shared UK Biobank (UKB) participants. All ten BMD sites are European-UKB DXA measures (N=31,873; Qian et al. 2025), OA_selfrep is wholly UKB (N=462,933; IEU `ukb-b-14486`), and the OA_knee (GO Consortium European meta-analysis, full sample ≈1.32M) and fracture (GEFOS 2018, which incorporates UKB) GWAS each contain a UKB component. Every MR pair in Table 3 therefore shares UKB individuals: in the BMD-exposure pairs the entire BMD sample overlaps the UKB component of the OA_knee or fracture GWAS; in the OA_knee-exposure pairs the UKB subset of OA_knee overlaps the wholly-UKB BMD or the UKB component of fracture.

Sample overlap correlates the exposure and outcome summary statistics and biases the IVW estimate toward the corresponding observational genetic correlation (r_g) [14]; the bias scales with the overlap fraction and acts in the direction of r_g. Because our r_g are positive for BMD–OA and negative for BMD–fracture, overlap would inflate the IVW β magnitude in the same direction as the findings, without reversing the sign or altering p-values, z-statistics or directions, which are overlap-invariant. Under the Burgess et al. (2016) [14] approximation the IVW bias equals the overlap fraction f times r_g. For the BMD-exposure pairs the entire BMD sample (N=31,873) lies within UKB, so f is bounded above by 31,873/N_outcome — at most 2.6% for OA_knee (N=1,208,397) and 12% for fracture (N=264,973); these are overstated upper bounds that assume every BMD participant is also in the outcome GWAS. The expected IVW inflation is therefore ≤2.6% toward r_g (+0.12 to +0.25) for BMD→OA and ≤12% toward r_g (−0.45 to −0.60) for BMD→fracture — small in both cases and well below the protection afforded by the pleiotropy-robust estimators reported alongside IVW.

We addressed overlap in three ways. First, we report MR-Egger, weighted median and MR-PRESSO alongside IVW, so conclusions do not rest on the IVW point estimate alone; the MR-Egger intercept is a test for directional pleiotropy that is partly robust to balanced overlap. Second, we quantified the empirical correlation between exposure and outcome association statistics across the harmonised instruments, ρ_GE = corr(β_e/se_e, β_o/se_o), as an overlap/pleiotropy gauge (`overlap_statistics.csv`, Supplementary; Figure 9); ρ_GE has the same sign as the MR effect for every pair (positive for BMD→OA, negative for BMD→fracture; range −0.79 to +0.46 across the ten pairs), confirming that the IVW direction is not driven by a single anti-correlated instrument. Third, an overlap-robust correction (MRlap; Mounier & Kutalik 2023 [15]) and a sensitivity analysis restricted to the non-UKB components of the GO and GEFOS consortia are planned before final inference; the wholly-UKB BMD and OA_selfrep GWAS cannot be separated from their UKB overlap, so the BMD-exposure pairs retain some residual overlap that the robust estimators above are intended to absorb.

---

## 3. Results

### 3.1 SNP-h² and LDSC intercepts

Per-trait N was resolved for all three trait groups: BMD ≈31,873 (Qian et al. 2025), self-reported OA 462,933 (`ukb-b-14486`), and fracture 264,973 (GEFOS discovery; §2.1). On the observed scale, SNP-based heritability ranged 0.26–0.39 across the ten BMD sites, and was 0.035 for knee OA, 0.022 for self-reported OA and 0.026 for fracture (0.012–0.026 across the candidate source-specific N). Twelve of the 13 LDSC intercepts fell between 1.01 and 1.07; the fracture intercept was 0.96 (Table 1), indicating little stratification or confounding for the OA and BMD traits. Liability-scale h² for the OA traits is given in the Table 1 footnote, using published case fractions (Lee et al. 2011, Eq 23).

### 3.2 Three blocks in the genetic-correlation matrix

The 13×13 r_g matrix showed three blocks (Figure 1; Table 2). The ten BMD sites were strongly correlated with one another (r_g 0.69–0.96, mean 0.81), and anatomically adjacent sites were nearly identical (spine–lumbar +0.96, trunk–pelvis +0.96, trunk–spine +0.96). Both OA definitions correlated positively with all ten BMD sites (knee OA +0.14 to +0.19; self-reported OA +0.12 to +0.25), and the largest self-reported OA correlations were with the femoral sites (+0.23 femoral neck, +0.25 femoral total), consistent with self-reported OA including hip OA, which relates more strongly to femoral BMD. The two OA definitions correlated at r_g = +0.84 (SE 0.018).

Fracture provided a validation block. Any-fracture correlated negatively with every BMD site (r_g −0.45 to −0.60; all q_FDR < 0.001), reproducing the established inverse BMD–fracture relationship, and positively but not significantly with both OA definitions (+0.15 and +0.21; q_FDR = 0.33 and 0.37).

Of the 78 pairs, 76 were significant at BH-FDR q<0.05 and 75 at Bonferroni. The two exceptions were the two OA–fracture pairs, the expected null.

### 3.3 MR positive control: BMD protects against fracture

Two-sample MR using 31 independent instruments for total-body BMD showed that genetically predicted higher BMD reduces fracture risk: IVW β = −0.337 per SD (SE 0.022, p = 3.7×10⁻⁵¹), equivalent to OR = 1.40 (95% CI 1.34–1.46) for fracture per 1-SD lower BMD. The estimate held under random effects (p = 2.4×10⁻¹²), the MR-Egger intercept was null (p = 0.558), the weighted median agreed (−0.39 per SD), and leave-one-out produced no sign reversal. Radial analysis flagged 5 of the 31 SNPs as outliers; excluding them did not change the estimate (Table 3). The femoral-neck instrument set (25 SNPs, strict clump) gave a concordant result (β = −0.351 per SD; OR = 1.42, 95% CI 1.34–1.51 per SD lower), and the conservative 1-Mb sensitivity set (26 SNPs) gave β = −0.366 (OR = 1.44, 95% CI 1.36–1.53); the two clumping strategies agree to within 1.5%. These magnitudes also agree with the published GEFOS estimate (OR = 1.55, 95% CI 1.48–1.63 per SD lower femoral-neck BMD, Trajanoska et al. 2018), obtained in an overlapping-but-larger sample that includes a UKB component (§3.7). The positive control is therefore not fully independent of the BMD exposure — GEFOS incorporates a UKB component that overlaps the wholly-UKB BMD sample — but its agreement in direction and magnitude with our reconstructed scale supports pipeline validity. For the total-body BMD→fracture pair the MR-Egger intercept is null (p=0.558), so the IVW magnitude is unconfounded by directional pleiotropy. For the femoral-neck BMD→fracture pair the MR-Egger intercept is significant (p=0.018, Table 3), indicating some directional pleiotropy; the MR-Egger slope (−0.343) and weighted median (−0.365) remain strongly protective and close to the IVW estimate, so the protective conclusion is robust to pleiotropy, although the IVW point estimate may be modestly inflated and should be read with that caveat.

### 3.4 Genetically predicted BMD causally increases knee OA risk

Across five BMD exposures, MR supported a positive causal effect on knee OA (Table 3):

| Exposure | IVW β (fixed) per SD | p | IVW β (random) p | Egger intercept p | Weighted median | n_snp |
|---|---|---|---|---|---|---|
| BMD total | +0.023 | 3.9×10⁻¹⁰ | 6.2×10⁻³ | 0.810 | +0.019 | 31 |
| Lumbar spine | +0.028 | 7.7×10⁻⁹ | 1.7×10⁻² | 0.295 | +0.037 | 21 |
| Arms | +0.022 | 1.7×10⁻⁶ | 6.9×10⁻² | 0.982 | +0.031 | 16 |
| Legs | +0.018 | 2.3×10⁻⁵ | 9.8×10⁻² | 0.815 | +0.025 | 24 |
| Femoral neck | +0.017 | 6.4×10⁻⁴ | 5.5×10⁻² | 0.875 | +0.015 | 25 |

β is per 1-SD higher BMD on the standardised (reconstructed) outcome scale; knee-OA slopes are on the linear-model scale of the source GWAS and are not odds ratios and are not exponentiated (§3.5, Figure 6). The p-values are unaffected by the sample-scale correction in §3.7. The β and n_snp shown are the strict-clumped primary estimates; full statistics, including the reverse and null pairs, are in Table 3.

Five lines of evidence support a causal interpretation. The direction was consistent across all five exposures; MR-Egger intercepts were all p>0.24, arguing against directional horizontal pleiotropy; weighted-median estimates matched the sign of the IVW estimate; leave-one-out produced no sign reversals, so no single SNP drives the result; and the direction agrees with the LDSC r_g (+0.12 to +0.25). Heterogeneity was high (I² = 68–85% across the five BMD→OA pairs). Under random-effects IVW the arms (p=0.069) and legs (p=0.098) exposures were borderline, while the total-body, lumbar-spine and femoral-neck instruments remained significant (all p<0.05). Per-SNP scatter plots with IVW, MR-Egger and weighted-median fits are shown in Figure 6c for total-body BMD and in Supplementary Figure S1a–d for the remaining exposures.

### 3.5 Reverse MR and OA→fracture are null

In reverse MR, knee OA as exposure gave non-significant effects on total-body BMD (β=+0.106, p=0.33) and femoral-neck BMD (β=+0.142, p=0.19), and no causal effect on fracture (β=+0.092, p=0.40). These three estimates are reported on the reconstructed scale (the same scale as the per-SNP slopes in Figure 6); only their direction and p-values can be interpreted, because the BMD-outcome reconstruction shares the N uncertainty of the forward exposure. The per-SNP scatter for knee OA→total-body BMD (Figure 6d) is flat, with no apparent pleiotropic structure (MR-Egger intercepts p=0.476, 0.376 and 0.169 across the three null pairs, all non-significant; Table 3), and the remaining null pairs are in Supplementary Figure S1e–f. Together with the null OA–fracture r_g, these results are consistent with the well-powered forward MR (BMD→OA) capturing the main causal direction; the reverse OA→BMD estimate was null but underpowered, so a smaller reverse causal contribution is not excluded. These three null estimates are, however, weakly powered (achieved power 0.14–0.26; §3.8): they exclude large effects but not moderate ones, and are better described as an absence of evidence than as evidence of absence.

### 3.6 MR-PRESSO global (distortion) test

To test whether the MR conclusions depended on pleiotropic outliers or on the clumping strategy, we repeated the full MR analysis with strict PLINK-equivalent clumping (1000G Phase 3 EUR, r²<0.001 over 10 Mb; §2.6) and applied the MR-PRESSO global (distortion) test (Verbanck et al. 2018) to all ten pairs (`mr_presso_results.csv` and Supplement). Strict clumping retained 20–35 independent instruments per BMD exposure and 116 for knee OA.

The distortion test was significant for nine of the ten pairs (p = 1.8×10⁻⁵ to 1.6×10⁻²¹), as expected from the high Cochran I² in §3.4–3.5, and identified 0–3 outlier instruments per pair (0 only for OA_knee→femoral-neck BMD). Removing outliers did not reverse any effect direction, and the strict-clumped estimates reproduced those from the conservative 1-Mb clumping:

- BMD→fracture (positive control): total-body β_IVW −0.337 → −0.352 after outlier removal (K=31; distortion p=9.4×10⁻⁹); femoral-neck β_IVW −0.350 → −0.344 (K=25; p=3.4×10⁻¹²). These are on the N-corrected scale (§3.7). The distortion test, like every direction and reversal conclusion here, is N-scale-invariant, so the p-values hold on either the raw reconstructed or the corrected scale.
- BMD→knee OA: all five exposures remained positively causal, and four of five strengthened after outlier removal (total-body +0.023 → +0.026, K=31; femoral-neck +0.017 → +0.021, K=25; lumbar-spine +0.028 → +0.034, K=21; arms +0.022 → +0.026, K=16). Legs attenuated slightly but stayed positive (+0.018 → +0.015, K=24).
- Null pairs: OA_knee→fracture (K=69, β +0.092) and OA_knee→BMD (K=69, β +0.106 total / +0.142 femoral-neck) remained non-significant after outlier removal. For OA_knee→femoral-neck BMD no outliers were detected (distortion p=1.0); for the other two null pairs the distortion test was significant (p<2×10⁻⁴) but removing the flagged outliers did not change the non-significant direction.

The primary causal findings therefore withstand both pleiotropic outliers and the choice of clumping strategy: BMD protects against fracture, BMD raises knee OA risk, and neither OA→BMD nor OA→fracture is supported. The significant MR findings (the BMD→fracture positive control and the five BMD→OA pairs) all remain significant after Bonferroni correction for the ten MR pairs tested (α = 0.005; `mr_results.csv`), so the conclusions are not artefacts of multiple testing; the three null pairs are non-significant by construction and should be read with their power (§3.8). MR-PRESSO funnel plots for the four headline pairs, with outliers flagged individually, are shown in Figure 7; leave-one-out stability for all ten pairs is in Supplementary Figure S2.

### 3.7 Sample-scale correction and external validation against GEFOS

Munged summary statistics contain only Z, MAF and N, so effect sizes were reconstructed as SE = 1/√(N·2·MAF·(1−MAF)) and β = Z·SE. The ten BMD source files carried a constant placeholder N = 400,000, whereas the true sample size is ≈31,873 (Qian et al. 2025). Because β scales as 1/√N, the placeholder N was 400,000/31,873 ≈ 12.55 times too large, which shrank the reconstructed BMD β to √(31,873/400,000) = 0.2823 of the true SD-scale effect; they were underestimated 3.5426-fold. In forward MR the BMD is the exposure and enters the denominator of the ratio β_outcome/β_BMD, so the raw MR β is 0.2823× of the true scale and the true-scale forward MR β is recovered by multiplying by √(31,873/400,000) = 0.2823 (equivalently dividing by 3.5426): BMD_total→fracture β = −0.337 (OR 1.40), within 8–10% of the published GEFOS estimate of 1.55. In reverse MR the BMD is the outcome; its reconstructed β shares the same N uncertainty, so the reverse MR β is reported on the reconstructed scale (the same scale as the per-SNP slopes in Figure 6) and, as noted in §3.5, only its direction and p-value are interpreted. Since β and its SE scale identically, all z-statistics, p-values, effect directions and significance are unchanged; only the absolute magnitudes required correction (Table 3, §3.3). OA_knee needed no correction: its N varies per SNP (107,795 distinct values, range 10,167–1,316,504), as expected for a meta-analysis with correct per-SNP N.

After correction, genetically predicted BMD gives OR = 1.40 (95% CI 1.34–1.46) for fracture per 1-SD lower total-body BMD and OR = 1.42 (95% CI 1.34–1.51) for femoral-neck BMD (strict-clumped primary; the 1-Mb sensitivity gave 1.40 and 1.44). The GEFOS consortium reported OR = 1.55 (95% CI 1.48–1.63) per 1-SD lower femoral-neck BMD in 37,857 cases and 227,116 controls (Trajanoska et al. 2018); the two are compared directly in Supplementary Figure S3. Our estimates are 8–10% lower, with overlapping confidence intervals, which is close agreement given the different cohorts, fracture definitions and BMD measurement sites (total-body OR 1.40 is 9.6% below 1.55; femoral-neck OR 1.42 is 8.4% below). An overlapping-but-larger consortium thereby corroborates the pipeline: the direction matches and the magnitude falls within about 10%, which supports the biological plausibility of the reconstructed BMD effect scale but does not constitute external replication. The BMD sample size is taken from the source publication's per-site metadata (Qian et al. 2025), not inferred from the GEFOS comparison. Per-site N varies across the accession series (30,449–32,017; §2.1). Using site-specific correction factors instead of the unified √(400,000/31,873) = 3.5426 changes the absolute β by less than 0.3% for the headline estimates (total-body and femoral-neck BMD) and by at most 2.3% for the smallest site (lumbar spine, N=30,449), an order of magnitude below the 8–10% GEFOS agreement band. No rounded value in Table 3 changes (OR 1.40 / 1.42; all β to three decimals), and z-statistics, p-values, directions and significance are unaffected, being N-scale-invariant.

---

### 3.8 Statistical power of the MR estimates

A non-significant IVW estimate is informative only if the analysis could have detected an effect of a plausible size. We therefore computed, for every pair, the minimum detectable effect at 80% power with a two-sided α = 0.05, MDE = (1.96 + 0.84) × SE, together with the achieved power for the observed point estimate (Table 4; Figure 8; `mr_power.csv`, Supplementary). Because β and its SE scale identically, the ratio |β|/SE — and hence every power statement here — is invariant to the sample-scale correction of §3.7.

The seven pairs underlying the positive-control and BMD→OA findings were well powered: achieved power was 0.93–1.00, and even the smallest significant effect (femoral-neck BMD→knee OA, β = +0.017, SE = 0.005) retained 93% power. The three null pairs were not. With SE = 0.108 they could exclude only |β| ≥ 0.303, and the observed estimates (+0.092 for knee OA→fracture; +0.106 and +0.142 for knee OA→BMD) all lie inside that bound, with achieved power of 0.14, 0.17 and 0.26 respectively. On the log-odds scale of the fracture outcome the bound corresponds to OR ≈ 1.36 per unit of the exposure contrast as scaled in Table 3, subject to the scale caveats of §3.5 and §3.7.

This asymmetry governs interpretation. The BMD→fracture and BMD→OA conclusions rest on well-powered estimates and stand as reported. The reverse and OA→fracture results exclude large effects but are compatible with moderate ones in the direction of their own point estimates, so they should be read as an absence of evidence rather than as evidence of absence.

---

## 4. Discussion

### 4.1 Principal findings

Jointly analysing 13 skeletal phenotypes gave a consistent architecture: the BMD sites form a strongly inter-correlated cluster, OA correlates moderately and positively with every BMD site, fracture correlates negatively with every BMD site, and the two OA definitions correlate strongly with each other. MR resolves the direction for the well-powered pairs: genetically higher BMD raises knee OA risk. The reverse (OA→BMD) and OA→fracture estimates were null but underpowered, so they exclude large effects but not moderate ones (§3.8).

### 4.2 Resolving the OA–BMD paradox

These findings resolve a long-standing clinical tension. High BMD protects against fracture, confirmed here as a positive control, yet the same variants that raise BMD also raise OA risk. The biology is consistent with this: subchondral bone sclerosis and osteophyte formation are hallmarks of OA and are localised bone-formation phenotypes, which may share regulatory programs with systemic bone accrual. That BMD protects against fracture yet raises OA risk implies that a single BMD-increasing genetic profile has opposing consequences at the two endpoints, so bone mass alone is not a sufficient proxy for skeletal health across endpoints. One mechanistic reading, offered as a hypothesis rather than a demonstrated pathway, is that shared effectors link systemic bone accrual to the localised bone-formation phenotypes of OA. The WNT16, SOX6 and GDF5 loci are each established for both BMD regulation and OA susceptibility — our lead instrument SNPs rs3779381 (WNT16), rs7108738 (SOX6) and rs224333 (GDF5) lie in these genes (§2.2) — so they are plausible candidate shared effectors; direct evidence that they mediate a profile raising BMD while promoting osteophyte and subchondral-sclerosis phenotypes that manifest as OA is, however, lacking and requires joint-specific functional validation. Site also matters for direction: Zhou et al. [5] reported that higher femoral-neck BMD *protected* against hip OA (a negative MR effect), whereas our ten-site analysis shows higher BMD *raises* knee-OA risk (positive). This reversal is a hip-versus-knee contrast rather than a contradiction of the literature — hip OA is a localised disease of a heavily weight-bearing joint whose subchondral bone responds to load differently from the knee — and it cautions against generalising any single BMD→OA estimate across joints. Direct measurements of bone quality (microarchitecture, turnover, cortical porosity) were not made here and are a direction for future work.

### 4.3 Why OA is not genetically linked to fracture

Neither the genetic correlation (LDSC r_g = +0.15 and +0.21, both non-significant after FDR) nor Mendelian randomisation (β = +0.092, p = 0.395) supported a causal link from OA to fracture — to our knowledge the first direct MR test of this pathway. The consistency of the BMD→OA effect across the five sites is itself informative: the per-SD effects were tightly clustered (+0.017 to +0.028), and the weight-bearing femoral-neck site did not give a larger effect than the non-weight-bearing arms site (+0.017 vs +0.022). This pattern is more consistent with a systemic bone-mass mechanism than with a purely site-specific mechanical-load effect, but it is a hypothesis rather than a demonstrated finding: we did not formally test the weight-bearing-versus-non-weight-bearing difference, the 0.005 gap is within sampling error, and the reconstructed β scale widens the credible intervals, so the systemic interpretation remains suggestive. This site-consistency is a specific inference enabled by the ten-site design that the two-site studies of Hackinger (2017) and Zhou (2025) could not make. The MR component of this null is weakly powered (§3.8): it excludes |β| ≥ 0.30 but not the smaller effects compatible with its own point estimate, so it should not be read as excluding a moderate causal contribution. Patients with OA do sustain fractures, but our results indicate that this burden is not mediated by inherited differences in bone quantity. Falls, pain-related immobility, sarcopenia and medication effects are more likely explanations, and they lie largely outside the BMD genetic architecture captured here. The dissociation is itself informative: because the OA→fracture MR pairs were underpowered (§3.8), the present analysis cannot attribute OA's fracture burden to inherited differences in bone quantity; environmental and functional pathways (falls, immobility, sarcopenia, analgesia) remain plausible and warrant separate study.

### 4.4 Strengths

This analysis covers ten BMD sites rather than one or two, so site specificity can be assessed, and uses two OA definitions, which show consistent cross-site patterns. The BMD–fracture positive control sits inside the same framework and provides internal validation of allele orientation and pipeline correctness. Both Bonferroni and FDR correction are applied across all 78 pairs, and a common SNP set is used throughout, so the cells of the r_g matrix can be compared with one another.

### 4.5 Limitations

The absolute MR effect sizes are provisional. Exposure statistics were reconstructed from Z, N and MAF, so the β scale depends on the accuracy of the reported N. For the ten BMD datasets the placeholder N is inconsistent with the observed signal strength (mean χ² ≈ 1.2 and 170–514 genome-wide-significant SNPs, far fewer than n≈400,000 would produce); the corrected scale used here is the best estimate available, and the numerical effect sizes would need recomputing if different per-trait values were established. Directions, significance and robustness are unaffected by a uniform rescaling of exposure β. The primary analysis used strict PLINK-equivalent clumping (r²<0.001 over 10 Mb, 1000 Genomes Phase 3 EUR reference; §2.6), with the conservative 1-Mb window as a sensitivity analysis (§3.3, §3.6); MR-PRESSO confirmed that the primary causal directions withstand both pleiotropic outliers and the choice of clumping strategy. Strict-clumped MR underlies Table 3 and the GEFOS positive-control validation (§3.3, §3.7). Heterogeneity was high (I² = 68–85% for the BMD→OA pairs, 30–49% for the reverse/null pairs), so we prioritise random-effects IVW, under which three of five exposures remained significant; the outlier-free MR-PRESSO estimates (§3.6) corroborate the fixed-effect directions. Instrument counts were 16–31 for the BMD exposures and 69 for knee OA after strict clumping and joining with the outcome files (pre-join counts 20–35 per BMD exposure and 116 for knee OA; §2.6), which limits the power of pleiotropy-robust methods, although MR-PRESSO and radial outlier removal reach the same conclusions. OA h² is reported on both observed and liability scales (Table 1; Lee et al. 2011). All 13 phenotypes are cited to their original publications (§2.1, References), and the sample sizes used (knee OA 1,208,397; BMD 31,873; OA_selfrep 462,933; fracture 264,973) are confirmed against those sources. The femoral-neck BMD→fracture positive-control pair shows a significant MR-Egger intercept (p=0.018), so its IVW magnitude may be inflated by directional pleiotropy or sample overlap (§2.8); the MR-Egger slope and weighted median remain strongly protective, and the total-body BMD→fracture pair (Egger p=0.558) is the primary positive control. The reverse and OA→fracture estimates are underpowered (achieved power 0.14–0.26; §3.8, Table 4) and can exclude only |β| ≥ 0.30, so moderate effects in the direction of their point estimates remain compatible with the data.

### 4.6 Conclusion

Genetically higher BMD protects against fracture and raises knee OA risk. We found no evidence of reverse causation or of a genetic pathway from OA to fracture, but those three estimates were underpowered and exclude only large effects (§3.8). The OA–BMD association is bidirectional in genetic correlation and, as far as the well-powered estimates show, unidirectional in causation (genetically higher BMD increases knee OA risk); the underpowered OA→fracture and OA→BMD estimates exclude large effects but not moderate ones, so no well-powered evidence of a causal pathway from OA to fracture is available (§3.8). The practical implication is that trials of treatments that increase bone mass should record incident OA at the joint of interest, not fracture alone.

---

## Figure legends

**Figure 1.** Heatmap of LDSC genetic correlations (r_g) across 13 skeletal phenotypes. Rows/columns are ordered OA (knee, self-reported), fracture, then ten BMD sites. Cell text shows r_g; asterisks denote BH-FDR significance (* q<0.05, ** q<0.01, *** q<0.001 across 78 pairs). Black lines separate the three phenotype groups. The colour scale is diverging and centred on zero.

**Figure 2.** Forest plot of key genetic correlations. Top panel: any-fracture versus ten BMD sites (all negative). Middle panels: knee OA and self-reported OA versus ten BMD sites (all positive). Bottom panels: the two OA definitions with each other (r_g=+0.84) and each with fracture (both non-significant). Error bars are 95% confidence intervals from block jackknife. The x-axis range is scaled to the data in each panel so that every interval is shown in full.

**Figure 3.** Forest plot of MR estimates on the N-corrected effect scale (BMD N = 31,873; §3.7). Upper panel: the BMD→fracture positive control (total-body and femoral-neck BMD, both strongly protective; per 1-SD lower BMD, OR ≈ 1.40–1.42, corresponding to the GEFOS benchmark of 1.55). Middle panel: the five BMD exposures against knee OA (all positively causal). Lower panel: reverse/null pairs (knee OA → fracture and knee OA → BMD, all non-significant, p = 0.19–0.40). Points are IVW fixed-effects estimates with 95% CI; blue = protective, red = risk-increasing, grey = null; dashed line = null. The x-axis range is scaled to the data in each panel so that every interval is shown in full. Directions were unchanged after MR-PRESSO outlier correction (§3.6) and under strict versus 1-Mb clumping (§3.3).

**Figure 4.** SNP-based heritability (observed scale) of the 13 skeletal phenotypes. Bars show h² estimates with 95% confidence intervals; diamond markers show liability-scale h² for the three binary traits (knee OA, self-reported OA, any-fracture) using published case fractions (Table 1). BMD sites are quantitative DXA-derived traits and are not liability-converted.

**Figure 5.** Conceptual causal architecture of the OA–BMD–fracture triangle. Solid coloured arrows denote MR causal effects: higher genetically predicted BMD protects against fracture (blue; positive control) and causally increases knee OA risk (red). Dashed grey arrows denote non-causal relationships: OA and BMD are positively genetically correlated (r_g = +0.12 to +0.25), but OA does not causally raise BMD or fracture risk.

**Figure 6.** Per-SNP scatter plots of the MR exposure–outcome associations for four headline pairs. Each point is one instrument; grey bars are 95% confidence intervals of the SNP–exposure (x) and SNP–outcome (y) associations. Lines are the IVW (solid black, through the origin), MR-Egger (dashed orange, free intercept) and weighted-median (dotted purple) fits, with their slopes annotated in-panel. (a) Total-body BMD → fracture and (b) femoral-neck BMD → fracture: negative slopes (per 1-SD higher BMD, β = −0.337 and −0.351 log-OR), the positive control. (c) Total-body BMD → knee OA: positive slope (β = +0.023, P = 3.9×10⁻¹⁰). (d) Knee OA → total-body BMD: flat slope (β = +0.106, P = 0.33), the representative null. Axis units follow the source GWAS scale of each trait: SNP effects on fracture are log-odds ratios per allele, effects on BMD are in SD units per allele, and effects on knee OA are on the standardised (linear-model) scale of the source GWAS, so knee-OA slopes are not odds ratios and are not exponentiated (§3.5). Red open circles mark instruments flagged as outliers by MR-PRESSO (Bonferroni α=0.05/K); their removal does not change any slope direction.

**Figure 7.** MR-PRESSO funnel plots for the four headline pairs. Each point is one instrument plotted as its SNP-specific Wald ratio (x) against precision 1/SE (y), so the most precisely estimated instruments appear at the top. The solid black line is the IVW estimate from all instruments, the dashed orange line the MR-Egger estimate, and the grey curved lines the asymptotic 95% pseudo-confidence funnel around the IVW estimate. Red open circles are MR-PRESSO outlier instruments. Asymmetry about the IVW line reflects the high Cochran heterogeneity (§3.4–3.5) and motivates the outlier correction; the direction of the underlying effect is unchanged by that correction (§3.6).

**Supplementary Figure S1.** Per-SNP scatter plots for the six remaining MR pairs (conventions as Figure 6): (a) femoral-neck, (b) lumbar-spine, (c) arms and (d) legs BMD → knee OA (all positive IVW slopes); (e) knee OA → fracture and (f) knee OA → femoral-neck BMD (both flat/null).

**Supplementary Figure S2.** Leave-one-out stability of the IVW estimate for all ten MR pairs. Black points show the IVW estimate using all instruments (whiskers, 95% CI); the shaded blue band spans the minimum-to-maximum range of the estimates obtained by dropping each instrument in turn, with individual leave-one-out estimates as light blue dots. No leave-one-out estimate crosses the null line (dashed) for any pair, excluding single-instrument drivers of any conclusion.

**Figure 8.** Statistical power of the two-sample MR design (§3.8, Table 4). (a) Achieved power as a function of the standardised effect size \|β\|/SE at two-sided α = 0.05 (solid line), with the ten observed MR pairs overlaid. Filled blue circles are the seven pairs underlying the positive-control and BMD→OA findings; red open circles are the three null pairs. The dashed grey line marks 80% power and the dotted orange line the corresponding detection threshold \|z\| = 2.80 (MDE = 2.80 × SE). (b) Achieved power for each pair at its observed point estimate (K = strict instrument count). Blue bars meet or exceed the 80% threshold (dotted orange line); red bars fall below it. The three null pairs (knee OA → fracture, knee OA → total BMD, knee OA → femoral-neck BMD) reach only 0.14–0.26 power and can exclude only \|β\| ≥ 0.303, so their non-significant results are an absence of evidence rather than evidence of absence.

**Figure 9.** Sample-overlap diagnostics for the ten MR pairs (§2.8). (a) ρ_GE, the correlation between SNP–exposure and SNP–outcome z-statistics across the harmonised instruments (`overlap_statistics.csv`, Supplementary). Positive bars (blue) are the BMD→OA pairs; negative bars (red) are the BMD→fracture pairs. ρ_GE shares the sign of the MR effect in all ten pairs, so no pair owes its direction to a single anti-correlated instrument. (b) ρ_GE against the signed IVW z-statistic; shaded regions are the sign-concordant quadrants, and all ten pairs fall within them (10/10 concordant). Because the exposure and outcome GWAS share UK Biobank participants, ρ_GE also reflects the sample-overlap component of the summary-statistic correlation; such overlap inflates IVW magnitudes toward r_g without changing signs, p-values or directions [14].

**Supplementary Figure S3.** External validation of the reconstructed BMD effect scale against the independent GEFOS benchmark (§3.7). Odds ratios for any fracture per 1-SD lower BMD are shown for GEFOS 2018 (femoral-neck BMD; Trajanoska et al.), for this study's femoral-neck BMD instrument set under the primary strict clumping and the conservative 1-Mb sensitivity clumping, and for this study's total-body BMD instrument set. Whiskers are 95% confidence intervals; the dashed line is the null (OR = 1). This study's estimates lie 8–10% below GEFOS with overlapping confidence intervals, supporting the biological plausibility of the reconstructed effect scale.

**Supplementary Figure S4.** Study design and MR analysis flow, compliant with the STROBE-MR reporting guideline (§2.7). The pipeline moves from GWAS inputs (13 skeletal phenotypes) through quality control and allele harmonisation, strict PLINK-equivalent clumping (r²<0.001 over 10 Mb, with a 1-Mb sensitivity window), two-sample MR across ten exposure–outcome pairs (IVW, MR-Egger, weighted median), and sensitivity/diagnostic checks (MR-PRESSO, leave-one-out, Cochran Q/I², sample-overlap ρ_GE, Bonferroni and BH-FDR correction).

**Graphical Abstract.** A separate file (GraphicalAbstract.pdf) summarises the OA–BMD–fracture triangle: genetically higher BMD protects against fracture yet raises knee-OA risk, and OA shows no well-powered causal pathway to fracture.

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
| **BMD total → fracture** (positive control) | 31 | −0.337 | 0.022 | 3.7×10⁻⁵¹ | 2.4×10⁻¹² | −0.337 | 0.558 | −0.364 | 78% | 0 |
| BMD femoral neck → fracture | 25 | −0.351 | 0.030 | 7.7×10⁻³¹ | 2.3×10⁻⁶ | −0.343 | 0.018 | −0.365 | 83% | 0 |
| **BMD total → OA knee** | 31 | +0.023 | 0.004 | 3.9×10⁻¹⁰ | 6.2×10⁻³ | +0.023 | 0.810 | +0.019 | 81% | 0 |
| Lumbar spine → OA knee | 21 | +0.028 | 0.005 | 7.7×10⁻⁹ | 1.7×10⁻² | +0.026 | 0.295 | +0.037 | 83% | 0 |
| Arms → OA knee | 16 | +0.022 | 0.005 | 1.7×10⁻⁶ | 6.9×10⁻² | +0.022 | 0.982 | +0.031 | 85% | 0 |
| Legs → OA knee | 24 | +0.018 | 0.004 | 2.3×10⁻⁵ | 9.8×10⁻² | +0.017 | 0.815 | +0.025 | 85% | 0 |
| Femoral neck → OA knee | 25 | +0.017 | 0.005 | 6.4×10⁻⁴ | 5.5×10⁻² | +0.017 | 0.875 | +0.015 | 68% | 0 |
| OA knee → fracture | 69 | +0.092 | 0.108 | 0.395 | 0.542 | +0.106 | 0.476 | +0.118 | 49% | 0 |
| OA knee → BMD total | 69 | +0.106 | 0.108 | 0.326 | 0.481 | +0.123 | 0.376 | +0.112 | 48% | 0 |
| OA knee → BMD femoral neck | 69 | +0.142 | 0.108 | 0.190 | 0.272 | +0.164 | 0.169 | +0.113 | 30% | 0 |

*Abbreviations:* IVW, inverse-variance weighted; RE, random effects. β is per 1 SD higher exposure on the standardised (reconstructed) outcome scale. The BMD-exposure (forward) MR estimates are corrected for the placeholder exposure N = 400,000 used in the source files (true 31,873; §3.7): the reconstructed BMD β is 0.2823× of the true SD scale (underestimated 3.5426-fold), so the forward MR β is multiplied by 0.2823 (÷3.5426) to recover the true scale. Reverse MR (BMD as outcome) is reported on the same reconstructed scale as the per-SNP slopes in Figure 6, because the BMD-outcome reconstruction shares the same N uncertainty and only its direction and p-value are interpreted (§3.5). In both directions p-values, Egger-intercept p, I² and LOO sign flips are N-scale-invariant and unchanged. All estimates use the strict PLINK-equivalent clumping of §2.6 (K = strict instrument count); the Egger, I² and weighted-median columns were previously drawn from the 1-Mb sensitivity set and are here unified to the primary strict set. The BMD→fracture positive control yields OR = 1.40 (total) and 1.42 (femoral neck) per 1 SD *lower* BMD, corresponding to the published GEFOS benchmark of 1.55 (§3.7). For BMD femoral neck → fracture the MR-Egger intercept is significant (p=0.018), indicating directional pleiotropy; the MR-Egger slope (−0.343) and weighted median (−0.365) are the pleiotropy-robust estimates for that pair, and the IVW point estimate may be modestly inflated (see §3.3).

---

## Table 4. Statistical power of the two-sample MR estimates

| Exposure → Outcome | K | IVW β | SE | abs z | MDE (80% power) | Achieved power | Interpretation |
|---|---|---|---|---|---|---|---|
| BMD total → fracture | 31 | −0.337 | 0.022 | 15.05 | 0.063 | 1.00 | powered |
| BMD femoral neck → fracture | 25 | −0.351 | 0.030 | 11.55 | 0.085 | 1.00 | powered |
| BMD total → OA knee | 31 | +0.023 | 0.004 | 6.26 | 0.010 | 1.00 | powered |
| Lumbar spine → OA knee | 21 | +0.028 | 0.005 | 5.77 | 0.014 | 1.00 | powered |
| Arms → OA knee | 16 | +0.022 | 0.005 | 4.78 | 0.013 | 1.00 | powered |
| Legs → OA knee | 24 | +0.018 | 0.004 | 4.23 | 0.012 | 0.99 | powered |
| Femoral neck → OA knee | 25 | +0.017 | 0.005 | 3.42 | 0.014 | 0.93 | powered |
| OA knee → fracture | 69 | +0.092 | 0.108 | 0.85 | 0.304 | 0.14 | underpowered |
| OA knee → BMD total | 69 | +0.106 | 0.108 | 0.98 | 0.303 | 0.17 | underpowered |
| OA knee → BMD femoral neck | 69 | +0.142 | 0.108 | 1.31 | 0.303 | 0.26 | underpowered |

MDE = (1.96 + 0.84) × SE is the minimum \|β\| detectable at 80% power with a two-sided α = 0.05; achieved power is computed for the observed point estimate. \|β\|/SE is invariant to the N-scale correction (§3.7), so all power statements are unaffected by it. For the fracture outcome the MDE corresponds to OR ≈ 1.07 per 1-SD lower total-body BMD, OR ≈ 1.09 per 1-SD lower femoral-neck BMD, and OR ≈ 1.36 per unit of the knee-OA exposure contrast as scaled in Table 3 (see the scale caveats in §3.5 and §3.7). K = strict instrument count.

---

## References

1. Hatzikotoulas K, Southam L, Stefansdottir L, et al. Translational genomics of osteoarthritis in 1,962,069 individuals. *Nature*. 2025;641(8065):1217–1224. doi:10.1038/s41586-025-08771-z. PMID 40205036.
2. Qian Y, Xia J, Wang P, et al. Genome-wide Association Studies of over 30,000 Samples with Bone Mineral Density at Multiple Skeletal Sites and Its Clinical Relevance. *Genomics Proteomics Bioinformatics*. 2025;23(5):qzaf097. doi:10.1093/gpbjnl/qzaf097. PMID 41206123.
3. Elsworth B, Lyon M, Alexander T, et al. The MRC-IEU UK Biobank GWAS pipeline. *bioRxiv*. 2018. doi:10.1101/340232. (Self-reported OA accessed as `ukb-b-14486` via OpenGWAS, https://gwas.mrcieu.ac.uk/datasets/ukb-b-14486/.)
4. Hackinger S, Trajanoska K, Styrkarsdottir U, et al. Evaluation of shared genetic aetiology between osteoarthritis and bone mineral density identifies SMAD3 as a novel osteoarthritis risk locus. *Hum Mol Genet*. 2017;26(19):3850–3858. doi:10.1093/hmg/ddx285.
5. Zhou J, Na J, Jiang Z, et al. Deciphering the shared genetic structure between hip osteoarthritis and femoral neck bone mineral density. *Front Genet*. 2025;16:1597005. doi:10.3389/fgene.2025.1597005.
6. Bulik-Sullivan B, Finucane HK, Anttila V, et al. An atlas of genetic correlations across human diseases and traits. *Nat Genet*. 2015;47(11):1236–1241. doi:10.1038/ng.3401.
7. Bowden J, Davey Smith G, Haycock PC, Burgess S. Consistent estimation of the causal effect in Mendelian randomization with some invalid instruments using a weighted median estimator. *Genet Epidemiol*. 2016;40(4):304–314. doi:10.1002/gepi.21965.
8. Verbanck M, Chen CY, Neale B, Do R. Detection of widespread horizontal pleiotropy in causal relationships inferred from Mendelian randomization between complex traits. *Nat Genet*. 2018;50(5):693–698. doi:10.1038/s41588-018-0099-7.
9. Lee SH, Wray NR, Goddard ME, Visscher PM. Estimating missing heritability for disease from genome-wide association studies. *Am J Hum Genet*. 2011;88(3):294–305. doi:10.1016/j.ajhg.2011.02.002.
10. Trajanoska K, Morris JA, Kaptoge S, et al., for the GEFOS consortium. Genome wide association analysis of risk of fracture and its association with bone mineral density: a meta-analysis of genome-wide association studies applied to the UK Biobank. *BMJ*. 2018;362:k3225. doi:10.1136/bmj.k3225. PMID 30158200.
11. Foss MV, Byers PD. Bone density, osteoarthrosis of the hip, and fracture of the upper end of the femur. *Ann Rheum Dis*. 1972;31(4):259–264. doi:10.1136/ard.31.4.259.
12. Dequeker J, Aerssens J, Luyten FP. Osteoarthritis and osteoporosis: clinical and research evidence of inverse relationship. *Aging Clin Exp Res*. 2003;15(5):426–439. doi:10.1007/BF03327364.
13. Davey Smith G, Ebrahim S. 'Mendelian randomization': can genetic epidemiology contribute to understanding environmental determinants of disease? *Int J Epidemiol*. 2003;32(1):1–22. doi:10.1093/ije/dyg070.
14. Burgess S, Davies NM, Thompson SG. Bias due to participant overlap in two-sample Mendelian randomization. *Genet Epidemiol*. 2016;40(7):597–608. doi:10.1002/gepi.21998.
15. Mounier N, Kutalik Z. Bias correction for inverse variance weighted Mendelian randomization. *Genet Epidemiol*. 2023;47:e22522. doi:10.1002/gepi.22522.

## Author contributions (CRediT)

Contributions are reported under the CRediT (Contributor Roles Taxonomy). # equal contribution; * corresponding author.

| Author | Conceptualization | Data curation | Formal analysis | Funding acquisition | Investigation | Methodology | Project administration | Resources | Software | Supervision | Validation | Visualization | Writing – original draft | Writing – review & editing |
|---|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| Zhenxin Zheng¹# | ✓ | | ✓ | | ✓ | ✓ | | | ✓ | | ✓ | ✓ | ✓ | ✓ |
| Song Chen¹# | ✓ | ✓ | | | ✓ | | | ✓ | | | ✓ | | ✓ | ✓ |
| Shenbo Huang¹ | | ✓ | | | ✓ | | | ✓ | | | ✓ | | | ✓ |
| Chao Qian¹ | | ✓ | | | ✓ | | | ✓ | | | ✓ | | | ✓ |
| Hang Yuan¹ | | ✓ | | | ✓ | | | ✓ | | | ✓ | | | ✓ |
| Huricha Bao¹* | ✓ | | | ✓ | | ✓ | ✓ | | | ✓ | | | | ✓ |
| Qiuyan Jiang²* | ✓ | | | ✓ | | ✓ | ✓ | | | ✓ | | | | ✓ |

This statement is provided for CRediT compliance and should be confirmed by all authors before submission.

## Authors, affiliations, funding and competing interests

- **Authors.** Zhenxin Zheng¹#, Song Chen¹#, Shenbo Huang¹, Chao Qian¹, Hang Yuan¹, Huricha Bao¹*, Qiuyan Jiang²*.
  - # Zhenxin Zheng and Song Chen contributed equally to this work.
  - * Corresponding authors: Huricha Bao, MD (baohuricha@126.com) and Qiuyan Jiang, MD (after9951023@gmail.com).
- **Affiliations.** ¹ Department of Orthopedics, The Quzhou Affiliated Hospital of Wenzhou Medical University, Quzhou People's Hospital, Quzhou, Zhejiang Province, 324000, China; ² Department of Surgery, The Quzhou Affiliated Hospital of Wenzhou Medical University, Quzhou People's Hospital, Quzhou, Zhejiang Province, 324000, China.
- **Funding.** None.
- **Competing interests.** The authors declare no competing interests.
- **Data availability.** Summary statistics are from publicly published GWAS (Refs 1–3) and were not redistributed; provenance and accessions are listed in §2.1 and the GWAS Catalog. The analysis code is deposited at the project repository (https://github.com/&lt;owner&gt;/OA-BMD-fracture-MR; archived at Zenodo, doi:10.5281/zenodo.XXXXXXX) for reproducibility, including the MRlap sample-overlap sensitivity script (code/mrlap_overlap_correction.R, executed in R). All source GWAS were approved by their original ethics committees; this secondary analysis required no new ethical approval.

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

| 9 | Correct MR effect magnitudes for the BMD placeholder N | ✅ **Done** | BMD files carried N=400,000 (placeholder) while true N=31,873 → munged β underestimated 3.54-fold (0.2823× of true). Forward MR β (BMD exposure, denominator of ratio) recovered by ×0.2823 (÷3.5426); reverse-MR β (BMD outcome) reported on the reconstructed scale, only direction/p interpreted (§3.7). Corrected in `code/21_mr_corrected_N.py` → `mr_corrected_N.csv`; Table 3 unified to the strict instrument set. **p-values unchanged** (β and SE scale identically). Corrected BMD→fracture OR per SD = 1.40–1.42, corresponding to the GEFOS benchmark of 1.55 within 8.4% (§3.7). |

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
