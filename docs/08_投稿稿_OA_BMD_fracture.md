Genetically higher bone mineral density raises knee osteoarthritis risk yet protects against fracture: genetic correlation and bidirectional Mendelian randomisation across 13 skeletal phenotypes

**Authors.** Zhenxin Zheng¹,#; Song Chen¹,#; Shenbo Huang¹; Chao Qian¹; Hang Yuan¹; Huricha Bao¹,*; Qiuyan Jiang²,*
¹ Department of Orthopedics, The Quzhou Affiliated Hospital of Wenzhou Medical University, Quzhou People's Hospital, Quzhou, Zhejiang 324000, China.
² Department of Surgery, The Quzhou Affiliated Hospital of Wenzhou Medical University, Quzhou People's Hospital, Quzhou, Zhejiang 324000, China.
Equal contribution; * joint corresponding authors.

**Running title.** BMD raises knee OA yet protects against fracture

## Highlights
- Higher genetically predicted BMD raises knee OA risk yet protects against fracture.
- BMD→OA was positive at all five analysed sites; 2 of 5 held under random effects.
- BMD instruments colocalise with fracture (9/31) not knee OA (0/31):endpoint-specific.
- Height bounds size pathway to 9.5–33.2%; 67–91% is neither height- nor BMI-mediated.
- A systemic bone-forming set-point remains a hypothesis consistent with these data.

**Keywords.** osteoarthritis; bone mineral density; fracture; linkage disequilibrium score regression; Mendelian randomisation; genetic correlation

---

## Abstract

**Objective.** Whether higher bone mineral density (BMD) shares aetiology with, or causes, osteoarthritis (OA) is unsettled. We quantified genetic sharing and causal direction across 13 skeletal phenotypes and tested whether one BMD-raising set-point has opposing fracture and joint consequences.

**Design.** Linkage-disequilibrium score regression (LDSC) and two-sample Mendelian randomisation (MR) were applied to 13 GWAS summary-statistic datasets (two OA definitions, ten regional BMD sites, any-fracture). LDSC covered all 78 phenotype pairs; MR ran in eight phenotypes with clumped instruments — five of ten prespecified BMD sites, both OA definitions, any-fracture — with BMD→fracture as positive control, multivariable MR conditioning on height and, separately, on BMI, MRlap correction, and reverse MR.

**Results.** BMD sites inter-correlated strongly (r_g 0.69–0.96) and negatively with fracture (−0.40 to −0.54); both OA definitions correlated positively with every BMD site (+0.12 to +0.25). Genetically higher BMD raised knee-OA risk at all five analysed sites — OR 1.060–1.099 per SD (95% CI 1.025–1.135) — and protected against fracture (OR 1.40, 95% CI 1.34–1.46 per SD lower BMD; p=3.7×10⁻⁵¹). Heterogeneity was high (I² 68–86%) but matched the positive control (78–83%), so the causal reading rests on pleiotropy-robust estimators rather than fixed-effect p alone; under random effects 2 of 5 remained significant. Conditioning on height attenuated the OA effect by 9.5–33.2% (a lower bound for size); conditioning separately on BMI left it positive and significant at all five sites (−3.7% to 46.0% attenuation), so the 67–91% height-independent residual is neither height- nor BMI-mediated — not identified as bone mass. Reverse estimates were non-significant and underpowered (power 0.14–0.26). MRlap moved all ten estimates further from the null (median 20.5%), reversing none.

**Conclusions.** One genetically indexed bone phenotype is protective at the whole-bone endpoint and harmful at the joint endpoint, 67–91% of the joint effect independent of height and not BMI-mediated. The joint effect is modest per SD (OR 1.06–1.10 areal BMD) yet attaches to a ubiquitous, modifiable exposure, measured at five sites only. These estimates do not support changing prescribing practice; excluding a clinically meaningful joint cost requires trials powered for incident OA recorded at the joint of interest alongside fracture.

---

## 1. Introduction

Osteoarthritis (OA) is a prevalent disorder of ageing [23]; areal bone mineral density (BMD), the clinical measure of skeletal strength, is measured in the same populations. Their epidemiological association has been debated since Foss and Byers reported raised femoral-head BMD in osteoarthritic hips in 1972 [11]. Observational studies, reviewed in [33], link higher areal BMD to increased OA prevalence and incidence, although osteoporosis and OA were long regarded as opposing conditions [12]. The direction matters clinically: if high BMD predisposes to OA, joint endpoints deserve monitoring in bone-anabolic trials, though whether such a cost is mediated by bone mass or skeletal size is itself unresolved.

Genome-wide data can separate these explanations. LDSC estimates genetic correlation (r_g) from summary statistics, tolerates sample overlap, and is unbiased by stratification when the intercept is near unity [6]. Two-sample MR uses genetic variants as instruments to interrogate causality [13]. Applied to OA and BMD, both have given partial answers [16,17]: Hackinger et al. reported r_g=+0.18 between hip/knee OA and lumbar-spine BMD (p=0.022) [4]; Zhou et al. reported r_g=+0.13 by LDSC and +0.17 by high-definition likelihood between hip OA and femoral-neck BMD, supporting a negative femoral-neck BMD→hip OA effect [5]; and a 2025 two-sample MR of BMD against OA progression reports the same site-by-joint dependence, femoral-neck BMD raising hip-replacement risk (OR 1.27, 95% CI 1.13–1.43) but not knee-replacement risk [31]. Our BMD→knee-OA analysis, by contrast, finds a positive effect contrasting with Zhou's negative femoral-neck estimate (discussed below).

The most direct predecessor is Hartley et al. [22], who used multivariable MR to estimate the causal effect of BMD on OA independently of BMI and reported a positive effect at both hip (OR 1.28, 1.05–1.57) and knee OA (OR 1.40, 1.20–1.63). Our study extends that work in four ways: (i) ten regional BMD sites rather than one composite, so site specificity is tested directly (five of ten entered MR; §2.3); (ii) conditioning on adult height and separately on BMI, separating size from adiposity; (iii) an internal BMD→fracture positive control, absent from Hartley et al.; and (iv) a full 13×13 LDSC matrix rather than pairwise estimates.

Three questions remain open. First, most analyses included at most two BMD sites, so site specificity is unknown despite OA being joint-localised and BMD site-regulated. Second, OA has usually been represented by a single definition, leaving their bone relationship uncertain. Third, the OA–BMD–fracture triangle has not been analysed jointly; modelling all three tests the pipeline and whether the OA–BMD association carries fracture relevance.

We addressed all three in 13 skeletal phenotypes using LDSC across all 78 phenotype pairs, two-sample MR with an internal positive control, and reverse MR.

## 2. Methods

This study followed the STROBE-MR reporting guideline [Skrivankova et al., BMJ 2021]. Analyses were performed in Python 3.13 with numpy, pandas, scipy and matplotlib. The sample-overlap correction used R 4.6.1 with MRlap v0.0.3.3. All code is in the project repository (Data availability).

**Data sources.** Thirteen GWAS summary-statistic datasets were analysed: radiographic knee OA (OA_knee; GO Consortium 2025, N≈1.21M) [1], self-reported OA (OA_selfrep; IEU OpenGWAS `ukb-b-14486`, N=462,933; 38,472 cases) [3], ten regional BMD sites (Qian et al. 2025, N≈1.873 each) [2], and any-fracture (GEFOS 2018, N=264,973; 37,857 cases) [10]; full accessions are in the Data availability statement. All four contribute sex-combined primary estimates, so every estimate is sex-averaged; a sex-stratified sensitivity analysis (F8) confirmed BMD→OA positivity in both sexes (§4). All LDSC used 1000 Genomes Phase 3 European LD scores; ancestry composition is in Table S14. The European-only reference against the 87.31%-European knee-OA GWAS means the near-1.0 intercept does not exclude population-stratification bias (primary limitation, §4).

**LDSC.** Genetic correlation (r_g) was estimated by LDSC [6] in Python (Bulik-Sullivan) on 690,150 HapMap3 SNPs with 1000 Genomes Phase 3 European LD scores, using delete-one-block jackknife (100 blocks) for SEs. Palindromic SNPs were excluded; alleles were harmonised to one reference so every matrix cell uses the same variants. Both Bonferroni (α=0.05/78) and Benjamini–Hochberg FDR were applied across 78 phenotype pairs.

**MR.** Instruments were selected at p<5×10⁻⁸ and clumped at r²<0.001 over 10 Mb against 1000 Genomes Phase 3 EUR (strict PLINK-equivalent clumping benchmarked against PLINK; 20–35 retained per BMD site, 116 for knee OA). The BMD exposure set followed a prespecified rule: of the 11 DXA-derived BMD sites in Qian et al. 2025, the cranial head site — the only non-load-bearing region, whose instruments carry non-skeletal pleiotropy (head BMD is genetically correlated with intracranial aneurysm [2]) — was excluded a priori, leaving ten prespecified canonical DXA regions; a site enters MR only if strict clumping completes. Clumping completed for five sites at submission (arms, femoral neck, lumbar spine, legs, total body), each ≥15 independent instruments; the remaining five showed comparable genome-wide signal (229–345 loci at p<5×10⁻⁸) but strict clumping was not completed (1000G EUR panel unavailable) and they were deferred rather than analysed with a non-prespecified set (Table S12). Instrument shortage drove no exclusion. Harmonisation against the outcome GWAS reduced the analysed sets to 16–31 instruments per forward pair and 69 per knee-OA-exposure pair; 47 of 116 knee-OA instruments (41%) were lost (Fig. S9), an attrition we report explicitly. The 69 retained instruments were identical across all three knee-OA outcomes; the dropped instruments were not effect-size selected (median −log10P 9.10 vs 9.89 retained; 43% stronger signal). Random-effects IVW used the multiplicative convention SE_random = SE_fixed × √[max(1, Q/df)] (TwoSampleMR), not additive DerSimonian–Laird τ²; where the convention changes an inference, both are reported (§3, Table 2 note).

**Sample overlap.** All pairs share UK Biobank participants, partially violating the two-sample assumption. Under the Burgess approximation [14], the IVW bias equals the overlap fraction × r_g, bounded above by 2.6% for BMD→OA and 12% for BMD→fracture, and acts in the direction of r_g; because r_g shares the sign of the MR effect in all ten pairs (as does the SNP-effect correlation ρ_GE, Fig. 9), this component acts away from the null. The two other IVW biases — weak-instrument (NOME) attenuation and winner's curse from selecting genome-wide-significant instruments — act towards the null. We report MR-Egger, weighted median and MR-PRESSO alongside IVW, and applied MRlap [15] (corrects weak-instrument bias and winner's curse while accounting for overlap) to all ten pairs. The joint correction moved every estimate further from the null (median +20.5%; observed-versus-corrected difference p<0.05 for 9/10 pairs), showing towards-null components dominate the away-from-null overlap component, reversing no direction or significance. The primary IVW estimates are conservative for the combined bias, not overlap alone; we report MRlap alongside IVW.

**Effect-scale reconstruction.** Two scale corrections were applied, each verified against the deposited data rather than assumed. (i) *BMD exposure.* The deposited files carry a placeholder N=400,000, but for every instrument the closed-form expectation SE = 1/√[2·N·f·(1−f)] matches N=31,873 to within 0.1%, so BMD effects enter on the verified per-SD scale. (ii) *Knee-OA outcome.* That file supplies only per-SNP N, Z and MAF, and reconstructing β=Z·SE with SE=1/√[2·N·f·(1−f)] returns the **linear-probability** effect of a binary outcome, not its log-odds; the log-odds scale requires the effective size N_eff=N·K(1−K). We did not assume K: of SNPs carrying both a deposited log-odds effect and the reconstructed Z, 39 matched, with |Z_deposited/Z_reconstructed| = 1.000000 (SD 5×10⁻¹⁴), confirming one association, and implied K(1−K)=N_eff/N median 0.0871 (IQR 0.0835–0.0881), i.e. K=0.0964; we applied 1/√[K(1−K)] = 3.389. The fracture outcome needed no correction because its file supplies log-odds β and SE directly — which is how the two binary outcomes came to be handled inconsistently. Both corrections scale β and SE by one constant, so z, p, direction, significance, Cochran Q, I², MR-Egger slope and intercept tests, weighted-median sign, PRESSO outlier calls and power are unchanged; only absolute effect sizes move — away from the null for BMD→OA (OR 1.060–1.099, not 1.017–1.028) and towards it for the reverse direction.

**Colocalisation.** Bayesian colocalisation (coloc.abf; Giambartolomei et al. 2014) tested the three a priori loci (WNT16, SOX6, GDF5) and all 31 BMD instrument loci (±500 kb; PP.H4≥0.8), with a BMD–fracture control.

**MVMR.** Multivariable MR conditioned the BMD→OA effect on adult height (Yengo et al. 2022, N=5.0M) [19] as a skeletal-size proxy, using the MVMR-IVW estimator [20] and MVMR-Egger [21], judged by the Sanderson–Smith conditional F (Fc≥10). The primary configuration requested up to n_h=100 height instruments (49–50 survived harmonisation; Yengo et al. 2022); a grid (n_h ∈ {0,100,300,1,000}) and both EUR and EUR_noUKB height sources are in the Supplementary Material. Height was the primary size proxy: BMI is a major knee-OA risk factor and part of the adiposity pathway through which a BMD-raising process could raise knee-OA risk, whereas height indexes skeletal size without adiposity. A second MVMR therefore conditioned the same effect on BMI (GIANT–UK Biobank meta-analysis, Pulit et al. 2019, N≈806,800) [38] with the same estimators, n_bmi=100 as the primary configuration (61–63 instruments survived harmonisation) and a grid (n_bmi ∈ {100,300,800}); conditional F statistics and all estimates are in Supplementary Table S16. Height and BMI were conditioned in separate analyses, so each proxy is adjusted one at a time rather than jointly.

## 3. Results

SNP-based heritability ranged 0.26–0.39 across the ten BMD sites and was 0.035 (knee OA), 0.022 (self-reported OA) and 0.026 (fracture); 12/13 LDSC intercepts fell between 1.009 and 1.069 (fracture 0.96). Because LDSC uses a European-only reference while the knee-OA GWAS is multi-ancestry (Table S14), the near-1.0 intercept does not exclude population-stratification or ancestry-mismatch bias (§4). The 13×13 r_g matrix showed three blocks (Figure 1, Table 2). The ten BMD sites inter-correlated strongly (r_g 0.69–0.96, mean 0.81). Both OA definitions correlated positively with all ten BMD sites (knee OA +0.14 to +0.19; self-reported OA +0.12 to +0.25; 20/20 significant) and with each other (r_g=+0.84). The two definitions also agreed at the variant level: across 1,980 genome-wide-significant knee-OA loci, 94.4% had concordant signed Z with self-reported OA, confirming they relate to bone in the same genetic direction. Any-fracture correlated negatively with every BMD site (r_g −0.40 to −0.54; all q_FDR<0.001). It correlated positively with knee OA (r_g=+0.119, q=3.3×10⁻⁹) but non-significantly with self-reported OA (+0.056, p=0.051); 77/78 pairs were significant at both BH-FDR and Bonferroni, the sole exception being self-reported-OA × fracture.

The positive control behaved as expected: two-sample MR using 31 independent total-body BMD instruments showed genetically predicted higher BMD reduces fracture risk (IVW β=−0.337 per SD, SE 0.022, p=3.7×10⁻⁵¹; OR 1.40, 95% CI 1.34–1.46, per SD lower BMD). The estimate held under random effects (p=2.4×10⁻¹²), the MR-Egger intercept was null (p=0.558), weighted median agreed (−0.364), and leave-one-out produced no sign reversal. The femoral-neck set (25 SNPs) gave a concordant result (β=−0.351; OR 1.42, 1.34–1.51), but its MR-Egger intercept was not null (p=0.018), a directional-pleiotropy signal we report rather than absorb; its pleiotropy-robust slope (−0.343) and weighted median (−0.365) agree with IVW, so the control calibrates the pipeline despite its own non-clean diagnostic. The conservative 1-Mb set gave β=−0.366 (OR 1.44); these agree with the published GEFOS estimate (OR 1.55, 1.48–1.63) [10].

MR supported a positive effect of BMD on knee OA across all five exposures (Table 2): total-body β=+0.078 (p=3.9×10⁻¹⁰), lumbar-spine +0.095 (p=7.7×10⁻⁹), arms +0.075 (p=1.7×10⁻⁶), legs +0.060 (p=2.3×10⁻⁵), femoral-neck +0.058 (p=6.4×10⁻⁴) — OR 1.060–1.099 per SD higher BMD, against OR 1.40 (1.20–1.63) per SD in the largest one-sample estimate [22], a 3.6–5.8-fold log-odds difference (adjudicated §4). MR-Egger intercepts were all p>0.24 and weighted medians matched the IVW sign at all five sites, with bootstrap 95% CIs excluding the null at three of five; with 16–31 instruments the intercept test has limited power, so a non-significant intercept excludes detectable unbalanced pleiotropy rather than pleiotropy itself. Leave-one-out produced no sign reversals and directions agreed with r_g (+0.12 to +0.25). Heterogeneity was high (I²=68–86%); under random-effects IVW only two of five exposures remained significant (total-body p=0.006, lumbar-spine p=0.017; a second random-effects convention additionally retains arms, p=0.033), so the count depends on the random-effects convention and we base the causal reading on estimator concordance rather than fixed-effect p-values. All five effects remained positively causal after MR-PRESSO outlier removal, with four of five strengthening. The positive control was comparably heterogeneous (I² 78% and 83%), indicating high I² is a property of the BMD instrument sets rather than an OA-specific IVW failure; a multi-estimator robustness overview of all ten pairs is in Supplementary Fig. S11.

Knee OA as exposure gave non-significant effects on total-body BMD (β=+0.031, p=0.33), femoral-neck BMD (β=+0.042, p=0.19) and fracture (β=+0.027, p=0.40); these were weakly powered (0.14–0.26; Table 3) and excluded only |β|≥0.09 in SD of BMD per unit OA liability. Removing the MR-PRESSO outliers flagged for these pairs moved two point estimates substantially (total-body BMD to +0.051, fracture to +0.059) without altering their non-significance, so the nulls reflect both low power and sensitivity to instrument selection.

Instrument strength and directionality. Every forward BMD instrument exceeded the conventional F>10 threshold (minimum per-SNP F 29.8–31.1; aggregate F 44–74; R²_exposure 0.034–0.063; Tables 3 and S13), so weak-instrument bias does not explain the forward effects. The Steiger directionality test supported BMD as the causally prior trait for every forward pair (z = 3.05–6.11, one-sided p = 4.8×10⁻⁴ to 5.0×10⁻¹⁰) and was non-significant for the three reverse pairs (|z|<0.1, p≈0.48–0.49). That non-result does not discriminate direction: with 69 instruments and a null reverse estimate the test lacks power to separate the models and is not read as evidence BMD rather than OA drives the association.

The MR-PRESSO distortion test was significant for 9/10 pairs (expected from the widespread horizontal pleiotropy common in complex-trait MR [8]) and identified 0–3 outliers per pair; removing outliers did not reverse any direction, and strict-clumped estimates reproduced the conservative 1-Mb clumping. All significant findings survived Bonferroni correction (α=0.005). The IVW estimate is therefore not read as a single clean causal parameter; the BMD→OA direction is taken from the concordance of pleiotropy-robust estimators (Table 2, Fig. S11) rather than the IVW point estimate alone.

On the verified per-SD scale, the BMD→fracture OR (1.42 femoral-neck) lay 8.4% below GEFOS (1.55) on the odds-ratio scale and 20.0% below on the log-odds scale, with overlapping CIs (we compare on the log scale). The total-body OR 1.40 was 9.6% below on the ratio and 23.1% on the log-odds scale (cross-site approximation). The MRlap correction moved all ten estimates further from the null (median +20.5%; difference p<0.05 for 9/10), the direction expected when the towards-null weak-instrument and winner's-curse components dominate the away-from-null overlap component, and it reversed no direction or significance conclusion; the cross-trait intercept of the correction was |λ| ≤ 0.037 for all ten pairs.

At WNT16, SOX6 and GDF5, loci nominated as shared effectors, colocalisation with knee OA was not supported (PP.H4 ≤0.005 at each; e.g. WNT16 knee-OA z=0.94, p=0.35, PP.H4=0.002; GDF5 knee-OA z=−13.7 but no BMD genome-wide signal within ±500 kb, PP.H4 0.003–0.017). Across all 31 BMD instrument loci, colocalisation was endpoint-split: with knee OA none reached PP.H4≥0.8 (median 0.005, 0/31), whereas with fracture 9/31 did (median 0.153; paired Wilcoxon one-sided p=0.0043, exact two-sided p=0.0086) — the same instrument set shares causal variants with bone but not joint; six loci favoured distinct variants (PP.H3≥0.8). One locus fell marginally short: rs11228269 (chr11:68.29 Mb, LRP5; BMD z=−8.65, knee-OA z=−4.06, concordant) gave PP.H4=0.799, 0.0006 below threshold, with residual mass on a BMD-only signal (PP.H1=0.156); we retain the pre-specified 0.8 cut-off and do not count it as colocalisation. The pipeline detected colocalisation where expected: at WNT16, BMD versus fracture gave PP.H4=0.999–1.000. The knee-OA null reflects power — under the MR estimate the expected knee-OA z at the strongest instrument is 0.85, and reaching genome-wide significance would require ~5×10⁷ participants (~42× the present sample).

How much of the BMD→OA effect is skeletal size? Conditioning the BMD→OA effect on height attenuated it by 9.5–33.2% across the two tabulated univariable baselines (13.5–33.2% against the primary Table 2 baselines and 9.5–21.2% against the MVMR-set baselines; primary n_h=100/EUR; Supplementary Tables S9 and S10 and Supplementary Fig. S10), leaving 67–91% independent of height. Because height indexes stature but not bone width or geometry [35], this residual is not identified as bone mass; the mass-versus-geometry split is resolved only by test F2. The positive control confirmed that the established BMD→fracture effect survives the same conditioning: BMD→fracture moved only −1.7% to −3.0% across the same two baselines (−2.8% against the MVMR-set baseline), the control being insensitive to baseline choice. Because conditioning on height captured a genuine, direct OA-raising height effect (θ_H=+0.145 to +0.165, p≈10⁻¹⁶ to 10⁻¹⁹), the attenuation reflects size confounding, not dilution by a null covariate. A parallel BMI-conditioned MVMR (Table S16) left θ_BMD positive and fixed-effect-significant at all five sites (total-body +0.069, p=2.6×10⁻⁸), attenuating −3.7% to 46.0% (the estimate rose at the femoral-neck site) — below the height attenuation at four of five sites; only the non-weight-bearing arms site exceeded it (46.0% vs 13.5%), losing significance under random effects. The MVMR-Egger intercept was nominally significant only at the lumbar site (p=0.027). The positive control moved 0.9%, and θ_BMI→fracture was null (p=0.61) while θ_BMI→knee-OA was strong (β≈+0.69): BMI raises joint- but not bone-endpoint risk, so the height-independent residual is not BMI-mediated. Under BMI conditioning the control's MVMR-Egger intercept was nominally significant at the femoral-neck site (p=2.1×10⁻⁴; total-body p=0.12), a pleiotropy signal reported rather than absorbed.

## 4. Discussion

Jointly analysing 13 skeletal phenotypes gave a consistent architecture: BMD sites form a strongly inter-correlated cluster, OA correlates positively with every BMD site, fracture negatively with every BMD site, and the two OA definitions correlate strongly. On pleiotropy-robust evidence, genetically higher BMD raises knee-OA risk; the reverse (OA→BMD) and OA→fracture estimates were non-significant and underpowered.

The opposing-endpoint pattern is an old clinical observation [11,12]; what is new is its demonstration for a single genetic instrument set within one calibrated pipeline. One non-exclusive explanation is partial mediation by skeletal size rather than bone mass [18,32]. Conditioning on height attenuated the effect by 9.5–33.2% (the fracture control moved only −1.7% to −3.0%) and height carried a direct OA-raising effect consistent with established epidemiology [24,25,26]; BMI conditioning likewise left the effect significant at all five sites (Table S16). The fracture control is height-transparent, and in the largest skeletal-geometry GWAS hip shape was causally implicated in fracture but not hip OA [35].

Zhou et al. [5] reported higher femoral-neck BMD protected against hip OA, whereas our ten-site analysis shows higher BMD raises knee-OA risk — a hip-versus-knee contrast cautioning against generalising any single BMD→OA estimate. Hartley et al. [22] reported positive BMD effects on both hip (OR 1.28, 1.05–1.57) and knee OA (OR 1.40, 1.20–1.63); Funck-Brentano et al. [36] likewise reported a positive overall effect (all-OA OR 1.14, 1.06–1.22); quantitative comparison is in Fig. S12. These contradict the directional form of P6 (motivated by Zhou's protective finding; P6 and F1 in Table S11). We adjudicate rather than cite selectively: the hip evidence disagrees with itself (protective [5], harmful [22,36]), and decisively for F1 the criterion is conjunctive — the hip effect must be positive AND of knee-like magnitude. Within the only same-pipeline comparison [22], hip and knee intervals overlap almost completely (1.05–1.57 vs 1.20–1.63) and on the log-odds scale the hip estimate lies 27% below the knee's (0.247 vs 0.336); Hartley establishes positivity, not equality. We therefore restate P6 as a magnitude claim (a hip effect positive and at least as large as the knee's would fail P6), deviating from the pre-registered directional form prompted by contradictory hip evidence, scoring only prospective hip estimates.

The directional agreement with Hartley et al. [22] is matched by a magnitude disagreement: their knee estimate is OR 1.40 (1.20–1.63) per SD higher eBMD against 0.058–0.094 here (OR 1.060–1.099 per SD higher BMD). Three features make this a difference of scale and design rather than biology. First, the exposures differ: Hartley instrumented heel-ultrasound eBMD (bone size) whereas we instrumented DXA areal BMD (area-normalised). Second, the designs differ: Hartley's is one-sample 2SLS (overlap 100%) whereas ours is bounded at ≤2.6% and MRlap moved every estimate a median 20.5% further from the null. Third, our pipeline is calibrated: the BMD→fracture control returns OR 1.40–1.42 against GEFOS 1.55. We read the directions as concordant and magnitudes unresolved, recording the disagreement as a limitation.

OA is a disorder of the osteochondral unit [27–30], in which subchondral bone remodelling feeds cartilage degradation. A systemic, polygenic bone-forming set-point readable as areal BMD would protect at the whole-bone endpoint and raise OA risk at the joint — the Systemic Osteogenic Set-point (SOS [34]), operationalised as areal BMD as a heuristic, not a fitted model. If the causal quantity is volumetric density, an areal-BMD estimate is diluted by skeletal size, so volumetric-BMD (or area-conditioned) instruments should yield a BMD→OA effect at least as large as the areal one — the quantitative form of F2, which areal data cannot test. Height captures neither bone width nor geometry [35], so the height-independent 67–91% residual is not established as bone mass (resolved only by F2). The reading generates six predictions (P1–P6) and ten falsifying tests (F1–F10), set out below and in Table S11.

**Pre-registered SOS predictions and falsifying tests.** *Predictions (derived).* **P1** opposing endpoints — supported (positive-control OR 1.40 per SD lower BMD; BMD→knee-OA 5/5 positive). **P2** systemic not local-mechanical — systemic clause supported (non-weight-bearing arms β = +0.075, p = 1.7×10⁻⁶), ordering clause under-powered (arms − femoral neck +0.017, one-sided p = 0.228). **P3** diffuse polygenic — 0/31 loci colocalise (PP.H4 ≥ 0.8); direction-consistent, non-discriminating. **P4** size a separable partial component — height attenuation 9.5–33.2% and BMI attenuation −3.7% to 46.0% (fracture controls −1.7% to −3.0% and 0.9%); 67–91% residual unresolved into mass versus geometry (F2). **P5** constitutional asymmetry — reverse null, consistent but not evidence. **P6** site-pathway dominance — weakened to a magnitude claim (hip evidence contradictory). *Falsifying tests.* **F1** hip-OA replication; **F2** volumetric/size decomposition; **F3** bone-versus-cartilage tissue enrichment (LDSC-SEG); **F4** eQTL colocalisation (cross-trait arm executed: fracture 9/31 versus knee-OA 0/31 loci at PP.H4≥0.8, paired Wilcoxon p=0.0043, §3; eQTLGen whole-blood arm executed, p=0.001156, but blood is neither bone nor cartilage, so the tissue-literal arm remains unexecuted); **F5** osteophyte-stratified outcome — F1–F3 and F5 prospective, requiring external data; F5 pre-registered here as the immediate next step and not performed. **F6** cross-site homogeneity (executed: Q = 6.92, 4 df, p = 0.14, τ² = 0) and **F7** joint-locus exclusion (executed: GDF5 not a BMD instrument at any site; excluding WNT16 or SOX6 changes no estimate's direction or significance) both passed. **F8** sex stratification (executed: positive in both sexes, §4). **F9** site-gradient ordering (arms−femoral neck +0.017, p=0.228). **F10** ancestry sensitivity (not executed: no public ancestry-stratified files). P4 is the discriminating prediction, resolved only by F2.

A direct MR test of OA→any-fracture was null (β=+0.027, p=0.40) despite a positive genetic correlation with fracture (r_g=+0.119, p=3.1×10⁻⁹). The two are reconcilable: the BMD-mediated component of that correlation runs oppositely (higher BMD raises knee-OA risk yet protects against fracture), so the positive r_g must arise from pathways outside BMD; genetic correlation is not causation. This is among the first direct MR tests of OA→any-fracture and excludes large but not moderate effects (power 0.14). Our BMD→OA finding is consistent with Alhassan et al. [16], who identified BMD among exposures with MR support for causal OA; the ten-site design extends that with site specificity and a quantified size decomposition.

**Strengths and limitations.** Strengths: ten BMD sites with two OA definitions, an internal BMD–fracture positive control calibrated against GEFOS (OR 1.40–1.42 vs 1.55), dual FDR/Bonferroni correction and directional concordance across IVW, MR-Egger, weighted median and MR-PRESSO. Limitations: the European-only BMD instruments and LDSC reference against an 87.31%-European knee-OA outcome mean forward estimates generalise most directly to European-ancestry populations (the near-1.0 intercept does not exclude ancestry-mismatch bias; Table S14). Five of ten prespecified BMD sites entered MR because strict clumping completed for them; the five deferred sites showed comparable genome-wide signal, so the analysed subset reflects compute availability rather than biology (Table S12). Areal BMD conflates size with mass [18,32,37]; height-conditioning bounds the size share at 9.5–33.2% (a lower bound), leaving the 67–91% height-independent residual unresolved into bone mass versus geometry (test F2); it is not BMI-mediated (Table S16), with height and BMI adjusted one at a time rather than jointly [35]. The primary MR is sex-averaged; the F8 sex-stratified MR (sex-specific outcome GWAS; ten-site random-effects meta on the verified per-SD scale) was positive in both sexes (KneeOA Female β=+0.090, Male β=+0.127; AllOA Female β=+0.097, Male β=+0.066; F-versus-M p=0.029 and 0.006), with no direction reversal; sex difference ran oppositely for the two definitions and the site-level interaction was null (Q p=0.55 and 0.27), so F8 is direction consistency rather than a resolved magnitude difference. Heterogeneity is high across all pairs including the positive control, so the causal claim rests on estimator concordance rather than fixed-effect p-values; the reverse and OA→fracture estimates are null but weakly powered (0.14–0.26), and residual overlap is bounded at ≤2.6% for BMD→OA but cannot be removed with summary data. 

**Conclusion.** Genetically higher BMD protects against fracture and, on the pleiotropy-robust evidence, raises knee-OA risk. The effect concerns the clinically measured areal-BMD phenotype at the five analysed sites; 67–91% of it is independent of height and not BMI-mediated, but not resolved into bone mass versus skeletal geometry; and reverse causation and an OA→fracture pathway were not supported in MR (the reverse estimates exclude only |β|≥0.09; the positive knee-OA × fracture r_g is genetic correlation, not causation). The joint-endpoint effect is modest in individual terms (OR 1.060–1.099 per SD higher BMD, 95% CI 1.025–1.135 across sites) but attaches to a ubiquitous, modifiable exposure; these estimates do not support changing prescribing practice, and excluding a clinically meaningful joint cost requires trials or pharmacogenetic designs powered for incident OA at the joint of interest alongside fracture.

---

## Figure legends

**Figure 1.** Heatmap of LDSC genetic correlations (r_g) across 13 skeletal phenotypes. Cell text shows r_g; asterisks denote BH-FDR significance. Black lines separate the three phenotype groups. The colour scale is diverging and centred on zero.


**Figure 2.** Forest plot of MR estimates on the N-corrected effect scale. Upper panel: BMD→fracture positive control (OR ≈1.40–1.42 vs GEFOS 1.55). Middle panel: five BMD exposures against knee OA (all positively causal). Lower panel: reverse/null pairs (all non-significant, p=0.19–0.40). Points are IVW fixed-effects estimates with 95% CI; blue = protective, red = risk-increasing, grey = null.


**Figure 3.** Conceptual causal architecture of the OA–BMD–fracture triangle. Solid arrows denote MR causal effects (BMD protects against fracture; BMD raises knee OA under pleiotropy-robust estimators); dashed grey arrows denote non-causal relationships.


**Figure 4.** Per-SNP scatter plots for four headline pairs (IVW, MR-Egger, weighted-median fits). (a,b) BMD→fracture (negative slopes, positive control). (c) Total-body BMD→knee OA (positive). (d) Knee OA→total-body BMD (flat/null). Red open circles mark MR-PRESSO outliers.


**Figure 5.** MRlap correction of the ten MR pairs (accounting for participant overlap while correcting weak-instrument bias and winner's curse); every pair moves away from the null (median +20.5%), showing that the towards-null components dominate the away-from-null overlap component, and the observed-versus-corrected difference is significant for 9/10 pairs (c).


**Supplementary Figure S1.** Per-SNP scatter plots for the six remaining MR pairs. **S2.** Leave-one-out stability of all ten pairs. **S3.** External validation of the reconstructed BMD scale against GEFOS. **S4.** Study design and MR flow (STROBE-MR). **S5.** Sensitivity of OA_knee h² to representative N (r_g is N-invariant). **S6.** Fracture sample-size decision. **S7.** Pleiotropy and heterogeneity diagnostics. **S8.** Colocalisation of BMD with knee OA and its power. **S9.** Instrument selection and attrition flow (STROBE-MR), tallied from the deposited clump and harmonised instrument tables. **S10.** Multivariable MR decomposition of the BMD→OA effect conditional on height: (a) univariable vs MVMR estimates; (b) attenuation percentages with the BMD→fracture positive control; the height-independent residual is not resolved into bone mass versus skeletal geometry. **S11.** Multi-estimator robustness grid for the ten MR pairs (IVW fixed/random, MR-Egger, weighted median, PRESSO outlier-removed, MRlap-corrected, MVMR). **S12.** Quantitative comparison with published MR estimates (Hartley 2022; Funck-Brentano 2019; Zhou 2025) on their native scales. **Supplementary Figure S13.** Forest plot of the selected genetic correlations: any-fracture versus ten BMD sites (all negative); knee OA and self-reported OA versus ten BMD sites (all positive); the two OA definitions with each other (r_g=+0.84), knee OA with fracture (r_g=+0.119, q=3.3×10⁻⁹) and self-reported OA with fracture (non-significant, +0.056). Error bars are 95% confidence intervals from block jackknife.
 **Supplementary Figure S14.** SNP-based heritability (observed scale) of the 13 skeletal phenotypes; diamond markers show liability-scale h² for the three binary traits.
 **Supplementary Figure S15.** MR-PRESSO funnel plots for the four headline pairs; red open circles are outlier instruments.
 **Supplementary Figure S16.** Statistical power of the two-sample MR design. Filled blue circles are the seven well-powered pairs; red open circles the three null pairs (power 0.14–0.26, excluding only |β|≥0.09).
 **Supplementary Figure S17.** Sample-overlap diagnostics (ρ_GE) for the ten MR pairs; ρ_GE shares the sign of the MR effect in all ten pairs.
 **Graphical Abstract.** The OA–BMD–fracture triangle: higher BMD protects against fracture yet raises knee-OA risk.


## Table 1. SNP-based heritabilities and LDSC intercepts

| Phenotype         | N used    | h² (observed) | SE     | Intercept |
| ----------------- | --------- | ------------- | ------ | --------- |
| OA_knee†          | 1,208,397 | 0.0354        | 0.0007 | 1.069     |
| OA_selfrep        | 462,933§  | 0.0223        | 0.0016 | 1.035     |
| fracture          | 264,973§§ | 0.0258        | 0.0014 | 0.960     |
| BMD total         | 32,000‡   | 0.391         | 0.012  | 1.025     |
| BMD femoral neck  | 32,000‡   | 0.257         | 0.009  | 1.027     |
| BMD femoral total | 32,000‡   | 0.317         | 0.011  | 1.031     |
| BMD lumbar spine  | 32,000‡   | 0.350         | 0.013  | 1.014     |
| BMD spine         | 32,000‡   | 0.357         | 0.014  | 1.016     |
| BMD arms          | 32,000‡   | 0.325         | 0.012  | 1.009     |
| BMD legs          | 32,000‡   | 0.325         | 0.012  | 1.025     |
| BMD pelvis         | 32,000‡   | 0.358         | 0.014  | 1.027     |
| BMD rib           | 32,000‡   | 0.303         | 0.011  | 1.012     |
| BMD trunk         | 32,000‡   | 0.366         | 0.013  | 1.020     |

§§ Fracture N adopted as 264,973 (37,857 cases / 227,116 controls; GEFOS discovery, Trajanoska et al. 2018), identified by the measured effective sample size N·K·(1−K)=31,605, which matches the GEFOS theoretical value of 32,449 within 2.6%. § OA_selfrep N confirmed at 462,933 (38,472 cases / 424,461 controls; IEU OpenGWAS `ukb-b-14486`). ‡ BMD h² was computed at a rounded representative N of 32,000 (the per-site GWAS N is 31,873); because h² ∝ 1/N the rounding makes the reported h² 0.40% lower than at the exact N, an effect far smaller than the reported SEs. † OA_knee h² is reported at the representative proband-subset N=1,208,397; h² is N-dependent, whereas r_g and its SE are N-invariant because r_g is a ratio of covariances estimated on the same SNP set. The LDSC regression uses the per-SNP N column; all r_g values are N-scale-invariant. ‡ BMD N verified per site (Qian et al. 2025; range 30,449–32,017); the unified N=31,873 is used. **Liability scale (Lee et al. 2011):** OA_knee h²_liab=0.094 (K=0.12); OA_selfrep h²_liab=0.073 (K=0.083); fracture h²_liab=0.062 (K=0.143). BMD is quantitative and not liability-converted. LDSC intercepts are near 1.0 (1.01–1.07) for 12 phenotypes; the fracture intercept is 0.96.

## Table 2. Two-sample Mendelian randomisation estimates

| Exposure → Outcome | n_snp | IVW β (fixed) | SE | p | IVW p (random) | OR (95% CI) | MR-Egger β (SE) | Egger p | Egger intercept p | Weighted median (95% CI) | I² | LOO sign flips |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **BMD total → fracture** (positive control) | 31 | −0.337 | 0.022 | 3.7×10⁻⁵¹ | 2.4×10⁻¹² | 1.401 (1.341–1.464) | −0.337 (0.049) | 4.0×10⁻¹² | 0.558 | −0.393 (−0.511 to −0.205) | 78% | 0 |
| BMD femoral neck → fracture | 25 | −0.351 | 0.030 | 7.7×10⁻³¹ | 2.3×10⁻⁶ | 1.420 (1.338–1.507) | −0.343 (0.068) | 4.9×10⁻⁷ | 0.018 | −0.365 (−0.414 to −0.129) | 83% | 0 |
| **BMD total → OA knee** | 31 | +0.078 | 0.013 | 3.9×10⁻¹⁰ | 0.006 | 1.081 (1.055–1.108) | +0.078 (0.029) | 0.007 | 0.810 | +0.065 (+0.027 to +0.177) | 81% | 0 |
| Lumbar spine → OA knee | 21 | +0.095 | 0.016 | 7.7×10⁻⁹ | 0.017 | 1.099 (1.064–1.135) | +0.089 (0.040) | 0.024 | 0.295 | +0.115 (+0.014 to +0.189) | 83% | 0 |
| Arms → OA knee | 16 | +0.075 | 0.016 | 1.7×10⁻⁶ | 0.070 | 1.078 (1.045–1.112) | +0.075 (0.043) | 0.080 | 0.982 | +0.026 (+0.019 to +0.217) | 86% | 0 |
| Legs → OA knee | 24 | +0.060 | 0.014 | 2.3×10⁻⁵ | 0.098 | 1.062 (1.033–1.092) | +0.058 (0.038) | 0.128 | 0.815 | +0.026 (−0.003 to +0.132) | 85% | 0 |
| Femoral neck → OA knee | 25 | +0.058 | 0.017 | 6.4×10⁻⁴ | 0.055 | 1.060 (1.025–1.095) | +0.058 (0.031) | 0.059 | 0.875 | +0.052 (−0.000 to +0.109) | 68% | 0 |
| OA knee → fracture | 69 | +0.027 | 0.032 | 0.395 | 0.542 | 0.973 (0.914–1.036) | +0.031 (0.045) | 0.490 | 0.476 | +0.019 (−0.112 to +0.126) | 49% | 0 |
| OA knee → BMD total | 69 | +0.031 | 0.032 | 0.326 | 0.481 | — | +0.036 (0.045) | 0.420 | 0.376 | +0.026 (−0.032 to +0.073) | 49% | 0 |
| OA knee → BMD femoral neck | 69 | +0.042 | 0.032 | 0.190 | 0.272 | — | +0.048 (0.038) | 0.206 | 0.169 | +0.023 (−0.063 to +0.122) | 30% | 0 |

Univariate instrument strength (F, R²_exposure) and the Steiger directionality test for every pair are in Table S13. β is the IVW estimate per 1 SD higher exposure; for the binary outcomes β is on the log-odds scale and OR = exp(|β|) with CI exp(|β| ± 1.96 SE), quoted per 1 SD *higher* BMD for the BMD→OA pairs and per 1 SD *lower* BMD for the BMD→fracture pairs, where higher BMD is protective; no OR is given for the reverse pairs, whose exposure is a binary liability rather than a 1-SD continuous exposure. MR-Egger p is a two-sided z-test on the slope. Weighted median is Bowden's estimator with weights (β_exposure/SE_outcome)², recomputed from the deposited harmonised instrument tables with 10,000-replicate bootstrap CIs (code/129); the IVW, MR-Egger, Cochran Q and I² columns reproduce from the same deposit to machine precision, and seven weighted-median values previously tabulated did not, so that column has been corrected (audit trail: code/output/genetics/wm_recomputed_vs_tabulated.csv). Forward BMD-exposure betas are per-SD on the verified N=31,873 scale (per-instrument closed-form SE check, within 0.1%); no rescaling was applied. Reverse MR is on the same verified SD scale. p-values, Egger-intercept p, I² and LOO sign flips are N-scale-invariant. Strict PLINK-equivalent clumping (K = strict instrument count). For BMD femoral neck → fracture the MR-Egger intercept is significant (p=0.018); the MR-Egger slope (−0.343) and weighted median (−0.365) are the pleiotropy-robust estimates. The random-effects p column uses the multiplicative convention SE_random = SE_fixed × √[max(1, Q/df)] (TwoSampleMR); an additive DerSimonian–Laird τ² convention would additionally retain the arms exposure under random effects (p=0.033, §3).

## Table 3. Statistical power of the two-sample MR estimates

| Exposure → Outcome          | K  | IVW β  | SE    | abs z | MDE (80% power) | Achieved power | Interpretation |
| --------------------------- | -- | ------ | ----- | ----- | --------------- | -------------- | -------------- |
| BMD total → fracture        | 31 | −0.337 | 0.022 | 15.05 | 0.063 | 1.00 | powered        |
| BMD femoral neck → fracture | 25 | −0.351 | 0.030 | 11.55 | 0.085 | 1.00 | powered        |
| BMD total → OA knee         | 31 | +0.078 | 0.013 | 6.26  | 0.035 | 1.00 | powered        |
| Lumbar spine → OA knee      | 21 | +0.095 | 0.016 | 5.77  | 0.046 | 1.00 | powered        |
| Arms → OA knee              | 16 | +0.075 | 0.016 | 4.78  | 0.044 | 1.00 | powered        |
| Legs → OA knee              | 24 | +0.060 | 0.014 | 4.23  | 0.040 | 0.99 | powered        |
| Femoral neck → OA knee      | 25 | +0.058 | 0.017 | 3.42  | 0.047 | 0.93 | powered        |
| OA knee → fracture          | 69 | +0.027 | 0.032 | 0.85  | 0.090 | 0.14 | underpowered   |
| OA knee → BMD total         | 69 | +0.031 | 0.032 | 0.98  | 0.089 | 0.17 | underpowered   |
| OA knee → BMD femoral neck  | 69 | +0.042 | 0.032 | 1.31  | 0.089 | 0.26 | underpowered   |

MDE = (1.96 + 0.84) × SE is the minimum |β| detectable at 80% power (two-sided α=0.05); achieved power is for the observed point estimate. |β|/SE is N-scale-invariant. K = strict instrument count.

## Supplementary Table S8. Colocalisation of total-body BMD with knee OA at all 31 instrument loci

Regions are ±500 kb around each strict-clumped total-body BMD instrument. PP.H0–PP.H4 are posterior probabilities under priors p1=p2=1×10⁻⁴, p12=5×10⁻⁶. No locus reaches PP.H4≥0.8. Six favour distinct causal variants (PP.H3≥0.8). The highest single-locus PP.H4 was 0.799 (rs11228269, chr11:68.29 Mb; BMD z=−8.65, knee-OA z=−4.06, direction-concordant; PP.H1=0.156), i.e. one locus fell 0.0006 short of the threshold, and the pre-specified 0.8 cut-off is retained. (Full locus-level table in the repository; see Data availability.)

## Supplementary Table S9. Multivariable MR of BMD → knee OA and BMD → fracture, conditioning on height

Primary configuration n_h=100 (EUR height; Yengo et al. 2022). n_H as tabulated (49–50) is the number that survived harmonisation, not the requested cap; n_BMD is 27/23/19/15/21 of 31/25/21/16/24 because six instruments are absent from the height GWAS. Attenuation = 100×(β_uni − β_MVMR)/β_uni. The two tabulated univariable baselines are β_uni from Table 2 and β_uni recomputed on the MVMR instrument set, and both are given across the full grid in Table S10. θ_H is the MVMR direct effect of height; Fc(BMD) and Fc(H) are Sanderson–Smith conditional F (≥10).

| Exposure → Outcome | n_BMD (of primary) | n_H used | β_uni (Table 2) | MVMR β_BMD | Atten. % vs Table 2 | θ_H | θ_H p | Fc(BMD) | Fc(H) | Egger p | I² |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| BMD total → OA knee | 27 (of 31) | 50 | +0.078 | +0.063 | 19.1 | +0.145 | 7.0×10⁻¹⁶ | 25.1 | 661.7 | 0.981 | 82.05 |
| Femoral-neck BMD → OA knee | 23 (of 25) | 49 | +0.058 | +0.049 | 14.6 | +0.165 | 1.4×10⁻¹⁹ | 15.8 | 561.8 | 0.855 | 79.71 |
| Lumbar-spine BMD → OA knee | 19 (of 21) | 50 | +0.095 | +0.068 | 27.6 | +0.153 | 2.5×10⁻¹⁷ | 16.8 | 823.5 | 0.798 | 82.73 |
| Arms BMD → OA knee | 15 (of 16) | 50 | +0.075 | +0.065 | 13.5 | +0.159 | 1.4×10⁻¹⁸ | 19.4 | 889.5 | 0.993 | 83.39 |
| Legs BMD → OA knee | 21 (of 24) | 49 | +0.060 | +0.040 | 33.2 | +0.152 | 4.7×10⁻¹⁷ | 20.7 | 448.7 | 0.988 | 83.40 |
| **BMD total → fracture (positive control)** | 27 (of 31) | 50 | −0.337 | −0.343 | −1.7 | +0.166 | 1.0×10⁻⁷ | 25.1 | 661.7 | 0.311 | 64.49 |

β per 1 SD on the N-corrected scale. Height-attributable fraction 9.5–33.2% across the two tabulated baselines (13.5–33.2% vs Table 2 and 9.5–21.2% vs the MVMR-set baseline), a lower bound for total size; the 67–91% height-independent residual is not resolved into mass versus geometry. Positive control moved −1.7% to −3.0% (control intact; −2.8% on the MVMR-set baseline).


## Supplementary Table S14. Population-ancestry composition of the 13 GWAS

| GWAS (phenotype) | Source | N (effective) | Ancestry composition | LDSC reference |
| --- | --- | --- | --- | --- |
| BMD (10 regional sites) | Qian et al. 2025 [2] | 31,873 each | 100% European (UK Biobank) | 1000G Phase 3 EUR |
| OA_knee (radiographic) | GO Consortium 2025 [1] | ≈1,208,397 | Multi-ancestry: EUR 87.31%, EAS 7.09%, AFR 3.08%, SAS 1.09%, HIS 0.91%, ADM 0.53% | 1000G Phase 3 EUR |
| OA_selfrep (self-reported) | IEU OpenGWAS `ukb-b-14486` [3] | 462,933 | European (UK Biobank) | 1000G Phase 3 EUR |
| Fracture (any) | GEFOS 2018 [10] | 264,973 | European | 1000G Phase 3 EUR |

The BMD exposure instruments are European-only, whereas the knee-OA outcome is 87.31% European (multi-ancestry meta-analysis); LDSC uses a European-only reference, so cross-population linkage-disequilibrium differences between the BMD instruments and the multi-ancestry knee-OA outcome can bias r_g and the MR estimates. The OA_knee percentages are taken from the GWAS Catalog entry and the originating Nature report of the multi-ancestry meta-analysis.


## Supplementary Table S15. Genetic correlations by phenotype block

| Block          | n pairs | r_g range        | r_g mean | Significant (q<0.05) |
| -------------- | ------- | ---------------- | -------- | -------------------- |
| OA × OA        | 1       | +0.838           | +0.838   | 1/1                  |
| OA × BMD       | 20      | +0.117 to +0.250 | +0.169   | 20/20                |
| Fracture × BMD | 10      | −0.535 to −0.401 | −0.460   | 10/10                |
| BMD × BMD      | 45      | +0.689 to +0.960 | +0.812   | 45/45                |
| OA × Fracture  | 2       | +0.056 to +0.119 | +0.087   | 1/2                  |


---

## References

1. Hatzikotoulas K, Southam L, Stefansdottir L, et al. Translational genomics of osteoarthritis in 1,962,069 individuals. *Nature*. 2025;641(8065):1217–1224. doi:10.1038/s41586-025-08771-z. PMID 40205036.
2. Qian Y, Xia J, Wang P, et al. Genome-wide Association Studies of over 30,000 Samples with Bone Mineral Density at Multiple Skeletal Sites and Its Clinical Relevance. *Genomics Proteomics Bioinformatics*. 2025;23(5):qzaf097. doi:10.1093/gpbjnl/qzaf097. PMID 41206123.
3. Elsworth B, Lyon M, Alexander T, et al. The MRC-IEU UK Biobank GWAS pipeline. *bioRxiv*. 2018. doi:10.1101/340232. (Self-reported OA accessed as `ukb-b-14486` via OpenGWAS, <https://gwas.mrcieu.ac.uk/datasets/ukb-b-14486/>.)
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
16. Alhassan E, Nguyen K, Hochberg MC, Mitchell BD. Causal Factors for Osteoarthritis: A Scoping Review of Mendelian Randomization Studies. *Arthritis Care Res (Hoboken)*. 2024;76(3):366–375. doi:10.1002/acr.25252. PMID 37846209.
17. Liu H, Zhu J, Jiang H, et al. Using Mendelian randomization analyses in osteoarthritis: state of art and future perspectives. *Ann Med*. 2025;57(1):2533428. doi:10.1080/07853890.2025.2533428. PMID 40697045.
18. Kröger H, Vainio P, Nieminen J, Kotaniemi A. Comparison of different models for interpreting bone mineral density measurements using DXA and MRI technology. *Bone*. 1995;17(2):157–159. doi:10.1016/s8756-3282(95)00162-x. PMID 8554924.
19. Yengo L, Vedantam S, Marouli E, et al. A saturated map of common genetic variants associates polygenic risk with human height. *Nature*. 2022;610(7933):704–712. doi:10.1038/s41586-022-05275-y.
20. Sanderson E, Davey Smith G, Windmeijer F, Bowden J. An examination of multivariable Mendelian randomization in the single-sample and two-sample summary data settings. *Int J Epidemiol*. 2018;48(3):713–727. doi:10.1093/ije/dyy262.
21. Burgess S, Jenkins DA, Davies NM, et al. Multivariable Mendelian randomization: the use of pleiotropic variants to estimate causal effects of multiple risk factors on a single outcome. *Am J Epidemiol*. 2015;181(4):251–260. doi:10.1093/aje/kwu283.
22. Hartley AE, Sanderson ECM, Granell R, Paternoster L, Zheng J, Davey Smith G, et al.; Genetics of Osteoarthritis Consortium; Gregson CL, Tobias JH. Using multivariable Mendelian randomization to estimate the causal effect of bone mineral density on osteoarthritis risk, independently of body mass index. *Int J Epidemiol*. 2022;51(4):1254–1267. doi:10.1093/ije/dyab251. PMID 34897459.
23. Dillon CF, Rasch EK, Gu Q, Hirsch R. Prevalence of knee osteoarthritis in the United States: arthritis data from the Third National Health and Nutrition Examination Survey 1991–94. *J Rheumatol*. 2006;33(11):2271–2279. PMID 17013996.
24. Welling M, Auvinen J, Lehenkari P, Männikkö M, Karppinen J, Eskola PJ. Association between height and osteoarthritis of the knee and hip: the Northern Finland Birth Cohort 1966 Study. *Int J Rheum Dis*. 2017;20(9):1095–1104. doi:10.1111/1756-185X.13059. PMID 28296269.
25. Hunter DJ, Niu J, Zhang Y, et al. Knee height, knee pain, and knee osteoarthritis: the Beijing Osteoarthritis Study. *Arthritis Rheum*. 2005;52(5):1418–1423. doi:10.1002/art.21017. PMID 15880346.
26. Sun Y, Li Y, Yu T, Zhang J. Causal associations of anthropometric measurements with osteoarthritis: a Mendelian randomization study. *PLoS One*. 2023;18(1):e0279198. doi:10.1371/journal.pone.0279198. PMID 36716300.
27. Goldring MB, Goldring SR. Articular cartilage and subchondral bone in the pathogenesis of osteoarthritis. *Nat Rev Rheumatol*. 2016;12(11):632–644. doi:10.1038/nrrheum.2016.148. PMID 27652499.
28. Burr DB, Gallant MA. Bone remodelling in osteoarthritis. *Nat Rev Rheumatol*. 2012;8(11):665–673. doi:10.1038/nrrheum.2012.130. PMID 22868925.
29. Zhen G, Wen C, Jia X, et al., Cao X. Inhibition of TGF-β signaling in mesenchymal stem cells of subchondral bone attenuates osteoarthritis. *Nat Med*. 2013;19(6):704–712. doi:10.1038/nm.3143. PMID 23685840.
30. Findlay DM, Kuliwaba JS. Bone–cartilage crosstalk: a conversation for understanding osteoarthritis. *Bone Res*. 2016;4:16028. doi:10.1038/boneres.2016.28. PMID 27672480.
31. Zhu R, Xing X, Bian J, Zhang X, Ge L, Cai G. Causal association between bone mineral density and the risk of joint replacement in patients with osteoarthritis: a Mendelian randomization study. *Clin Rheumatol*. 2025;44(2):823–830. doi:10.1007/s10067-024-07289-5. PMID 39738847.
32. Looker AC, Beck TJ, Orwoll ES. Does body size account for gender differences in femur bone density and geometry? *J Bone Miner Res*. 2001;16(7):1291–1299. doi:10.1359/jbmr.2001.16.7.1291.
33. Nevitt MC, Felson DT. High bone density and radiographic osteoarthritis: questions answered and unanswered. *Osteoarthritis Cartilage*. 2020;28(9):1151–1153. doi:10.1016/j.joca.2020.05.006. PMID 32417559.
34. Hardcastle SA, Dieppe P, Gregson CL, et al. Osteophytes, enthesophytes, and high bone mass: a bone-forming triad with potential relevance in osteoarthritis. *Arthritis Rheumatol*. 2014;66(9):2429–2439. doi:10.1002/art.38729. PMID 24910132.
35. Faber BG, Frysz M, Zheng J, et al. The genetic architecture of hip shape and its role in the development of hip osteoarthritis and fracture. *Hum Mol Genet*. 2025;34(3):207–217. doi:10.1093/hmg/ddae169.
36. Funck-Brentano T, Nethander M, Moverare Skrtic S, Richette P, Ohlsson C. Causal factors for knee, hip, and hand osteoarthritis: a Mendelian randomization study in the UK Biobank. *Arthritis Rheumatol*. 2019;71(10):1634–1641. doi:10.1002/art.40928. PMID 31099188.
37. Estrada K, Styrkarsdottir U, Evangelou E, et al. Genome-wide meta-analysis identifies 56 bone mineral density loci, provides insights into the genetic architecture of osteoporosis, and reveals a role for RSPO3. *PLoS Genet*. 2012;8(1):e1001217. doi:10.1371/journal.pgen.1001217. PMID 22253714.
38. Pulit SL, Stoneman C, Morris AP, et al. Meta-analysis of genome-wide association studies for body fat distribution in 694,649 individuals of European ancestry. *Hum Mol Genet*. 2019;28(1):166–174. doi:10.1093/hmg/ddy327. PMID 30239722. (BMI summary statistics from the GIANT–UK Biobank meta-analysis: Zenodo, https://doi.org/10.5281/zenodo.1251813.)

## Acknowledgements and declarations

**Acknowledgements.** The authors thank the investigators and participants of the source GWAS consortia for making summary statistics publicly available. No funding was received for this work.

**Declaration of generative AI in scientific writing.** An AI writing assistant was used during preparation for language editing, formatting and internal-consistency checking; all analyses, results, interpretations and conclusions are the authors' own and were verified against the deposited code and outputs. The authors take full responsibility for the content of the article.

## Author contributions (CRediT)

Contributions are reported under the CRediT (Contributor Roles Taxonomy). # equal contribution; * corresponding author.

| Author          | Conceptualization | Data curation | Formal analysis | Funding acquisition | Investigation | Methodology | Project administration | Resources | Software | Supervision | Validation | Visualization | Writing – original draft | Writing – review & editing |
| --------------- | :---------------: | :-----------: | :-------------: | :-----------------: | :-----------: | :---------: | :--------------------: | :-------: | :------: | :---------: | :--------: | :-----------: | :----------------------: | :------------------------: |
| Zhenxin Zheng¹# |         ✓         |               |        ✓        |                     |       ✓       |      ✓      |                        |           |     ✓    |             |      ✓     |       ✓       |             ✓            |              ✓             |
| Song Chen¹#     |         ✓         |       ✓       |                 |                     |       ✓       |             |                        |     ✓     |          |             |      ✓     |               |             ✓            |              ✓             |
| Shenbo Huang¹   |                   |       ✓       |                 |                     |       ✓       |             |                        |     ✓     |          |             |      ✓     |               |                          |              ✓             |
| Chao Qian¹      |                   |       ✓       |                 |                     |       ✓       |             |                        |     ✓     |          |             |      ✓     |               |                          |              ✓             |
| Hang Yuan¹      |                   |       ✓       |                 |                     |       ✓       |             |                        |     ✓     |          |             |      ✓     |               |                          |              ✓             |
| Huricha Bao¹*   |         ✓         |               |                 |          ✓          |               |      ✓      |            ✓           |           |          |      ✓      |            |               |                          |              ✓             |
| Qiuyan Jiang²*  |         ✓         |               |                 |          ✓          |               |      ✓      |            ✓           |           |          |      ✓      |            |               |                          |              ✓             |

## Authors, affiliations, funding and competing interests

- **Funding.** None.
- **Competing interests.** The authors declare no competing interests.
- **Data availability.** Summary statistics are from publicly published GWAS (Refs 1–3,10,19,38) and were not redistributed. Provenance and accessions are listed above. The analysis code and all reproduction material have been deposited in a public GitHub repository (https://github.com/Fisherzzx001/OA-BMD-fracture-MR) and archived with a Zenodo DOI via the GitHub–Zenodo integration (https://doi.org/10.5281/zenodo.23400860; this concept DOI always resolves to the most recent version and lists every version of the deposit). The deposit includes every script under `code/` together with the result tables and source data under `code/output/` required to reproduce every table and figure (including the MRlap sample-overlap correction); bulk intermediate files that are regenerated deterministically from the source summary statistics by the deposited scripts (for example PLINK clumping inputs and harmonised files) are not archived, together with the supplementary files enumerated in the Supplementary Material section (STROBE-MR checklist, Table S10 MVMR grid, Table S11 SOS predictions and falsifying tests). All source GWAS were approved by their original ethics committees. This secondary analysis required no new ethical approval. An independent preregistration of this study — distinct from the archived comorbidity-network protocol of the OA project — is deposited in the repository (`docs/118_OSF预注册_OA_BMD_fracture.md`) and archived via the Zenodo DOI above; it declares P1–P6 as derived predictions and F1–F10 as registered tests (F4 cross-trait and eQTLGen arms, and F6–F9, executed with decision rules fixed before execution; F1–F3, F5 and F10 prospective, requiring external data). An OSF mirror is pending.
- **Declaration of competing interest.** The authors declare that they have no known competing financial interests or personal relationships that could have appeared to influence the work reported in this paper.

---

## Supplementary Material

The following supplementary files accompany this article and are referenced in the text. All are provided as machine-readable Markdown/CSV in the submission repository (`docs/Supplementary/`).

- **Figure S1.** Per-SNP scatter plots for the six remaining MR pairs.
- **Figure S2.** Leave-one-out stability of all ten pairs.
- **Figure S3.** External validation of the reconstructed BMD scale against GEFOS.
- **Figure S4.** Study design and MR flow (STROBE-MR).
- **Figure S5.** Sensitivity of OA_knee h² to representative N (r_g is N-invariant).
- **Figure S6.** Fracture sample-size decision.
- **Figure S7.** Pleiotropy and heterogeneity diagnostics.
- **Figure S8.** Colocalisation of BMD with knee OA and its power.
- **Figure S9.** Instrument selection and attrition flow (STROBE-MR).
- **Figure S10.** Multivariable MR decomposition of the BMD→OA effect conditional on height.
- **Figure S11.** Multi-estimator robustness grid for the ten MR pairs.
- **Figure S12.** Quantitative comparison with published MR estimates.
- **Figure S13.** Forest plot of the selected genetic correlations (any-fracture, knee OA and self-reported OA versus the ten BMD sites; the two OA definitions with each other and with fracture).
- **Figure S14.** SNP-based heritability (observed scale) of the 13 skeletal phenotypes, with liability-scale h² for the three binary traits.
- **Figure S15.** MR-PRESSO funnel plots for the four headline pairs.
- **Figure S16.** Statistical power of the two-sample MR design.
- **Figure S17.** Sample-overlap diagnostics (ρ_GE) for the ten MR pairs.
- **Table S8.** Colocalisation of total-body BMD with knee OA at all 31 instrument loci.
- **Table S9.** Multivariable MR of BMD→knee OA and BMD→fracture, conditioning on height (primary n_h=100 / EUR configuration).
- **Table S10.** Multivariable MR grid: full sensitivity grid (EUR and EUR_noUKB; n_h ∈ {0,100,300,1,000}).
- **Table S11.** Systemic Osteogenic Set-point (SOS) theory: predictions (P1–P6) and falsifying tests (F1–F7).
- **Table S12.** Exposure (BMD-site) selection and instrument counts: the 11-site source, the a priori cranial-head exclusion, genome-wide-significant locus counts per site, and which sites entered MR (`variant_selection_table.csv` in `code/output/genetics/`).
- **Table S13.** Univariate instrument strength and Steiger directionality for the ten MR pairs: per-SNP and aggregate F, R²_exposure, Steiger z and one-sided p, and direction supported (`instrument_strength_steiger.csv` in `code/output/genetics/`).
- **Table S14.** Population-ancestry composition of the 13 GWAS and the LDSC reference used for each.
- **Table S15.** Genetic correlations by phenotype block (n pairs, r_g range, r_g mean, proportion significant at q<0.05).
- **Table S16.** Multivariable MR of BMD→knee OA and BMD→fracture, conditioning on BMI (primary n_bmi=100 / GIANT–UK Biobank configuration).
- **STROBE-MR checklist.** 22-item reporting-compliance table (Skrivankova et al., BMJ 2021).
- **Graphical Abstract.** The OA–BMD–fracture triangle: higher BMD protects against fracture yet raises knee-OA risk.

Supplementary figures S1–S17 and the Graphical Abstract are provided as high-resolution PNG in the submission repository (`Figures/`); supplementary tables S8–S16 and the STROBE-MR checklist are provided as Markdown/CSV.
