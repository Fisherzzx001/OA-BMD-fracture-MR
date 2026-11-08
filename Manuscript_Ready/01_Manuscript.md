Genetically higher BMD raises knee osteoarthritis risk yet protects against fracture on pleiotropy-robust evidence

**Authors.** Zhenxin Zheng¹,#; Song Chen¹,#; Shenbo Huang¹; Chao Qian¹; Hang Yuan¹; Huricha Bao¹,*; Qiuyan Jiang²,*
¹ Department of Orthopedics, The Quzhou Affiliated Hospital of Wenzhou Medical University, Quzhou People's Hospital, Quzhou, Zhejiang 324000, China.
² Department of Surgery, The Quzhou Affiliated Hospital of Wenzhou Medical University, Quzhou People's Hospital, Quzhou, Zhejiang 324000, China.
Equal contribution; * joint corresponding authors.

**Running title.** BMD raises knee OA yet protects against fracture

## Highlights
- Higher genetically predicted BMD raises knee OA risk yet protects against fracture.
- BMD→OA was positive at all five analysed sites; 2 of 5 held under random effects.
- BMD variants colocalise with fracture (9/31), not knee OA (0/31): endpoint-specific.
- Height bounds size pathway to 9.5–33.2%; 67–91% is neither height- nor BMI-mediated.
- A systemic bone-forming set-point is a testable hypothesis consistent with data.

**Keywords.** osteoarthritis; bone mineral density; fracture; linkage disequilibrium score regression; Mendelian randomisation; genetic correlation

---

## Abstract

**Objective.** Whether higher bone mineral density (BMD) shares aetiology with, or causes, osteoarthritis (OA) is unsettled. We quantified genetic sharing and causal direction across 13 skeletal phenotypes and tested whether one BMD-raising set-point has opposing fracture and joint consequences.

**Design.** LDSC and two-sample MR were applied to 13 GWAS summary-statistic datasets (two OA definitions, ten regional BMD sites, any-fracture). LDSC covered all 78 phenotype pairs; MR ran in eight phenotypes with clumped instruments — five of ten prespecified BMD sites, both OA definitions and any-fracture — using BMD→fracture as the positive control. Further configurations conditioned on height and separately on BMI, applied MRlap sample-overlap correction, and reverse MR.

**Results.** BMD sites inter-correlated strongly (r_g 0.69–0.96) and negatively with fracture (−0.40 to −0.54); both OA definitions correlated positively with every BMD site (+0.12 to +0.25). Genetically higher BMD raised knee-OA risk at all five analysed sites — OR 1.060–1.099 per SD (95% CI 1.025–1.135) — and protected against fracture (OR 1.40, 95% CI 1.34–1.46 per SD lower BMD; p=3.7×10⁻⁵¹). Heterogeneity was high (I² 68–86%) but matched the positive control (78–83%), so the causal reading rests on pleiotropy-robust estimators; under random effects 2 of 5 remained significant. Conditioning on height attenuated the OA effect by 9.5–33.2% (a lower bound for size); conditioning separately on BMI left it positive and significant at all five sites (−3.7% to 46.0% attenuation), so the 67–91% height-independent residual is neither height- nor BMI-mediated — not identified as bone mass. Reverse estimates were non-significant and underpowered (power 0.14–0.26), so the bidirectional MR provided no directional evidence. MRlap moved all ten estimates further from the null (median 20.5%), reversing none.

**Conclusions.** One genetically indexed bone phenotype is protective at the whole-bone endpoint and harmful at the joint endpoint, 67–91% of the joint effect independent of height and not BMI-mediated. The joint effect is modest per SD (OR 1.06–1.10 areal BMD) yet attaches to a ubiquitous, modifiable exposure. These estimates do not support changing prescribing practice; one completed knee-OA structural-progression trial randomised a BMD-raising drug against placebo and reported benefit, yet three years is not lifelong exposure.

---

## 1. Introduction

Osteoarthritis (OA) is a prevalent disorder of ageing [1]; areal bone mineral density (BMD), the clinical measure of skeletal strength, is measured in the same populations. Their epidemiological association has been debated since Foss and Byers reported raised femoral-head BMD in osteoarthritic hips in 1972 [2]. Observational studies, reviewed in [3], link higher areal BMD to increased OA prevalence and incidence, although osteoporosis and OA were long regarded as opposing conditions [4]. The direction matters clinically: if high BMD predisposes to OA, joint endpoints deserve monitoring in bone-anabolic trials, though whether such a cost is mediated by bone mass or skeletal size is itself unresolved.

Genome-wide data can separate these explanations. LDSC estimates genetic correlation (r_g) from summary statistics, tolerates sample overlap, and is unbiased by stratification when the intercept is near unity [5]. Two-sample MR uses genetic variants as instruments to interrogate causality [6]. Applied to OA and BMD, scoping and methodological reviews have given partial answers [7, 8]: Hackinger et al. reported r_g=+0.18 between hip/knee OA and lumbar-spine BMD (p=0.022) [9]; Zhou et al. reported r_g=+0.13 by LDSC and +0.17 by high-definition likelihood between hip OA and femoral-neck BMD, supporting a negative femoral-neck BMD→hip OA effect [10]; and a 2025 two-sample MR reported total BMD as a significant positive mediator of knee-OA risk (mediation proportion 3.59%) [11]. Our BMD→knee-OA analysis, by contrast, finds a positive effect contrasting with Zhou's negative femoral-neck estimate (discussed below).

The most direct predecessor is Hartley et al. [12], who used multivariable MR to estimate the causal effect of BMD on OA independently of BMI and reported a positive effect at both hip (OR 1.28, 1.05–1.57) and knee OA (OR 1.40, 1.20–1.63). Our study extends that work in four ways: (i) ten regional BMD sites rather than one composite, so site specificity is tested directly (five of ten entered MR; §2.3); (ii) conditioning on adult height and separately on BMI, separating size from adiposity; (iii) an internal BMD→fracture positive control; and (iv) a full 13×13 LDSC matrix rather than pairwise estimates.

Three questions remain open. First, most analyses included at most two BMD sites, so site specificity is unknown despite OA being joint-localised and BMD site-regulated. Second, OA has usually been represented by a single definition, leaving their bone relationship uncertain. Third, the OA–BMD–fracture triangle has not been analysed jointly; modelling all three tests the pipeline and whether the OA–BMD association carries fracture relevance.

We addressed all three in 13 skeletal phenotypes using LDSC across all 78 phenotype pairs, two-sample MR with an internal positive control, and reverse MR.

## 2. Methods

This study followed the STROBE-MR reporting guideline [13]. Analyses were performed in Python 3.13 with numpy, pandas, scipy and matplotlib. The sample-overlap correction used R 4.6.1 with MRlap v0.0.3.3.

**Data sources.** Thirteen GWAS summary-statistic datasets were analysed: radiographic knee OA (OA_knee; GO Consortium 2025, N≈1.21M) [14], self-reported OA (OA_selfrep; IEU OpenGWAS `ukb-b-14486`, N=462,933; 38,472 cases) [15], ten regional BMD sites (Qian et al. 2025, N≈1.873 each) [16], and any-fracture (GEFOS 2018, N=264,973; 37,857 cases) [17]; full accessions are in the Data availability statement. All four sources contribute sex-combined primary estimates, so every estimate is sex-averaged; a sex-stratified sensitivity analysis (F8) confirmed BMD→OA positivity in both sexes (§4). All LDSC used 1000 Genomes Phase 3 European LD scores; ancestry composition is in Table S14. The European-only reference against the 87.31%-European knee-OA GWAS means the near-1.0 intercept does not exclude population-stratification bias (primary limitation, §4).

**LDSC.** Genetic correlation (r_g) was estimated by LDSC [5] in Python (Bulik-Sullivan) on 690,150 HapMap3 SNPs with 1000 Genomes Phase 3 European LD scores, using delete-one-block jackknife (100 blocks) for SEs. Palindromic SNPs were excluded; alleles were harmonised to one reference so every matrix cell uses the same variants. Both Bonferroni (α=0.05/78) and Benjamini–Hochberg FDR were applied across 78 phenotype pairs. A spectral decomposition of the resulting 13×13 r_g matrix was performed as an exploratory probe of latent structure (§3, Fig. S22).

**MR.** Valid instruments must satisfy three assumptions: relevance; independence from exposure–outcome confounders; and exclusion restriction — no effect on the outcome except through the exposure. We address each below. Instruments were selected at p<5×10⁻⁸ and clumped at r²<0.001 over 10 Mb against 1000 Genomes Phase 3 EUR (strict PLINK-equivalent clumping benchmarked against PLINK; 20–35 retained per BMD site, 116 for knee OA). The BMD exposure set followed a prespecified rule: of the 11 DXA-derived BMD sites in Qian et al. 2025, the cranial head site — the only non-load-bearing region, whose instruments carry non-skeletal pleiotropy (head BMD is genetically correlated with intracranial aneurysm [16]) — was excluded a priori, leaving ten canonical DXA regions; a site enters MR only if strict clumping completes. Clumping completed for five sites (arms, femoral neck, lumbar spine, legs, total body), each ≥15 independent instruments; the remaining five showed comparable genome-wide signal (229–345 loci at p<5×10⁻⁸) but strict clumping was not completed (1000G EUR panel unavailable) and they were deferred rather than analysed with a non-prespecified set (Table S12). Instrument shortage drove no exclusion. Harmonisation against the outcome GWAS reduced the analysed sets to 16–31 instruments per forward pair and 69 per knee-OA-exposure pair; 47 of 116 knee-OA instruments (41%) were lost (Fig. S9), an attrition we report. The 69 retained instruments were identical across all three knee-OA outcomes; the dropped instruments were not effect-size selected (median −log10P 9.10 vs 9.89 retained; 43% stronger signal). Random-effects IVW used the multiplicative convention SE_random = SE_fixed × √[max(1, Q/df)] (TwoSampleMR), not additive DerSimonian–Laird τ²; where the convention changes an inference, both are reported (§3, Table 2 note).

**Sample overlap.** All pairs share UK Biobank participants, partially violating the two-sample assumption. Under the Burgess approximation [18], the IVW bias equals the overlap fraction × r_g and acts in the direction of r_g (pairwise UKB overlaps and bias bounds in Table S14); because r_g shares the sign of the MR effect in all ten pairs (as does the SNP-effect correlation ρ_GE, Fig. S17), it acts away from the null. The two other IVW biases — weak-instrument (NOME) attenuation and winner's curse from instrument selection — act towards the null. We report MR-Egger [19], weighted median [20] and MR-PRESSO alongside IVW, and applied MRlap [21] (corrects weak-instrument bias and winner's curse while accounting for overlap) to all ten pairs. The joint correction moved every estimate further from the null (median +20.5%; observed-versus-corrected difference p<0.05 for 9/10 pairs), showing towards-null components dominate the away-from-null overlap component, reversing no direction or significance. Primary IVW estimates are conservative for the combined bias, not overlap alone; the internally estimated per-pair cross-trait intercepts and the default-setting sensitivity are in Table S14. Worst-case Burgess biases from the tabulated overlaps are ≤0.7% of the BMD→knee-OA estimate and ≤4.5% of the positive-control estimate (Table S14); MRlap retained all ten directions (median +20.5% away from the null), so the positive joint-endpoint effect survives overlap correction.

**Effect-scale reconstruction.** Two verified scale corrections were applied. (i) *BMD exposure*: the deposited placeholder N=400,000 gives SE = 1/√[2·N·f·(1−f)] matching the verified N=31,873 to within 0.1%, so no rescaling is needed. (ii) *Knee-OA outcome*: the file supplies only N, Z and MAF; naive β=Z·SE is a linear-probability effect, not log-odds, which requires N_eff=N·K(1−K). K was estimated, not assumed: at 39 SNPs with deposited GO log-odds SE, N_eff=1/[SE²·2f(1−f)] over munged N_total gives K(1−K) (median 0.0871), i.e. K=0.0964, so we applied 1/√[K(1−K)] = 3.389 (Table S17). The Z ratio is unity by construction and cannot identify K; it is a same-source check only. Fracture needed no correction (log-odds supplied directly). Each correction multiplies β and SE by one constant, preserving the per-SNP z; for the femoral-neck, lumbar and arms BMD→OA pairs it also removed one instrument each (26→25, 22→21, 17→16), so their pooled z and p use the final set. Table 2 reports the final reconstructed estimates; absolute BMD→OA effects move away from the null (OR 1.060–1.099, not 1.017–1.028) and the reverse effects towards it.

**Colocalisation.** Bayesian colocalisation (coloc.abf; Giambartolomei et al. 2014 [22]) tested the three a priori loci (WNT16, SOX6, GDF5) and all 31 BMD instrument loci (±500 kb; PP.H4≥0.8), with a BMD–fracture control.

**MVMR.** Multivariable MR conditioned the BMD→OA effect on adult height (Yengo et al. 2022, N=5.0M) [23] as a skeletal-size proxy, using the MVMR-IVW estimator [24] and MVMR-Egger [25], judged by the Sanderson–Smith conditional F (Fc≥10). The primary configuration requested up to n_h=100 height instruments (49–50 survived harmonisation); a grid (n_h ∈ {0,100,300,1,000}) and both EUR and EUR_noUKB sources are in the Supplementary Material. Height was the primary size proxy: BMI is a major knee-OA risk factor and part of the adiposity pathway through which a BMD-raising process could raise knee-OA risk, whereas height indexes skeletal size alone. A second MVMR conditioned the same effect on BMI (GIANT–UK Biobank, Pulit et al. 2019, N≈806,800) [26] with the same estimators, n_bmi=100 primary (61–63 instruments survived) and a grid (n_bmi ∈ {100,300,800}); conditional F statistics and estimates are in Supplementary Table S16. Height and BMI were conditioned in separate analyses, so each proxy is adjusted one at a time rather than jointly.

## 3. Results

SNP-based heritability ranged 0.26–0.39 across the ten BMD sites and was 0.035 (knee OA), 0.022 (self-reported OA) and 0.026 (fracture); 12/13 LDSC intercepts fell between 1.009 and 1.069 (fracture 0.96). Because LDSC uses a European-only reference while the knee-OA GWAS is multi-ancestry (Table S14), the near-1.0 intercept does not exclude population-stratification or ancestry-mismatch bias (§4). The 13×13 r_g matrix showed three blocks (Figure 1, Table 2). The ten BMD sites inter-correlated strongly (r_g 0.69–0.96, mean 0.81). Both OA definitions correlated positively with all ten BMD sites (knee OA +0.14 to +0.19; self-reported OA +0.12 to +0.25; 20/20 significant) and with each other (r_g=+0.84). The two definitions also agreed at the variant level: across 1,980 genome-wide-significant knee-OA loci, 94.4% had concordant signed Z with self-reported OA. Any-fracture correlated negatively with every BMD site (r_g −0.40 to −0.54; all q_FDR<0.001). It correlated positively with knee OA (r_g=+0.119, q=3.3×10⁻⁹) but non-significantly with self-reported OA (+0.056, p=0.051); 77/78 pairs were significant at both BH-FDR and Bonferroni, the sole exception being self-reported-OA × fracture. A spectral decomposition resolved two orthogonal axes — a bone-mass/fracture axis (PC1 66.7% of variance; all ten BMD sites loading +0.30, fracture −0.18, both OA definitions ≈0) and an OA axis (PC2 14.1%; OA_knee +0.69, OA_selfrep +0.68, BMD ≈0) — so BMD and OA occupy separate genetic dimensions (Fig. S22; exploratory, OA h² 0.035).

The positive control behaved as expected: two-sample MR using 31 independent total-body BMD instruments showed genetically predicted higher BMD reduces fracture risk (IVW β=−0.337 per SD, SE 0.022, p=3.7×10⁻⁵¹; OR 1.40, 95% CI 1.34–1.46, per SD lower BMD). The estimate held under random effects (p=2.4×10⁻¹²), the MR-Egger intercept was null (p=0.558), weighted median agreed (−0.364), and leave-one-out produced no sign reversal. The femoral-neck set (25 SNPs) gave a concordant result (β=−0.351; OR 1.42, 1.34–1.51), but its MR-Egger intercept was not null (p=0.018), a directional-pleiotropy signal we report rather than absorb; its pleiotropy-robust slope (−0.343) and weighted median (−0.365) agree with IVW, so the control calibrates the pipeline despite its own non-clean diagnostic. The conservative 1-Mb set gave β=−0.366 (OR 1.44); these agree with the published GEFOS estimate (OR 1.55, 1.48–1.63) [17].

MR supported a positive effect of BMD on knee OA across all five exposures (Table 2): total-body β=+0.078 (p=3.9×10⁻¹⁰), lumbar-spine +0.095 (p=7.7×10⁻⁹), arms +0.075 (p=1.7×10⁻⁶), legs +0.060 (p=2.3×10⁻⁵), femoral-neck +0.058 (p=6.4×10⁻⁴) — OR 1.060–1.099 per SD higher BMD, against OR 1.40 (1.20–1.63) per SD in the largest one-sample estimate [12], a 3.6–5.8-fold log-odds difference (adjudicated §4). MR-Egger intercepts [19] were all p>0.24 and weighted medians matched the IVW sign at all five sites, with bootstrap 95% CIs excluding the null at three of five; with 16–31 instruments the intercept test has limited power, so a non-significant intercept excludes detectable unbalanced pleiotropy rather than pleiotropy itself. Leave-one-out produced no sign reversals and directions agreed with r_g (+0.12 to +0.25). Heterogeneity was high (I²=68–86%); under random-effects IVW only two of five exposures remained significant (total-body p=0.006, lumbar-spine p=0.017; a second random-effects convention additionally retains arms, p=0.033), so the count depends on the random-effects convention (Table S17). All five effects remained positively causal after MR-PRESSO outlier removal, with four of five strengthening (Table S17). The positive control was comparably heterogeneous (I² 78% and 83%), indicating high I² is a property of the BMD instrument sets rather than an OA-specific IVW failure.

Knee OA as exposure gave non-significant effects on total-body BMD (β=+0.031, p=0.33), femoral-neck BMD (β=+0.042, p=0.19) and fracture (β=+0.027, p=0.40); these were weakly powered (0.14–0.26; Table 3) and excluded only |β|≥0.09 in SD of BMD per unit OA liability. Removing the MR-PRESSO outliers flagged for these pairs moved two point estimates substantially (total-body BMD to +0.051, fracture to +0.059) without altering their non-significance.

Instrument strength and directionality. Every forward BMD instrument exceeded the conventional F>10 threshold (minimum per-SNP F 29.8–31.1; aggregate F 44–74; R²_exposure 0.034–0.063; Tables 3 and S13), so weak-instrument bias does not explain the forward effects. The Steiger directionality test supported BMD as the causally prior trait for every forward pair (z = 3.05–6.11, one-sided p = 4.8×10⁻⁴ to 5.0×10⁻¹⁰) and was non-significant for the three reverse pairs (|z|<0.1, p≈0.48–0.49). That non-result does not discriminate direction: with 69 instruments and a null reverse estimate the test lacks power to separate the models and is not read as evidence BMD rather than OA drives the association.

The MR-PRESSO distortion test was significant for 9/10 pairs (expected from the widespread horizontal pleiotropy common in complex-trait MR [27]) and identified 0–3 outliers per pair; removing outliers did not reverse any direction, and strict-clumped estimates reproduced the conservative 1-Mb clumping. All significant findings survived Bonferroni correction (α=0.005). The IVW estimate is therefore not read as a single clean causal parameter (Table 2, Fig. S11, S20).

On the verified per-SD scale, the BMD→fracture OR (1.42 femoral-neck) lay 8.4% below GEFOS (1.55) on the odds-ratio scale and 20.0% below on the log-odds scale, with overlapping CIs (we compare on the log scale). The total-body OR 1.40 was 9.6% below on the ratio and 23.1% on the log-odds scale (cross-site approximation).

At WNT16, SOX6 and GDF5, loci nominated as shared effectors, colocalisation with knee OA was not supported (PP.H4 ≤0.005 at each; e.g. WNT16 knee-OA z=0.94, p=0.35, PP.H4=0.002; GDF5 knee-OA z=−13.7 but no BMD genome-wide signal within ±500 kb, PP.H4 0.003–0.017). Across all 31 BMD instrument loci, colocalisation was endpoint-split: with knee OA none reached PP.H4≥0.8 (median 0.005, 0/31), whereas with fracture 9/31 did (median 0.153; paired Wilcoxon one-sided p=0.0043, exact two-sided p=0.0086) — the same instrument set shares causal variants with bone but not joint; six loci favoured distinct variants (PP.H3≥0.8) (Fig. S18, S19). One locus fell marginally short: rs11228269 (chr11:68.29 Mb, LRP5; BMD z=−8.65, knee-OA z=−4.06, concordant) gave PP.H4=0.799, 0.0006 below threshold, with residual mass on a BMD-only signal (PP.H1=0.156); we retain the pre-specified 0.8 cut-off and do not count it as colocalisation (1/31 even at a 0.7 threshold; Table S17). The pipeline detected colocalisation where expected: at WNT16, BMD versus fracture gave PP.H4=0.999–1.000. The knee-OA null reflects power — under the MR estimate the expected knee-OA z at the strongest instrument is 0.85, and reaching genome-wide significance would require ~5×10⁷ participants (~42× the present sample).

How much of the BMD→OA effect is skeletal size? Conditioning the BMD→OA effect on height attenuated it by 9.5–33.2% across the two tabulated univariable baselines (13.5–33.2% against the primary Table 2 baselines and 9.5–21.2% against the MVMR-set baselines; primary n_h=100/EUR; Supplementary Tables S9 and S10 and Supplementary Fig. S10), leaving 67–91% independent of height. Because height indexes stature but not bone width or geometry [28], this residual is not identified as bone mass; test F2 addresses the geometry component (§4, Table S28). The positive control confirmed that the established BMD→fracture effect survives the same conditioning: BMD→fracture moved only −1.7% to −3.0% across the same two baselines (−2.8% against the MVMR-set baseline), the control being insensitive to baseline choice. Because conditioning on height captured a genuine, direct OA-raising height effect (θ_H=+0.145 to +0.165, p≈10⁻¹⁶ to 10⁻¹⁹), the attenuation reflects size confounding, not dilution by a null covariate. A parallel BMI-conditioned MVMR (Table S16) left θ_BMD positive and fixed-effect-significant at all five sites (total-body +0.069, p=2.6×10⁻⁸), attenuating −3.7% to 46.0% (the estimate rose at the femoral-neck site) — below the height attenuation at four of five sites; only the non-weight-bearing arms site exceeded it (46.0% vs 13.5%), losing significance under random effects. The MVMR-Egger intercept was nominally significant only at the lumbar site (p=0.027). The positive control moved 0.9%, and θ_BMI→fracture was null (p=0.61) while θ_BMI→knee-OA was strong (β≈+0.69): BMI raises joint- but not bone-endpoint risk, so the height-independent residual is not BMI-mediated. Under BMI conditioning the control's MVMR-Egger intercept was nominally significant at the femoral-neck site (p=2.1×10⁻⁴; total-body p=0.12), a pleiotropy signal reported rather than absorbed.

## 4. Discussion

The architecture is consistent across 13 skeletal phenotypes (Fig. 1, Table 2): on pleiotropy-robust evidence, genetically higher BMD raises knee-OA risk — 2/5 sites significant under random-effects IVW — and protects against fracture.

The opposing-endpoint pattern is old [2, 4]; what is new is its demonstration for one instrument set in one calibrated pipeline. One non-exclusive explanation is partial mediation by skeletal size rather than bone mass [29–30]. Conditioning on height attenuated the effect by 9.5–33.2% (fracture control −1.7% to −3.0%) and height carried a direct OA-raising effect [31–33]; the control is height-transparent, and hip shape was causally implicated in fracture but not hip OA [28].

Zhou et al. [10] reported higher femoral-neck BMD protected against hip OA, whereas our ten-site analysis shows higher BMD raises knee-OA risk. Hartley et al. [12] reported positive BMD effects on both hip (OR 1.28, 1.05–1.57) and knee OA (OR 1.40, 1.20–1.63); Funck-Brentano et al. [34] reported a positive overall effect (all-OA OR 1.14, 1.06–1.22) (Fig. S12). The three studies instrument different phenotypes — femoral-neck BMD, heel-ultrasound eBMD and DXA areal BMD at multiple sites — so the hip evidence is not on one measurement scale; we adjudicate rather than cite selectively, and F1 demands a positive, knee-like hip effect. In the only same-pipeline comparison [12], hip and knee intervals overlap (1.05–1.57 vs 1.20–1.63), the hip log-odds 27% below the knee's (0.247 vs 0.336) — positivity, not equality. Our prospective hip estimate resolves F1: the same 31 instruments give BMD on hip OA IVW β=+0.097 (OR 1.102, 1.068–1.136 per SD), positive and knee-like (hip minus knee z=0.94, p=0.35; Table S27), restoring P6's pre-registered direction.

The directional agreement with Hartley et al. [12] is matched by a magnitude disagreement: their knee estimate is OR 1.40 (1.20–1.63) per SD higher eBMD against 0.058–0.094 here (OR 1.060–1.099 per SD higher BMD). This is a difference of scale and design, not biology: Hartley instrumented heel-ultrasound eBMD (bone size), we DXA areal BMD (area-normalised); Hartley's one-sample 2SLS (overlap 100%) contrasts with ours bounded at ≤2.6%, MRlap moving every estimate a median 20.5% further from the null. Directions concord, magnitudes unresolved.

OA is an osteochondral-unit disorder [35–38]; the set-point we instrument acts upstream of subchondral remodelling. A systemic, polygenic bone-forming set-point readable as areal BMD would protect at the whole-bone endpoint and raise OA risk at the joint — the Systemic Osteogenic Set-point (SOS) (Figs. S21, S23), which we pose as a testable hypothesis distinct in kind from the bone-forming triad (osteophytes, enthesophytes and high bone mass; OR 3.00, 95% CI 1.96–4.58 [39]) and from the broader osteoporosis–osteoarthritis 'bone–cartilage paradox' [3, 4]. Subchondral Wnt/RSPO3–NFATC1 signalling is a candidate effector axis (S30, note vi); anti-resorptives also raise areal BMD, yet bisphosphonate users had a lower knee-replacement rate (RR 0.71, 0.66–0.77) [40], implicating formation rather than density. Because areal BMD fuses mass, size and geometry, F2 asks whether geometry transmits the effect. Its geometry arm is a formal conditional test: a BMD-plus-hip-shape MVMR (59 shape variants with GEFOS BMD effects) leaves the knee-OA conditional BMD effect positive in the better-powered modes, so geometry does not absorb the signal; hip-OA estimates are heterogeneous and sign-flipping, and the volumetric arm remains untested. Unlike the descriptive bone–cartilage paradox [3, 4] or the phenotype-level bone-forming triad [39], SOS is a falsifiable causal assertion, not a restatement of the pattern. Its discriminating kernel is one latent set-point: the same inherited variants must fix opposite signs at whole-bone and joint endpoints, and distinct variant sets or concordant signs would refute it (Fig. S23). The arms result (β=+0.075, p=1.7×10⁻⁶) supports its systemic clause, yet P2's ordering clause was negative (p=0.228) and P2 is partially supported. Supported is the directional pattern, not the osteogenic mechanism: formation-side evidence is indirect and the set-point core is untested (F2 volume arm, F10). The reading generates P1–P6 and F1–F10 (Table S11); survival counts executed tests only, not register completeness (Fig. S23).

**Pre-registered SOS predictions and falsifying tests.** *Predictions (derived).* **P1** opposing endpoints — supported (positive-control OR 1.40 per SD lower BMD; BMD→knee-OA 5/5 positive). **P2** systemic not local-mechanical — partially supported: systemic clause supported (non-weight-bearing arms β = +0.075, p = 1.7×10⁻⁶), ordering clause tested and negative (arms − femoral neck +0.017, one-sided p = 0.228; arms − mean weight-bearing +0.002, p = 0.448). **P3** diffuse polygenic — 0/31 loci colocalise (PP.H4 ≥ 0.8); direction-consistent, non-discriminating. **P4** size a separable partial component — height attenuation 9.5–33.2% and BMI attenuation −3.7% to 46.0% (fracture controls −1.7% to −3.0% and 0.9%); 67–91% residual non-size, non-adiposity; mass versus geometry: geometry excluded in-pipeline, volume arm untested (F2). **P5** constitutional asymmetry — reverse null, direction-consistent but underpowered (0.14–0.26), a hypothesis under test. **P6** site-pathway dominance — restored to its directional form on the prospective hip estimate (F1); cross-site homogeneity (F6, τ²=0) favours systemic over site-specific. *Falsifying tests.* **F1** hip-OA replication (executed: β=+0.097, OR 1.102, 1.068–1.136; positive and knee-like; Table S27); **F2** mass–geometry decomposition (now a formal conditional multivariable MR; knee-OA conditional effect positive, hip-OA heterogeneous; volume arm untested); **F3** bone-versus-cartilage tissue enrichment (gene-set, executed this revision); **F4** eQTL colocalisation (cross-trait: fracture 9/31 vs knee-OA 0/31 at PP.H4≥0.8, p=0.0043; eQTLGen blood p=0.001156; cartilage and skeletal-muscle (bone-proxy) arms (FunGen; Wilson 2025) 0–1/31 loci — Table S26); **F5** osteophyte-stratified outcome — executed as public-data surrogate: site-stratified MR (spine/hand osteophyte-dominant, knee cartilage-dominant) supports formation-side, aligned with osteophyte-rich OA endotypes (S30, note vi); hip-osteophyte GWAS remains unpublished; full in S29. **F6** cross-site homogeneity (executed: Q = 6.92, 4 df, p = 0.14, τ² = 0) and **F7** joint-locus exclusion (executed: GDF5 not a BMD instrument at any site; excluding WNT16 or SOX6 changes no estimate's direction or significance) both passed. **F8** sex stratification (executed: positive in both sexes, §4). **F9** site-gradient ordering (arms−femoral neck +0.017, p=0.228). **F10** ancestry sensitivity (registered cross-ancestry replication unexecuted for want of public files; an ancestry-composition proxy on the mixed outcome in S32 gives 7/7 positive IVW estimates, Q_between=0.29, p=0.592). P4 is the discriminating prediction; F2, its decisive test, remains unexecuted.

A direct MR test of OA→any-fracture was null (β=+0.027, p=0.40) despite a positive genetic correlation with fracture (r_g=+0.119, p=3.1×10⁻⁹). The BMD-mediated component runs oppositely (higher BMD raises knee-OA risk yet protects against fracture). Conditioning the OA→fracture effect on BMD_total in a multivariable MR (69 OA instruments; Table S31) left it unchanged (β=+0.034, 95% CI −0.029 to 0.097), confirming BMD is not a confounder of the null and that the positive r_g is non-causal; Egger-MVMR showed no horizontal pleiotropy (intercept p=0.20). Higher BMD in OA has not meant fewer fractures, the residual association attributed to falls [41]. This is among the first direct MR tests of OA→any-fracture, excluding large but not moderate effects (power 0.14). Our BMD→OA finding is consistent with Alhassan et al. [7]; the reverse-null reading of P5 is therefore provisional (see the companion OA↔bone funnel, Supplementary Table S30, note vi).

**Strengths and limitations.** Strengths: ten BMD sites with two OA definitions, an internal BMD–fracture positive control calibrated against GEFOS (OR 1.40–1.42 vs 1.55), dual FDR/Bonferroni correction and directional concordance across IVW, MR-Egger, weighted median and MR-PRESSO. Limitations: the European-only BMD instruments and LDSC reference against an 87.31%-European knee-OA outcome mean forward estimates generalise most directly to European-ancestry populations (Table S14). Five of ten prespecified sites entered MR (Table S12); the deferred five showed comparable signal, so the subset reflects compute availability not biology. Areal BMD conflates size with mass [29–30, 42]; F2 excludes the geometry axis in-pipeline but the volume arm stays untested. The sex-stratified F8 MR was positive in both sexes (KneeOA Female β=+0.090, Male β=+0.127; AllOA Female β=+0.097, Male β=+0.066; F-versus-M p=0.029 and 0.006), with no reversal — direction consistency, not a resolved magnitude difference. All pairs are heterogeneous; reverse and OA→fracture estimates are null but weakly powered (0.14–0.26), and residual overlap is bounded at ≤2.6% and cannot be removed with summary data. 

**Conclusion.** Genetically higher BMD protects against fracture and, on the pleiotropy-robust evidence, raises knee-OA risk (2/5 sites significant under random-effects IVW). The effect concerns the clinically measured areal-BMD phenotype at the five analysed sites; 67–91% of it is height- and BMI-independent, and not transmitted by bone geometry (F2), albeit on weak conditional instruments (conditional F below 4 in all but one mode), so this is direction-indicative rather than magnitude-precise (Table S28, Part B); and reverse causation and an OA→fracture pathway were not supported in MR (the reverse estimates exclude only |β|≥0.09; the positive knee-OA × fracture r_g is genetic correlation, not causation). The joint-endpoint effect is modest in individual terms (OR 1.060–1.099 per SD higher BMD, 95% CI 1.025–1.135 across sites) but attaches to a ubiquitous, modifiable exposure; these estimates do not support changing prescribing practice; one trial of a BMD-raising drug reported structural benefit [43], but strontium inflates DXA areal BMD as a measurement artefact and is restricted for cardiovascular risk, and a second trial (TREASURE) found no osteophyte progression despite greater BMD gain [44]; three years is not lifelong exposure.

---

## Figure legends

**Figure 1.** Heatmap of LDSC genetic correlations (r_g) across 13 skeletal phenotypes. Cell text shows r_g; asterisks denote BH-FDR significance. Black lines separate the three phenotype groups. The colour scale is diverging and centred on zero.


**Figure 2.** Forest plot of MR estimates on the N-corrected effect scale. Upper panel: BMD→fracture positive control (OR ≈1.40–1.42 vs GEFOS 1.55). Middle panel: five BMD exposures against knee OA (all positively causal). Lower panel: reverse/null pairs (all non-significant, p=0.19–0.40). Points are IVW fixed-effects estimates with 95% CI; blue = protective, red = risk-increasing, grey = null.


**Figure 3.** Conceptual causal architecture of the OA–BMD–fracture triangle. Solid arrows denote MR causal effects (BMD protects against fracture; BMD raises knee OA under pleiotropy-robust estimators); dashed grey arrows denote non-causal relationships.


**Figure 4.** Per-SNP scatter plots for four headline pairs (IVW, MR-Egger, weighted-median fits). (a,b) BMD→fracture (negative slopes, positive control). (c) Total-body BMD→knee OA (positive). (d) Knee OA→total-body BMD (flat/null). Red open circles mark MR-PRESSO outliers.


**Figure 5.** MRlap correction of the ten MR pairs (accounting for participant overlap while correcting weak-instrument bias and winner's curse); every pair moves away from the null (median +20.5%), showing that the towards-null components dominate the away-from-null overlap component, and the observed-versus-corrected difference is significant for 9/10 pairs (c).


**Supplementary Figure S1.** Per-SNP scatter plots for the six remaining MR pairs. **S2.** Leave-one-out stability of all ten pairs. **S3.** External validation of the reconstructed BMD scale against GEFOS. **S4.** Study design and MR flow (STROBE-MR). **S5.** Sensitivity of OA_knee h² to representative N (r_g is N-invariant). **S6.** Fracture sample-size decision. **S7.** Pleiotropy and heterogeneity diagnostics. **S8.** Colocalisation of BMD with knee OA and its power. **S9.** Instrument selection and attrition flow (STROBE-MR), tallied from the deposited clump and harmonised instrument tables. Instruments are lost through palindromic or strand-ambiguous status, or absence from the outcome GWAS, and not through effect-size selection (§2). **S10.** Multivariable MR decomposition of the BMD→OA effect conditional on height: (a) univariable vs MVMR estimates; (b) attenuation percentages with the BMD→fracture positive control; the height-independent residual is not resolved into bone mass versus skeletal geometry. **S11.** Multi-estimator robustness grid for the ten MR pairs (IVW fixed/random, MR-Egger, weighted median, PRESSO outlier-removed, MRlap-corrected, MVMR). **S12.** Quantitative comparison with published MR estimates (Hartley 2022; Funck-Brentano 2019; Zhou 2025) on their native scales. **Supplementary Figure S13.** Forest plot of the selected genetic correlations: any-fracture versus ten BMD sites (all negative); knee OA and self-reported OA versus ten BMD sites (all positive); the two OA definitions with each other (r_g=+0.84), knee OA with fracture (r_g=+0.119, q=3.3×10⁻⁹) and self-reported OA with fracture (non-significant, +0.056). Error bars are 95% confidence intervals from block jackknife.
 **Supplementary Figure S14.** SNP-based heritability (observed scale) of the 13 skeletal phenotypes; diamond markers show liability-scale h² for the three binary traits. The liability-scale diamonds use the GWAS case-fraction K (OA_knee h²_liab = 0.106 at K = 0.096; OA_selfrep 0.073 at K = 0.083; fracture 0.062 at K = 0.143); Table 1 reports the population-prevalence convention (OA_knee 0.094 at K = 0.12). Both are valid conversions of the same observed h² under different K assumptions.
 **Supplementary Figure S15.** MR-PRESSO funnel plots for the four headline pairs; red open circles are outlier instruments.
 **Supplementary Figure S16.** Statistical power of the two-sample MR design. Filled blue circles are the seven well-powered pairs; red open circles the three null pairs (power 0.14–0.26, excluding only |β|≥0.09).
 **Supplementary Figure S17.** Sample-overlap diagnostics (ρ_GE) for the ten MR pairs; ρ_GE shares the sign of the MR effect in all ten pairs.
 **Supplementary Figure S18.** Endpoint-specific colocalisation at the 31 BMD instrument loci. Split horizontal bars show the posterior probability of a shared causal variant (PP.H4) for BMD↔knee OA (blue) versus BMD↔fracture (orange) at each of the 31 strict-clumped total-body BMD instruments. No locus reached PP.H4≥0.8 with knee OA (0/31; median 0.005), whereas 9/31 did with fracture (median 0.153). The dashed line marks the PP.H4=0.8 colocalisation threshold; one locus (rs11228269, LRP5) fell 0.0006 short (PP.H4=0.799) and is not counted. A paired Wilcoxon signed-rank test (exact) gave one-sided p=0.0043 and two-sided p=0.0086, confirming the instrument set shares causal variants with bone but not joint.
 **Supplementary Figure S19.** Colocalisation at three nominated effector loci across five BMD sites. Bars show PP.H4 for BMD↔knee OA (blue) and BMD↔fracture (orange) at WNT16, SOX6 and GDF5, each tested against the five BMD sites that entered MR (total body, lumbar spine, femoral neck, legs, arms). BMD and fracture colocalise strongly at WNT16 (PP.H4 0.999–1.000) and moderately at SOX6 (0.40–0.42 at four sites), whereas BMD and knee OA do not (PP.H4 ≤0.0017 and ≤0.0012). The exception is GDF5 at the lumbar-spine site: BMD↔knee OA shows PP.H4=0.725 while BMD↔fracture shows 0.036 — the single locus–site where the endpoint-specific pattern reverses; the remaining GDF5 sites show no colocalisation with either outcome (PP.H4 ≤0.017).
 **Supplementary Figure S20.** Multi-estimator robustness of the ten MR pairs. Direction–sign concordance across five estimators (IVW fixed-effects, MR-Egger, weighted median, leave-one-out, MR-PRESSO outlier-removed), grouped as BMD→knee OA (A, five pairs), BMD→fracture (B, two pairs) and reverse / OA→outcome (C, three pairs). Cell arrows give each estimator's direction (↑ risk-increasing for the exposure→outcome reading, ↓ protective); an asterisk marks p<0.05 for IVW and MR-Egger; the MR-PRESSO column tags the heterogeneity call (het = significant distortion, clean = none). Colour: blue = direction matching the IVW reading with p<0.05; light blue = matching direction but not significant; vermillion = direction opposite the IVW reading; amber = leave-one-out unstable or heterogeneity detected. All ten pairs retain a consistent direction across every estimator; the three reverse/null pairs (C) stay null and heterogeneous under all five.
 **Supplementary Figure S21.** The Systemic Osteogenic Set-point (SOS) hypothesis. A polygenic osteogenic set-point, operationalised as areal BMD, protects at the whole-bone endpoint (lower fracture risk; OR 1.40 per SD lower BMD, IVW β=−0.337) yet raises knee-OA risk at the joint (non-weight-bearing arms β=+0.075, p=1.7×10⁻⁶). The joint effect is read as subchondral remodelling (formation-side) rather than density: the bone-forming triad (osteophytes, enthesophytes, high bone mass) carries an OA OR of 3.00 [1.96–4.58] (Hardcastle 2014 [39]), and bisphosphonate users had a lower knee-replacement rate (RR 0.71 [0.66–0.77]) (Cui 2022 [40]). SOS is a testable hypothesis (F1–F10), not a demonstrated mechanism; its geometry arm is now a formal conditional BMD-plus-hip-shape multivariable MR (knee-OA conditional BMD effect positive in all estimable modes, hip-OA heterogeneous), and is partly informed externally by purely geometric hip-shape modes causing fracture (HSM2 OR 1.27 [1.12–1.44]) but not hip OA (Faber 2025 [28]).
 **Graphical Abstract.** The OA–BMD–fracture triangle: higher BMD protects against fracture yet raises knee-OA risk.

**Supplementary Figure S22.** Latent structure of the 13-phenotype genetic-correlation matrix. PCA of the 13×13 LDSC r_g matrix; points are phenotype loadings on PC1 (bone-mass/fracture axis) versus PC2 (OA axis), coloured by group. BMD sites (green) separate from OA (orange) and fracture (blue); the two axes are orthogonal, so BMD and OA occupy distinct genetic dimensions.

**Supplementary Figure S23.** SOS hypothesis–prediction–falsification framework. The central Systemic Osteogenic Set-point (SOS) hypothesis — one polygenic osteogenic set-point, operationalised as areal BMD, which protects the whole bone yet raises knee-OA risk at the joint — generates six derived predictions (P1–P6, middle row) and is tested by ten falsifying tests (F1–F10, lower two rows; the full register is in Table S11). Connectors link each prediction to the test(s) that target it; test boxes are colour-coded by outcome (green = supported/pass, amber = partial or with an unexecuted arm, grey = null or non-discriminating, red = not executed). All statistics are taken verbatim from §4 and Tables S25–S31. The verdict counts executed tests only: F2's volumetric arm is unexecuted, F10 has an executed ancestry-composition robustness proxy but no registered cross-ancestry replication (S32), and F5 rests on a public-data surrogate, so "9/10" denotes survival among executed tests, not completeness of the register; the osteogenic set-point mechanism is supported indirectly, not demonstrated (§4). Theory boundary. The discriminating kernel is not a label for the observed pattern; it is non-trivial because BMD, OA and fracture occupy distinct genetic dimensions (S22: the PC1 bone-mass/fracture axis is orthogonal to the PC2 OA axis), so a single set-point cannot be reduced to a generic bone–OA correlation. One latent set-point, read through a shared set of inherited variants, must fix opposite signs at whole-bone and joint endpoints, with sign set by the set-point rather than by endpoint-specific pathways. The descriptive bone–cartilage paradox [3, 4] and the phenotype-level bone-forming triad [39] are rivals because neither requires one variant set nor sign opposition; what SOS adds is the single-set-point discriminator — one variant set fixing sign-opposing causal effects — rather than a restatement of the descriptive pattern. The kernel is falsifiable: endpoint-specific instrument sets, loss of sign opposition under volumetric BMD instruments, or ancestry-specific reversal would refute it. F4 (cross-trait colocalisation), F6 (cross-site homogeneity) and F7 (joint-locus exclusion) probe it only indirectly; the F2 volumetric arm (true vBMD instruments) and the F10 cross-ancestry test are the priority subsequent tests that could close the two remaining gaps. The F10 gap is attenuated, not removed: the outcome GWAS [14] performed GWAS meta-analyses within four ancestry groups (East Asian, African, South Asian, Hispanic) and detected no ancestry-specific study-wide significant associations, indicating the osteoarthritis genetic architecture is not ancestry-specific at the GWAS level and supporting the cross-ancestry transportability of the European-instrument estimates, while direct non-European replication stays blocked by the absence of non-European BMD exposure summary statistics; the ancestry-composition robustness proxy in S32 is a within-outcome consistency check and does not substitute for that replication.

**Supplementary Figure S24.** Site-stratified genetic-correlation panel (data in Table S30). **A**, BMD×OA genetic correlation (r_g) by joint site (S30a): hip 0.131 (95% CI 0.057–0.205), knee 0.152 (0.088–0.217), spine 0.196 (0.104–0.287), all P < 1×10⁻³; hand 0.071 (−0.039 to 0.181, P = 0.20, n.s.). **B**, OA site-to-site r_g with knee as reference (S30b): knee–hip 0.569, knee–spine 0.620, knee–hand 0.591, all P < 10⁻⁵⁰, showing OA sites share most of their genetic architecture. **C**, SNP heritability on the observed scale (S30c): total-body BMD h² = 0.348 sits an order of magnitude above the OA sites (0.009–0.033). The r_g values quantify polygenic sharing and are distinct from the causal MR estimates in Table 2 and S29; the site gradient of r_g does not track the MR gradient in S29.

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

§§ Fracture N adopted as 264,973 (37,857 cases / 227,116 controls; GEFOS discovery, Trajanoska et al. 2018), identified by the measured effective sample size N·K·(1−K)=31,605, which matches the GEFOS theoretical value of 32,449 within 2.6%. § OA_selfrep N confirmed at 462,933 (38,472 cases / 424,461 controls; IEU OpenGWAS `ukb-b-14486`). ‡ BMD h² was computed at a rounded representative N of 32,000 (the per-site GWAS N is 31,873); because h² ∝ 1/N the rounding makes the reported h² 0.40% lower than at the exact N, an effect far smaller than the reported SEs. † OA_knee h² is reported at the representative proband-subset N=1,208,397; h² is N-dependent, whereas r_g and its SE are N-invariant because r_g is a ratio of covariances estimated on the same SNP set. The LDSC regression uses the per-SNP N column; all r_g values are N-scale-invariant. ‡‡ BMD N verified per site (Qian et al. 2025; range 30,449–32,017); the unified N=31,873 is used. **Liability scale (Lee et al. 2011) [45]:** OA_knee h²_liab=0.094 (K=0.12); OA_selfrep h²_liab=0.073 (K=0.083); fracture h²_liab=0.062 (K=0.143). BMD is quantitative and not liability-converted. LDSC intercepts are near 1.0 (1.01–1.07) for 12 phenotypes; the fracture intercept is 0.96.

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
| OA knee → fracture | 69 | +0.027 | 0.032 | 0.395 | 0.542 | 1.028 (0.965–1.094) | +0.031 (0.045) | 0.490 | 0.476 | +0.019 (−0.112 to +0.126) | 49% | 0 |
| OA knee → BMD total | 69 | +0.031 | 0.032 | 0.326 | 0.481 | — | +0.036 (0.045) | 0.420 | 0.376 | +0.026 (−0.032 to +0.073) | 49% | 0 |
| OA knee → BMD femoral neck | 69 | +0.042 | 0.032 | 0.190 | 0.272 | — | +0.048 (0.038) | 0.206 | 0.169 | +0.023 (−0.063 to +0.122) | 30% | 0 |

Univariate instrument strength (F, R²_exposure) and the Steiger directionality test for every pair are in Table S13. β is the IVW estimate per 1 SD higher exposure; for the binary outcomes β is on the log-odds scale and OR = exp(|β|) with CI exp(|β| ± 1.96 SE), quoted per 1 SD *higher* BMD for the BMD→OA pairs and per 1 SD *lower* BMD for the BMD→fracture pairs, where higher BMD is protective; no OR is given for the reverse pairs, whose exposure is a binary liability rather than a 1-SD continuous exposure. MR-Egger p is a two-sided z-test on the slope. Weighted median is Bowden's estimator with weights (β_exposure/SE_outcome)², recomputed from the deposited harmonised instrument tables with 10,000-replicate bootstrap CIs (code/129); the IVW, MR-Egger, Cochran Q and I² columns reproduce from the same deposit to machine precision, and seven weighted-median values previously tabulated did not, so that column was recomputed from the deposited harmonised instrument tables (code/129) and corrected. Forward BMD-exposure betas are per-SD on the verified N=31,873 scale (per-instrument closed-form SE check, within 0.1%); no rescaling was applied. Reverse MR is on the same verified SD scale. p-values, Egger-intercept p, I² and LOO sign flips are N-scale-invariant. Strict PLINK-equivalent clumping (K = strict instrument count). For BMD femoral neck → fracture the MR-Egger intercept is significant (p=0.018); the MR-Egger slope (−0.343) and weighted median (−0.365) are the pleiotropy-robust estimates. The random-effects p column uses the multiplicative convention SE_random = SE_fixed × √[max(1, Q/df)] (TwoSampleMR); an additive DerSimonian–Laird τ² convention would additionally retain the arms exposure under random effects (p=0.033, §3).

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

MDE = (1.96 + 0.84) × SE is the minimum |β| detectable at 80% power (two-sided α=0.05); achieved power is for the observed point estimate. |β|/SE is N-scale-invariant. K = strict instrument count. Achieved power refers to fixed-effect IVW; under the primary random-effects IVW only two of five forward BMD→knee-OA estimates remain significant (Table S17-b; §3).

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

| GWAS (phenotype) | Source | N (effective) | UK Biobank participants | Ancestry composition | LDSC reference |
| --- | --- | --- | --- | --- | --- |
| BMD (10 regional sites) | Qian et al. 2025 [16] | 31,873 each | 31,873 (100%) | 100% European (UK Biobank) | 1000G Phase 3 EUR |
| OA_knee (radiographic) | GO Consortium 2025 [14] | ≈1,208,397 | includes UKB* | Multi-ancestry: EUR 87.31%, EAS 7.09%, AFR 3.08%, SAS 1.09%, HIS 0.91%, ADM 0.53% | 1000G Phase 3 EUR |
| OA_selfrep (self-reported) | IEU OpenGWAS `ukb-b-14486` [15] | 462,933 | 462,933 (100%) | European (UK Biobank) | 1000G Phase 3 EUR |
| Fracture (any) | GEFOS 2018 [17] | 264,973 | 22,330 (1,832 cases + 20,498 controls) | European | 1000G Phase 3 EUR |

The BMD exposure instruments are European-only, whereas the knee-OA outcome is 87.31% European (multi-ancestry meta-analysis); LDSC uses a European-only reference, so cross-population linkage-disequilibrium differences between the BMD instruments and the multi-ancestry knee-OA outcome can bias r_g and the MR estimates. The OA_knee percentages are taken from the GWAS Catalog entry and the originating Nature report of the multi-ancestry meta-analysis.

*OA_knee is multi-ancestry with UKB as a recruitment country, but the GO Consortium did not report the UKB count separately, so overlap with BMD is an upper bound of 31,873 / 1,208,397 = 2.6%. Pairwise UKB overlap and Burgess bias (see §2): BMD is 100% UKB (31,873); OA_selfrep 100% UKB (462,933), overlap with BMD = 31,873 (6.9% of OA_selfrep); fracture (GEFOS 2018) includes 22,330 UKB participants (1,832 cases + 20,498 controls; 8.4% of its 264,973), so BMD↔fracture overlap = 22,330 (8.4% of fracture, 70% of BMD). Burgess bias = overlap fraction × r_g: BMD→OA_knee ≤ ~0.3–0.7% (overlap ≤2.6% × r_g 0.12–0.25); BMD→OA_selfrep ~0.8–1.7% (6.9% × 0.12–0.25); BMD→fracture ~3.4–4.5% on the conservative 8.4%-of-fracture basis (up to ~12% on the 70%-of-BMD basis); all act in the direction of r_g. MRlap (§2) accounts for overlap directly. *MRlap input and internally estimated overlap.* No overlap proportion was supplied to MRlap (`code/run_mrlap.R` L83–85 passes only `ld`, `hm3`, `do_pruning`, `MR_pruning_LD` and `MR_threshold`); it ran on its default setting. MRlap does not return a scalar sample-overlap fraction: under the default, participant overlap is accommodated through the cross-trait LDSC intercept λ (`int_crosstrait`), estimated internally from the summary statistics. The deposited per-pair values (`code/output/mrlap_overlap_results.csv`) are BMD_total→fracture −0.0361, BMD_femur_neck→fracture −0.0245, BMD_total→OA_knee −0.0051, BMD_femur_neck→OA_knee −0.0007, BMD_lumbar→OA_knee +0.0010, BMD_arms→OA_knee −0.0101, BMD_leg→OA_knee −0.0056, OA_knee→fracture +0.0153, OA_knee→BMD_total −0.0051 and OA_knee→BMD_femur_neck −0.0007 (range −0.0361 to +0.0153; |λ| ≤ 0.037 in all ten pairs, §3). These are intercepts, not overlap fractions: the fractions we report are the tabulated UKB counts above, derived independently, and none was fed to MRlap. *Sensitivity to the default.* All ten |λ| are small and of mixed sign, so the internally estimated overlap signal is minor relative to the reported effects. Because MRlap's default does not take a user-specified overlap fraction, we assess sensitivity against the independently tabulated counts: the Burgess bounds above (≤2.6–8.4% overlap; biases ≤0.7–4.5%) bound the away-from-null component, and every corrected estimate retained its direction (10/10; median +20.5% away from the null; `code/output/mrlap_direction_verification.txt`), so no headline conclusion depends on the default.


## Supplementary Table S15. Genetic correlations by phenotype block

| Block          | n pairs | r_g range        | r_g mean | Significant (q<0.05) |
| -------------- | ------- | ---------------- | -------- | -------------------- |
| OA × OA        | 1       | +0.838           | +0.838   | 1/1                  |
| OA × BMD       | 20      | +0.117 to +0.250 | +0.169   | 20/20                |
| Fracture × BMD | 10      | −0.535 to −0.401 | −0.460   | 10/10                |
| BMD × BMD      | 45      | +0.689 to +0.960 | +0.812   | 45/45                |
| OA × Fracture  | 2       | +0.056 to +0.119 | +0.087   | 1/2

## Supplementary Table S31. Multivariable MR of OA_knee and BMD_total on any-fracture (conditioning on BMD)

Exposures [OA_knee, BMD_total], outcome any-fracture, 69 OA instruments (code/264). IVW-MVMR and Egger-MVMR. The OA→fracture effect is unchanged after conditioning on BMD (θ_OA|BMD β=+0.034, 95% CI −0.029 to 0.097; OR 1.03), confirming BMD is not a confounder of the null; Egger-MVMR intercept p=0.20 (no horizontal pleiotropy). BMD retains a protective effect (θ_BMD|OA β=−0.223, 95% CI −0.394 to −0.051; OR 0.80 per SD higher BMD), matching the univariate positive control.

| Effect | IVW β | 95% CI | OR (95% CI) | Egger β | Egger 95% CI | Egger OR (95% CI) |
| --- | --- | --- | --- | --- | --- | --- |
| θ_OA\|BMD (OA→fracture, BMD-conditioned) | +0.034 | −0.029 to 0.097 | 1.035 (0.971–1.102) | +0.039 | −0.024 to 0.103 | 1.040 (0.976–1.108) |
| θ_BMD\|OA (BMD→fracture, OA-conditioned) | −0.223 | −0.394 to −0.051 | 0.800 (0.674–0.950) | −0.234 | −0.407 to −0.062 | 0.791 (0.666–0.940) |                  |


---

## Supplementary Table S25. Bone-versus-cartilage tissue enrichment — competitive gene-set test (F3, this revision)

Tissue specificity of the 13-phenotype genetic architecture was tested by a competitive gene-set enrichment (MAGMA-style): each GWAS SNP was assigned to the nearest gene body within ±50 kb; gene-level association = the minimum p across assigned SNPs; each GO gene set was tested against the remaining genes by a one-sided Mann–Whitney U on −log10 p. Positive enrichment Δ means the set's genes are more associated than background. Bone-side sets: bone development, ossification, skeletal-system development. Cartilage-side sets: cartilage development, chondrocyte differentiation, extracellular-matrix organisation, skeletal-joint morphogenesis. Source GWAS: OA knee, total-body BMD, any-fracture (sumstats in OA_FLAGSHIP/07_reference/ldsc; code/171).

| Trait | Gene set | Tissue | n (set) | mean −log10p set | mean −log10p bg | enrichment Δ | MWW p (1-sided) |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| OA_knee | bone_development | bone | 83 | 3.415 | 2.610 | 0.805 | 2.21e-04 |
| OA_knee | cartilage_development | cartilage | 72 | 5.036 | 2.606 | 2.430 | 1.45e-08 |
| OA_knee | ossification | bone | 71 | 3.741 | 2.609 | 1.131 | 2.78e-03 |
| OA_knee | chondrocyte_differentiation | cartilage | 81 | 4.932 | 2.605 | 2.326 | 1.58e-08 |
| OA_knee | bone_morphogenesis | other | 9 | 2.084 | 2.613 | -0.528 | 6.07e-01 |
| OA_knee | extracellular_matrix_organization | cartilage | 37 | 3.084 | 2.612 | 0.472 | 3.27e-02 |
| OA_knee | skeletal_system_development | bone | 42 | 5.354 | 2.608 | 2.746 | 9.96e-07 |
| OA_knee | skeletal_joint_morphogenesis | cartilage | 82 | 4.057 | 2.608 | 1.449 | 4.58e-10 |
| OA_knee | bone_remodeling | other | 1 | 4.853 | 2.612 | 2.241 | NA |
| BMD_total | bone_development | bone | 83 | 2.121 | 1.671 | 0.450 | 3.07e-03 |
| BMD_total | cartilage_development | cartilage | 72 | 2.198 | 1.671 | 0.526 | 1.35e-02 |
| BMD_total | ossification | bone | 71 | 2.285 | 1.671 | 0.614 | 1.00e-01 |
| BMD_total | chondrocyte_differentiation | cartilage | 81 | 2.277 | 1.671 | 0.606 | 2.11e-02 |
| BMD_total | bone_morphogenesis | other | 9 | 2.124 | 1.673 | 0.452 | 1.74e-02 |
| BMD_total | extracellular_matrix_organization | cartilage | 37 | 1.903 | 1.673 | 0.230 | 1.52e-02 |
| BMD_total | skeletal_system_development | bone | 42 | 3.403 | 1.670 | 1.733 | 1.26e-04 |
| BMD_total | skeletal_joint_morphogenesis | cartilage | 82 | 2.530 | 1.670 | 0.860 | 3.16e-07 |
| BMD_total | bone_remodeling | other | 1 | 3.645 | 1.673 | 1.972 | NA |
| fracture | bone_development | bone | 83 | 1.599 | 1.455 | 0.144 | 3.14e-02 |
| fracture | cartilage_development | cartilage | 72 | 1.895 | 1.454 | 0.441 | 2.16e-04 |
| fracture | ossification | bone | 71 | 1.538 | 1.455 | 0.083 | 3.52e-01 |
| fracture | chondrocyte_differentiation | cartilage | 81 | 1.812 | 1.454 | 0.357 | 8.38e-04 |
| fracture | bone_morphogenesis | other | 9 | 1.590 | 1.455 | 0.135 | 1.81e-01 |
| fracture | extracellular_matrix_organization | cartilage | 37 | 1.696 | 1.455 | 0.241 | 5.33e-04 |
| fracture | skeletal_system_development | bone | 42 | 1.825 | 1.455 | 0.370 | 5.42e-04 |
| fracture | skeletal_joint_morphogenesis | cartilage | 82 | 1.902 | 1.454 | 0.448 | 1.36e-06 |
| fracture | bone_remodeling | other | 1 | 1.445 | 1.455 | -0.010 | NA |

## Supplementary Table S26. eQTL colocalisation — tissue-literal (cartilage and skeletal-muscle bone-proxy) arms of F4 (this revision)

At each of the 31 total-body BMD instrument loci (±500 kb), tissue eQTL summary statistics were colocalised with OA knee, total-body BMD and any-fracture by coloc.abf (priors p1=p2=1e-4, p12=5e-6; eQTL variance prior SD=0.15). Two tissues were tested as the genetic counterparts of the two components of the osteochondral unit: (i) cartilage, from FunGen low-grade osteoarthritis cartilage (Southam et al.; n≈200 donors); and (ii) skeletal muscle, used as a bone-proximal proxy from the Wilson et al. (AJHG 2025) meta-analysis (GTEx+FUSION, n=1,002). PP.H4 = posterior probability of a shared causal variant. Neither tissue reached PP.H4≥0.8 for BMD or fracture, and only one muscle (bone-proxy) locus reached it for OA knee (PP.H4=0.887, a single borderline signal). The BMD instruments therefore do not colocalise with cartilage or skeletal-muscle regulatory variation at genome-wide-significant density; the strongest cartilage signal was with total-body BMD (top locus PP.H4=0.475). This is consistent with a formation-side (bone/subchondral) rather than cartilage-matrix or muscle-mediated interpretation of the BMD→OA effect, and parallels the cross-trait arm in which the same instruments colocalise with fracture (bone) but not OA (joint). Per-locus PP.H4 values are in f4_literal_cartilage_*.csv and f4_literal_muscle_*.csv (code/output/genetics/).

| GWAS trait | Cartilage: n loci | Cartilage median PP.H4 | Cartilage n PP.H4≥0.8 | Muscle: n loci | Muscle median PP.H4 | Muscle n PP.H4≥0.8 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| OA knee | 31 | 0.009 | 0 | 31 | 0.002 | 1 |
| BMD total | 31 | 0.061 | 0 | 31 | 0.000 | 0 |
| fracture | 31 | 0.002 | 0 | 31 | 0.004 | 0 |

---

## Supplementary Table S27. Prospective hip-osteoarthritis replication (falsifying test F1): two-sample MR of total-body BMD on hip OA, and the conjunctive adjudication against the knee-OA estimate

F1 was pre-registered as the prospective replication of the BMD->OA finding at the hip, with a conjunctive criterion: the hip effect must be positive AND of knee-like magnitude. The hip outcome is the Hatzikotoulas 2025 mixed-ancestry hip-OA GWAS (GCST90566798; 95,395 European cases, 1,011,691 European controls plus East Asian, African, South Asian, Hispanic and admixed groups; maximum per-variant 97,328 cases and 1,055,375 controls), and the knee outcome is the same study's knee-OA GWAS (GCST90566800; maximum per-variant 172,256 cases and 1,144,248 controls) — the file from which the manuscript's own knee-OA estimate was derived, so that the two rows below are produced by an identical pipeline on an identical instrument set (31 total-body BMD instruments, Qian 2025).

| Outcome | n_snp | IVW beta (fixed) | SE | p | IVW p (random) | OR per SD (95% CI) | MR-Egger beta (SE) | Egger p | Egger intercept p | Weighted median (95% CI) | Cochran Q | I2 | LOO sign flips |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| hip OA (GCST90566798) | 31 | +0.097 | 0.016 | 6.6e-10 | 5.8e-04 | 1.102 (1.068-1.136) | +0.097 (0.016) | 6.1e-10 | 0.29 | +0.096 (+0.047 to +0.181) | 96.6 (30 df) | 69% | 0 |
| knee OA (GCST90566800) | 31 | +0.078 | 0.012 | 3.3e-10 | 6.0e-03 | 1.081 (1.055-1.108) | +0.078 (0.012) | 3.4e-10 | 0.63 | +0.066 (+0.027 to +0.174) | 156.6 (30 df) | 81% | 0 |

**Adjudication.** Difference in IVW log-odds (hip minus knee) = +0.0188; SE under the conservative assumption of independent outcomes = 0.0200; z = +0.94, two-sided p = 0.35. Ratio of hip to knee log-odds = 1.24. The 95% intervals overlap from 1.0684 to 1.1079. Instrument strength: mean F = 807, minimum F = 168; SNP-level R2 for the exposure 29.4 versus 15.4 for the outcome, Steiger-type directionality +0.312. Leaving out the 8 variants with absolute studentised residual above 2 leaves the estimate positive and strengthened (beta = +0.1300, SE 0.0180, p = 5.0e-13); leave-one-out produced no sign reversal (range +0.0867 to +0.1113).

**Verification against the published knee-OA row.** Re-deriving the knee-OA estimate from the raw GWAS Catalog file reproduces the manuscript's Table 2 values: IVW beta +0.078 (this table +0.078), OR 1.081 (1.055-1.108), I2 81% (80.8%) and weighted median +0.065 (+0.027 to +0.177) (+0.066, +0.027 to +0.174). The MR-Egger intercept is non-significant in both (this table p = 0.63; manuscript p = 0.810).

**Decision.** F1 passes its conjunctive criterion. The prospective same-pipeline hip estimate is positive (p = 6.6e-10) and statistically indistinguishable from the knee estimate (z = 0.94, p = 0.35); the hip estimate lies 24% above the knee estimate, within the magnitude band spanned by the only prior same-pipeline comparison (Hartley et al., hip 27% below knee). P6 is therefore restored to its pre-registered directional form.

**Per-instrument values (hip OA outcome; the 31 total-body BMD instruments).** beta_e and se_e are the total-body BMD exposure effects per allele; hip_lnOR is the hip-OA log-odds on the BMD-increasing allele with its standard error.

| SNP | beta_e (BMD) | se_e | hip ln(OR) | se | hip OR | N (cases+controls) |
| --- | --- | --- | --- | --- | --- | --- |
| rs10048146 | -0.0601 | 0.0101 | -0.0028 | 0.0075 | 0.997 | 793,934 |
| rs10808100 | -0.0876 | 0.0084 | -0.0165 | 0.0068 | 0.984 | 764,809 |
| rs10931982 | -0.0547 | 0.0095 | -0.0218 | 0.0073 | 0.978 | 1,068,445 |
| rs11228269 | -0.0770 | 0.0089 | -0.0198 | 0.0069 | 0.980 | 779,271 |
| rs12568930 | +0.1014 | 0.0103 | +0.0160 | 0.0083 | 1.016 | 794,837 |
| rs12600549 | +0.0457 | 0.0081 | +0.0148 | 0.0054 | 1.015 | 1,151,848 |
| rs1262476 | -0.0525 | 0.0091 | -0.0141 | 0.0064 | 0.986 | 1,107,952 |
| rs17150506 | -0.0696 | 0.0096 | -0.0051 | 0.0064 | 0.995 | 1,149,812 |
| rs17212655 | -0.0761 | 0.0106 | +0.0135 | 0.0069 | 1.014 | 1,151,490 |
| rs2414095 | -0.0476 | 0.0083 | +0.0042 | 0.0066 | 1.004 | 769,928 |
| rs2566755 | +0.0826 | 0.0096 | +0.0152 | 0.0074 | 1.015 | 798,833 |
| rs3020333 | +0.0759 | 0.0080 | +0.0137 | 0.0053 | 1.014 | 1,135,963 |
| rs3755955 | -0.0958 | 0.0111 | -0.0085 | 0.0073 | 0.992 | 1,149,668 |
| rs3779381 | +0.1652 | 0.0090 | +0.0158 | 0.0061 | 1.016 | 1,150,501 |
| rs3929244 | -0.0563 | 0.0090 | +0.0084 | 0.0061 | 1.008 | 1,140,961 |
| rs4846607 | -0.0496 | 0.0081 | -0.0016 | 0.0053 | 0.998 | 1,137,241 |
| rs4974056 | +0.0485 | 0.0079 | +0.0088 | 0.0059 | 1.009 | 760,772 |
| rs6040061 | -0.0487 | 0.0079 | -0.0114 | 0.0052 | 0.989 | 1,143,852 |
| rs6679677 | +0.0748 | 0.0134 | -0.0200 | 0.0093 | 0.980 | 791,684 |
| rs6846509 | +0.0761 | 0.0088 | +0.0010 | 0.0061 | 1.001 | 1,117,295 |
| rs6956946 | -0.0474 | 0.0083 | -0.0229 | 0.0068 | 0.977 | 757,349 |
| rs7010267 | +0.0940 | 0.0080 | -0.0070 | 0.0053 | 0.993 | 1,150,501 |
| rs7277076 | +0.0457 | 0.0080 | +0.0074 | 0.0056 | 1.007 | 1,116,974 |
| rs777355 | -0.0719 | 0.0079 | -0.0070 | 0.0052 | 0.993 | 1,151,169 |
| rs7938055 | +0.0795 | 0.0090 | -0.0134 | 0.0060 | 0.987 | 1,151,848 |
| rs876701 | +0.0675 | 0.0105 | +0.0115 | 0.0081 | 1.012 | 804,117 |
| rs884205 | -0.0644 | 0.0092 | -0.0312 | 0.0063 | 0.969 | 1,116,927 |
| rs9293511 | -0.0518 | 0.0081 | +0.0028 | 0.0061 | 1.003 | 795,156 |
| rs930900 | +0.0556 | 0.0086 | +0.0139 | 0.0068 | 1.014 | 925,856 |
| rs9515218 | -0.0464 | 0.0081 | -0.0044 | 0.0056 | 0.996 | 1,108,931 |
| rs9594759 | +0.0585 | 0.0080 | +0.0123 | 0.0056 | 1.012 | 1,099,288 |

## Supplementary Table S28. Geometry arm of the mass-versus-geometry decomposition (falsifying test F2): two-sample MR of ten hip-shape modes on hip OA and knee OA, plus a conditional BMD-plus-hip-shape multivariable MR

Areal BMD fuses bone mass, size and geometry, so the BMD->OA effect could in principle be transmitted by bone geometry rather than by mass or formation. Height indexes stature but not bone width or geometry, leaving this split to F2. The ten statistical hip-shape modes (HSM1-HSM10) of Faber et al. 2025 (PMID 39574169; deposited in the GWAS Catalog as GCST90482704-GCST90482713) are derived from DXA imaging and encode geometry only, with no density component; they are used here as pure-geometry instruments in the same pipeline and with the same hip and knee outcomes as Table S27. Instruments are the genome-wide-significant lead variants of each mode (P < 5e-9, 193 distinct positions); exposure effects are the authors' reported betas with direction, and standard errors are derived from the reported P values. Hip-shape mode 5 yielded only one instrument and is not analysable. A conditional BMD-plus-hip-shape multivariable MR extending this screen is reported in the second part of this table.

### Part A. Single-variable two-sample MR of each hip-shape mode on hip OA and knee OA

| Mode | Outcome | n_snp | IVW beta | SE | p (fixed) | OR (95% CI) | p (random) | Cochran Q | I2 | MR-Egger beta | Egger intercept p | Weighted median | Mean F |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HSM1 | hip OA | 21 | +0.0775 | 0.0259 | 0.0028 | 1.081 (1.027-1.137) | 0.496 | 385.6 (20 df) | 95% | +0.1102 | 0.000192 | -0.1183 | 49 |
| HSM1 | knee OA | 21 | -0.0074 | 0.0206 | 0.721 | 0.993 (0.953-1.034) | 0.916 | 229.6 (20 df) | 91% | -0.0285 | 0.00225 | +0.0280 | 49 |
| HSM2 | hip OA | 26 | +0.0879 | 0.0204 | 1.59e-05 | 1.092 (1.049-1.136) | 0.253 | 355.7 (25 df) | 93% | +0.0903 | 0.744 | +0.0238 | 50 |
| HSM2 | knee OA | 24 | -0.0293 | 0.0172 | 0.0883 | 0.971 (0.939-1.004) | 0.531 | 170.0 (23 df) | 86% | -0.0167 | 0.0482 | -0.0038 | 50 |
| HSM3 | hip OA | 19 | -0.0205 | 0.0266 | 0.442 | 0.980 (0.930-1.032) | 0.834 | 243.2 (18 df) | 93% | -0.0307 | 0.576 | -0.0578 | 47 |
| HSM3 | knee OA | 19 | -0.0704 | 0.0211 | 0.000842 | 0.932 (0.894-0.971) | 0.128 | 86.5 (18 df) | 79% | -0.0858 | 0.292 | -0.0440 | 47 |
| HSM4 | hip OA | 26 | +0.0021 | 0.0220 | 0.926 | 1.002 (0.960-1.046) | 0.981 | 374.4 (25 df) | 93% | +0.0069 | 0.0445 | -0.0537 | 51 |
| HSM4 | knee OA | 24 | +0.0481 | 0.0176 | 0.00628 | 1.049 (1.014-1.086) | 0.175 | 93.4 (23 df) | 75% | +0.0488 | 0.822 | +0.0500 | 51 |
| HSM5 | hip OA | 1 | not analysable (insufficient instruments) | | | | | | | | | | |
| HSM5 | knee OA | 1 | not analysable (insufficient instruments) | | | | | | | | | | |
| HSM6 | hip OA | 17 | -0.1697 | 0.0296 | 9.5e-09 | 0.844 (0.796-0.894) | 0.14 | 241.5 (16 df) | 93% | -0.1631 | 7.93e-06 | -0.0220 | 48 |
| HSM6 | knee OA | 16 | -0.1184 | 0.0240 | 8.27e-07 | 0.888 (0.847-0.931) | 0.198 | 219.9 (15 df) | 93% | -0.1356 | 8.1e-09 | -0.1850 | 48 |
| HSM7 | hip OA | 7 | -0.0882 | 0.0493 | 0.0737 | 0.916 (0.831-1.008) | 0.55 | 53.8 (6 df) | 89% | -0.0712 | 0.0483 | -0.1200 | 39 |
| HSM7 | knee OA | 7 | -0.0604 | 0.0391 | 0.122 | 0.941 (0.872-1.016) | 0.54 | 38.2 (6 df) | 84% | -0.0749 | 0.0266 | -0.0700 | 39 |
| HSM8 | hip OA | 19 | +0.1183 | 0.0261 | 5.92e-06 | 1.126 (1.069-1.185) | 0.173 | 198.7 (18 df) | 91% | +0.1679 | 1.05e-09 | +0.0720 | 48 |
| HSM8 | knee OA | 19 | +0.0000 | 0.0209 | 0.998 | 1.000 (0.960-1.042) | 0.999 | 98.4 (18 df) | 82% | +0.0303 | 5.61e-06 | +0.0040 | 48 |
| HSM9 | hip OA | 13 | +0.0118 | 0.0346 | 0.734 | 1.012 (0.946-1.083) | 0.934 | 199.7 (12 df) | 94% | +0.0847 | 8.03e-13 | +0.0457 | 49 |
| HSM9 | knee OA | 12 | -0.0221 | 0.0284 | 0.437 | 0.978 (0.925-1.034) | 0.749 | 64.6 (11 df) | 83% | -0.0444 | 0.00238 | +0.0357 | 49 |
| HSM10 | hip OA | 15 | -0.0197 | 0.0287 | 0.493 | 0.981 (0.927-1.037) | 0.822 | 129.4 (14 df) | 89% | -0.0286 | 3.82e-07 | -0.0500 | 52 |
| HSM10 | knee OA | 15 | +0.0746 | 0.0227 | 0.00104 | 1.077 (1.030-1.126) | 0.25 | 113.8 (14 df) | 88% | +0.0781 | 0.0309 | +0.0070 | 52 |

**Pooled across modes.** Hip OA: pooled beta +0.0186 (SE 0.0090), OR 1.019 (1.001-1.037); across-mode heterogeneity Q = 81.1 (8 df), I2 = 90%. Knee OA: pooled beta -0.0130 (SE 0.0072), OR 0.987 (0.973-1.001); Q = 56.5 (8 df), I2 = 86%. The largest single-mode estimate is 0.170 (hip) and 0.118 (knee) in absolute log-odds, against 0.097 and 0.078 for total-body BMD on the same outcomes; the smallest per-mode P value is 9.5e-09 (hip) and 8.3e-07 (knee), and under the random-effects convention used throughout the manuscript no mode reaches the Bonferroni threshold of 0.005 across the ten modes (hip False, knee False); 4 of ten modes (hip) and 3 of ten (knee) would qualify on fixed-effect P values alone, which is why the random-effects convention governs.

**Interpretation of Part A.** The pooled geometry arm is an order of magnitude smaller than the BMD arm and carries no consistent direction (hip and knee pooled estimates have opposite signs), so the BMD->OA effect is not accounted for by the hip-geometry axes tested. Individual modes, however, are heterogeneous (I2 75-95%) and several show significant MR-Egger intercepts, so the per-mode estimates are pleiotropy-affected and are not read causally; the single-variable design is informative rather than a clean falsification, and it is reported as such. The conditional BMD-plus-hip-shape MVMR in Part B now supplies the formal test and reaches the same qualitative conclusion for knee OA (and is inconclusive for hip OA). The literal volumetric-instrument form of F2 (areal BMD replaced by a volume- or area-conditioned BMD instrument) remains untested because no public summary-statistics file for a volumetric BMD trait was found in the GWAS Catalog; the manuscript's height conditioning (9.5-33.2% attenuation) and BMI conditioning already cover the stature and adiposity axes.

### Part B. Conditional multivariable MR (BMD + hip-shape): formal geometry-conditioning of falsifying test F2

**Design.** A shared-instrument MVMR is infeasible because the 31 total-body BMD instruments and the HSM lead variants are disjoint (zero overlap; see Table S17), so we condition within the HSM-instrument frame: each HSM lead variant is used as the instrument, and GEFOS total-body BMD beta/SE (tidy_femur_total.tsv, 754,224 SNPs) are retrieved on the same variants. Of the 193 HSM lead variants, 59 carried a usable GEFOS BMD effect and entered the conditional MR; n_snp in Part B is therefore the GEFOS-overlapping subset of each mode's lead variants and is smaller than n_snp in Part A. IVW-MVMR follows Sanderson & Windmeijer (2015) and estimates theta_BMD|HSM (BMD's direct effect on OA holding hip-shape geometry constant) and theta_HSM|BMD (geometry's direct effect holding BMD constant), with only outcome SEs used for weighting. Egger-MVMR gives the horizontal-pleiotropy intercept; the Sanderson-Smith conditional F assesses instrument strength for the added BMD coefficient. Modes with <4 overlapping instruments are reported as insufficient (need >=4 for stable identification). The hip-shape geometry axes tested here are purely geometric and, externally, cause fracture but not hip OA (HSM2 OR 1.27, 1.12–1.44 [28]); a null or heterogeneous result is therefore informative about the BMD->OA path specifically.

| Mode | Outcome | n_snp | theta_BMD|HSM (SE) | p_BMD | theta_HSM|BMD (SE) | p_HSM | cond. F | Egger-int p |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HSM1 | hip OA | 5 | -1.288 (0.270) | 1.9e-06 | +0.274 (0.060) | 4.1e-06 | 0.11 | <1e-15 |
| HSM2 | hip OA | 7 | +0.267 (0.107) | 0.013 | +0.134 (0.038) | 4.8e-04 | 0.41 | 8.3e-08 |
| HSM3 | hip OA | 3 | insufficient (need >=4) | | | | | |
| HSM4 | hip OA | 9 | +0.045 (0.151) | 0.77 | -0.080 (0.042) | 0.058 | 3.52 | 0.0013 |
| HSM5 | hip OA | 0 | insufficient (need >=4) | | | | | |
| HSM6 | hip OA | 4 | -0.921 (0.248) | 2.1e-04 | -0.754 (0.068) | <1e-15 | 0.74 | 0.34 |
| HSM7 | hip OA | 3 | insufficient (need >=4) | | | | | |
| HSM8 | hip OA | 5 | -0.411 (0.416) | 0.32 | +0.036 (0.066) | 0.58 | 3.49 | 0.17 |
| HSM9 | hip OA | 1 | insufficient (need >=4) | | | | | |
| HSM10 | hip OA | 4 | +0.162 (0.146) | 0.27 | +0.041 (0.076) | 0.59 | 11.87 | 0.090 |
| HSM1 | knee OA | 5 | -0.127 (0.214) | 0.55 | +0.090 (0.047) | 0.054 | 0.11 | 0.73 |
| HSM2 | knee OA | 7 | +0.133 (0.085) | 0.12 | -0.074 (0.031) | 0.015 | 0.41 | 0.099 |
| HSM3 | knee OA | 3 | insufficient (need >=4) | | | | | |
| HSM4 | knee OA | 8 | +0.283 (0.150) | 0.059 | -0.018 (0.034) | 0.59 | 2.62 | 0.83 |
| HSM5 | knee OA | 0 | insufficient (need >=4) | | | | | |
| HSM6 | knee OA | 3 | insufficient (need >=4) | | | | | |
| HSM7 | knee OA | 3 | insufficient (need >=4) | | | | | |
| HSM8 | knee OA | 5 | +1.065 (0.331) | 0.0013 | +0.177 (0.052) | 6.2e-04 | 3.49 | 0.064 |
| HSM9 | knee OA | 1 | insufficient (need >=4) | | | | | |
| HSM10 | knee OA | 4 | +0.539 (0.116) | 3.4e-06 | +0.248 (0.061) | 4.2e-05 | 11.87 | 0.096 |

**Reading of Part B.** In the knee-OA models — where the marginal BMD->OA effect is robustly positive (Table 2) — the conditional BMD effect given geometry (theta_BMD|HSM) is positive in four of the five estimable modes (HSM2 +0.133, HSM4 +0.283, HSM8 +1.065, HSM10 +0.539) and negative but non-significant in HSM1 (-0.127, p=0.55); the two best-powered modes (HSM8, HSM10, n=5 and 4 with the largest conditional F) are both clearly positive and significant (p=0.0013 and 3.4e-06). Because conditioning on hip-shape geometry leaves the BMD->knee-OA effect intact and positive, bone geometry does not absorb (mediate) the BMD signal — this formally upgrades F2 from a single-variable null screen to a conditional test. In the hip-OA models the conditional estimates are heterogeneous and sign-flipping (HSM6: marginal BMD->OA +0.19 -> conditional -0.92; HSM1 and HSM6 negative, HSM2/HSM4/HSM10 positive), giving no consistent direction; this is consistent with the generally weaker hip-OA evidence and leaves the hip-OA geometry conclusion inconclusive rather than confirmed.

**Why the magnitudes are not quoted as a single coefficient.** The HSM lead variants are strong instruments for geometry but weak for the added BMD coefficient: the conditional F is far below the conventional 10 (range 0.11-11.87; only HSM10 reaches ~11.9, all others <4 and most <1), per-mode n is 4-9, I2 is 48-98%, and Egger-MVMR intercepts are significant in several modes (e.g. HSM1/HSM2 hip, p<1e-7). The conditional BMD coefficients are therefore direction-indicative, not magnitude-precise; cross-mode pooling is not attempted because the per-mode estimates are direction-heterogeneous and the instrument sets are not independent samples. theta_HSM|BMD (geometry's direct OA effect given BMD) is significant in several modes but direction-heterogeneous (e.g. HSM2 hip +0.134, p=4.8e-4; HSM6 hip -0.754, p<1e-15; HSM8 knee +0.177, p=6.2e-4), corroborating the Part A finding that no hip-shape axis reproducibly drives OA.

**Conclusion for F2.** The formal conditional MVMR corroborates the single-variable geometry arm: for knee OA the BMD effect survives geometry conditioning, so bone geometry is not the transmission path; for hip OA the conditional test is underpowered and inconclusive. The literal volumetric-instrument form of F2 (areal BMD replaced by a volume- or area-conditioned BMD instrument) remains untested because no public summary-statistics file for a volumetric BMD trait was found in the GWAS Catalog; the manuscript's height conditioning (9.5-33.2% attenuation) and BMI conditioning already cover the stature and adiposity axes.

## Supplementary Table S29. Osteophyte-stratified public-data surrogate for falsifying test F5: site-stratified MR of BMD on osteophyte-dominant versus cartilage-dominant OA, with F4-style Bayesian colocalisation

F5 was pre-registered as a hip-osteophyte-stratified outcome. The hip-osteophyte GWAS exists only as a conference abstract without public summary statistics, so the registered test could not be run. As a public-data surrogate we substitute site-stratified OA outcomes in which spinal and hand OA are osteophyte-dominant (bone-forming) while knee OA is cartilage-loss-dominant, and ask whether the BMD→OA effect is stronger where osteophyte formation dominates — the formation-side prediction of the Systemic Osteogenic Set-point (SOS). The same 31 BMD_total instruments and the same pipeline as Table 2 are used; ORs are per SD higher genetically determined areal BMD (total femur).

**Part A. Site-stratified two-sample MR (BMD_total → OA).**

| OA site | Dominance | IVW OR per SD (95% CI) | p | Weighted-median OR per SD | MR-Egger OR (slope) | Egger intercept (SE) | Egger intercept p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Knee | Cartilage-dominant | 1.0812 (1.0552–1.1079) | 3.3×10⁻¹⁰ | 1.0470 | 1.0811 | +0.00043 (0.00089) | 0.629 |
| Spine | Osteophyte-dominant | 1.1291 (1.0900–1.1696) | 1.4×10⁻¹¹ | 1.1132 | 1.1295 | −0.00087 (0.00130) | 0.502 |
| Hand | Osteophyte-dominant | 1.1236 (1.0423–1.2113) | 2.4×10⁻³ | 1.0820 | 1.1232 | +0.00171 (0.00272) | 0.529 |

*Notes.* n = 31 BMD_total instruments per site. The spine and hand (osteophyte-dominant) ORs are numerically larger than the knee (cartilage-dominant) OR by ≈4–4.5% on the odds scale (spine−knee log-OR difference = +0.043; hand−knee = +0.039). All three MR-Egger intercepts are non-significant, so the gradient is not driven by directional pleiotropy; confidence intervals overlap broadly and the differences are not individually significant.

**Part B. F4-style Bayesian colocalisation extended to osteophyte-dominant sites (PP.H4 ≥ 0.8; ±500 kb window; priors p1=p2=1×10⁻⁴, p12=5×10⁻⁶, sd_prior quant=0.15 / cc=0.20).**

| OA site | Loci analysed (≥20 shared SNPs) | Loci PP.H4 ≥ 0.8 | Interpretation |
| --- | --- | --- | --- |
| Knee | 31 | 0 | Reproduces the published cross-trait F4 result (knee OA 0/31) |
| Spine | 31 | 0 | No shared causal variant at any instrument locus |
| Hand | 31 | 0 | No shared causal variant at any instrument locus |

*Notes.* Colocalisation is 0/31 at every site, identical to the published knee-OA cross-trait arm of F4 (Tables S8 and S26). The diffuse, polygenic architecture of the BMD→OA effect is therefore site-invariant: colocalisation does not separate osteophyte-dominant from cartilage-dominant sites. The only site-discriminating signal is the MR magnitude gradient in Part A, which is modest and not significant on its own.

**Synthesis for F5.** The magnitude gradient (osteophyte-dominant spine and hand > cartilage-dominant knee) is weakly supportive of the formation-side SOS prediction but not decisive: confidence intervals overlap and no site difference reaches significance. The colocalisation arm provides no discriminating evidence (0/31 everywhere). We therefore report F5 as partially addressed by a public-data surrogate rather than fully executed: the registered hip-osteophyte test remains unperformed for want of public data, and this surrogate substitutes site dominance (spine/hand osteophyte-dominant versus knee cartilage-dominant) as a proxy. The limitation is stated in the pre-registration-status note and the Table S11 ledger.

## Supplementary Table S30. Cross-trait genetic correlation between total-body BMD and site-stratified osteoarthritis (LD Score regression)

Genetic correlations (r_g) between total-body DXA BMD (GCST90568455; verified N = 31,986) and the four site-stratified OA outcomes — hip, knee, spine and hand (Hatzikotoulas et al. 2025, GWAS Catalog GCST90566798/6800/6801/6803, GRCh37) — were estimated with LD Score regression (ldsc v2.0.1, Bulik-Sullivan et al. 2015 [5]) using European-ancestry 1000 Genomes Phase 3 HapMap3 LD Scores (eur_w_ld_chr_1–22; 1,290,028 regression SNPs; the HLA region excluded). All traits were munged against the same w_hm3.snplist (1,217,311 SNPs), so the regression SNP set is identical across traits (685,394–685,477 SNPs per BMD–OA pair). Because r_g = b_gencov / √(b_h2,1 · b_h2,2) is a ratio of LD Score regression slopes, a global rescaling of N cancels exactly; the verified N is therefore required only before quoting h², not r_g.

**Table S30a. Genetic correlation of total-body BMD with site-stratified OA (verified N = 31,986).**

| OA site | r_g | SE | z | P | h²_obs of OA (SE) | genetic-covariance intercept (SE) |
| --- | --- | --- | --- | --- | --- | --- |
| Hip | 0.1309 | 0.0378 | 3.4639 | 5.32×10⁻⁴ | 0.0284 (0.0016) | 0.0063 (0.0066) |
| Knee | 0.1523 | 0.0328 | 4.6443 | 3.41×10⁻⁶ | 0.0319 (0.0014) | −0.0041 (0.0067) |
| Spine | 0.1957 | 0.0468 | 4.1840 | 2.86×10⁻⁵ | 0.0130 (0.0009) | −0.0023 (0.0064) |
| Hand | 0.0713 | 0.0561 | 1.2708 | 0.204 | 0.0088 (0.0010) | −0.0056 (0.0057) |

*Notes.* Three of four BMD–OA correlations are significantly positive but small (0.13–0.20), consistent with a shared polygenic background for bone formation and OA (Fig. S24A); they are not evidence of a causal BMD→OA effect (the causal estimates are the two-sample MR results in Table 2 and S29). All four genetic-covariance intercepts are close to zero, indicating no substantial sample overlap or uncorrected stratification. The h²_obs column is estimated inside each cross-trait regression on the intersection SNP set, so it differs marginally from the single-trait h² in Table S30c; both are observed-scale diagnostics. Repeating the analysis at the placeholder N = 400,000 changed r_g by ≤0.0002 (maximum, knee) and at N = 31,873 gave identical values to three decimals, confirming N-invariance empirically. BMD×OA regressions used 691,715 intersecting SNPs (685,462 / 685,469 / 685,477 / 685,394 SNPs with valid alleles for hip / knee / spine / hand — identical to the pre-correction run).

**Table S30b. Genetic correlation between OA sites (knee as reference; OA inputs unchanged).**

| Pair | r_g | SE | z | P |
| --- | --- | --- | --- | --- |
| Knee–hip | 0.5692 | 0.0191 | 29.76 | 1.41×10⁻¹⁹⁴ |
| Knee–spine | 0.6200 | 0.0233 | 26.63 | 3.43×10⁻¹⁵⁶ |
| Knee–hand | 0.5905 | 0.0392 | 15.08 | 2.11×10⁻⁵¹ |

*Notes.* OA sites share most of their genetic architecture (r_g 0.569–0.620, all P < 10⁻⁵⁰), consistent with the diffuse-polygenic reading of the 0/31 colocalisation result (F4). Small between-site differences in r_g should not be over-interpreted against this shared background.

**Table S30c. Single-trait SNP heritability and LDSC intercept (diagnostic).**

| Trait | h² observed (SE) | Intercept (SE) | λGC | Mean χ² |
| --- | --- | --- | --- | --- |
| Total-body BMD, N = 31,986 | 0.3479 (0.0355) | 1.0293 (0.0094) | 1.1619 | 1.2486 |
| Total-body BMD, N = 31,873 (sensitivity) | 0.3491 (0.0356) | 1.0293 (0.0094) | 1.1619 | 1.2486 |
| OA hip | 0.0283 (0.0013) | 1.0559 (0.0113) | 1.4316 | 1.6424 |
| OA knee | 0.0326 (0.0011) | 1.0532 (0.0111) | 1.5845 | 1.8227 |
| OA spine | 0.0131 (0.0008) | 1.0398 (0.0074) | 1.2136 | 1.2503 |
| OA hand | 0.0087 (0.0009) | 1.0073 (0.0069) | 1.105 | 1.1299 |

*Notes.* Correcting the placeholder N raised BMD h² from 0.0278 to 0.3479 (×12.51; theoretical ×12.51), within the published DXA-BMD SNP-h² range (≈0.25–0.40), while the intercept was unchanged — the correction acts on scale only. For binary OA outcomes, observed-scale h² is diagnostic; a liability-scale h² would additionally require population prevalence.

**Interpretation and honest boundaries.** (i) Direction, not magnitude, is the finding: the positive BMD–OA r_g is a polygenic-sharing signal, distinct from the MR causal estimates. (ii) The site gradient of r_g does not track the MR gradient in S29: LDSC ranks spine highest and hand lowest and non-significant, whereas S29 shows osteophyte-dominant spine and hand ORs above cartilage-dominant knee; the two analyses measure different quantities and are not mutually corroborating. (iii) Hand OA has the lowest h² and mean χ² (1.13), so its null r_g is plausibly a power limitation rather than a biological exception. (iv) LD Scores are European-only while the OA GWAS are mixed-ancestry (≈87.3% European for knee OA); a near-unity intercept does not exclude ancestry-mismatch bias, the same limitation that applies to the main-text LDSC. (v) No cross-study corroboration is claimed: the r_g in supplementary Table S15 derive from a different OA GWAS set and a different 13-phenotype design; numerical similarity between S15 and S30 is descriptive, not an independent replication. (vi) The formation-side SOS reading gains endotype and mechanistic support. Genetically distinct hypertrophic (osteophyte-rich) OA endotypes [46] corroborate an osteophyte-dominant joint pathway, and the RSPO3–NFATC1 subchondral-bone switch [47] — whereby genetic liability to low bone mass protects against OA by increasing subchondral compliance — supplies, by reverse inference, the mechanistic anchor for a high-bone set-point as OA-promoting. A companion five-layer OA–eBMD–fracture funnel [48] (the authors' own unpublished deposit) reports a reverse OA↔bone signal under eBMD instruments; this bounds rather than refutes the constitutional-asymmetry reading of P5, which here concerns DXA areal BMD and fracture under low power (0.14–0.26). (vii) Lineage of the bone–cartilage paradox. SOS stands on a descriptive literature, not a new observation: population surveys found higher axial and appendicular bone mass in generalised OA, with greater bone width obscuring forearm BMD [49], and iliac-crest and subchondral-bone studies led to the proposal that primary OA may initially be a bone disease, stiffer denser bone transmitting more force to cartilage [50]. What SOS adds is a falsifiable causal assertion — one polygenic set-point read as areal BMD, tested by P1–P6/F1–F10 with one instrument set — rather than a restatement of that association; the 2025 MR evidence [11, 12] and the osteophyte-rich endotypes [46] are consistent with, but do not demonstrate, its formation-side core (§4). SOS therefore differs from the bone-forming triad [39] in kind, not degree: the triad describes co-occurring phenotypes, whereas SOS requires one variant set with sign-opposing causal effects at two endpoints — a constraint neither the triad nor the descriptive paradox [3, 4] entails (Fig. S23).

## Supplementary Table S32. Ancestry-composition robustness proxy for falsifying test F10 on the mixed-ancestry knee-OA outcome (code/314)

**Status.** The registered F10 cross-ancestry replication is **not executed**: Hatzikotoulas 2025 released site-stratified (hip/knee/spine/hand) but no ancestry-stratified OA summary statistics, and the previously staged `GCST90566811` is the X-chromosome knee-OA GWAS (0/31 autosomal BMD instruments matched), not a European-restricted autosomal outcome. This table reports a narrower, executable question: is the pooled BMD→knee-OA estimate driven by the variants whose summary statistics were contributed by the most ancestry-diverse set of cohorts? The proxy is a robustness/consistency check **within one mixed-ancestry GWAS**, not a replication of the effect in a non-European population.

**Proxy.** In a multi-ancestry meta-analysis, per-variant `NCASES+NCONTROLS` is the sum over contributing cohorts; across the 31 instruments it ranges 817,610–1,315,578 (median 1,267,074). We map it linearly onto a relative European-purity weight anchored at the study-level 87.31% European composition (Table S14): `t = (N−Nmin)/(Nmax−Nmin)`, `purity = 1−t`, `f_EUR,proxy = 1 − 0.1269·t`. Per-variant N also reflects imputation quality and array coverage, so this is a proxy, not a measured ancestry fraction.

**Estimator.** All rows use the identical 31 European-only total-body BMD instruments and the identical IVW engine (code/13); outcome effects are oriented to the BMD effect allele and their SEs derived from the deposited 95% CIs. P0 is a validation and reproduces the locked headline OR/SD 1.0812 exactly. No previously reported β/SE/p/OR/CI is changed.

| Analysis | n | IVW β (SE) | OR/SD (95% CI) | P | Cochran Q / df (P) | Egger intercept P |
| --- | --- | --- | --- | --- | --- | --- |
| P0 All 31 instruments (reference, locked) | 31 | +0.0781 (0.0124) | 1.0812 (1.0552–1.1079) | 3.31e-10 | 156.6 / 30 (p<0.001) | 0.835 |
| P1 European-core stratum (N ≤ median) | 16 | +0.0687 (0.0214) | 1.0711 (1.0271–1.1171) | 1.34e-03 | 112.4 / 15 (p<0.001) | 0.389 |
| P2 Multi-ancestry-enriched stratum (N > median) | 15 | +0.0828 (0.0152) | 1.0863 (1.0543–1.1193) | 5.63e-08 | 44.0 / 14 (p<0.001) | 0.396 |
| P4 Relative-European-purity-weighted IVW | 31 | +0.0788 (0.0216) | 1.0820 (1.0370–1.1288) | 2.71e-04 | — | — |
| P5 Delete 10 most multi-ancestry instruments | 21 | +0.0877 (0.0183) | 1.0917 (1.0531–1.1316) | 1.72e-06 | 124.0 / 20 (p<0.001) | 0.737 |
| P6 Delete 10 most European-core instruments | 21 | +0.0773 (0.0141) | 1.0803 (1.0508–1.1106) | 4.48e-08 | 71.0 / 20 (p<0.001) | 0.573 |
| P7 1000G-EUR-EAF-concordant subset (abs dEAF ≤ median) | 16 | +0.0793 (0.0156) | 1.0826 (1.0499–1.1162) | 3.81e-07 | 72.8 / 15 (p<0.001) | 0.680 |

**P3 between-stratum test.** The European-core and multi-ancestry-enriched strata give β=+0.0687 (SE 0.0214) and β=+0.0828 (SE 0.0152); Q_between=0.29 on 1 df, P=0.592 — no detectable ancestry-composition dependence.

**P7 orthogonal reference check.** Against 503 1000 Genomes Phase 3 European samples, outcome effect-allele frequencies matched at 31/31 loci (r=0.9825, max |ΔEAF|=0.0446); restricting to the EAF-concordant half gives OR/SD 1.0826 (P=3.81×10⁻⁷). Allele-frequency concordance does not by itself establish ancestry transportability.

**P8 mixture bound.** Under a simple sample-size-weighted two-component mixture (an approximation: meta-analytic weights are inverse-variance, not exactly N-weighted), if the European component effect were exactly zero the 12.69% non-European component would have to carry β=+0.6152 (OR/SD 1.850), i.e. 7.88× the pooled estimate, to generate the observed pooled effect — a magnitude incompatible with the ancestry-stratified null reported in [14] (no ancestry-specific study-wide significant associations across EAS/AFR/SAS/HIS meta-analyses).

**Interpretation and honest boundaries.** (i) Seven of seven estimates are positive and nominally significant, and all remain within a narrow OR range under stratification, continuous weighting, deletion and allele-frequency restriction: the pooled BMD→knee-OA estimate is not attributable to either ancestry-composition stratum alone. (ii) This is robustness to differential cohort contribution within one mixed-ancestry GWAS, not a test of whether the effect replicates in non-European populations; the BMD instruments and the LDSC reference remain European-only. (iii) The per-variant N proxy is confounded by imputation/array coverage, so the strata are not clean ancestry strata and the between-stratum P value is a diagnostic, not a formal ancestry-by-effect interaction test. (iv) F10 therefore moves from "no executed analysis" to "an executed robustness proxy with the registered cross-ancestry replication still unexecuted"; the S23 verdict of 9/10 executed tests is unchanged because a proxy is not a registered-test execution. (v) Closing F10 still requires non-European BMD exposure summary statistics or ancestry-stratified OA GWAS; F2's volumetric arm likewise awaits volumetric BMD/hip-osteophyte data from the Faber/Tobias releases. Full numerical output: `code/output/genetics/L7_F10_ancestry_proxy.json` (`python3 code/314_f10_ancestry_proxy.py`).

## References
1. Dillon CF, Rasch EK, Gu Q, Hirsch R. Prevalence of knee osteoarthritis in the United States: arthritis data from the Third National Health and Nutrition Examination Survey 1991–94. *J Rheumatol*. 2006;33(11):2271–2279. PMID 17013996.
2. Foss MV, Byers PD. Bone density, osteoarthrosis of the hip, and fracture of the upper end of the femur. *Ann Rheum Dis*. 1972;31(4):259–264. doi:10.1136/ard.31.4.259.
3. Nevitt MC, Felson DT. High bone density and radiographic osteoarthritis: questions answered and unanswered. *Osteoarthritis Cartilage*. 2020;28(9):1151–1153. doi:10.1016/j.joca.2020.05.006. PMID 32417559.
4. Dequeker J, Aerssens J, Luyten FP. Osteoarthritis and osteoporosis: clinical and research evidence of inverse relationship. *Aging Clin Exp Res*. 2003;15(5):426–439. doi:10.1007/BF03327364.
5. Bulik-Sullivan B, Finucane HK, Anttila V, et al. An atlas of genetic correlations across human diseases and traits. *Nat Genet*. 2015;47(11):1236–1241. doi:10.1038/ng.3401.
6. Davey Smith G, Ebrahim S. 'Mendelian randomization': can genetic epidemiology contribute to understanding environmental determinants of disease? *Int J Epidemiol*. 2003;32(1):1–22. doi:10.1093/ije/dyg070.
7. Alhassan E, Nguyen K, Hochberg MC, Mitchell BD. Causal Factors for Osteoarthritis: A Scoping Review of Mendelian Randomization Studies. *Arthritis Care Res (Hoboken)*. 2024;76(3):366–375. doi:10.1002/acr.25252. PMID 37846209.
8. Liu H, Zhu J, Jiang H, et al. Using Mendelian randomization analyses in osteoarthritis: state of art and future perspectives. *Ann Med*. 2025;57(1):2533428. doi:10.1080/07853890.2025.2533428. PMID 40697045.
9. Hackinger S, Trajanoska K, Styrkarsdottir U, et al. Evaluation of shared genetic aetiology between osteoarthritis and bone mineral density identifies SMAD3 as a novel osteoarthritis risk locus. *Hum Mol Genet*. 2017;26(19):3850–3858. doi:10.1093/hmg/ddx285.
10. Zhou J, Na J, Jiang Z, et al. Deciphering the shared genetic structure between hip osteoarthritis and femoral neck bone mineral density. *Front Genet*. 2025;16:1597005. doi:10.3389/fgene.2025.1597005.
11. Yu H, Wang J, Xu X, Li H, Guo J. Revealing the mediating mechanisms between BMI and osteoarthritis: a Mendelian randomization and mediation analysis. *Aging Clin Exp Res*. 2025;37(1):119. doi:10.1007/s40520-025-03035-2.
12. Hartley AE, Sanderson ECM, Granell R, Paternoster L, Zheng J, Davey Smith G, et al.; Genetics of Osteoarthritis Consortium; Gregson CL, Tobias JH. Using multivariable Mendelian randomization to estimate the causal effect of bone mineral density on osteoarthritis risk, independently of body mass index. *Int J Epidemiol*. 2022;51(4):1254–1267. doi:10.1093/ije/dyab251. PMID 34897459.
13. Skrivankova VW, Richmond RC, Woolf BAR, Yarmolinsky J, Davies NM, Swanson SA, VanderWeele TJ, Higgins JPT, Timpson NJ, Dimou N, Langenberg C, Golub RM, Loder EW, Gallo V, Tybjaerg-Hansen A, Davey Smith G, Egger M, Richards JB. Strengthening the Reporting of Observational Studies in Epidemiology Using Mendelian Randomization: the STROBE-MR statement. *JAMA*. 2021;326(16):1614–1621. doi:10.1001/jama.2021.18236. PMID 34698778.
14. Hatzikotoulas K, Southam L, Stefansdottir L, et al. Translational genomics of osteoarthritis in 1,962,069 individuals. *Nature*. 2025;641(8065):1217–1224. doi:10.1038/s41586-025-08771-z. PMID 40205036.
15. Elsworth B, Lyon M, Alexander T, et al. The MRC-IEU UK Biobank GWAS pipeline. *bioRxiv*. 2018. doi:10.1101/340232. (Self-reported OA accessed as `ukb-b-14486` via OpenGWAS, <https://gwas.mrcieu.ac.uk/datasets/ukb-b-14486/>.)
16. Qian Y, Xia J, Wang P, et al. Genome-wide Association Studies of over 30,000 Samples with Bone Mineral Density at Multiple Skeletal Sites and Its Clinical Relevance. *Genomics Proteomics Bioinformatics*. 2025;23(5):qzaf097. doi:10.1093/gpbjnl/qzaf097. PMID 41206123.
17. Trajanoska K, Morris JA, Kaptoge S, et al., for the GEFOS consortium. Genome wide association analysis of risk of fracture and its association with bone mineral density: a meta-analysis of genome-wide association studies applied to the UK Biobank. *BMJ*. 2018;362:k3225. doi:10.1136/bmj.k3225. PMID 30158200.
18. Burgess S, Davies NM, Thompson SG. Bias due to participant overlap in two-sample Mendelian randomization. *Genet Epidemiol*. 2016;40(7):597–608. doi:10.1002/gepi.21998.
19. Bowden J, Davey Smith G, Burgess S. Mendelian randomization with invalid instruments: effect estimation and bias detection through Egger regression. *Int J Epidemiol*. 2015;44(2):512–525. doi:10.1093/ije/dyv080. PMID 26050253.
20. Bowden J, Davey Smith G, Haycock PC, Burgess S. Consistent estimation of the causal effect in Mendelian randomization with some invalid instruments using a weighted median estimator. *Genet Epidemiol*. 2016;40(4):304–314. doi:10.1002/gepi.21965.
21. Mounier N, Kutalik Z. Bias correction for inverse variance weighted Mendelian randomization. *Genet Epidemiol*. 2023;47:e22522. doi:10.1002/gepi.22522.
22. Giambartolomei C, Vukcevic D, Schadt EE, Franke L, Hingorani AD, Wallace C, Plagnol V. Bayesian test for colocalisation between pairs of genetic association studies using summary statistics. *PLoS Genet*. 2014;10(5):e1004383. doi:10.1371/journal.pgen.1004383. PMID 24830394.
23. Yengo L, Vedantam S, Marouli E, et al. A saturated map of common genetic variants associates polygenic risk with human height. *Nature*. 2022;610(7933):704–712. doi:10.1038/s41586-022-05275-y.
24. Sanderson E, Davey Smith G, Windmeijer F, Bowden J. An examination of multivariable Mendelian randomization in the single-sample and two-sample summary data settings. *Int J Epidemiol*. 2018;48(3):713–727. doi:10.1093/ije/dyy262.
25. Burgess S, Jenkins DA, Davies NM, et al. Multivariable Mendelian randomization: the use of pleiotropic variants to estimate causal effects of multiple risk factors on a single outcome. *Am J Epidemiol*. 2015;181(4):251–260. doi:10.1093/aje/kwu283.
26. Pulit SL, Stoneman C, Morris AP, et al. Meta-analysis of genome-wide association studies for body fat distribution in 694,649 individuals of European ancestry. *Hum Mol Genet*. 2019;28(1):166–174. doi:10.1093/hmg/ddy327. PMID 30239722. (BMI summary statistics from the GIANT–UK Biobank meta-analysis: Zenodo, https://doi.org/10.5281/zenodo.1251813.)
27. Verbanck M, Chen CY, Neale B, Do R. Detection of widespread horizontal pleiotropy in causal relationships inferred from Mendelian randomization between complex traits. *Nat Genet*. 2018;50(5):693–698. doi:10.1038/s41588-018-0099-7.
28. Faber BG, Frysz M, Zheng J, et al. The genetic architecture of hip shape and its role in the development of hip osteoarthritis and fracture. *Hum Mol Genet*. 2025;34(3):207–217. doi:10.1093/hmg/ddae169.
29. Kröger H, Vainio P, Nieminen J, Kotaniemi A. Comparison of different models for interpreting bone mineral density measurements using DXA and MRI technology. *Bone*. 1995;17(2):157–159. doi:10.1016/s8756-3282(95)00162-x. PMID 8554924.
30. Looker AC, Beck TJ, Orwoll ES. Does body size account for gender differences in femur bone density and geometry? *J Bone Miner Res*. 2001;16(7):1291–1299. doi:10.1359/jbmr.2001.16.7.1291.
31. Welling M, Auvinen J, Lehenkari P, Männikkö M, Karppinen J, Eskola PJ. Association between height and osteoarthritis of the knee and hip: the Northern Finland Birth Cohort 1966 Study. *Int J Rheum Dis*. 2017;20(9):1095–1104. doi:10.1111/1756-185X.13059. PMID 28296269.
32. Hunter DJ, Niu J, Zhang Y, et al. Knee height, knee pain, and knee osteoarthritis: the Beijing Osteoarthritis Study. *Arthritis Rheum*. 2005;52(5):1418–1423. doi:10.1002/art.21017. PMID 15880346.
33. Sun Y, Li Y, Yu T, Zhang J. Causal associations of anthropometric measurements with osteoarthritis: a Mendelian randomization study. *PLoS One*. 2023;18(1):e0279198. doi:10.1371/journal.pone.0279198. PMID 36716300.
34. Funck-Brentano T, Nethander M, Moverare Skrtic S, Richette P, Ohlsson C. Causal factors for knee, hip, and hand osteoarthritis: a Mendelian randomization study in the UK Biobank. *Arthritis Rheumatol*. 2019;71(10):1634–1641. doi:10.1002/art.40928. PMID 31099188.
35. Goldring MB, Goldring SR. Articular cartilage and subchondral bone in the pathogenesis of osteoarthritis. *Nat Rev Rheumatol*. 2016;12(11):632–644. doi:10.1038/nrrheum.2016.148. PMID 27652499.
36. Burr DB, Gallant MA. Bone remodelling in osteoarthritis. *Nat Rev Rheumatol*. 2012;8(11):665–673. doi:10.1038/nrrheum.2012.130. PMID 22868925.
37. Zhen G, Wen C, Jia X, et al., Cao X. Inhibition of TGF-β signaling in mesenchymal stem cells of subchondral bone attenuates osteoarthritis. *Nat Med*. 2013;19(6):704–712. doi:10.1038/nm.3143. PMID 23685840.
38. Findlay DM, Kuliwaba JS. Bone–cartilage crosstalk: a conversation for understanding osteoarthritis. *Bone Res*. 2016;4:16028. doi:10.1038/boneres.2016.28. PMID 27672480.
39. Hardcastle SA, Dieppe P, Gregson CL, et al. Osteophytes, enthesophytes, and high bone mass: a bone-forming triad with potential relevance in osteoarthritis. *Arthritis Rheumatol*. 2014;66(9):2429–2439. doi:10.1002/art.38729. PMID 24910132.
40. Cui B, Chen Y, Tian Y, Liu H, Huang Y, Yin G, Xie Q. Effects of medications on incidence and risk of knee and hip joint replacement in patients with osteoarthritis: a systematic review and meta-analysis. *Adv Rheumatol*. 2022;62:22. doi:10.1186/s42358-022-00253-4. PMID 35761320.
41. Hardcastle SA, Dieppe P, Gregson CL, Davey Smith G, Tobias JH. Osteoarthritis and bone mineral density: are strong bones bad for joints? *BoneKEy Rep*. 2015;4:624. doi:10.1038/bonekey.2014.119. PMID 25628884.
42. Estrada K, Styrkarsdottir U, Evangelou E, et al. Genome-wide meta-analysis identifies 56 bone mineral density loci, provides insights into the genetic architecture of osteoporosis, and reveals a role for RSPO3. *PLoS Genet*. 2012;8(1):e1001217. doi:10.1371/journal.pgen.1001217. PMID 22253714.
43. Reginster JY, Badurski J, Bellamy N, Bensen W, Chapurlat R, Chevalier X, Christiansen C, Genant H, Navarro F, Nasonov E, Sambrook PN, Spector TD, Cooper C. Efficacy and safety of strontium ranelate in the treatment of knee osteoarthritis: results of a double-blind, randomised placebo-controlled trial. *Ann Rheum Dis*. 2013;72(2):179–186. doi:10.1136/annrheumdis-2012-202231. PMID 23117245. (SEKOIA; ISRCTN41323372.)
44. Sobue Y, Kosugiyama H, Asai S, et al. A randomized controlled trial comparing romosozumab and denosumab in elderly women with primary osteoporosis and knee osteoarthritis. *Sci Rep*. 2025;15:22441. doi:10.1038/s41598-025-05187-7. PMID 40594305.
45. Lee SH, Wray NR, Goddard ME, Visscher PM. Estimating missing heritability for disease from genome-wide association studies. *Am J Hum Genet*. 2011;88(3):294–305. doi:10.1016/j.ajhg.2011.02.002.
46. Faber BG, Jung M, Ebsim R, et al. Structural phenotypes of osteoarthritis are clinically and genetically distinct: findings from 59,539 UK Biobank participants. medRxiv. 2026. doi:10.64898/2026.02.08.26345686.
47. Peng L, Wei J, Yang H, et al. Deciphering the RSPO3-NFATC1 switch in the bone-cartilage paradox: genetic liability to osteoporosis protects against osteoarthritis by modulating subchondral bone compliance. Osteoarthritis Cartilage. 2026 May 30:S1063-4584(26)00921-0. doi:10.1016/j.joca.2026.05.013. PMID 42218990.
48. Zheng Z. Shared genetic architecture of osteoarthritis and bone traits (OA × eBMD × fracture): a five-layer cross-trait multi-omics funnel. Zenodo. 2026. doi:10.5281/zenodo.23123786. (Authors' own unpublished deposit.)
49. Dequeker J. Inverse relationship of interface between osteoporosis and osteoarthritis. *J Rheumatol*. 1997;24(4):795–798. PMID 9101521.
50. Dequeker J, Mokassa L, Aerssens J. Bone density and osteoarthritis. *J Rheumatol Suppl*. 1995;43:98–100. PMID 7752151.

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
- **Data availability.** Summary statistics are from publicly published GWAS (Refs 14–17, 23, 26) and were not redistributed; all source summary statistics were accessed on 17 July 2026. Provenance and accessions are listed above. The analysis code and all reproduction material have been deposited in a public GitHub repository (https://github.com/Fisherzzx001/OA-BMD-fracture-MR) and archived with a Zenodo DOI via the GitHub–Zenodo integration (https://doi.org/10.5281/zenodo.23400860; this concept DOI always resolves to the most recent version and lists every version of the deposit). The deposit includes every script under `code/` together with the result tables and source data under `code/output/` required to reproduce every table and figure (including the MRlap sample-overlap correction); bulk intermediate files that are regenerated deterministically from the source summary statistics by the deposited scripts (for example PLINK clumping inputs and harmonised files) are not archived, together with the supplementary files enumerated in the Supplementary Material section (STROBE-MR checklist, Table S10 MVMR grid, Table S11 SOS predictions and falsifying tests). All source GWAS were approved by their original ethics committees. This secondary analysis required no new ethical approval. An independent preregistration of this study by the authors (our own protocol) — distinct from the archived comorbidity-network protocol of the OA project — is deposited in the repository (`docs/118_OSF预注册_OA_BMD_fracture.md`) and archived via the Zenodo DOI above; it declares P1–P6 as derived predictions and F1–F10 as registered tests (F4 cross-trait and eQTLGen arms, and F6–F9, executed with decision rules fixed before execution; F1 executed prospectively against a pre-registered conjunctive criterion and passed; F2 executed in its geometry arm, the volume arm remaining untested for want of public data; F3 executed; F10's registered cross-ancestry replication not performed for want of public ancestry-stratified data, with an ancestry-composition robustness proxy reported in Table S32; F5, declared not performed at preregistration for want of a hip-osteophyte GWAS, is executed here as a public-data surrogate (site-stratified MR and F4-style colocalisation, Table S29), because site-stratified OA outcomes — spine and hand osteophyte-dominant, knee cartilage-dominant — are publicly available). An OSF mirror is being deposited.
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
- **Figure S18.** Endpoint-specific colocalisation at the 31 BMD instrument loci (PP.H4, BMD↔knee OA versus BMD↔fracture).
- **Figure S19.** Colocalisation at WNT16, SOX6 and GDF5 across the five BMD sites that entered MR.
- **Figure S20.** Multi-estimator robustness grid of the ten MR pairs (IVW, MR-Egger, weighted median, leave-one-out, MR-PRESSO).
- **Figure S21.** The Systemic Osteogenic Set-point (SOS) hypothesis schematic.
- **Figure S22.** Latent structure of the 13-phenotype genetic-correlation matrix (PCA of the 13×13 r_g matrix; BMD, OA and fracture load on separate axes).
- **Table S8.** Colocalisation of total-body BMD with knee OA at all 31 instrument loci.
- **Table S9.** Multivariable MR of BMD→knee OA and BMD→fracture, conditioning on height (primary n_h=100 / EUR configuration).
- **Table S10.** Multivariable MR grid: full sensitivity grid (EUR and EUR_noUKB; n_h ∈ {0,100,300,1,000}).
- **Table S11.** Systemic Osteogenic Set-point (SOS) theory: predictions (P1–P6) and falsifying tests (F1–F10).
- **Table S12.** Exposure (BMD-site) selection and instrument counts: the 11-site source, the a priori cranial-head exclusion, genome-wide-significant locus counts per site, and which sites entered MR (`variant_selection_table.csv` in `code/output/genetics/`).
- **Table S13.** Univariate instrument strength and Steiger directionality for the ten MR pairs: per-SNP and aggregate F, R²_exposure, Steiger z and one-sided p, and direction supported (`instrument_strength_steiger.csv` in `code/output/genetics/`).
- **Table S14.** Population-ancestry composition of the 13 GWAS and the LDSC reference used for each.
- **Table S15.** Genetic correlations by phenotype block (n pairs, r_g range, r_g mean, proportion significant at q<0.05).
- **Table S16.** Multivariable MR of BMD→knee OA and BMD→fracture, conditioning on BMI (primary n_bmi=100 / GIANT–UK Biobank configuration).
- **Table S17.** Sensitivity analyses requested at review: (a) uncertainty of the effect-scale reconstruction constant K; (b) significance counts under three random-effects conventions; (c) change in the IVW estimate on removal of MR-PRESSO outliers for all ten pairs; (d) colocalisation counts as a function of the PP.H4 threshold.
- **Table S31.** Multivariable MR of OA_knee and BMD_total on any-fracture, conditioning on BMD (code/264): the OA→fracture effect is unchanged after conditioning on BMD, confirming BMD is not a confounder of the null.
- **Table S25.** Bone-versus-cartilage tissue enrichment — competitive gene-set test (F3, this revision): OA knee enriches cartilage-development and chondrocyte-differentiation gene sets more strongly than bone-development sets; total-body BMD enriches skeletal-system-development most strongly (code/171).
- **Table S26.** eQTL colocalisation — tissue-literal (cartilage and skeletal-muscle bone-proxy) arms of F4 (this revision): at the 31 BMD instrument loci, cartilage eQTL colocalisation with OA knee, BMD total or fracture is 0/31 at PP.H4≥0.8, and skeletal-muscle (bone-proxy) eQTL is 0–1/31 (one borderline OA-knee locus, PP.H4=0.887) — neither tissue colocalises with BMD/fracture, consistent with a formation-side rather than cartilage-matrix or muscle-mediated interpretation (code/172a, 172b, 173).
- **Table S27.** Prospective hip-osteoarthritis replication (falsifying test F1): two-sample MR of total-body BMD on hip OA, and the conjunctive adjudication against the knee-OA estimate.
- **Table S28.** Geometry arm of the mass-versus-geometry decomposition (falsifying test F2): two-sample MR of ten hip-shape modes on hip OA and knee OA, plus a conditional BMD-plus-hip-shape multivariable MR formalising the geometry test.
- **Table S29.** Osteophyte-stratified public-data surrogate for falsifying test F5: site-stratified two-sample MR of BMD_total on spine, hand and knee OA (osteophyte-dominant versus cartilage-dominant) with F4-style Bayesian colocalisation at the 31 instrument loci (code/200, code/201).
- **Table S30.** Cross-trait genetic correlation between total-body BMD and site-stratified OA (hip, knee, spine, hand) from LD Score regression (ldsc v2.0.1, European HapMap3 LD Scores; verified BMD N = 31,986): Tables S30a (BMD×OA r_g), S30b (OA site-to-site r_g), S30c (single-trait h² and LDSC intercept) with interpretation and honest boundaries (docs/Supplementary/Table_S29_ldsc_site_rg.md; code/207).
- **Table S32.** Ancestry-composition robustness proxy for falsifying test F10: per-variant-N-stratified, purity-weighted, deletion and 1000G-EUR-EAF-concordance IVW sensitivities on the mixed-ancestry knee-OA outcome, with the registered cross-ancestry replication explicitly remaining unexecuted (`docs/Supplementary/Table_S32_f10_ancestry_proxy.md`; code/314).
- **STROBE-MR checklist.** Item-by-item reporting-compliance table for all 20 STROBE-MR items, with locations in the manuscript (Skrivankova et al. 2021 [13]).
- **Graphical Abstract.** The OA–BMD–fracture triangle: higher BMD protects against fracture yet raises knee-OA risk.

Supplementary figures S1–S24 and the Graphical Abstract are provided as high-resolution PNG in the submission repository (`Figures/`); supplementary tables S8–S17, S25–S32 and the STROBE-MR checklist are provided as Markdown/CSV.
