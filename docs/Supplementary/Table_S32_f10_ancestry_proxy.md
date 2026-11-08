# Table S32. F10 ancestry-composition robustness proxy on the mixed-ancestry knee-OA outcome

**Status.** The registered F10 cross-ancestry replication is **not executed**: Hatzikotoulas 2025 released site-stratified (hip/knee/spine/hand) but no ancestry-stratified OA summary statistics, and the previously staged `GCST90566811` is the X-chromosome knee-OA GWAS (0/31 autosomal BMD instruments matched), not a European-restricted autosomal outcome. This table reports a narrower, **executable** question: is the pooled BMD→knee-OA estimate driven by the variants whose summary statistics were contributed by the most ancestry-diverse set of cohorts?

**Proxy.** In a multi-ancestry meta-analysis, per-variant `NCASES+NCONTROLS` is the sum over contributing cohorts; across the 31 instruments it ranges 817,610–1,315,578 (median 1,267,074). We map it linearly onto a relative European-purity weight anchored at the study-level 87.31% European composition (Table S14): `t = (N−Nmin)/(Nmax−Nmin)`, `purity = 1−t`, `f_EUR,proxy = 1 − 0.1269·t`. Per-variant N also reflects imputation quality and array coverage, so this is a **proxy, not a measured ancestry fraction**; every row below is a robustness/consistency check, **not** a replication of the effect in a non-European population.

**Estimator.** All rows use the identical 31 European-only total-body BMD instruments and the identical IVW engine (`code/13`); outcome effects are oriented to the BMD effect allele and their SEs derived from the deposited 95% CIs. P0 is a validation and must reproduce the locked headline OR/SD 1.0812 — it does, exactly. No previously reported β/SE/p/OR/CI is changed.

| Analysis | n | IVW β (SE) | OR/SD (95% CI) | P | Cochran Q / df (P) | Egger intercept P |
|---|---|---|---|---|---|---|
| P0 All 31 instruments (reference, locked) | 31 | +0.0781 (0.0124) | 1.0812 (1.0552–1.1079) | 3.31e-10 | 156.6 / 30 (p=0.000) | 0.835 |
| P1 European-core stratum (N <= median, n=16) | 16 | +0.0687 (0.0214) | 1.0711 (1.0271–1.1171) | 1.34e-03 | 112.4 / 15 (p=0.000) | 0.389 |
| P2 Multi-ancestry-enriched stratum (N > median, n=15) | 15 | +0.0828 (0.0152) | 1.0863 (1.0543–1.1193) | 5.63e-08 | 44.0 / 14 (p=0.000) | 0.396 |
| P4 relative-European-purity-weighted IVW | 31 | +0.0788 (0.0216) | 1.0820 (1.0370–1.1288) | 2.71e-04 | — | — |
| P5 Delete 10 most multi-ancestry instruments (n=21) | 21 | +0.0877 (0.0183) | 1.0917 (1.0531–1.1316) | 1.72e-06 | 124.0 / 20 (p=0.000) | 0.737 |
| P6 Delete 10 most European-core instruments (n=21) | 21 | +0.0773 (0.0141) | 1.0803 (1.0508–1.1106) | 4.48e-08 | 71.0 / 20 (p=0.000) | 0.573 |
| P7 1000G-EUR-EAF-concordant subset (abs dEAF <= median, n=16) | 16 | +0.0793 (0.0156) | 1.0826 (1.0499–1.1162) | 3.81e-07 | 72.8 / 15 (p=0.000) | 0.680 |

**P3 between-stratum test.** The European-core and multi-ancestry-enriched strata give β=+0.0687 (SE 0.0214) and β=+0.0828 (SE 0.0152); Q_between=0.29 on 1 df, P=0.592 — no detectable ancestry-composition dependence.

**P7 orthogonal reference check.** Against 503 1000 Genomes Phase 3 European samples, outcome effect-allele frequencies matched at 31/31 loci (r=0.9825, max |ΔEAF|=0.0446); restricting to the EAF-concordant half gives OR/SD 1.0826 (P=3.81e-07). Allele-frequency concordance does not by itself establish ancestry transportability.

**P8 mixture bound.** Under a simple sample-size-weighted two-component mixture (an approximation: meta-analytic weights are inverse-variance, not exactly N-weighted), if the European component effect were exactly zero the 12.69% non-European component would have to carry β=+0.6152 (OR/SD 1.850), i.e. 7.88× the pooled estimate, to generate the observed pooled effect — a magnitude incompatible with the ancestry-stratified null reported in [14] (no ancestry-specific study-wide significant associations across EAS/AFR/SAS/HIS meta-analyses).

**Interpretation and honest boundaries.** 7/7 estimates are positive and nominally significant, and all remain within a narrow OR range under stratification, continuous weighting, deletion and allele-frequency restriction: the pooled BMD→knee-OA estimate is not attributable to either ancestry-composition stratum alone. (i) This is robustness to differential cohort contribution *within one mixed-ancestry GWAS*, not a test of whether the effect replicates in non-European populations; the BMD instruments and the LDSC reference remain European-only. (ii) The per-variant N proxy is confounded by imputation/array coverage, so the strata are not clean ancestry strata; the between-stratum P value should be read as a diagnostic, not as a formal ancestry-by-effect interaction test. (iii) F10 therefore changes from 'no executed analysis' to 'an executed robustness proxy with the registered cross-ancestry replication still unexecuted'; the S23 verdict of 9/10 executed tests is unchanged because a proxy is not a registered-test execution. (iv) Closing F10 still requires non-European BMD exposure summary statistics or ancestry-stratified OA GWAS; F2's volumetric arm likewise awaits volumetric BMD/hip-osteophyte data from the Faber/Tobias releases.

Reproduce: `python3 code/314_f10_ancestry_proxy.py` → `code/output/genetics/L7_F10_ancestry_proxy.json`.
