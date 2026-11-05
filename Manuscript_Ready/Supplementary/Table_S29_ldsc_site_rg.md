# Supplementary Table S29 · Cross-trait genetic correlation between total-body BMD and site-stratified osteoarthritis

This file is the standalone, submission-ready version of manuscript supplementary **Table S29** (embedded in `docs/08_投稿稿_OA_BMD_fracture.md`). It reports cross-trait LD Score regression genetic correlations (r_g) between total-body DXA BMD and the four site-stratified OA outcomes, and the single-trait h² / LDSC intercepts. No reported β, SE, p, OR or CI in the manuscript is altered by this table; values are verified by `code/207_bmd_verified_ldsc.py` (`COMPARISON.md`), which is a re-run, not a restatement.

## Supplementary Methods

**Cross-trait LD Score regression.** Genetic correlations (r_g) between total-body DXA-derived bone mineral density and four site-stratified osteoarthritis outcomes (hip, knee, spine, hand) were estimated with LD Score regression (ldsc v2.0.1) using European-ancestry LD Scores computed from 1000 Genomes Phase 3 HapMap3 SNPs (`eur_w_ld_chr_1–22`, 1,290,028 regression SNPs; M = 1,290,050, M_5_50 = 1,173,569; the HLA region excluded), supplied as both `--ref-ld-chr` and `--w-ld-chr`. OA summary statistics were from Hatzikotoulas et al. 2025 (GWAS Catalog accessions GCST90566798, GCST90566800, GCST90566801, GCST90566803), reported on GRCh37; position-encoded variant identifiers were mapped to rsIDs via the 1000 Genomes Phase 3 EUR panel, and non-biallelic, strand-ambiguous and unmapped variants were removed. Because these are multi-cohort meta-analyses, the per-variant sample size N = NCASES + NCONTROLS was retained for every row (`--n-min 1`); applying munge's default 90th-percentile N filter would have systematically removed SNPs contributed by smaller cohorts. Z-scores were reconstructed from the reported P-values, with the odds ratio used only to set the sign, because the standard-error field was missing throughout. All five traits were munged against the same `w_hm3.snplist` (1,217,311 SNPs, SNP/A1/A2), so effect alleles and the regression SNP set are identical across traits; SNPs absent from a GWAS were written with missing Z and dropped by ldsc.

**Sample size of the BMD exposure.** The legacy BMD file carried a constant placeholder N = 400,000. The source trait is total-body DXA BMD (GCST90568455, UK Biobank field 23236, Qian et al. 2025), curated at N = 31,986 European-ancestry participants, a value independently supported by reverse-engineering the reported SE at the lead locus (≈29,300–32,000). We therefore re-munged BMD with the verified constant N (via `--ignore N --N 31986`, since munge applies `--N` only when no N column is present) and repeated all analyses at N = 31,873, the femur-total N used for the manuscript's MR N-correction, as a sensitivity analysis. Because r_g = b_gencov/√(b_h2,1·b_h2,2) is a ratio of LD Score regression slopes, a global rescaling of N cancels exactly; r_g is N-invariant, whereas single-trait h² scales as 1/N. The N correction is therefore required only before quoting h², not before quoting r_g.

**Reporting scale and limitations.** Cross-trait LDSC intercepts are reported to assess sample overlap and stratification. Single-trait h² is on the observed scale; for binary OA outcomes a liability-scale h² would additionally require population prevalence, so OA h² values are diagnostic only. BMD was assigned a single constant N because per-variant effective sample sizes are not distributed with the source statistics. All estimates are European-reference and sex-combined. LDSC r_g, two-sample MR and Bayesian colocalisation estimate different quantities — genome-wide polygenic sharing, the effect of 31 instrument variables, and single-locus shared causal variants respectively — and are not mutually corroborating.

## Results

**Table S29a. Genetic correlation of total-body BMD with site-stratified OA (verified N = 31,986).**

| OA site | r_g | SE | z | P | h²_obs of OA (SE) | genetic-covariance intercept (SE) |
| --- | --- | --- | --- | --- | --- | --- |
| Hip | 0.1309 | 0.0378 | 3.4639 | 5.32 × 10⁻⁴ | 0.0284 (0.0016) | 0.0063 (0.0066) |
| Knee | 0.1523 | 0.0328 | 4.6443 | 3.41 × 10⁻⁶ | 0.0319 (0.0014) | −0.0041 (0.0067) |
| Spine | 0.1957 | 0.0468 | 4.1840 | 2.86 × 10⁻⁵ | 0.0130 (0.0009) | −0.0023 (0.0064) |
| Hand | 0.0713 | 0.0561 | 1.2708 | 0.204 | 0.0088 (0.0010) | −0.0056 (0.0057) |

BMD×OA regressions used 691,715 intersecting SNPs, with 685,462 / 685,469 / 685,477 / 685,394 SNPs carrying valid alleles for hip / knee / spine / hand. These counts are **identical to the pre-correction run**, confirming that re-munging BMD with `--merge-alleles` and the verified N changed neither the regression SNP set nor the allele alignment (the source BMD file already contained 691,715 HapMap3 SNPs with zero allele mismatches). The h²_obs column is estimated inside each cross-trait regression on the intersection SNP set and its weights, so it differs slightly from the single-trait h² in Table S29c; both are observed-scale diagnostics. Repeating the analysis at the placeholder N = 400,000 changed r_g by ≤0.0002 (maximum, knee) and at N = 31,873 gave identical values to three decimals, confirming N-invariance empirically. All four genetic-covariance intercepts are close to zero.

**Table S29b. Genetic correlation between OA sites (knee as reference; OA inputs unchanged).**

| Pair | r_g | SE | z | P |
| --- | --- | --- | --- | --- |
| Knee–hip | 0.5692 | 0.0191 | 29.76 | 1.41 × 10⁻¹⁹⁴ |
| Knee–spine | 0.6200 | 0.0233 | 26.63 | 3.43 × 10⁻¹⁵⁶ |
| Knee–hand | 0.5905 | 0.0392 | 15.08 | 2.11 × 10⁻⁵¹ |

**Table S29c. Single-trait SNP heritability and LDSC intercept (diagnostic).**

| Trait | h² observed (SE) | Intercept (SE) | λGC | Mean χ² |
| --- | --- | --- | --- | --- |
| Total-body BMD, N = 31,986 | 0.3479 (0.0355) | 1.0293 (0.0094) | 1.1619 | 1.2486 |
| Total-body BMD, N = 31,873 (sensitivity) | 0.3491 (0.0356) | 1.0293 (0.0094) | 1.1619 | 1.2486 |
| OA hip | 0.0283 (0.0013) | 1.0559 (0.0113) | 1.4316 | 1.6424 |
| OA knee | 0.0326 (0.0011) | 1.0532 (0.0111) | 1.5845 | 1.8227 |
| OA spine | 0.0131 (0.0008) | 1.0398 (0.0074) | 1.2136 | 1.2503 |
| OA hand | 0.0087 (0.0009) | 1.0073 (0.0069) | 1.105 | 1.1299 |

Correcting the placeholder N raised BMD h² from 0.0278 to 0.3479 (×12.514; theoretical ×12.505), within the published DXA-BMD SNP-h² range (≈0.25–0.40), while the intercept was unchanged — the correction acts on scale only. The 31,986 / 31,873 discrepancy (0.35%) does not affect any rounded value.

## Interpretation and honest boundaries

1. **Direction, not magnitude, is the finding.** Three of four BMD–OA correlations are significantly positive but small (0.13–0.20), consistent with a shared polygenic background for bone formation and OA. These correlations are *not* evidence of a causal BMD→OA effect; the causal estimates are the two-sample MR results (main text Table 2; Table S28).
2. **The site gradient does not track the MR gradient.** In the site-stratified MR (Table S28) the osteophyte-dominant spine and hand sites have larger odds ratios than cartilage-dominant knee OA (1.129, 1.124 vs 1.081). Here r_g is highest for spine, lowest and non-significant for hand, and intermediate for hip. The "osteophyte-dominant ⇒ stronger BMD correlation" conjecture is therefore **not supported** by LDSC; the two analyses measure different quantities. Hand OA has the lowest h² and mean χ² (1.13), so its null r_g is more plausibly a power limitation than a biological exception.
3. **OA sites share most of their genetic architecture** (r_g 0.569–0.620, all P < 10⁻⁵⁰), consistent with the diffuse-polygenic reading of the 0/31 colocalisation result (F4). Small between-site differences in r_g should not be over-interpreted against this shared background.
4. **No cross-study corroboration is claimed.** The r_g values in supplementary Table S15 derive from a different set of OA GWAS (GO Consortium 2025 knee OA, self-reported OA) and a different 13-phenotype design. Numerical similarity between S15 and S29 is descriptive, not an independent replication.
5. **Ancestry.** LD Scores are European-only while the OA GWAS are mixed-ancestry (≈87.3% European for knee OA); a near-unity intercept does not exclude ancestry-mismatch bias, the same limitation that applies to the main-text LDSC.

## Reproduction

`code/204_prep_oa_gwas_for_munge.py` → `code/205_run_ldsc_munge_rg.py --weights standard` (exploratory baseline, placeholder BMD N) → `code/206_summarise_ldsc_rg.py` → `code/207_bmd_verified_ldsc.py` (verified N + `--merge-alleles`, this table). Machine-readable values: `code/output/genetics/ldsc_rg_verified/verified_results.csv`; full comparison and boundary notes: `COMPARISON.md` in the same directory.
