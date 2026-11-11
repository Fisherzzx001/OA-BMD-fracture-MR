<!-- Standalone supplementary table extracted verbatim from
     docs/08_投稿稿_OA_BMD_fracture.md by code/332_extract_inline_supplementary.py.
     Source slice sha256: 26923b5a762cde3fb4d12bf66890aaf522d2b779ba32f6336ecc66488e9e4e80 . No statistic was rewritten. -->

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

## Reproduction

`code/204_prep_oa_gwas_for_munge.py` → `code/205_run_ldsc_munge_rg.py --weights standard` (exploratory baseline, placeholder BMD N) → `code/206_summarise_ldsc_rg.py` → `code/207_bmd_verified_ldsc.py` (verified N + `--merge-alleles`, this table). Machine-readable values: `code/output/genetics/ldsc_rg_verified/verified_results.csv`; full comparison and boundary notes: `COMPARISON.md` in the same directory.
