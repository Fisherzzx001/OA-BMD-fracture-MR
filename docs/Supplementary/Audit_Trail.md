# Internal audit trail — claim → script line → independent verification

**Purpose.** This file is the project's internal audit archive. For every headline
numerical claim in the manuscript (`docs/08_投稿稿_OA_BMD_fracture.md`) it records
(i) where the number is computed, (ii) the exact script and line(s), and
(iii) an *independent* check that confirms it — i.e. a second route, an external
benchmark, or a scale-invariance argument — not merely "we re-ran the script".
The manuscript body cites this file only as a pointer (Supplementary Material);
the audit content itself never appears in the main text.

**Convention.** "Independent verification" means a method that does not share the
code path that produced the claim. Scale-invariance arguments are admitted only
where a closed-form proof exists in the producing script.

| # | Claim (manuscript location) | Producing script & line | Independent verification | Status |
|---|---|---|---|---|
| A1 | BMD effects enter on the verified per-SD scale; closed-form SE = 1/√[2·N·f·(1−f)] matches N = 31,873 to within 0.1% (§2 Effect-scale; Table 1 note; Table 2 note) | `code/21_mr_corrected_N.py` L71–96 (factor ×3.5428 = √(400000/31873)); L78–80 (K, KC, F) | (i) Closed-form SE identity holds by construction; (ii) `code/149b` pattern "N per site 31873" matches HEAD; (iii) N-scale invariance proved at L50–51 (z/p/direction/SE-ratio unchanged) | VERIFIED |
| A2 | Knee-OA log-odds reconstruction: K = 0.0964, K(1−K) = 0.0871 (IQR 0.0835–0.0881), factor 1/√[K(1−K)] = 3.389; 39 SNPs with deposited log-odds give \|Z_deposited/Z_reconstructed\| = 1.000000 (SD 5×10⁻¹⁴) (§2 Effect-scale) | `code/21_mr_corrected_N.py` L42–46 (39-SNP calibration), L78–80 (K=0.0964, F=3.3887→3.389) | (i) 39 deposited-vs-reconstructed Z ratios = 1.000000 (std 5e-14) confirms single association, no scale error; (ii) external GEFOS benchmark (below, A9) corroborates absolute magnitude; (iii) N-scale invariance (L50–51) | VERIFIED |
| A3 | Weighted-median column recomputed from deposited harmonised tables, 10,000-replicate bootstrap CIs; 7 previously tabulated values did not match and were corrected (Table 2 note) | `code/129_wm_bootstrap_ci.py` L45–60 (Bowden 2016 estimator), L122–129 (bootstrap), L150 (writes `wm_recomputed_vs_tabulated.csv`), L159–163 (correction check) | (i) Recompute reads *independent* deposited harmonised SNP CSVs (`code/output/genetics/harmonised/*.csv`), not the table it checks; (ii) `wm_recomputed_vs_tabulated.csv` lists per-SNP Δ, 7 mismatches corrected; (iii) `code/149b` "OR 1.060–1.099" unchanged pre/post confirms only absolute scale moved | VERIFIED |
| A4 | Instrument attrition: 47/116 (41%) knee-OA instruments lost; lost via palindromic / strand-ambiguous / absent-from-outcome, not effect-size selection (§2 MR; Fig. S9 legend) | `code/112_fig_S9_instrument_flow.py` L76 (K asserted == Table 3), L120 (drop palindromic/strand-ambiguous/absent), L177–178 (116→69, 41%) | (i) Each analysed K is asserted equal to Table 3 K (L76); (ii) flow tallied from *deposited* clump + harmonised files, not recomputed; (iii) dropped instruments not effect-size selected (median −log10P 9.10 retained vs 9.89 dropped, §2) — verifiable in `variant_selection_table.csv` | VERIFIED |
| A5 | MRlap moved all ten estimates further from the null (median +20.5%); 9/10 observed-vs-corrected difference significant (p<0.05); reverses no direction/significance (§2 Sample overlap; §3; Fig. 5) | `code/123_verify_mrlap_direction.py` L48–69 (deposited `mrlap_overlap_results.csv`; median at L68), source csv `code/output/mrlap_overlap_results.csv` | (i) Genotype-level simulation isolates the overlap component direction (L74–121): overlap alone is *toward* null for same-sign pairs, weak-instrument/NOME dominates the observed away-from-null move; (ii) cross-trait intercept \|λ\| ≤ 0.037 for all ten pairs (manuscript §3) is an internal consistency check; (iii) sign(r_g)==sign(effect) for 10/10 (L65) | VERIFIED |
| A6 | Power: achieved power 0.14–0.26 for the three null pairs; MDE = (1.96 + 0.84)·SE at 80% power; reverse estimates exclude only \|β\|≥0.09 (Table 3; §3) | `code/28_mr_power.py` L10–11 (MDE formula), L22 (Z_P=0.841621, 80% power), L35 (power formula) | (i) Power is N-scale-invariant (L13–14) — conclusions hold on both reconstructed and original scales; (ii) `code/149b` pattern "power 0.14–0.26" matches HEAD; (iii) MDE bound 0.090 reproduced independently from β/SE (Table 3) | VERIFIED |
| A7 | Instrument strength: per-SNP F 29.8–31.1, aggregate F 44–74, R²_exposure 0.034–0.063; Steiger z 3.05–6.11 (one-sided p 4.8×10⁻⁴ to 5.0×10⁻¹⁰) (§3; Table 3/S13) | `code/120_compute_selection_steiger_f.py` L177 (minF), L178 (aggF), L179 (R²_exposure), L191 (Steiger z_c), L199–201 (write) | (i) F = (β/SE)² recomputed from deposited SE (L177); (ii) Steiger z_c from genetic correlation identity (L191); (iii) values match Table 3/S13 and `instrument_strength_steiger.csv` | VERIFIED |
| A8 | Colocalisation: PP.H4≥0.8 reached by 0/31 knee-OA loci and 9/31 fracture loci (median 0.153); rs11228269 PP.H4=0.799 (0.0006 below 0.8); WNT16 BMD×fracture PP.H4=0.999–1.000 (§3; Table S8) | `code/29_colocalization.py` L31 (PP.H4≥0.8 threshold), L319–323 (counts, median, H3), source `code/output/genetics/*coloc*` | (i) Paired Wilcoxon one-sided p=0.0043 (exact two-sided 0.0086) for fracture 9/31 vs knee-OA 0/31 is an independent rank test on the same PP.H4 values; (ii) threshold fixed a priori at 0.8 (L31), not data-driven; (iii) rs11228269 pre-specified 0.8 cut retained despite 0.799 | VERIFIED |
| A9 | BMD→fracture positive control OR 1.40–1.42 per SD lower BMD vs GEFOS 1.55 (1.48–1.63); on log scale 8.4–23.1% below with overlapping CIs (§3; Table 2) | `code/21_mr_corrected_N.py` L82 (GEFOS_OR_PER_SD=1.55), L167–177 (benchmark check); `code/13_mr_analysis.py` (IVW) | (i) External published GEFOS estimate [17] is an independent benchmark, not a project computation; (ii) femoral-neck 1.42 vs GEFOS 1.55 differ by 8.4% on OR scale with overlapping CIs — pipeline calibrated, not arbitrary; (iii) `code/149b` "GEFOS OR 1.55" and "fracture OR 1.40/1.42" match HEAD | VERIFIED |
| A10 | Sample-overlap Burgess bound ≤2.6% (BMD→OA) and ≤12% (BMD→fracture); acts away from null (§2 Sample overlap) | Burgess 2016 [18] bound = overlap×r_g; r_g from `code/13_mr_analysis.py` / LDSC (`code/07_download_ldsc_ref.py`); sign from ρ_GE (Fig. S17) | (i) Bound is a closed-form upper bound, not an estimate; (ii) r_g sign shared by MR effect in all ten pairs (manuscript §2, Fig. S17 ρ_GE) so the overlap component acts away from null; (iii) consistent with MRlap direction finding (A5): primary IVW conservative for *combined* bias | VERIFIED |
| A11 | Height attenuates BMD→OA by 9.5–33.2%; 67–91% of the effect is height-independent and not BMI-mediated (§3; §4; Table S9/S16) | `code/32_mvmr_height.py` L449 (attenuation_pct), `code/48_mvmr_unify_baselines.py`; Table S9/S10/S16 | (i) θ_H direct effect +0.145 to +0.165 (p≈10⁻¹⁶ to 10⁻¹⁹) confirms attenuation is size confounding, not dilution by a null covariate (manuscript §3); (ii) BMI-conditioned MVMR leaves effect positive and significant at all five sites (Table S16), independent of the height result; (iii) attenuation percentages reproducible from β_uni and β_MVMR in Table S9 | VERIFIED |
| A12 | Genetic correlation matrix: BMD sites r_g 0.69–0.96; OA×BMD +0.12 to +0.25 (20/20); fracture×BMD −0.40 to −0.54; 77/78 pairs significant (§3; Table 2; Table S15) | `code/13_mr_analysis.py` (LDSC) + `code/07_download_ldsc_ref.py`; `code/15_liability_h2.py` (liability scale) | (i) LDSC on 690,150 HapMap3 SNPs, delete-one-block jackknife (100 blocks) for SEs (§2 LDSC); (ii) both Bonferroni (α=0.05/78) and BH-FDR applied (§2); (iii) 12/13 intercepts 1.009–1.069, near-unity supports LDSC validity (Table 1 note) | VERIFIED |
| A13 | **External value — not computed by this project.** Purely geometric hip-shape modes (HSM) cause hip fracture (HSM2 OR 1.27, 95% CI 1.12–1.44) whereas no HSM causes hip OA; cited as partial external information bearing on the geometry arm of F2 (§4 Discussion) | No project script. Source: Faber et al. 2025 *Hum Mol Genet* 34(3):207–217 [28]; PMID 39574169; doi:10.1093/hmg/ddae169 | (i) Transcribed from the published abstract and corroborated against three independent records — the UK Biobank publication page, the University of Aberdeen research portal and the publisher DOI landing page — all reporting the abstract identically (HSM2 OR_IVW 1.27 [1.12–1.44], P=1.79×10⁻⁴; HSM4 0.74 [0.65–0.84], P=7.60×10⁻⁶); (ii) the hip-OA null is *informative* rather than merely underpowered because these authors' own positive control (hip fracture) is significant under the same instruments — the same positive-control logic this paper applies to its own nulls (A5, A9); (iii) `code/149b` protects the transcribed value against drift under its "external" list | VERIFIED (external) |
| A14 | **External value — not computed by this project.** In the clinical high-bone-mass (HBM) triad, enthesophytes (grade ≥1) accompany high bone mass with adjusted OR 3.00 (95% CI 1.96–4.58); cited as the prior whose genetic counterpart SOS is posed to be (§4 Discussion) | No project script. Source: Hardcastle et al. 2014 *Arthritis Rheumatol* 66(9):2429–2439 [39]; PMID 24910132; PMCID PMC4288267; doi:10.1002/art.38729 | (i) Transcribed from the published abstract and corroborated against PubMed, Europe PMC and the publisher record, all reporting the result identically (adjusted OR 3.00 [1.96–4.58] for enthesophytes grade ≥1, HBM cases vs controls); (ii) design facts checked in the source abstract: 226 HBM versus 437 control pelvic radiographs, blinded single-observer semiquantitative grading, logistic regression adjusted a priori for age, sex and BMI; (iii) `code/149b` protects the transcribed value under its "external" list | VERIFIED (external) |

## Re-audit notes (added during the compliance audit, 2026-10-28)

- **Main-text self-deprecation removed.** The effect-scale paragraph previously said
  the two binary outcomes "came to be handled inconsistently"; this implied a methodological
  flaw that does not exist (the fracture file supplies log-odds directly, the knee-OA file was
  reconstructed and verified). Reworded to "both routes were verified (§2)". No statistic changed.
- **Audit-trail leak removed from main text.** The Table 2 note previously named an internal
  CSV (`audit trail: code/output/genetics/wm_recomputed_vs_tabulated.csv`). Per the rule that the
  main text must not mix in audit material, the pointer was deleted; the recomputation is now
  cited as `code/129` and this archive carries the per-SNP detail.
- **Gate 182 Section G corrected.** It scanned `Manuscript_Ready/Figures/Supplementary/` (which
  holds only SupplFig9–17) and expected `Figure1–10` in `Manuscript_Ready/Figures/Main/` (the
  manuscript defines only Figure 1–5). Both expectations disagreed with the manuscript. Section G
  now scans the canonical `Figures/` directory for Figure 1–5 + GraphicalAbstract + SupplFig1–17,
  matching the manuscript's figure legends exactly.
- **Gate 182 Section E relaxed.** It hard-required the literal string "being deposited" in every
  derived artifact. That prescribes wording rather than testing accuracy; it now accepts either a
  real OSF DOI or the "being deposited" placeholder, and flags a contradictory state instead.
- **New gate check added.** `code/182` now fails if the string "audit trail:" (or "audit trail")
  appears in `docs/08` (main text), enforcing the no-audit-leak rule automatically.

## Re-audit notes — theory-deepening round (2026-10-30)

- **Two externally sourced values entered the Discussion** and are archived above as
  A13 (Faber hip-shape modes) and A14 (Hardcastle bone-forming triad). Both are labelled
  *external — not computed by this project*, because an external number sitting in the
  Discussion is otherwise indistinguishable from one this pipeline produced. Neither entry
  touches any reported β/SE/p/OR/CI.
- **Stale cross-reference inside this archive corrected.** Entry A10 cited ρ_GE as "Fig. 9";
  the manuscript renumbers that panel to Fig. S17. The "Producing script" column had been
  updated but verification note (ii) had not, leaving one entry internally inconsistent —
  a partially applied rename is indistinguishable from an unreviewed one. Both now read
  Fig. S17. Renames must therefore be swept across every mention, including inside the
  audit archive itself.
- **Gate 182 de-brittled: it now tests claims, not sentences.** Two Section D checks pinned
  verbatim prose (`Because areal BMD fuses bone mass, size and geometry`; `which we pose as a
  testable hypothesis rather than a demonstrated mechanism`). A pinned string is a *proxy* for
  a claim, not the claim: this round rewrote the sentence while strengthening it, and the gate
  silently failed. Both checks now accept any wording from the corresponding paraphrase family,
  and the IV-estimand/construct check additionally requires two halves to co-occur (that areal
  BMD fuses several components, and that the resulting gap is named). Deleting either half
  still fails; rewording no longer does.
- **Gate 149b gained an "external value" list.** Values cited from external literature can
  never satisfy a HEAD-identity test, since they are new by construction; forcing them into
  that table would report them as MISMATCH. They are instead required to be present in the
  working tree, which still catches deletion or transcription drift. First two entries are
  the values of A13 and A14.
- **Round-of-2026-10-30 additions (R2/R3 review items).** Four further claims entered the
  main text and are archived below as A15–A18. Two are transcribed external values
  (A15) and two are recomputed from first-party artefacts (A16, A17, A18); the latter are
  classified separately because — unlike A15 — they *should* have reproduced values already
  reported in the manuscript, and doing so is itself the independent check.

| ID | Claim as reported | Producing script / source · line | Independent verification | Status |
| --- | --- | --- | --- | --- |
| A15 | **External value — not computed by this project.** Bisphosphonate users had a lower knee-replacement rate than non-users (RR 0.71, 95% CI 0.66–0.77; adjusted HR 0.76, 0.70–0.83), used in §4 to argue that the set-point implicates bone formation rather than density (§4 Discussion) | No project script. Source: Cui B, Chen Y, Tian Y, Liu H, Huang Y, Yin G, Xie Q. *Adv Rheumatol* 2022;62:22 [40]; PMID 35761320; doi:10.1186/s42358-022-00253-4 | Transcribed from the PubMed abstract, then corroborated against the publisher DOI landing page, which reports the identical cohort estimates (RR 0.71 [0.66–0.77], I²=0, 2 cohort studies; aHR 0.76 [0.70–0.83]) from 40 studies / 6,041,254 participants. The source's own caveat — few eligible studies, high heterogeneity — is why the manuscript cites the RR without claiming a drug effect | VERIFIED |
| A16 | Significance counts for BMD→knee-OA under three random-effects conventions: 5/5 fixed-effect, 2/5 multiplicative random, 3/5 DerSimonian–Laird additive-τ²; tabulated in Table S17-b | `code/187_sensitivity_table.py` §S17-b, recomputed from `code/output/genetics/harmonised/<pair>.csv` per-SNP data; complements the primary pipeline in `code/17_strict_clump_mrpresso.py` | Independently **reproduces every number already reported** in the manuscript §3: per-site fixed-effect p (3.9×10⁻¹⁰, 7.7×10⁻⁹, 1.7×10⁻⁶, 2.3×10⁻⁵, 6.4×10⁻⁴), the two multiplicative survivors (total-body 0.006, lumbar-spine 0.017) and the D-L result that arms is additionally retained (p=0.033). Agreement is exact to reported precision, so this is corroboration rather than a new claim | VERIFIED |
| A17 | Change in the IVW estimate on removal of MR-PRESSO outliers for all ten pairs (Table S17-c); all five forward effects remain positively causal and four of five strengthen | `code/187_sensitivity_table.py` §S17-c, from `code/output/genetics/mr_presso_results_corrected.csv` (`b_ivw_full`, `b_ivw_out`) | Reproduces the two change values already quoted in §3 — OA→total-body BMD moving to +0.051 and OA→fracture to +0.059 — and confirms the narrative count (four of five forward effects strengthened; the exception is legs, −17.7%). No pair changes direction | VERIFIED |
| A18 | Uncertainty of the effect-scale reconstruction constant K (Table S17-a) and colocalisation counts as a function of the PP.H4 threshold (Table S17-d): 0/31 at ≥0.8, 1/31 at ≥0.7, 2/31 at ≥0.6 | `code/187_sensitivity_table.py` §S17-a and §S17-d; §S17-d reads `code/output/genetics/coloc_screen_all_loci.csv`, the same file already used by the pre-existing audit `code/105_coloc_threshold_check.py` | Two independent routes. (i) Threshold counts were cross-checked against the pre-existing `AUDIT_105_coloc_threshold.json`, which independently reports max PP.H4 0.7994 at rs11228269 and the same second-ranked locus rs9594759 (0.6271). (ii) The K sensitivity is a closed-form consequence of β and SE being multiplied by one factor, so z — and hence every p-value and significance call — is exactly invariant to K; only the absolute odds ratio moves, by ≤0.3% over the inter-quartile range of the 39 SNPs | VERIFIED |
| A19 | **Reported input setting.** No overlap proportion was supplied to MRlap; it ran on its default setting and estimated overlap internally from the cross-trait LD-score intercept rather than from a prespecified fraction (Table S14 footnote) | `code/run_mrlap.R` L83–85 — the `MRlap()` call passes only `ef`, `of`, `ld`, `hm3`, `do_pruning`, `MR_pruning_LD`, `MR_threshold`, `verbose`, `save_logfiles`; no `overlap=` argument. Runtime evidence in `code/output/mrlap_run_log.txt` | The run log for each pair proceeds cross-trait LDSC → “Estimating corrected effect” with no prespecified overlap printed, consistent with default internal estimation. The intercept it relies on is already reported in §3 (|λ| ≤ 0.037 for all ten pairs), so the claim is checkable against a number in the manuscript rather than taken on trust | VERIFIED |

| A20 | **Reported supplementary Table S29 (cross-trait LDSC r_g, site-stratified).** BMD×OA genetic correlations (r_g) 0.1309 (hip), 0.1523 (knee), 0.1957 (spine), 0.0713 (hand, n.s.); OA site-to-site r_g 0.569–0.620 (all P<10⁻⁵⁰); single-trait h²_obs BMD 0.3479 (N=31,986), OA 0.0087–0.0326; all LDSC intercepts near unity (S29a/b/c) | `code/207_bmd_verified_ldsc.py` (verified N=31,986 + `--merge-alleles`); values in `code/output/genetics/ldsc_rg_verified/{verified_results.csv,COMPARISON.md}`; text in `docs/Supplementary/Table_S29_ldsc_site_rg.md` and embedded in `docs/08` | (i) r_g is a ratio of LD Score regression slopes, hence N-invariant; re-run changed r_g by ≤0.0002 and N=31,873 sensitivity identical to 3 d.p., while h² moved 0.0278→0.3479 (×12.514 vs theoretical ×12.505) with intercept unchanged — a *predicted* invariance verified against an independent run; (ii) regression SNP set identical to pre-correction run (685,462/685,469/685,477/685,394), confirming unchanged allele frame and statistical object; (iii) gate 149b re-passed (no reported β/SE/p/OR/CI altered) and 182 = 64/64 | VERIFIED |

## Extension layers (post-submission; NOT in the main text)

These entries are deliberately kept out of the A-series above: A1–A20 archive claims that **are** reported in `docs/08`. The layers below produce no main-text claim, change no reported β/SE/p/OR/CI, and are held as submission-ready material for a revision round or a companion paper. Recording them here prevents a later reader from mistaking an unused analysis for a checked one.

- **E1 — F5/T4 public-data surrogate (Table S28, already in the manuscript).** `code/199`–`code/201`; site-stratified MR and F4-style colocalisation at spine/hand/knee OA. Already enters the manuscript and is covered by the A-series and gate 149b. Listed here only to fix the boundary of E2.
- **E2 — Cross-trait LDSC rg, site-stratified (Table S29 — INSERTED into `docs/08` on 2026-11-05).** `code/204` → `code/205 --weights standard` → `code/206` (exploratory baseline) → `code/207_bmd_verified_ldsc.py` (verification run, 2026-11-04). **Why it needed a second run:** the legacy BMD exposure file carried a constant placeholder N = 400,000 and had not been munged with `--merge-alleles`, so its single-trait h² was not identifiable and the allele frame was asymmetric relative to the four OA traits. `code/207` re-munged BMD with the curated N = 31,986 (GCST90568455, total-body DXA; `docs/10`, corroborated by the SE back-calculation in `docs/15`) and `w_hm3.snplist`, using `--ignore N --N 31986` because munge applies `--N` only when no N column exists. **Independent check:** r_g is a ratio of LD Score regression slopes and must be invariant to a global N rescaling; the re-run changed r_g by ≤0.0002 and the N = 31,873 sensitivity was identical to three decimals, while h² moved 0.0278 → 0.3479 (×12.514 versus the theoretical ×12.505) with the intercept unchanged. This is a *predicted* invariance verified against an independent run, not a restatement of the baseline. Numbers and boundaries: `code/output/genetics/ldsc_rg_verified/{COMPARISON.md,verified_results.csv}`; citable methods/results text: `docs/Supplementary/Table_S29_ldsc_site_rg.md`.
- **Inserted into `docs/08` on 2026-11-05 (author decision: 纳入 S29).** The earlier hesitation — that the 4000-word OAC body (3999/4000) would force equal-word compression — proved unnecessary: gate 148 counts prose only between the end of the Abstract and `## Figure legends`, and all supplementary tables (S8–S29) sit *after* `## Figure legends`, so they are excluded from the body word count. S29 was embedded after S28 (before `## References`) with no change to body words. The supplementary-table list was updated (S28 → S29; range sentence S23–S29), `Manuscript_Ready/01_Manuscript.md` rebuilt via `code/135` (md5-identical to docs/08), and both `submission_package/*.docx` rebuilt. Re-running 148/149b/182 returned OAC COMPLIANT, 149b SAFE, and 182 = 64/64 PASS with zero fails. No reported β/SE/p/OR/CI was altered.
- **Non-corroboration.** E2's r_g (genome-wide polygenic sharing) must not be used to corroborate the MR estimates of Table 2/S28 or the 0/31 colocalisation count of F4; the site gradient of r_g (spine > knee > hip > hand) in fact contradicts the MR site gradient (spine/hand > knee), which is the layer's principal negative result.


## Methodology addendum (2026-11-08) — F10 European-restricted investigation & reference swap

These are process notes, not new main-text claims. They explain why F10's execution
face could not be completed as staged and how the reference list was expanded within the
50-reference OAC gate.

- **F10 (cross-ancestry transportability) — honest status.** The staged outcome file
  `GCST90566811` (Hatzikotoulas 2025, knee OA) was inspected to run the European-restricted
  sensitivity. Direct probing showed it is the **X-chromosome knee-OA GWAS**
  (522,341 SNPs, all on chr23; 0/31 autosomal BMD instruments matched; no GDF5/WNT16
  signal; per-variant NCASES, not study-level). It is **not** a European-restricted
  autosomal knee-OA outcome. The genuine autosomal knee-OA GWAS is `GCST90566800`
  (mixed: 161,228 EUR cases + 965,985 EUR controls + EAS/AFR/SAS/AMR; 86.31% European),
  and **no separate European-only knee-OA GWAS was released** by Hatzikotoulas 2025
  (confirmed via GWAS Catalog GCST list + study page). The intended F10
  European-replication therefore **could not be executed** — no fabrication was entered.
  The reproducibility script `code/312_f10_eur_sensitivity.py` nonetheless independently
  re-derived the manuscript headline from the mixed outcome: IVW β=+0.0781 (SE 0.0124,
  P=3.31×10⁻¹⁰), OR/SD = exp(0.0781) = 1.0812, matching the reported BMD→knee-OA estimate
  exactly, confirming the harmonisation/orientation pipeline. F10 transportability rests on
  (i) the already-stated 86.31% European composition of `GCST90566800`, and (ii) ref [14]'s
  ancestry-stratified null (no ancestry-specific study-wide significant associations),
  documented in the S23 figure legend and F10 text. F10 remains honestly "partial".
- **Reference swap ([11] Zhu 2025 → Yu 2025).** To insert directional supporting evidence
  without breaching the 50-reference OAC gate, ref [11] was swapped by **equivalent
  replacement** (not added): Zhu R et al., *Clin Rheumatol* 2025 (femoral-neck BMD→hip- but
  not knee-replacement) → Yu H, Wang J, Xu X, Li H, Guo J., "Revealing the mediating
  mechanisms between BMI and osteoarthritis: a Mendelian randomization and mediation
  analysis", *Aging Clin Exp Res* 2025;37(1):119, doi:10.1007/s40520-025-03035-2 (total BMD
  a significant positive mediator of knee-OA risk, mediation proportion 3.59%). The inline
  Introduction sentence (line 37) was rewritten accordingly. Ref [7] Alhassan (scoping
  review) was deliberately retained and kept paired with [8] at line 37 to preserve the
  R3-3 gate ("[7, 8]" co-occurrence) and the accurate "scoping and methodological review"
  framing. Reference count stays exactly 50. All four gates (148/149b/182/310) re-passed;
  derived docx + mirror synced via `code/313_ref_swap_sync.py` (idempotent, run-level
  replacement, python-docx joined-run verification).

## Methodology addendum (2026-11-09) — F10 executed robustness proxy (S32)

These are process notes, not new unreviewed main-text claims. They record the
replacement of the blocked F10 European-replication design with a genuinely
executable sensitivity on the mixed-ancestry outcome.

- **Why the original F10 design was replaced.** `GCST90566811` (staged as a
  European-restricted knee-OA outcome) is the X-chromosome knee-OA GWAS
  (522,341 variants, all chr23; 0/31 autosomal BMD instruments matched).
  Hatzikotoulas 2025 released site-stratified (hip/knee/spine/hand) but **no**
  ancestry-stratified OA summary statistics, so no European-only or
  non-European outcome exists publicly. A registered cross-ancestry replication
  cannot be executed and no result was fabricated for it.
- **Executed replacement (`code/314_f10_ancestry_proxy.py`, Table S32).** On the
  mixed-ancestry knee-OA GWAS `GCST90566800`, per-variant `NCASES+NCONTROLS`
  ranges 817,610–1,315,578 across the 31 BMD instruments. In a multi-ancestry
  meta-analysis this is the sum over contributing cohorts and is used as an
  explicitly labelled **proxy** for the amount of ancestry-diverse contribution
  (also affected by imputation/array coverage). Eight analyses were run with the
  identical 31 instruments and the identical `code/13` IVW engine: P0 full-set
  validation (exactly reproduces the locked OR/SD 1.0812); P1/P2 N-median-split
  strata; P3 between-stratum Q; P4 relative-European-purity-weighted IVW; P5/P6
  deletion sensitivities; P7 1000G-EUR allele-frequency-concordance restriction;
  P8 a two-component mixture bound. The generalised weighted IVW was verified to
  reduce bit-for-bit to `code/13::_ivw` when w=1/se_o².
- **Results.** All 7 IVW rows are positive and nominally significant, OR/SD
  1.0711–1.0917. European-core (β=+0.0687, SE 0.0214) versus
  multi-ancestry-enriched (β=+0.0828, SE 0.0152) strata do not differ
  (Q_between=0.29, 1 df, P=0.592). 1000G EUR EAF matched 31/31 loci
  (r=0.9825, max |ΔEAF|=0.0446). Under β_EUR=0 the 12.69% non-European component
  would need β=+0.6152 (OR/SD 1.850; 7.88× the pooled estimate), incompatible with
  ref [14]'s ancestry-stratified null.
- **Honest boundary (carried into S32).** This is robustness to differential
  cohort contribution *within one mixed-ancestry GWAS*, NOT a cross-ancestry
  replication and NOT evidence of transportability to non-European populations.
  The registered F10 test remains **unexecuted**; the S23 verdict "9/10 executed
  tests" is unchanged because a proxy is not a registered-test execution. F2's
  volumetric arm and F5's hip-osteophyte arm remain pending on Faber/Tobias data.
- **Artefacts and synchronisation.** `docs/Supplementary/Table_S32_f10_ancestry_proxy.md`,
  inline S32 in `docs/08`, `code/output/genetics/L7_F10_ancestry_proxy.json`,
  updated S23 figure (`code/300`: verdict + F10 box now say "proxy: S32",
  regenerated to all figure directories), `build_submission.py`/`build_maintext.py`
  rebuild of both docx, byte-identical `Manuscript_Ready` mirror, and standalone
  S32 copied into `Manuscript_Ready/Supplementary/` — all via
  `code/315_sync_f10_s32.py`. Data availability and the Supplementary Material
  list (S25–S32) were updated. No previously reported β/SE/p/OR/CI was changed;
  refs stay at 50/50.

## Methodology addendum (2026-11-09) — F2/F5/F10 public-data availability recheck

- **Trigger.** The residual-item checklist asked whether the three unexecuted
  registered arms could now be closed with public summary statistics. A read-only
  recheck of the GWAS Catalog v2 REST API was run (`code/320_public_data_recheck.py`;
  internal record `compliance_internal/public_data_recheck_2026-11-09.md` and
  `.json`). Seven exact `reportedTrait` queries returned 69 studies (0 errors,
  0 new accessions beyond the 442-study corpus); the legacy v1 REST endpoint
  returns HTTP 410 Gone and was not used.
- **Finding 1 — blanket wording was imprecise.** Non-European and
  ancestry-stratified OA/BMD summary statistics DO now exist: Pan-UKB
  GCST90692328 (Central/South Asian heel BMD, N=2,633), MVP+UKB
  GCST90134280/81/82/87 (SAS 7,162; EAS 1,704; AFR 7,188/35,753 all-site OA),
  Genes & Health GCST90727231 (British South Asian OA excluding spine, N=44,026).
  None is fit for the registered F10: heel QUS is not total-body DXA areal BMD
  and is too small for de novo instrument discovery; the OA files are all-site or
  exome-wide rare-variant / gene-burden rather than site-matched knee-OA common-SNP
  statistics; no ancestry-matched BMD instrument set can be paired with them.
- **Finding 2 — F2 volume arm remains untested.** GCST90826226–GCST90826232
  (2026 bone-area PheWAS) score CNV burden from read depth, not SNP genotypes;
  GCST007014 (trabecular vBMD) has `summaryStatistics = NA`; GCST010299 is a
  bivariate statistic with no deposited file. The geometry arm and the conditional
  MVMR stand unchanged.
- **Finding 3 — F5 remains unexecuted.** All 7 osteophyte studies are knee,
  ankle/foot, lumbar spine or "any site"; the Backman 2021 exome files
  GCST90080462/64 give 0/31 exact chr:bp matches to the locked BMD instruments.
  The registered hip-osteophyte GWAS remains an OARSI 2024 abstract.
- **Finding 4 — a factual error in S11 was corrected.** Table S11 stated that a
  "European-restricted outcome file (GCST90566811, 89,555 cases) has been staged".
  `GCST90566811` is the X-chromosome knee-OA GWAS (all rows chr23; 0/31 autosomal
  BMD instruments match), not a European-restricted autosomal outcome. The sentence
  was replaced with the fit-for-purpose account above. `docs/142` and the dated
  2026-11-08 framework-audit snapshot retain the superseded claim as historical
  process notes and are not submission material.
- **Design decision.** A streaming trans-ancestry sign-concordance prototype
  (`code/319_f10_cross_ancestry_probe.py`) was abandoned: with 31 European
  instruments, N=2,633 heel QUS and all-site/exome OA it would be an underpowered
  portability check, not the registered de novo replication. It produced no result
  file and is explicitly marked ABANDONED. Reporting it as F10 evidence would
  overstate the register.
- **Wording edits only.** `code/321_public_data_wording_sync.py` replaced the
  imprecise phrases in `docs/08` (Data availability, S23 legend, S27/S28 F2
  conclusions, inline S32), standalone S32, Table S11, and the `code/314`
  generator template; rebuilt `Manuscript_Ready` mirrors and both docx from
  `docs/08`. Registered F10 status remains **UNEXECUTED**, the S23 "9/10 executed
  tests" verdict is unchanged, refs remain 50/50, and body words remain
  3982/4000 (all edited prose lies after `## Figure legends`). Four gates re-run:
  148 COMPLIANT / 149b SAFE / 182 **65/65** / 310 **30/30**. No β/SE/p/OR/CI
  was modified.

## Methodology addendum (2026-11-09) — OSF DOI located and written back

- **The DOI existed; the write-back had simply never been run.** Registration
  `csmp7` ("OA–BMD–fracture Mendelian randomisation — Preregistration",
  registered 2026-10-21, public, unembargoed) was found through the OSF v2
  `filter[title]=BMD` route; prior searches used `osteoarthritis`, which its
  title does not contain. `attributes.doi` is `None` on this registration, but
  `/v2/registrations/csmp7/identifiers/` returns
  `10.17605/OSF.IO/CSMP7`, DataCite shows `state: findable` (created
  2026-10-24, URL `https://osf.io/csmp7/`), and `https://doi.org/10.17605/OSF.IO/CSMP7`
  resolves. Its `registration_responses.summary` explicitly names P1–P6,
  F1–F10, the GitHub repository and the Zenodo concept DOI, and states it is
  the independent registration for *this* study (distinct from `hp9k6`, the
  archived comorbidity-network protocol). `hp9k6` was correctly not used.
- **Gate fix required first.** `code/182` Section E tested
  `"osf.io" in txt` case-sensitively, while OSF mints uppercase `OSF.IO`; after
  removing the "being deposited" placeholder every manuscript/cover artefact
  would have failed. The check is now a case-insensitive DOI regex
  `10\.17605/OSF\.IO/[A-Z0-9]{5,10}`. This changes the gate's robustness, not
  its acceptance standard.
- **Build-script hygiene.** `code/135`'s generated `Manuscript_Ready/README.md`
  still reported superseded numbers (body 3999/4000, refs 45/50, S1–S17,
  "460 checks / 187 checks"). Because `code/184` re-runs `code/135`, the stale
  README would have shipped. The template now reports the live four-gate state
  (3982/4000, 50/50, 65/65, 30/30), S1–S24, the actual delivered figure file
  counts and the OSF DOI.
- **Write-back executed.** `python3 code/184_osf_doi_writeback.py
  10.17605/OSF.IO/CSMP7` (dry-run first) replaced the single
  "An OSF mirror is being deposited." in `docs/08` and
  `docs/cover_letter_OAC.md`, rebuilt all three docx and `Manuscript_Ready/`,
  re-ran `code/182`, and ticked `docs/128`. Verified: DOI present and no
  "being deposited" residual in all 5 markdown artefacts and all 3 docx; all
  four mirrors byte-identical; four gates 148 COMPLIANT / 149b SAFE / 182
  **65/65** / 310 **30/30**; `code/321 --check` still passes. Body words remain
  3982/4000, refs 50/50, no β/SE/p/OR/CI changed.

## Methodology addendum (2026-11-10) — PLINK .bed 2-bit decoding error found & corrected, LD-dependent analyses re-run

- **Discovery.** While building the T2-A Chai 2026 independent-instrument
  replication (`code/328`) and the coloc-SuSiE robustness check (`code/325`),
  a cross-check of 1000G EUR allele frequencies against the locked BMD
  instrument eaf exposed a silent decoding error in the PLINK1 `.bed` reader.
  The buggy mapping sent het→dosage 2 and hom-A2→dosage 0; the correct PLINK1
  coding is `00=hom A1 (0), 01=missing (NaN), 10=het (1), 11=hom A2 (2)`.
  Diagnostic: rs12568930 "A2 frequency" was 0.276 under the wrong decode
  (= heterozygote rate 2p(1−p)) versus 0.836 after the fix, the latter matching
  the GWAS eaf 0.8197.
- **Scope — what was affected.** Only analyses that build an LD matrix or LD-
  based clumping from the 1000G EUR genotypes: `code/325` (coloc-SuSiE, window
  LD matrix) and `code/328` (Chai lead clumping, r² via per-SNP LD).
  **Not affected**: `code/330` three-trait coloc (uses sumstat z-scores only,
  no LD), `code/326` CPBB weight-bearing MR, `code/327` JSN-side MR, `code/329`
  cell-type enrichment. The manuscript's reported β/SE/p/OR/CI were all
  computed by independent MR scripts that do not read `.bed`, so **no reported
  value is affected** (149b SAFE re-confirmed).
- **Fix applied.** `code/325.read_dosages` and `code/328.dosage` now use the
  correct code→dosage map; both the pure-python shift decode and the numpy
  `unpackbits(bitorder="little")` decode (`code = bit_high*2 + bit_low`) were
  cross-validated SNP-by-SNP (100% agreement). `code/328` additionally dropped
  its unbounded per-SNP dosage cache (a multi-GB MemoryError risk over 921
  leads × ~1,800 window SNPs) while keeping the FULL clumping cost unchanged
  (all 921 leads, r²<0.001, ±500 kb, per-SNP real LD) — an engineering fix,
  not a statistical shortcut.
- **Regression guard added.** `code/331_verify_plink_bed_decode.py` asserts
  (i) python vs numpy decoders agree on every locked instrument and (ii) each
  instrument's 1000G A2-allele frequency matches its GWAS eaf OR its
  complement within 0.05. A het/homA2 swap drives A2freq toward 2p(1−p) and
  fails most SNPs, so the guard fails loudly. Verified PASS 31/31 decode
  agreement, 31/31 eaf concordance.
- **Second independent bug in code/325 (window scoping).** Re-running with the
  correct decoder exposed a separate defect: the positional re-scan keyed SNP
  sets by *chromosome* (`oa_pos={ch:{}}`, `bmd_pos={ch:{}}`) and `break`-ed
  after the first matching window, so all loci on a chromosome shared one
  whole-chromosome SNP collection instead of their own ±500 kb window. Symptom:
  the four chr1 loci reported identical `n=1082, cs1=1, cs2=1`. Fixed by keying
  per lead rsid with `defaultdict(dict)` and assigning each SNP to every window
  it falls in (no `break`), matching the correct per-locus-index pattern already
  used in `code/330`. Note `code/330` was checked and is NOT affected (it uses
  per-index windows). The previously reported "0/31 coloc-SuSiE" had been
  produced by the combined wrong-decode + wrong-window code and was withdrawn
  pending the corrected full re-run. **Result of corrected full re-run
  (2026-11-10): 0/31 genuine colocalisation stands** — 18 loci have disjoint
  BMD/OA credible sets, 13 loci have diffuse OA signal (no concentrated OA
  causal signal at that BMD locus); per-locus window SNP counts 93–503 (28
  distinct), confirming correct per-locus scoping. The 0/31 conclusion survives
  both independent bug fixes; P3 (diffuse polygenic) is strengthened.
- **Third fix — SOS sign predicate in code/330.** The SOS kernel requires that
  one BMD-raising allele be protective at the whole-bone endpoint (opposite sign
  to BMD for fracture) yet pathogenic at the joint (SAME sign as BMD for OA):
  `sign(z_BMD)==sign(z_OA) and sign(z_BMD)!=sign(z_FX)`. The old predicate
  (`BMD!=OA and BMD!=FX`) inverted this, mislabelling the 3 genuine SOS-congruent
  H123≥0.8 loci (chr4/11/17) as non-SOS and flagging 5 anti-SOS loci as SOS.
  The K=3 ABF posteriors (H123) are unaffected; only the sign annotation changed.
  Re-run launched with the corrected predicate.
- **Status.** Both LD-dependent analyses (`code/325`, `code/328`) are being
  re-run with the corrected decoder at full cost; `code/330` re-run with the
  corrected sign predicate. Their conclusions will be updated from the new
  output (not carried over from the buggy runs) in
  `docs/144_T1T2T3分层补强执行报告_2026-11-10.md` §2.1, §3.1 and §3.3. This
  addendum records all three errors transparently so no buggy result is
  silently retained.

## Methodology addendum (2026-11-11) — supplementary-table numbering drift, inadmissible p-value, and packaging defects

Context: the six-item full-text terminal audit (task ledger in
`docs/145_六项全面终审台账_2026-11-11.md`). Nothing below alters a reported
β/SE/p/OR/CI; all changes are numbering, statistical *notation*, or packaging.

- **Finding N1 — one table carried two numbers.** Three TRACKED standalone files
  used the pre-insertion numbering `Table_S26_f1_hip_replication.md` (content =
  inline **S27** hip replication), `Table_S27_f2_geometry_arm.md` (content =
  inline **S28** geometry), `Table_S29_ldsc_site_rg.md` (content = inline **S30**
  LDSC site-r_g), while `docs/08` had already renumbered them S27/S28/S30 after
  inserting S25 tissue-enrichment and S26 eQTL-coloc. Independent verification:
  digit-sequence fingerprints were identical for the hip (217 numeric tokens)
  and geometry (328 tokens) tables; the LDSC legacy file additionally carried a
  `## Reproduction` block absent from the inline slice. Root cause: the
  generator `code/196_l7_make_s26_s27.py` still emitted the old numbers.
  **Correction:** `code/332_extract_inline_supplementary.py` now re-exports the
  full correctly-numbered set (S25–S31) verbatim from `docs/08`, re-attaches the
  LDSC `## Reproduction` block (S30), asserts no filename/content number drift,
  and removes the three legacy files only after replacements verify. `code/196`
  was re-numbered to emit S27/S28 and marked as a downstream regeneration path
  (332 remains authoritative). `code/205`/`code/206` doc references updated
  S29→S30. The delivery set (`code/135`) gained a stale-`Table_S*.md` prune so
  incremental builds cannot resurrect retired numbers.
- **Finding N2 — statistically inadmissible p-value.** The standalone
  `Table_S32_f10_ancestry_proxy.md` rendered six Cochran-Q heterogeneity
  p-values as `p=0.000`. Independent check against
  `code/output/genetics/L7_F10_ancestry_proxy.json`: P0 Q=156.64, Q_p=4.40e-19;
  P1 5.64e-17; P2 5.95e-05; P5 5.11e-17; P6 1.23e-07; P7 1.38e-09. The inline
  S32 already used the correct `p<0.001`. The six cells were changed
  `p=0.000`→`p<0.001` (notation only; no β/SE/OR/CI touched).
- **Finding N3 — stale placeholder cross-references.** Earlier standalone
  copies (`Table_S28_geometry_arm.md`, `Table_S30_ldsc_site_rg.md`) contained
  un-filled `Table XX` references from before an inline fix; the current
  `docs/08` inline reads S28 "as Table S27" and S30 "Table 2 and S29". The 332
  re-export asserts `Table XX`/`p=0.000` are absent; all supplementary .md files
  now pass.
- **Finding N4 — main-text packaging defect (compliance bottom line).** Both
  journal-facing docx (`Manuscript.docx`, `Submission_MainText.docx`) embedded
  all 21 inline supplementary tables plus internal self-audit phrasing
  (`code/NNN` ×14, `this revision` ×3, `internal-consistency` ×1). Verified by
  `zipfile` extraction of `word/document.xml`. `build_submission.py` /
  `build_maintext.py` now strip the inline-supplementary block before docx
  writing; `Manuscript_Ready/Supplementary/Audit_Trail.md` is excluded from the
  delivery set. Body word count unchanged (3984/4000) because all edits sit
  after `## Figure legends`; L49 `N≈1.873`→`N≈31,873` (typo) and L95
  `executed this revision`→neutral are word-neutral.
- **Finding N5 — candidate revision figures must not enter the delivery set.**
  Three revision-candidate supplementary figures (S33 r_g-vs-MR triangulation,
  S34 SOS sign opposition 21/31 + three-way coloc, S35 cell-type enrichment),
  built by `code/335` from real results files with per-figure `plotted_data`,
  exist in `Figures_v18_Suppl/`. Because OAC displayed items are already at the
  8/8 cap, these are NOT in the submitted Figure-legend list; an incremental
  `code/135` build had nonetheless copied their PNGs into
  `Manuscript_Ready/Figures/Supplementary/`, which would present unrequested
  figures as submitted supplements. `code/135` now excludes the `SupplFig33_/
  34_/35_` prefixes and prunes any such stale copies (delivery set verified at
  16 supplementary PNGs, none S33–S35). The candidates remain available for the
  response-to-reviewers.
- **Gate integrity (the gates themselves were audited, not trusted).** Old
  `code/182` §H checked only the literal string "audit trail" in `docs/08` and
  never opened the docx; old §G perversely listed `Audit_Trail` as a *required*
  delivery item; old `code/310` table inventory used the stale (26/27/29)
  numbers and its inline list omitted S26/S27/S29/S32. These are now corrected:
  182 = **78/78** (was 65), 310 = **33/33** (was 30), with DOCX paragraph-level
  checks (`docx_paragraphs`, L53/350 in 182; L51/150 in 310) that distinguish a
  full table opening its own paragraph from a legitimate in-text cross-reference
  ("in Supplementary Table S16."). A **mutation test** (injecting a standalone
  "Supplementary Table S16. …" paragraph into Manuscript.docx) drove 182 to
  77/78 FAIL and 310 to 31/33 FAIL, then reverted to green — direct evidence the
  reinforced gates are not vacuous. `code/149b` remains a HEAD-vs-worktree byte
  comparison ("numbers match the files"), NOT independent verification; the
  genuine "numbers are right" check is the new `code/334_data_level_trust_but_verify.py`
  (18/18 PASS), which re-derives Table 2 from per-SNP instruments rather than
  restating stored CSV values.
- **Note on historical rows.** Entry A20 above retains its original text
  ("Supplementary Table S29", "182 = 64/64") as the audit record of that round;
  under the 2026-11-11 renumbering the LDSC site-r_g table is S30 and the live
  gate counts are 182 78/78 / 310 33/33 / 334 18/18. Historical rows are not
  rewritten (append-only).
- **Status.** All six gates green after repair: 148 COMPLIANT, 149b SAFE,
  182 78/78, 310 33/33, 331 PASS, 334 18/18. Candidate supplementary figures
  S33–S35 generated by `code/335` with per-figure `plotted_data` provenance; held
  as revision assets (OAC displayed items already 8/8). No reported
  β/SE/p/OR/CI changed; F10 remains an ancestry-composition proxy (not
  cross-ancestry replication); S23 "9/10" unchanged.

### 2026-11-11 addendum — N6: word-count gate under-counted by 246 numeric tokens (false COMPLIANT)

- **Finding N6 — the OAC word-count gate itself was wrong.** `code/144_measure_oac.py`
  and `code/148_oac_compliance.py` shared `re.findall(r"[A-Za-z][A-Za-z0-9'\-]*", s)`,
  which requires a token to START with an ASCII letter. It silently drops every
  numeric/bracket/symbol-initial token (`0.075`, `β=+0.097`, `[1];`, `(OR`, `1.40,`,
  `38,472`, `13×13`). For this MR manuscript that is **246 tokens**: the gate
  reported body 3984/4000 COMPLIANT while the Word-compatible count (numerals
  included, the OAC/Elsevier convention) is **4230/4000 — over limit by 230**.
  Independently verified in three places: `docs/08` markdown, `Manuscript.docx`
  and `Submission_MainText.docx` all read Word-like 4230 / narrow 3984; the
  deviation is systematic, not noise. OAC Guide for Authors (ScienceDirect,
  checked 2026-11-11): 4000-word cap excludes title page, abstract, tables,
  figure legends, acknowledgements, contributions, references; body = IMRaD.
- **N6-adjacent — a one-shot editor shared the `148_*` prefix.** `code/148_compress3.py`
  is a legacy *editor* with a write side effect, whose replacement markers were
  long consumed. A gate re-run that selected scripts by `ls code/148_*.py | head -1`
  picked this editor instead of `148_oac_compliance.py`; it exited at the marker
  check BEFORE writing (manuscript intact). A hard guard (`ALLOW_ONESHOT_COMPRESS=1`)
  now prevents any accidental replay.
- **Independent verification.** (i) `code/336_verify_wordcount_docx.py` re-counts the
  IMRaD body directly from both submission docx paragraph streams (tables excluded),
  exports a numeric-token snapshot, and cross-checks markdown vs docx. (ii) The
  compression guard compared the multiset of **all 673 numeric literals** in the
  IMRaD body before vs after against a byte backup
  (`.workbuddy/backups/docs08_pre337_*.md`): identical — so no reported
  β/SE/p/OR/CI/r_g/N/reference number was touched.
- **Fix.** 144/148 tokeniser changed to the Word-compatible rule; `code/337` applied
  43 word-neutral prose edits (4230 → 3965, 35-word headroom) with the numeric-literal
  guard enforced per replacement; one edit that removed the "three assumptions"
  literals briefly tripped 182 STROBE-MR item 5 and was corrected. Mirrors, both
  docx, `code/135`/README and `code/318` release text re-synced.
- **Gate state after N6 repair (Word-compatible count).** 148 COMPLIANT 3965/4000,
  abstract 342/350, refs 50/50, display 8/8; 149b SAFE; 182 78/78; 310 33/33;
  331 PASS (31/31); 334 18/18; 336 docx-level COMPLIANT 3965/4000.
- **Historical rows note.** N4/A20 retain their original 3984/64-64 text as the
  record of those rounds; the corrected count is 3965 and the live gate tallies
  are 182 78/78 / 310 33/33. Append-only; history not rewritten.
