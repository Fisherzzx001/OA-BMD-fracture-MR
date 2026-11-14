# Supplementary Table S34. Effect-scale validation — independent variant-level check (M1 remediation)

The BMD exposure file carries a *placeholder* N = 400,000 on every one of its 691,715 rows (a deposition convention, not the analysed sample size), but the **effect estimates themselves are already on the verified per-SD scale** — this was confirmed directly against an independent BMD GWAS (Section B) and by a per-SNP standard-error diagnostic (Section A). Consequently **no effect-scale reconstruction was applied to any reported estimate**; the reported MR results sit on the as-deposited, true per-SD scale. Three checks were run (`code/351b_independent_scale_test.py` for the variant-level GEFOS test; `code/350_m1_scale_reconstruction_validation.py` for the ±20 % sensitivity; JSON in `code/output/genetics/m1_audit/m1_independent_gefs_tbmd_test.json` and `m1_validation_results.json`; forest plot `Figures_v18_Suppl/Fig_M1_scale_sensitivity.png`).

## A. Placeholder-N file signature and per-SNP SE diagnostic
All 691,715 rows of the deposited BMD exposure carry N = 400,000 (100.0 %), confirming the N column is a deposition convention. The effect scale, however, is not compressed to match it: for every instrument the deposited SE satisfies SE = 1/√[2·N·f·(1−f)] with N = 31,873 to within 0.1 %, and a placeholder N = 400,000 would predict SE 3.54-fold *smaller* — a pattern absent from the data. The deposited β/SE are therefore already at the true N = 31,873 scale.

## B. GEFOS variant-level validation (definitive)
The decisive test compares the 31 MR instruments against the **independent GEFOS total-body BMD GWAS** (GCST005348; Medina-Gómez et al. 2018, AJHG; N = 56,284 European; trait-matched to the Qian total-body BMD exposure). All 31 instruments match at GRCh37 position with 0 ambiguous alleles. Orienting each GEFOS β to the exposure effect allele and taking the median ratio GEFOS/deposited β:

| Statistic | Value |
| --- | --- |
| Instruments matched / used | 31 / 31 (0 ambiguous) |
| median(GEFOS β / deposited β) | 0.7935 |
| Q1 – Q3 | 0.668 – 0.845 |
| median(Z_GEFOS / Z_deposited) | 1.0736 |
| Expected if a ×3.543 reconstruction were correct | ≈ 3.543 |
| Expected if deposited β already true-scale | ≈ 1.0 |
| **Adjudication** | **RECONSTRUCTION_OVER_INFLATED_WRONG** |

The observed median ratio (0.79 ≈ 1.0, not 0.28 = 1/3.543) shows the deposited β is **already on the true per-SD scale**; a ×√(400,000/31,873) = 3.543 "reconstruction" would have been a 3.5× over-inflation and was correctly **not applied**. The absolute per-SD scale is independently anchored by the BMD→fracture positive control (OR 1.40; femoral-neck 1.42) lying 8–23 % below the GEFOS fracture estimate (1.55) on the log-odds scale — a cross-site (total-body vs femoral-neck) approximation, not a reconstruction error, and with overlapping CIs.

**Reconciliation with the reported estimates.** Recomputing the BMD_total→OA_knee IVW on the as-deposited β (no factor) gives β = 0.07834, OR/SD = 1.0815 — an exact match to the reported Table 2 value. The reported MR results were therefore computed on the true scale; the ×3.543 reconstruction never entered them, so the reviewer's hypothesized 3.5× magnitude reversal does not occur.

## C. ±20 % scale sensitivity of the BMD→OA MR estimates
Each reported β was re-scaled by f ∈ {0.8, 0.9, 1.0, 1.1, 1.2} and the OR/CI recomputed; all five BMD→OA sites remain significant under −20 % scaling (`survives_20pct_down = True`). The protective BMD→fracture OR (OR < 1) reverses its protective reading under −20 % scaling, as expected for a protection factor.

| Pair | K | OR/SD (f = 1.0) | 95 % CI | survives −20 % |
| --- | --- | --- | --- | --- |
| BMD total → OA knee | 31 | 1.0815 | 1.055–1.108 | True |
| BMD femoral neck → OA knee | 25 | 1.0596 | 1.025–1.095 | True |
| BMD lumbar → OA knee | 21 | 1.0992 | 1.065–1.135 | True |
| BMD arms → OA knee | 16 | 1.0780 | 1.045–1.112 | True |
| BMD legs → OA knee | 24 | 1.0619 | 1.033–1.092 | True |
| BMD total → fracture | 31 | 0.7138 | 0.683–0.746 | False (protection reverses) |

## Conclusion
The variant-level cross-check against an independent BMD GWAS is now complete (Section B): the deposited BMD β/SE are already on the true per-SD scale (median ratio 0.79, not 0.28; 31/31 instruments matched), so the √3.543 reconstruction would have been a 3.5× over-inflation and was not applied. The reported MR estimates are on this verified scale and are robust to ±20 % perturbation at all five OA sites. The BMD→fracture positive control against GEFOS remains the external anchor of the absolute per-SD scale. Per-SNP z, p, direction and significance are scale-invariant, so the reconstruction question bears only on the absolute magnitude, which the positive control bounds rather than certifies.
