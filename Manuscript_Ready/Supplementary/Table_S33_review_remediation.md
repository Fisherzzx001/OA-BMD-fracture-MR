# Supplementary Table S33. Section-2 review remediation audit

Companion numeric record for the P1-2/3/4 and P2-2 remediations. All values are
recomputed by read-only scripts `code/346`, `code/347` and `code/348` from
first-party artefacts; no reported beta, SE, p-value or odds ratio is altered.

## Part A — Burgess overlap-bias denominators (BMD↔fracture)

22,330 UKB participants are shared (GEFOS fracture; 1,832 cases + 20,498 controls).
BMD=31,873 (100% UKB); fracture=264,973; r_g(BMD_total,fracture)=−0.4868,
r_g(BMD_femur_neck,fracture)=−0.5086. Bias ≈ overlap fraction × |r_g|.

| Overlap denominator | Fraction | Bias, total-body | Bias, femoral-neck |
| --- | ---: | ---: | ---: |
| 8.4% of fracture sample | 0.084 | 4.10% | 4.29% |
| Geometric mean √(0.701×0.084) | 0.243 | 11.83% | 12.36% |
| 70% of BMD sample (exposure-referenced) | 0.701 | 34.10% | 35.64% |

**Label correction.** The prior S14 note wrote "~12% on the 70%-of-BMD basis";
12% is in fact the geometric-mean-overlap figure (24.3% × |r_g|). Literal
70%-of-BMD × |r_g| is ~34–36%. Burgess's exact formula (olsbias × overlap.prop ×
1/F) carries a 1/F attenuation that the overlap × r_g shortcut omits, so all
three entries are conservative ceilings and even 34–36% is a worst-case bound.
Reproduce: `python3 code/346_p1p2_overlap_dual_denominator.py --apply`.

## Part B — MRlap-corrected positive control vs GEFOS

MRlap re-clumps to 48 (total-body) / 34 (femoral-neck) instruments on its own
exposure scale, so its raw corrected β is NOT the Table-2 OR. We transfer its
scale-free relative displacement corrected/observed onto the Table-2 IVW β
(assumes the relative overlap/weak-instrument displacement is approximately
instrument-set invariant) — a sensitivity, not a substitute point estimate.

| Control | MRlap observed→corrected | Relative shift | Table-2 IVW OR | Transferred OR | Gap to GEFOS 1.55, log-odds |
| --- | --- | ---: | ---: | ---: | ---: |
| BMD total→fracture | −0.1222→−0.1409 | +15.3% | 1.401 | 1.475 | 23.1% → 11.4% |
| BMD femoral neck→fracture | −0.1032→−0.1232 | +19.4% | 1.420 | 1.521 | 19.9% → 4.4% |

Under ~70% sharing of the BMD sample the BMD→fracture control is not a strictly
two-sample comparison. MRlap moves the control toward GEFOS, so residual
shortfall is not inflated by overlap; the reported primary OR remains 1.40/1.42.

## Part C — Ancestry mismatch: EUR-only replication (P1-3)

A EUR-only subset of the *same* multi-ancestry knee-OA GWAS (GCST90566800) is not
publicly reconstructible: the file has per-SNP NCASES/NCONTROLS heterogeneity and
no per-SNP ancestry label. We therefore use an independent EUR-only knee-OA GWAS
(FinnGen; same 1000G Phase-3 EUR LD reference) as a directional replication, and
report the GO↔FinnGen knee-OA r_g to show definition/ancestry heterogeneity.

| BMD site | GO multi-ancestry r_g | FinnGen EUR-only r_g | Sign agrees |
| --- | ---: | ---: | :---: |
| Arms | +0.163 | +0.073 | yes |
| Femoral neck | +0.173 | +0.051 | yes |
| Femur total | +0.153 | +0.069 | yes |
| Lumbar spine | +0.143 | +0.027 | yes |
| Legs | +0.188 | +0.069 | yes |
| Total body | +0.154 | +0.054 | yes |

GO vs FinnGen knee-OA r_g=+0.137 (SE 0.016; intercept 0.258), far below 1, so
the two knee-OA definitions are genetically heterogeneous and magnitude gaps
cannot be attributed to ancestry alone. All five MR-entered BMD sites retain a
positive r_g in the EUR-only data, so ancestry mismatch is unlikely to have
created the *sign* of the BMD–OA correlation; it may still affect magnitude, so
the stratification limitation stands. Liability h² sensitivity: OA_knee
h²_liab=0.106 at K=0.096 (GWAS case-fraction) vs 0.094 at K=0.12 (population
prevalence). Reproduce: `python3 code/347_p1p3_ancestry_ldsc_sensitivity.py --apply`.

## Part D — Reverse-MR instrument-set attrition and Egger flip (P2-2)

OA_knee→BMD_total: the pre-reconstruction fit (placeholder N=400,000 scale,
before the knee-OA log-odds rescaling) used 79 harmonised instruments with a
nominally significant MR-Egger intercept (p=0.0392); after reconstruction
(×1/√[K(1−K)]=3.389) and re-harmonisation the delivered authoritative table
(`mr_results_corrected.csv`) uses 69 instruments with a null intercept (p=0.376).
The 116→69 knee-OA attrition is reported in Fig. S9; the additional 79→69
contraction belongs to the earlier uncorrected-scale pipeline, excluded from the
deliverable under the single-authoritative-table rule (P0-2). The 69-instrument
corrected set is therefore the reported analysis, and the intercept flip shows
the directional-pleiotropy diagnosis is instrument-set sensitive: we do not read
the reverse MR as evidence either for or against direction. Raw-scale values are
retained only in the internal audit trail, not distributed.
