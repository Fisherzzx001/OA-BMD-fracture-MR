<!-- Standalone supplementary table S35. Source: docs/08_投稿稿_OA_BMD_fracture.md
     (## Supplementary Table S35., lines 508–520), maintained by hand and delivered
     as an independent file because build_submission strips the inline block from the
     docx. Reproducible via code/run_eur_only_mr.py. No statistic rewritten. -->

## Supplementary Table S35. European-only knee-osteoarthritis MR replication of BMD → knee OA (code/run_eur_only_mr.py)

**Design.** An ancestry-restricted replication of the primary BMD→knee-OA MR. The European-only BMD instruments (Qian et al. 2025) were taken against an independent European-only knee-OA GWAS (FinnGen), harmonised to the reference A1 allele and analysed by inverse-variance-weighted (IVW), MR-Egger and weighted-median MR. This is a within-European sensitivity of the primary effect, not a cross-ancestry replication; no value reported in the main text is altered.

| BMD exposure site | n | IVW β (SE) | OR per SD (95% CI) | p | MR-Egger intercept p | Weighted-median β |
| --- | --- | --- | --- | --- | --- | --- |
| Total body | 13 | +0.0651 (0.0268) | 1.067 (1.013–1.125) | 1.5×10⁻² | 0.40 | +0.0652 |
| Femoral neck | 12 | +0.0754 (0.0352) | 1.078 (1.006–1.155) | 3.3×10⁻² | 0.82 | +0.1217 |
| Lumbar spine | 9 | +0.1103 (0.0381) | 1.117 (1.036–1.203) | 3.8×10⁻³ | 0.88 | +0.1041 |
| Arms | 6 | +0.0883 (0.0320) | 1.092 (1.026–1.163) | 5.8×10⁻³ | 0.30 | +0.0611 |
| Legs | 10 | +0.0606 (0.0295) | 1.063 (1.003–1.126) | 4.0×10⁻² | 3.5×10⁻³ | +0.0621 |

**Interpretation.** All five BMD sites show a positive, nominally significant effect of higher BMD on knee-OA risk in the European-only analysis (OR/SD 1.06–1.12), matching the direction and magnitude of the multi-ancestry GO knee-OA result (OR/SD 1.060–1.099). The replication argues against ancestry mismatch as the source of the positive BMD→knee-OA signal. MR-Egger intercepts are non-significant at four of five sites; the legs intercept (p=3.5×10⁻³) is a small-sample artefact with 10 instruments and is not accompanied by a sign change (weighted median +0.062). Because FinnGen knee OA and the GO knee-OA GWAS are independent European samples, this table is a within-ancestry consistency check; cross-ancestry replication remains unexecuted (§S32, F10).
