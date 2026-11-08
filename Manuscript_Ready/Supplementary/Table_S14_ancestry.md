# Supplementary Table S14. Population-ancestry composition of the 13 GWAS

The BMD exposure instruments are European-only, whereas the knee-OA outcome is 87.31% European (multi-ancestry meta-analysis); LDSC uses a European-only reference, so cross-population LD differences between the BMD instruments and the multi-ancestry knee-OA outcome can bias r_g and the MR estimates. This ancestry mismatch is a primary limitation (§4). The OA_knee percentages are taken from the GWAS Catalog entry and the originating Nature report of the multi-ancestry meta-analysis.

| GWAS (phenotype) | Source | N (effective) | Ancestry composition | LDSC reference |
| --- | --- | --- | --- | --- |
| BMD (10 regional sites) | Qian et al. 2025 [2] | 31,873 each | 100% European (UK Biobank) | 1000G Phase 3 EUR |
| OA_knee (radiographic) | GO Consortium 2025 [1] | ≈1,208,397 | Multi-ancestry: EUR 87.31%, EAS 7.09%, AFR 3.08%, SAS 1.09%, HIS 0.91%, ADM 0.53% | 1000G Phase 3 EUR |
| OA_selfrep (self-reported) | IEU OpenGWAS `ukb-b-14486` [3] | 462,933 | European (UK Biobank) | 1000G Phase 3 EUR |
| Fracture (any) | GEFOS 2018 [10] | 264,973 | European | 1000G Phase 3 EUR |

**MRlap overlap correction** (full per-pair values in the manuscript's Table S14 note). No overlap proportion was supplied; MRlap v0.0.3.3 ran on its default and does not return a scalar overlap fraction — overlap is accommodated through the internally estimated cross-trait LDSC intercept λ (`int_crosstrait` in `code/output/mrlap_overlap_results.csv`), which ranged from −0.0361 to +0.0153 (|λ| ≤ 0.037) across the ten pairs. The UKB overlap counts tabulated above are reported independently and were not fed to MRlap; all ten corrected estimates retained their direction (median +20.5% away from the null), so no conclusion depends on the default.

**Burgess overlap-bias bounds (bias = overlap fraction × r_g; §2 of the manuscript).** The tabulated counts are BMD 31,873 participants (100% UKB), OA_selfrep 462,933 (100% UKB), fracture 22,330 UKB participants (8.4% of 264,973; 70% of BMD), and OA_knee UKB participation not separately reported (overlap upper bound 31,873 / 1,208,397 = 2.6%). Hence BMD→OA_knee bias ≤ ~0.3–0.7% (2.6% × r_g 0.12–0.25), BMD→OA_selfrep ~0.8–1.7% (6.9% × 0.12–0.25) and BMD→fracture ~3.4–4.5% on the conservative 8.4%-of-fracture basis (up to ~12% on the 70%-of-BMD basis); all act in the direction of r_g. MRlap accounts for overlap directly and retained all ten directions (median +20.5% away from the null), so the positive BMD→knee-OA effect survives the worst-case overlap correction.
