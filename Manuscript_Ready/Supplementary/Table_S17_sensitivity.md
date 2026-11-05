## Supplementary Table S17. Sensitivity analyses requested at review

All values below are computed from first-party artefacts in the submission repository (`code/output/genetics/`); script `code/187_sensitivity_table.py` regenerates this file from those artefacts. No reported beta, SE, p-value or odds ratio is altered anywhere in this table.

### S17-a. Uncertainty of the effect-scale reconstruction constant K

| Scenario | K(1−K) | Implied K | Scale factor 1/√[K(1−K)] | Factor vs reported (%) | OR at the upper end of the reported range | OR at the lower end |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Reported (median) | 0.0871 | 0.0964 | 3.388 | +0.00 | 1.099 | 1.060 |
| Lower quartile of the 39-SNP distribution | 0.0835 | 0.0920 | 3.461 | +2.13 | 1.101 | 1.061 |
| Upper quartile of the 39-SNP distribution | 0.0881 | 0.0976 | 3.369 | -0.57 | 1.098 | 1.060 |
| K deliberately over-specified +10% | 0.0948 | 0.1060 | 3.248 | -4.14 | 1.095 | 1.057 |
| K deliberately under-specified −10% | 0.0792 | 0.0868 | 3.553 | +4.85 | 1.104 | 1.063 |

Reported values sit in the row 'Reported (median)' (factor 3.389, OR 1.060–1.099). Because the reconstruction multiplies beta and SE by the *same* constant, z and therefore every p-value, significance call, Cochran Q, I², Egger-intercept test, PRESSO call and power value are exactly invariant to K; only the absolute effect size moves. Across the inter-quartile range of the 39 SNPs carrying both a deposited log-odds effect and the reconstructed Z, the scale factor varies by +2.13% to -0.57%, which shifts the headline odds ratios by less than 0.3% in absolute terms. Even a deliberately pessimistic ±10% mis-specification of K moves the odds ratios by no more than about ±0.5%, well inside their stated confidence intervals, and changes no inference.

### S17-b. Number of statistically significant BMD→knee-OA exposures under three random-effects conventions

| BMD site | IVW fixed-effect p | Multiplicative random (SE×√max(1,Q/df)) p | DerSimonian–Laird additive-τ² p | τ² | I² (%) |
| --- | ---: | ---: | ---: | ---: | ---: |
| total | 3.86e-10 | 0.0062 | 0.00481 | 0.0212 | 81 |
| lumbar | 7.7e-09 | 0.0166 | 0.0495 | 0.0273 | 83 |
| arms | 1.71e-06 | 0.0699 | 0.0333 | 0.0257 | 86 |
| leg | 2.29e-05 | 0.0976 | 0.0735 | 0.0283 | 85 |
| femur neck | 0.000636 | 0.0547 | 0.0538 | 0.0157 | 68 |

Significant at p<0.05: **5/5** under the fixed-effect model, **2/5** under the multiplicative convention used as primary, and **3/5** under an additive DerSimonian–Laird τ² implementation. All five point estimates are positive and all five agree in direction under every convention; only the significance count depends on the convention, which is why the causal reading rests on estimator concordance rather than on fixed-effect p-values.

### S17-c. Change in the IVW estimate on removal of MR-PRESSO outliers (all ten pairs)

| Pair | Instruments | Outliers removed | IVW beta (full) | IVW beta (outlier-removed) | Absolute change | Relative change (%) | Direction changed? |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| total → fracture | 31 | 2 | -0.3371 | -0.3520 | -0.0148 | -4.4 | no |
| femur_neck → fracture | 25 | 3 | -0.3507 | -0.3447 | +0.0060 | +1.7 | no |
| total → knee | 31 | 3 | +0.0783 | +0.0884 | +0.0101 | +12.9 | no |
| femur_neck → knee | 25 | 3 | +0.0579 | +0.0708 | +0.0129 | +22.3 | no |
| lumbar → knee | 21 | 3 | +0.0946 | +0.1133 | +0.0187 | +19.8 | no |
| arms → knee | 16 | 2 | +0.0751 | +0.0866 | +0.0115 | +15.3 | no |
| leg → knee | 24 | 3 | +0.0600 | +0.0494 | -0.0106 | -17.7 | no |
| OA_knee → fracture | 69 | 2 | +0.0272 | +0.0591 | +0.0319 | +117.2 | no |
| OA_knee → BMD_total | 69 | 1 | +0.0313 | +0.0506 | +0.0193 | +61.5 | no |
| OA_knee → BMD_femur_neck | 69 | 0 | +0.0419 | +0.0419 | +0.0000 | +0.0 | no |

No pair changed direction on outlier removal; its effect moved away from the null for the majority of pairs, so the outlier correction is not load-bearing for any headline conclusion.

### S17-d. Colocalisation counts across the knee-OA endpoint as a function of the PP.H4 threshold

| Threshold | Loci reaching PP.H4 ≥ threshold (of 31 BMD instrument loci) | Comment |
| --- | ---: | --- |
| ≥0.90 | 0/31 |  |
| ≥0.80 | 0/31 | pre-specified cut-off; the applied analysis |
| ≥0.75 | 1/31 |  |
| ≥0.70 | 1/31 | relaxed: admits the single marginal locus rs11228269 (PP.H4 0.7994) |
| ≥0.60 | 2/31 | admits a second locus (rs9594759, PP.H4 0.6271) |
| ≥0.50 | 3/31 |  |

The five highest PP.H4 values against knee OA are rs11228269 (0.7994), rs9594759 (0.6271), rs12600549 (0.5530), rs9293511 (0.4518), rs6040061 (0.2941). Against the fracture endpoint 9/31 loci reach PP.H4 ≥0.8, so the endpoint-split contrast reported in the manuscript is preserved at every threshold from 0.5 to 0.9: relaxing the cut-off to 0.7 adds exactly one locus and leaves the interpretation unchanged.

