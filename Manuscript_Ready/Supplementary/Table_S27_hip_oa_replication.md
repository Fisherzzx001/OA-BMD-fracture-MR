<!-- Standalone supplementary table extracted verbatim from
     docs/08_投稿稿_OA_BMD_fracture.md by code/332_extract_inline_supplementary.py.
     Source slice sha256: 87d979b5d7f95dc59f1a436184ba9bbf2f2a155c7e65e37bf418a7cedd3bb034 . No statistic was rewritten. -->

## Supplementary Table S27. Within-study cross-site hip-osteoarthritis concordance (falsifying test F1): two-sample MR of total-body BMD on hip OA, and the conjunctive adjudication against the knee-OA estimate

F1 was pre-registered as a prospective cross-site test of the BMD->OA finding at the hip (within the same GWAS, not an independent cohort), with a conjunctive criterion: the hip effect must be positive AND of knee-like magnitude. The hip outcome is the Hatzikotoulas 2025 mixed-ancestry hip-OA GWAS (GCST90566798; 95,395 European cases, 1,011,691 European controls plus East Asian, African, South Asian, Hispanic and admixed groups; maximum per-variant 97,328 cases and 1,055,375 controls), and the knee outcome is the same study's knee-OA GWAS (GCST90566800; maximum per-variant 172,256 cases and 1,144,248 controls) — the file from which the manuscript's own knee-OA estimate was derived, so that the two rows below are produced by an identical pipeline on an identical instrument set (31 total-body BMD instruments, Qian 2025).

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
