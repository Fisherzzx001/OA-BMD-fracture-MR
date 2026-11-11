<!-- Standalone supplementary table extracted verbatim from
     docs/08_投稿稿_OA_BMD_fracture.md by code/332_extract_inline_supplementary.py.
     Source slice sha256: 6da705048e19c195ae774bc7afd082b45fe8d9b0821c6160f253b08b7d98c9cc . No statistic was rewritten. -->

## Supplementary Table S29. Osteophyte-stratified public-data surrogate for falsifying test F5: site-stratified MR of BMD on osteophyte-dominant versus cartilage-dominant OA, with F4-style Bayesian colocalisation

F5 was pre-registered as a hip-osteophyte-stratified outcome. The hip-osteophyte GWAS exists only as a conference abstract without public summary statistics, so the registered test could not be run. As a public-data surrogate we substitute site-stratified OA outcomes in which spinal and hand OA are osteophyte-dominant (bone-forming) while knee OA is cartilage-loss-dominant, and ask whether the BMD→OA effect is stronger where osteophyte formation dominates — the formation-side prediction of the Systemic Osteogenic Set-point (SOS). The same 31 BMD_total instruments and the same pipeline as Table 2 are used; ORs are per SD higher genetically determined areal BMD (total femur).

**Part A. Site-stratified two-sample MR (BMD_total → OA).**

| OA site | Dominance | IVW OR per SD (95% CI) | p | Weighted-median OR per SD | MR-Egger OR (slope) | Egger intercept (SE) | Egger intercept p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Knee | Cartilage-dominant | 1.0812 (1.0552–1.1079) | 3.3×10⁻¹⁰ | 1.0470 | 1.0811 | +0.00043 (0.00089) | 0.629 |
| Spine | Osteophyte-dominant | 1.1291 (1.0900–1.1696) | 1.4×10⁻¹¹ | 1.1132 | 1.1295 | −0.00087 (0.00130) | 0.502 |
| Hand | Osteophyte-dominant | 1.1236 (1.0423–1.2113) | 2.4×10⁻³ | 1.0820 | 1.1232 | +0.00171 (0.00272) | 0.529 |

*Notes.* n = 31 BMD_total instruments per site. The spine and hand (osteophyte-dominant) ORs are numerically larger than the knee (cartilage-dominant) OR by ≈4–4.5% on the odds scale (spine−knee log-OR difference = +0.043; hand−knee = +0.039). All three MR-Egger intercepts are non-significant, so the gradient is not driven by directional pleiotropy; confidence intervals overlap broadly and the differences are not individually significant.

**Part B. F4-style Bayesian colocalisation extended to osteophyte-dominant sites (PP.H4 ≥ 0.8; ±500 kb window; priors p1=p2=1×10⁻⁴, p12=5×10⁻⁶, sd_prior quant=0.15 / cc=0.20).**

| OA site | Loci analysed (≥20 shared SNPs) | Loci PP.H4 ≥ 0.8 | Interpretation |
| --- | --- | --- | --- |
| Knee | 31 | 0 | Reproduces the published cross-trait F4 result (knee OA 0/31) |
| Spine | 31 | 0 | No shared causal variant at any instrument locus |
| Hand | 31 | 0 | No shared causal variant at any instrument locus |

*Notes.* Colocalisation is 0/31 at every site, identical to the published knee-OA cross-trait arm of F4 (Tables S8 and S26). The diffuse, polygenic architecture of the BMD→OA effect is therefore site-invariant: colocalisation does not separate osteophyte-dominant from cartilage-dominant sites. The only site-discriminating signal is the MR magnitude gradient in Part A, which is modest and not significant on its own.

**Synthesis for F5.** The magnitude gradient (osteophyte-dominant spine and hand > cartilage-dominant knee) is weakly supportive of the formation-side SOS prediction but not decisive: confidence intervals overlap and no site difference reaches significance. The colocalisation arm provides no discriminating evidence (0/31 everywhere). We therefore report F5 as partially addressed by a public-data surrogate rather than fully executed: the registered hip-osteophyte test remains unperformed for want of public data, and this surrogate substitutes site dominance (spine/hand osteophyte-dominant versus knee cartilage-dominant) as a proxy. The limitation is stated in the pre-registration-status note and the Table S11 ledger.
