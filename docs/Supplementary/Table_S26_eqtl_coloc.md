<!-- Standalone supplementary table extracted verbatim from
     docs/08_投稿稿_OA_BMD_fracture.md by code/332_extract_inline_supplementary.py.
     Source slice sha256: e57813ffe2de7bb7631824e0e87895606962a94503226d4b3452a76b3e8bccf0 . No statistic was rewritten. -->

## Supplementary Table S26. eQTL colocalisation — tissue-literal (cartilage and skeletal-muscle bone-proxy) arms of F4 (this revision)

At each of the 31 total-body BMD instrument loci (±500 kb), tissue eQTL summary statistics were colocalised with OA knee, total-body BMD and any-fracture by coloc.abf (priors p1=p2=1e-4, p12=5e-6; eQTL variance prior SD=0.15). Two tissues were tested as the genetic counterparts of the two components of the osteochondral unit: (i) cartilage, from FunGen low-grade osteoarthritis cartilage (Southam et al.; n≈200 donors); and (ii) skeletal muscle, used as a bone-proximal proxy from the Wilson et al. (AJHG 2025) meta-analysis (GTEx+FUSION, n=1,002). PP.H4 = posterior probability of a shared causal variant. Neither tissue reached PP.H4≥0.8 for BMD or fracture, and only one muscle (bone-proxy) locus reached it for OA knee (PP.H4=0.887, a single borderline signal). The BMD instruments therefore do not colocalise with cartilage or skeletal-muscle regulatory variation at genome-wide-significant density; the strongest cartilage signal was with total-body BMD (top locus PP.H4=0.475). This is consistent with a formation-side (bone/subchondral) rather than cartilage-matrix or muscle-mediated interpretation of the BMD→OA effect, and parallels the cross-trait arm in which the same instruments colocalise with fracture (bone) but not OA (joint). Per-locus PP.H4 values are in f4_literal_cartilage_*.csv and f4_literal_muscle_*.csv (code/output/genetics/).

| GWAS trait | Cartilage: n loci | Cartilage median PP.H4 | Cartilage n PP.H4≥0.8 | Muscle: n loci | Muscle median PP.H4 | Muscle n PP.H4≥0.8 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| OA knee | 31 | 0.009 | 0 | 31 | 0.002 | 1 |
| BMD total | 31 | 0.061 | 0 | 31 | 0.000 | 0 |
| fracture | 31 | 0.002 | 0 | 31 | 0.004 | 0 |

---
