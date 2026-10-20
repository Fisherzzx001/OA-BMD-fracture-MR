# STROBE-MR Checklist (Skrivankova et al., BMJ 2021;372:n223)

**Study:** Genetically higher bone mineral density raises knee osteoarthritis but protects against fracture — Mendelian randomisation of 13 skeletal phenotypes.
**Reporting guideline:** STROBE-MR (formerly STROBE-MR extension). Each item is marked **Reported** (with location in the main manuscript `docs/08`) or **N/A**.

| # | STROBE-MR item | Status | Location in manuscript |
|---|---|---|---|
| 1 | **Title** identifies the report as a Mendelian randomisation study | Reported | Title: "…Mendelian randomisation of 13 skeletal phenotypes" |
| 2 | **Abstract** follows STROBE-MR abstract structure (background, objectives, methods, results, conclusions) | Reported | Abstract (background/objectives/methods/results/conclusions present; unstructured per OAC) |
| 3 | **Introduction – Background:** rationale for the MR investigation | Reported | §1 Introduction (Foss & Byers 1972 → present debate; why MR) |
| 4 | **Introduction – Objectives:** specific objectives/hypotheses | Reported | §1 final paragraph (three open questions; 13-phenotype design) |
| 5 | **Methods – Data sources:** GWAS provenance, search, eligibility, date accessed | Reported | §2 Data sources (13 GWAS; Refs 1–3,10; accessions) |
| 6 | **Methods – GWAS description:** phenotype definition, participants, QC, effect allele, units, population | Reported | §2 Data sources; Table 1 (h², N, intercepts); §3 heritability |
| 7 | **Methods – Exposure instrument:** selection, number, pruning/clumping, population, proxies, palindromes | Reported | §2 MR (p<5×10⁻⁸; r²<0.001/10 Mb; 20–35 per BMD; 116 knee OA; HapMap3 EUR; palindromes excluded in LDSC) |
| 8 | **Methods – Outcome instrument / outcome GWAS:** same detail | Reported | §2 Data sources (knee OA GO 2025 N≈1.21M; self-rep IEU ukb-b-14486; fracture GEFOS) |
| 9 | **Methods – Assumption 1 (relevance):** instrument–exposure association strength (F, R²) | Reported | §2 MR (instrument counts; Table 3; S9 conditional F ≥10); §2 Effect-scale reconstruction |
| 10 | **Methods – Assumption 2 (independence):** pleiotropy assessment plan | Reported | §2 MR (MR-Egger intercept, MR-PRESSO, leave-one-out); §3 robustness |
| 11 | **Methods – Assumption 3 (exclusion restriction):** no confounders of instrument–outcome | Reported | §2 Sample overlap (Burgess bound); §2 MRlap; §4 limitations (residual overlap) |
| 12 | **Methods – Statistical analysis:** MR estimators, heterogeneity, sensitivity, software | Reported | §2 MR / LDSC / Colocalisation / MVMR (IVW, MR-Egger, weighted median, MR-PRESSO, LOO; Python 3.13; MRlap R 4.6.1) |
| 13 | **Methods – Sample overlap:** assessment and handling | Reported | §2 Sample overlap (Burgess approximation; bound 2.6% BMD→OA, 12% BMD→fracture; MRlap correction) |
| 14 | **Results – GWAS descriptives:** N, cases, h², intercept, population | Reported | Table 1; §3 Genetic-correlation architecture (intercepts 1.01–1.07) |
| 15 | **Results – Exposure instrument:** SNPs, F-statistic, R² | Reported | §2 MR (counts); Table 3 (n_snp); Table S13 (univariate F_min, F_agg, R²_exposure); S9 (MVMR conditional F ≥10) |
| 16 | **Results – Outcome instrument:** same | Reported | Table 3 (knee OA 69 SNPs); §2 Data sources |
| 17 | **Results – Assumption testing:** heterogeneity (Cochran Q/I²), pleiotropy (MR-Egger, MR-PRESSO), sensitivity | Reported | §3 (I² 68–86%; MR-Egger intercept p>0.24; MR-PRESSO 9/10 significant, 0–3 outliers); Table 3 |
| 18 | **Results – Main MR estimates:** with CIs, p-values, direction | Reported | Table 3 (IVW β, SE, p, random-effects p, Egger, weighted median, I², LOO); §3 |
| 19 | **Results – Supplementary/stratified analyses:** MVMR, coloc, reverse, MRlap | Reported | §3 (positive control, reverse MR, coloc, MVMR S9, MRlap Fig 10); S8, S9, S10 |
| 20 | **Results – Sample overlap results:** empirical overlap diagnostics | Reported | §3 (Fig 9 ρ_GE; §2 MRlap median +20.5% away from null) |
| 21 | **Discussion – Key results:** interpreted in context | Reported | §4 Principal findings; Resolving the OA–BMD paradox; Site and direction |
| 22 | **Discussion – Limitations, interpretation, generalisability, conclusions, funding, ethics** | Reported | §4 Strengths and limitations; Conclusion; §2 Data sources (ethics: secondary analysis, no new approval); Declarations (Funding: None; Competing interests: none); Data availability |

**Overall:** all 22 STROBE-MR items are addressed in the main manuscript and accompanying supplementary files. No item is Not Applicable.
