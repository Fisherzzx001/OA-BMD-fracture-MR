# UK Biobank 申请材料（可直接粘贴到 AMS 申请系统）

> 项目标题（≤ 80 字符）  
> **Harmonised comorbidity networks in osteoarthritis: a binational comparison with CHARLS**

> 关键：申请时**务必以学生身份**（Student status）提交 + 导师/合作者在 PI 栏内，
> 锁定 **£500 学生费**而非 £3,000–9,000 全额。这是 B 档能否在 ¥1 万预算内闭环的决定性一步。
> 本项目已对应 OSF 预注册（链接在最终提交时填入），UKB 申请号（95180）作为参考。
> 参考申请 95180 的方法学（Ding 2025 Arch Public Health）已在 OSF 公开文档中归档。

---

## 1. Lay summary（300–500 词）

Osteoarthritis (OA) is the most common joint disease in older adults and rarely occurs alone—people with OA typically have two or more other chronic conditions such as high blood pressure, diabetes, or back pain. These "comorbidities" drive most of the disability and healthcare costs in OA, but we still do not know whether their patterns differ between countries. A clear answer matters because it would tell doctors whether OA management needs to be tailored by country (e.g., focusing on metabolic diseases in Chinese patients and on bone/muscle diseases in British patients) or whether the same global guideline is enough.

We will answer this by analysing two large, well-characterised population cohorts: the UK Biobank (UK, ~500,000 adults) and the China Health and Retirement Longitudinal Study (CHARLS, ~17,000 adults aged 45+). In each cohort we will identify people with OA, build a "comorbidity network" showing which diseases co-occur, and compare the structure between the two countries. We will use a recently-developed harmonisation method (the "comparable network" framework) to make sure that differences we see reflect real biology, not just differences in how diseases are recorded in each database.

This study is low risk (no direct participant contact; only secondary analysis of de-identified data) and has clear public-health benefit. The methods and code will be made open-source for other researchers to use.

---

## 2. Background and rationale

OA is the leading cause of years lived with disability in adults over 50 globally, and its systemic impact is increasingly recognised. Comorbidities—not the joint disease itself—drive most of the excess mortality and cost. Two recent publications illustrate the state of the field: (i) Ding et al. (Arch Public Health 2025; PMID 40604842) quantified age-stratified single-comorbidity hazard ratios in OA using UK Biobank and CHARLS but did not analyse network structure; (ii) the present proposal extends their approach by applying and extending the "comorbidity network" framework recently validated for fracture multimorbidity (Meng et al., Osteoporos Int 2026; PMID 42342916). No published study has yet compared the **structural** comorbidity network of OA patients between a Chinese and a European population with statistical inference, leaving a clear gap that prevents precision public-health recommendations.

**Why UK Biobank is essential for this question**: UKB uniquely provides (a) nationwide coverage with both self-reported and ICD-10-coded diagnoses; (b) deep linkage to primary care and hospital records; (c) ~100,000 OA cases for a high-powered 1:1 matched comparison. No other open-access resource offers this combination.

---

## 3. Aims and objectives

**Primary aim**: To determine whether the comorbidity network structure of adults with OA differs between the UK and China, using harmonised methods.

**Specific objectives**:
1. Build harmonised comorbidity networks (16 chronic conditions) in OA cases and matched non-OA controls in UKB and CHARLS.
2. Quantify network complexity with five indicators: average degree, density, modularity Q, centralisation, and clustering coefficient (the first two are pre-registered primary; the latter three are secondary).
3. Test within-cohort excess (OA vs control, H1) using paired-label permutation (10,000 permutations).
4. Test cross-country module-level differences (H2/H4) using two-sample bootstrap (5,000 resamples) and an aetiological-module pre-specified aggregation score (mean within-module O/E).
5. Compare hub diseases (degree centrality) between cohorts (H3).

**Pre-registered hypotheses** (locked on OSF before data extraction; see `docs/01_*`):
- H1: OA patients in both cohorts have a denser comorbidity network than matched controls.
- H2: Chinese OA networks show stronger metabolic-module aggregation; UK OA networks show stronger musculoskeletal-module aggregation.
- H3: Top-5 hub diseases in Chinese OA network are metabolic; in UK they are bone/musculoskeletal.
- H4 (insurance): even if H1 is not confirmed, module-level structure differs between countries.

---

## 4. Study design and statistical analysis plan

Type: dual-cohort, population-based, matched case–control study with cross-cohort network comparison.

**Sample**: UK Biobank participants with a first OA diagnosis (ICD-10 M15–M19 from hospital records **or** self-reported osteoarthritis in Field 20002) and 1:1 age-, sex-, and assessment-centre-matched non-OA controls.

**Variables**:
- Exposure: OA status (binary), with optional restriction to knee (M17) and hip (M16) OA in secondary analyses.
- Outcomes: 16 chronic conditions harmonised with CHARLS (hypertension, diabetes, dyslipidaemia, heart disease, stroke, asthma, COPD, liver disease, kidney disease, digestive disease, cancer, depression, dementia, chronic back pain, neck/shoulder pain, multi-site joint pain), each coded as binary.
- Covariates: age, sex, BMI, smoking, alcohol, education, Townsend deprivation index.

**Analytic pipeline** (all steps pre-registered; see `code/lib_network.py` for implementation):
1. Disease-specific ICD-10 and self-report code mapping (`code/fields.py`, `code/02_harmonise.py`).
2. 1:1 matching on age (±1 year), sex, and assessment centre (`code/01_build_cohort.py`).
3. Network construction per cohort: nodes = diseases (prevalence > 1%), edges = significant excess co-occurrence (observed/expected ratio > 1, BH-FDR < 0.05).
4. Equal-density sparsification for cross-country comparison (K = min of OA-network edge counts in each cohort).
5. Five network indicators (average degree, density, modularity, centralisation, clustering).
6. H1 test: paired-label permutation (10,000 per cohort) of OA vs control.
7. H2/H4 test: two-sample bootstrap (5,000 resamples) of cross-country differences in module aggregation score (mean within-module O/E).
8. H3 test: degree-centrality comparison with BH-FDR correction.
9. Eight pre-specified sensitivity analyses (definition breadth, missing-data imputation, weighting, sparsification threshold, incidence-only OA, temporal network, two-year washout, sex interaction).

**Software and reproducibility**: open-source Python pipeline (MIT) + R translation available; GitHub release with Zenodo DOI; UKB analysis conducted on the UK Biobank Research Analysis Platform (RAP).

---

## 5. Data access and security

All work will be performed within the UKB Research Analysis Platform (RAP), with no local download of individual-level data. Only aggregate results and code will leave the RAP.

---

## 6. Ethical considerations

UK Biobank has received ethics approval from the North West Multi-centre Research Ethics Committee (REC reference 11/NW/0382). This study is a secondary analysis of de-identified data; no re-identification is attempted. The PI's host institution will provide a local institutional approval letter as part of the application.

The CHARLS arm (covered under Peking University IRB 00001052-11015) will use the publicly available de-identified dataset and does not require additional UKB-side ethics review.

---

## 7. Research team

| Role | Name | Affiliation | Contribution |
|---|---|---|---|
| PI / applicant | (your name, student status) | Hospital & University | Project lead, data extraction, analysis, writing |
| Co-PI / guarantor | (supervisor, faculty) | (institution) | UKB institutional sponsorship, clinical interpretation, manuscript guarantor |
| Statistical co-PI | (epidemiologist) | School of Public Health | Sampling design, network inference |
| UK collaborator | (to be confirmed, see §9) | UK-based institution | UKB application sponsorship, manuscript review |

---

## 8. Timeline and outputs

| Month | Milestone |
|---|---|
| 0 | UKB application submitted; OSF pre-registration |
| 1–2 | UKB approval; data extraction via RAP |
| 3–4 | Harmonisation with CHARLS; matched cohorts |
| 4–8 | Primary analysis (H1–H4) |
| 8–10 | Sensitivity analyses (8) |
| 10–12 | Manuscript drafting; preprint on medRxiv |
| 12–15 | Submission to *Osteoarthritis and Cartilage* (primary target) |
| 15–18 | Peer review and revision |

Open-access publications and code (no IP or commercial restrictions).

---

## 9. Resources to confirm before submission

- [ ] UKB institutional sponsorship letter from a UK-based PI (or a co-investigator at a UK institution; this is **non-negotiable** and is the most common cause of delay). Strategies:
  1. Identify an existing OA / comorbidity / rheumatology researcher in the UK with active UKB access and offer co-authorship + methodological contribution.
  2. Use your home institution's UK collaborating department.
- [ ] Supervisor / PI of record (must hold a permanent or temporary academic position).
- [ ] OSF project created and the pre-registration PDF uploaded.
- [ ] GitHub repository (public, MIT license) with the analysis code at the application time.
- [ ] Conflict of interest disclosure: none.

---

## 10. Estimated cost

| Item | Cost |
|---|---|
| UK Biobank access (student fee) | **£500** |
| UK Biobank RAP compute (default project allowance) | £0–1,000 |
| Open access fee (waivers available for many journals; if needed) | up to £2,000 |
| **Total** | **< £3,500** |

> Submitting under a non-student / non-UKB-collaborator pathway would cost £3,000–9,000 in access fees plus standard compute, which is why §9.1 is critical.
