# OA ↔ BMD ↔ Fracture 跨性状遗传分析

> **Bone Mineral Density and Osteoarthritis Share Genetic Architecture but Not Causal Pathways to Fracture: A Cross-Trait Analysis of 13 Skeletal Phenotypes**
>
> 代码与中间结果仓库（投稿前整改 #6）。对应稿件：`docs/08_投稿稿_OA_BMD_fracture.md`

## 研究概述
对 13 个骨骼表型（2 个 OA 定义、10 个部位 BMD、任何骨折）进行：
1. **LDSC** 跨 78 对表型的遗传相关（rg）矩阵（原生 Python 3 实现，非 py2 移植）；
2. **两样本 MR**（IVW 固定/随机、MR-Egger、加权中位数、Cochran Q/I²、radial 离群、留一法）含阳性对照（BMD→骨折）；
3. **MR-PRESSO** 全局（distortion）检验与严格 PLINK clumping（r²<0.001, 10 Mb）。

## 目录结构
```
code/
  01_build_cohort.py ... 07_download_ldsc_ref.py    # 队列 / 网络分析（B 档观察性部分）
  08_genetic_correlation.py   # 原生 LDSC：h² / gencov / rg
  09_rg_matrix.py             # 13 表型 rg 矩阵 + 等位基因 harmonisation（含 SIGN_FLIP）
  10_fracture_row.py          # 骨折 GWAS 载入
  11_rg_table.py 12_final_figures.py
  13_mr_analysis.py           # 两样本 MR（1-Mb 近似 clump）
  14_figure3_mr_forest.py     # Figure 3 MR 森林图
  15_liability_h2.py          # OA h² 转 liability scale（Lee 2011 Eq 23）
  16_recompute_h2_truen.py    # 按真实 N 重算 h²（BMD≈32k）
  16b_fracture_h2.py          # 骨折 h² 对 N 的敏感性
  17_strict_clump_mrpresso.py # 严格 PLINK clump + MR-PRESSO（本整改新增）
  output/genetics/            # 中间结果：rg_*, mr_results.csv, mr_presso_results.csv,
                              #           Figure1-3.png, *_h2*.csv, tophits_*, QC_*
docs/
  08_投稿稿_OA_BMD_fracture.md     # 稿件
  09_数据下载实操指南.md           # UKB/CHARLS 下载
  10_参考文献与原始GWAS来源.md     # 原始 GWAS 溯源 + 参考文献（整改 #2）
  11_补充材料_OA_selfrep_Z校正.md # Z 方向校正证据（整改 #8）
  12_补充材料_MR_PRESSO.md         # MR-PRESSO 逐对结果（整改 #3）
```

## 环境与依赖
- Python 3.13
- `numpy`, `pandas`, `scipy`, `matplotlib`
- 安装：`python -m venv env && env/Scripts/pip install numpy pandas scipy matplotlib`

## 运行顺序（遗传学部分）
主分析管线（Python 3.13，依赖见根目录 `requirements.txt`）：
```bash
python 07_download_ldsc_ref.py      # 仅首次：1000G EUR LD 参考（写入 code/ref/，已 git-ignore）
python 09_rg_matrix.py              # 13 表型 rg 矩阵（S1/Figure 1–2）
python 13_mr_analysis.py            # MR（1-Mb 近似 clump；旧 Table 3 来源）
python 17_strict_clump_mrpresso.py  # 严格 PLINK 等效 clump（r²<0.001 / 10 Mb）+ MR-PRESSO（主分析）
python 18_plink_clump_pure.py       # 纯 PLINK clump 复现校验
python 19_N_sensitivity_decision.py # N-尺度敏感性决策（严格 vs 1-Mb）
python 21_mr_corrected_N.py         # 按重建 N 尺度校正 MR 效应（Table 3 主分析）
python 22_correct_figure3.py        # Figure 3 森林图效应尺度校正
python 23_audit_rederive.py         # 审计：从源重算校验
python 24_rebuild_table3_corrected.py # Table 3 严格集重算
python 25_audit2_rederive.py        # 第二轮审计重算
python 26_export_harmonised_snp.py  # 导出 10 对每-SNP 协调文件（harmonised/*.csv）
python 27_overlap_statistics.py     # 样本重叠 ρ_GE（Figure 9 / S）
python 28_mr_power.py               # MR 统计效能（Table 4 / Figure 8）
# 文章图（Figure 1–9 + S1–S4）：
python build_figures.py             # Fig 1–5, 8, 9, S3（rg 热图 / 森林 / 因果网络 / h² / 效能 / 重叠）
python build_figures_v9.py          # Fig S1–S3（MR 散点 / 留一法 / GEFOS 验证）
python build_figures_mr.py          # Fig 6, 7, S1, S2（MR 散点/漏斗/留一法）
python build_figures_v14.py         # 顶刊美化重绘（Okabe-Ito 色盲安全 + PDF 矢量），含 GraphicalAbstract / S4
# 稿件：
python build_submission.py          # 由 docs/08 生成 Manuscript.docx（OUTDIR 可指定输出目录）
python build_maintext.py            # 生成盲法 Submission_MainText.docx
```
审计脚本（`90_audit_mr_recompute.py` / `91_audit_crosscheck.py` / `92_audit_figures_v9.py`）用于投稿前独立复盘，输出 `AUDIT_*.csv`，详见 `docs/15`–`docs/21`。

### 样本重叠敏感性分析（MRlap，需 R）
```r
# 在本地 R (>= 4.1) 环境运行；本构建沙箱无 R
# install.packages(c("data.table","dplyr","readr"))
# remotes::install_github("MRCIEU/MRlap")
Rscript code/mrlap_overlap_correction.R
```
脚本读取 `output/genetics/harmonised/<EXP>__<OUT>.csv`，对 10 对做 MRlap 重叠校正，
输出 `output/genetics/mrlap_overlap_results.csv`（MRlap β / se / p 与估计重叠相关 ρ，对比朴素 IVW）。
结果用于核验稿件 §2.8 的 Burgess 2016 上界（BMD→OA ≤2.6%、BMD→fracture ≤12%）是否保守。


## 数据可用性声明（Data availability）
本研究所用汇总统计均来自公开发表的 GWAS（来源与登录号见 `docs/10_参考文献与原始GWAS来源.md`），
受各原联盟数据使用条款约束，**本仓库不重新分发原始汇总统计文件**。本仓库提供：
1. 全部重格式化 / harmonisation / 分析代码；
2. 中间结果（`output/genetics/` 下 rg 矩阵、MR 估计、MR-PRESSO、Figures、各表型 h² 与 top-hits）；
3. 严格 clumping 输入文件（`*_clump_input.txt`）。

原始汇总统计请按以下来源申请 / 下载：
- OA_knee：GO Consortium 2025（*Nature* 641:1217；DOI 10.1038/s41586-025-08771-z；GWAS Catalog）
- BMD ×10：Qian et al. 2025（*Genomics Proteomics Bioinformatics*；DOI 10.1093/gpbjnl/qzaf097；WBBC wbbc.westlake.edu.cn / GWAS Catalog）
- OA_selfrep、fracture：Neale Lab UK Biobank Round 2（OpenGWAS `ukb-d-*`；nealelab.is/uk-biobank）
- LD 参考：1000 Genomes Phase 3 EUR（`eur_w_ld_chr`）

## 引用
分析代码：本仓库。方法学见稿件 References（Bulik-Sullivan 2015; Bowden 2016; Verbanck 2018; Lee 2011）。

## License
代码：MIT。汇总统计数据归属原联盟，依其条款使用。

---

## 推送到 GitHub + 归档到 Zenodo（本机无 gh CLI / 无 R，按步骤手动执行）
> 本地仓库已在项目根目录 `OA_B_PROJECT/` 初始化（已 `git add -A` 并提交，`.gitignore` 已排除 `code/ref/` 等第三方/大文件）。

### 1) GitHub
```bash
cd "C:/Users/Administrator/WorkBuddy/2026-09-08-22-17-40/OA_B_PROJECT"
# 在 GitHub 新建空仓库（建议名 OA-BMD-fracture-MR，公开）后：
git remote add origin https://github.com/<USER>/OA-BMD-fracture-MR.git
git branch -M main
git push -u origin main
```

### 2) Zenodo 归档（生成 DOI，供稿件 Data availability 引用）
- 在 Zenodo 用 GitHub 账号授权并开启 "GitHub" 同步，连接上述仓库；
- 打一个带版本标签的 release：`git tag v15 && git push origin v15`；
- Zenodo 自动抓取并归档，生成 `doi:10.5281/zenodo.XXXXXXX`；
- 将真实 `XXXXXXX` 与 GitHub 仓库 URL 回填到根目录 `.zenodo.json` 的 `related_identifiers`、
  `docs/08` 的 Data availability 段，以及各作者 ORCID（当前为 `0000-0000-0000-0000` 占位）。

> 注意：本分析**不重新分发原始 GWAS 汇总统计**（受原联盟条款约束），仅发布代码与每-SNP 协调文件，符合公开仓储政策。


