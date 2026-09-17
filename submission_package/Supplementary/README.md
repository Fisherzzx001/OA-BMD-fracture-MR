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
```bash
python 07_download_ldsc_ref.py      # 仅首次：1000G EUR LD 参考
python 09_rg_matrix.py              # rg 矩阵
python 13_mr_analysis.py            # MR（1-Mb 近似 clump）
python 17_strict_clump_mrpresso.py  # 严格 clump 命令 + MR-PRESSO
python 14_figure3_mr_forest.py      # Figure 3
python 16_recompute_h2_truen.py     # 真实 N 重算 h²
python 15_liability_h2.py            # OA liability h²
```

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

## 推送到 GitHub（本机无 gh CLI，用 git 直接推）
> 本地仓库已在项目根目录 `OA_B_PROJECT/` 初始化（已 `git add -A` 并提交，`.gitignore` 已排除 `code/ref/` 等第三方/大文件）。

```bash
cd "C:/Users/Administrator/WorkBuddy/2026-09-08-22-17-40/OA_B_PROJECT"
# 在 GitHub 新建空仓库后：
git remote add origin https://github.com/<USER>/<REPO>.git
git branch -M main
git push -u origin main
```

