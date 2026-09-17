# 投稿材料包 · Submission Package — OA / BMD / Fracture cross-trait MR

**论文题目：** Bone Mineral Density and Osteoarthritis Share Genetic Architecture but Not Causal Pathways to Fracture: A Cross-Trait Analysis of 13 Skeletal Phenotypes
**目标期刊（建议）：** *Bone* 或 *Osteoporosis International*（首选）；*Osteoarthritis and Cartilage* 为冲刺目标。
**状态：** 稿件科学/统计/作者信息全部就绪，可投稿。
**打包日期：** 2026-09-17

---

## 1. 包内文件清单

| 文件 / 文件夹 | 内容 | 用途 |
|---|---|---|
| `Manuscript.docx` | 最终投稿稿（标题页 + 摘要 + 正文 + 图注 + 表 + 参考文献，单文件完整版） | 主稿件上传（单文件版） |
| `Submission_TitlePage.docx` | 单独标题页（作者/单位/基金/COI/Running title） | 多数期刊要求标题页与正文分传 |
| `Submission_MainText.docx` | 盲法正文（仅标题 + 摘要 + 正文，不含作者信息） | 双盲审稿用正文 |
| `Cover_Letter.docx` | 投稿信 | 投稿系统 Cover Letter 上传 |
| `Figures/` | Figure1_rg_heatmap_13.png、Figure2_forest_key_results.png、Figure3_MR_forest.png | 三张主图（含图注见稿件） |
| `Tables/` | `rg_matrix_13.csv`、`mr_corrected_N.csv`、`mr_presso_results.csv`、`oa_h2_liability.csv`、`strict_clump_summary.csv`、`rg_h2_trueN.csv` | 支持性数据（rg 矩阵、MR、MR-PRESSO、h²、严格 clump 汇总、h²_trueN） |
| `Supplementary/` | `10_参考文献与原始GWAS来源.md`、`11_补充材料_OA_selfrep_Z校正.md`、`12_补充材料_MR_PRESSO.md`、`README.md`（代码仓库说明） | 补充材料 / 方法学溯源 |

> 说明：稿件中的 Table 1–3（h²、rg 模块、MR 估计）已内嵌于 `Manuscript.docx`；`Tables/` 提供其底层数据以便审稿人与编辑复核。

## 2. 投稿前最后一项待办（用户操作）

- [ ] **在 GWAS Catalog 中核实 BMD 原始汇总数据的 accession：`GCST90568457`（Qian et al. 2025）**，确认 N=31,873 与来源论文一致。这是稿件中唯一依赖外部核实的编号；正文已注明"should be confirmed directly in the GWAS Catalog before submission"。

## 3. 引用占位符处理说明（已自动修复）

稿件初稿中的行内 `[REF]` 占位符已全部解析：
- Foss & Byers 1972 → 新增为 **参考文献 11**（骨密度与髋 OA / 骨折的经典文献，*Ann Rheum Dis* 1972;31(3):259–264）。
- LDSC 截距稳健性 → 参考文献 6（Bulik-Sullivan 2015）。
- 两样本 MR 方法学 → 参考文献 7（Bowden weighted median）。
- Hackinger 2017 → 参考文献 4；Zhou 2025 → 参考文献 5。

## 4. 期刊格式提示

- **Bone / Osteoporosis International** 通常要求：标题页（作者+单位+基金+COI，本包已并入 `Manuscript.docx` 顶部）+ 摘要（结构式，已含 Objective/Methods/Results/Conclusion）+ Keywords（已含）+ 正文 + 参考文献 + 图注 + 表。可按期刊系统拆分标题页与盲法正文。
- 数学公式以行内纯文本给出（β = Z·SE 等），初次投稿可接受；若期刊要求 MathML/图片，可后续替换。
- 图表为 300 dpi PNG，已满足多数期刊要求。

## 5. 版本与来源

- 稿件源文件：`docs/08_投稿稿_OA_BMD_fracture.md`（git 已提交，commit `794b2f9`）。
- 统计重算与审计：两次独立复盘审计（`docs/15`、`docs/16`）全部条目已修复并提交。
- 代码仓库：见 `Supplementary/README.md`（含严格 clump 与 MR-PRESSO 复现命令）。
