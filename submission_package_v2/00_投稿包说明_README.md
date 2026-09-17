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
| `Figures/` | **美化/新增后的最终图（5 张，PNG + 矢量 PDF）**：`Figure1_rg_heatmap_13`、`Figure2_forest_key_results`、`Figure3_MR_forest`、`Figure4_h2_bar`、`Figure5_causal_network` | 投稿用最终图（图注见稿件 §Figure legends） |
| `Figures_Original/` | **原始图（保留，未改动）**：原 `Figure1_rg_heatmap_13.png`、`Figure2_forest_key_results.png`、`Figure3_MR_forest.png` | 版本对照 / 保留原始输出 |
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
- 图表为 300 dpi PNG（另附矢量 PDF 供排版与后续编辑），已满足多数期刊要求。

## 5. 版本与来源

- 稿件源文件：`docs/08_投稿稿_OA_BMD_fracture.md`（git 已提交）。
- 统计重算与审计：两次独立复盘审计（`docs/15`、`docs/16`）全部条目已修复并提交。
- 代码仓库：见 `Supplementary/README.md`（含严格 clump 与 MR-PRESSO 复现命令）。

---

## 6. 图片更新说明（本次新增 / 美化）

### 6.1 做了什么

依据论文数据与统计结果，对原有 3 张主图做了**重绘美化**，并**新增 2 张图**，以强化"三模块遗传结构 + 因果方向"这一核心叙事：

| 图 | 状态 | 内容与改进 |
|---|---|---|
| **Figure 1** | 重绘美化 | 13×13 遗传相关热力图。改用发散色阶（RdBu_r，以 0 为中心）、单元格内标注 r_g 数值与 BH-FDR 显著性星号（\* q<0.05 / \*\* q<0.01 / \*\*\* q<0.001，在 78 个表型对上计算）、黑线分隔 OA / fracture / BMD 三个模块、附色标。 |
| **Figure 2** | 重绘美化 | 关键遗传相关森林图。改为 4 面板（a 骨折×10 个 BMD、b 膝 OA×BMD、c 自报 OA×BMD、d OA–OA 与 OA–骨折），误差棒为 block jackknife 95% CI，显著用蓝色、不显著（两个 OA–骨折对）用灰色，直观体现"预期中的零结果"。 |
| **Figure 3** | 重绘美化 | MR 森林图。改为 3 面板（a BMD→骨折阳性对照、b BMD→膝 OA、c 反向/零结果），β 为每 1 SD 暴露的效应量并标注 95% CI；蓝=保护、红=风险升高、灰=无效应；各面板独立缩放，b 面板放大后 CI 清晰可见。 |
| **Figure 4** | **新增** | SNP 遗传力（h²，观测尺度）横向条形图，13 个表型按分组配色（OA 紫 / 骨折 青 / BMD 蓝），含 95% CI；菱形标记二分类表型（膝 OA、自报 OA、骨折）的**易感性尺度 h²**。把原 Table 1 的数值可视化，强化 §3.1。 |
| **Figure 5** | **新增** | 概念性因果网络图（OA–BMD–骨折三角）。实线箭头为 MR 因果效应（蓝：BMD→骨折保护；红：BMD→OA 因果升高风险），虚线为"非因果"（OA 与 BMD 遗传相关，但 OA 不因果性升高 BMD 或骨折风险）。适合作为 Graphical Abstract / TOC 图。 |

### 6.2 使用的绘图规范（Skill）

- **nature-figure**（Nature 系列规范）：Arial 无衬线字体、≥7 pt、色盲安全配色（蓝 #4878CF / 红 #D65F5F / 青 #82C6E2 / 紫 #956CB4 / 橙 #E69F00）、去除顶部与右侧边框、小写粗体面板标签（a/b/c/d）、300 dpi、同时输出可编辑矢量 PDF（`pdf.fonttype=42`）。Figure 1–3、Figure 5 按此规范绘制。
- **science-figures**（scienceplots 期刊风格）：用于 Figure 4 柱状图的原型生成与期刊样式套用（`gen_figure.py --type bar --style nature`），随后为获得横向布局、分组配色、精确数值标注与易感性尺度标记而按 nature-figure 规范重绘为终版。

### 6.3 建议但本次未生成的图（可选补充）

- **MR 散点图**（每个关键 MR 对的 per-SNP β 散点 + IVW 拟合线）是 MR 论文常见的强化图。但工作区内仅保留了 clump 输入（`SNP\tP`，无 β）与聚合后的 MR 结果，**未保留逐 SNP 的 exposure–outcome 合并数据集**，因此无法在不重跑 harmonisation 的前提下真实复现。若需要，可重跑 `code/` 中的 harmonisation 流程并导出逐 SNP 数据，即可补出 3–4 张 MR 散点图（建议放 Supplementary）。
- **MR-PRESSO 漏斗图**（离群 SNP 可视化）同样依赖上述逐 SNP 数据，可作为后续补充。
