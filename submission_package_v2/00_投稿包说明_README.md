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
| `Figures/` | **美化/新增后的最终图（7 张主图 + 2 张补充图，均 PNG + 矢量 PDF）**：`Figure1_rg_heatmap_13`、`Figure2_forest_key_results`、`Figure3_MR_forest`、`Figure4_h2_bar`、`Figure5_causal_network`、**`Figure6_MR_scatter`、`Figure7_MR_funnel`（新增）**、**`SupplFig1_MR_scatter_remaining`、`SupplFig2_leave_one_out`（新增补充图）** | 投稿用最终图（图注见稿件 §Figure legends） |
| `Figures_Original/` | **原始图（保留，未改动）**：原 `Figure1_rg_heatmap_13.png`、`Figure2_forest_key_results.png`、`Figure3_MR_forest.png` | 版本对照 / 保留原始输出 |
| `Tables/` | `rg_matrix_13.csv`、`mr_corrected_N.csv`、`mr_presso_results.csv`、`oa_h2_liability.csv`、`strict_clump_summary.csv`、`rg_h2_trueN.csv` | 支持性数据（rg 矩阵、MR、MR-PRESSO、h²、严格 clump 汇总、h²_trueN） |
| `Tables/Harmonised_SNP/` | **逐 SNP harmonised 数据（新增）**：10 个 `<exposure>__<outcome>.csv`（每行一个工具变量：SNP/CHR/BP、exposure 与 outcome 的 β/SE/P、EAF、等位基因翻转标记、Wald 比值、MR-PRESSO 逐 SNP 离群 P 与标记、留一法估计）+ `_validation.csv`（重跑 vs 已发表 mr_corrected_N.csv 的一致性校验，Δβ≈1e-17）+ `_loo.csv` | MR 散点图 / 漏斗图 / 留一法图的底层数据，供审稿人复核 |
| `Supplementary/` | `10_参考文献与原始GWAS来源.md`、`11_补充材料_OA_selfrep_Z校正.md`、`12_补充材料_MR_PRESSO.md`、**`26_export_harmonised_snp.py`（新增，逐 SNP 导出与校验脚本）**、`README.md`（代码仓库说明） | 补充材料 / 方法学溯源 |

> 说明：稿件中的 Table 1–3（h²、rg 模块、MR 估计）已内嵌于 `Manuscript.docx`；`Tables/` 提供其底层数据以便审稿人与编辑复核。

## 2. 投稿前数据核实（✅ 已完成，2026-09-17）

- [x] ~~在 GWAS Catalog 中核实 BMD 原始汇总数据的 accession~~ — **已核实**：GCST90568457 = DXA-BMD femur total（UKB field 23291），N=31,873 欧洲血统 UKB（Qian et al. 2025，PMID 41206123，summary statistics 公开）。同系列逐部位 accession/N 已一并核验并写入稿件 §2.1（lumbar GCST90568448 N=30,449；total body GCST90568455 N=31,986；femoral neck GCST90568458 N=32,017；legs/pelvis/femur-total N=31,873）。逐部位 N 对校正因子的影响已量化：headline OR（1.40/1.42）变化 <0.3%，最大偏差 lumbar 2.3%，所有舍入报告值与全部 z/p/方向不变（§3.7）。

**当前状态：无遗留投稿前待办。**

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

### 6.3 MR 散点图 / MR-PRESSO 漏斗图（本次已补齐，Figure 6 / Figure 7 + 2 张补充图）

上一版说明中提到"缺逐 SNP 数据、无法真实绘制 MR 散点图"。**本次已通过重跑 harmonisation 补齐**：

1. **重跑 harmonisation（`code/26_export_harmonised_snp.py`，已收入 `Supplementary/`）**
   - 直接从磁盘上的原始 GWAS 汇总统计重建：12 个 munge 后 LDSC sumstats（Z/N/MAF → β/SE）+ 骨折原始 TSV（真实 log-OR β/SE），**不使用任何模拟或插值**；
   - 仪器集 = 严格 PLINK 等价 clump（r²<0.001 / 10 Mb，与稿件 §2.6 主口径一致）；
   - 效应标度 = N 校正因子（BMD ×3.5426，OA_selfrep ×0.8695），与稿件 §3.7 完全一致；
   - **一致性校验**：重跑的 IVW β / SE / P / K 与已发表的 `mr_corrected_N.csv` 逐位比对，10/10 对 Δβ ≈ 1e-17（浮点噪声级），全部 ✅ 后才允许作图（见 `Tables/Harmonised_SNP/_validation.csv`）。

2. **新增图（Nature 风格，PNG + 矢量 PDF）**
   - **Figure 6 — MR 散点图（2×2）**：total-BMD→fracture（阳性对照）、femoral-neck-BMD→fracture（阳性对照）、total-BMD→knee OA（主发现）、knee-OA→total-BMD（代表性零结果）。每点一个工具变量（含 x/y 95% CI 误差棒），IVW（过原点实线）、MR-Egger（虚线）、加权中位数（点线）三条拟合线；MR-PRESSO 离群 SNP 用红色空心圆标出。
   - **Figure 7 — MR-PRESSO 漏斗图（2×2）**：Wald 比值 × 1/SE 精度，附 IVW 竖线、MR-Egger 竖线与 95% pseudo-CI 漏斗线；离群 SNP 同样红色标出——直观呈现 §3.6 的异质性与离群结构。
   - **Supplementary Figure S1 — 其余 6 对散点图（2×3）**：4 个 BMD→膝 OA 正向对 + 2 个零结果对。
   - **Supplementary Figure S2 — 留一法稳定性（全部 10 对）**：每行一对，黑点为全 SNP IVW 估计（含 95% CI），蓝色带为逐一剔除每个 SNP 后的估计范围；10 对全部无符号翻转。

3. **稿件同步更新**：§3.4/§3.5/§3.6 正文加入 Figure 6 / Figure 7 / SupplFig S1–S2 交叉引用；图注区新增 Figure 6、Figure 7 与补充图 S1、S2 的完整图例；`Manuscript.docx`、`Submission_MainText.docx` 已重新生成（盲法正文校验无作者信息）。

---

## 7. 语言润色（去 AI 痕迹，v5）

对包内**全部英文文字**（稿件正文、图注、表注、投稿信、medRxiv/OSF 英文声明）做了一次去 AI 痕迹的语言润色，按 SCI 顶刊写作习惯执行。

**处理原则**：只改措辞与节奏，**不动任何数据、统计量、引用与结论**。已用脚本逐条替换并做数字多重集比对校验：除下文"顺带修正"一项外，全部数字出现次数与旧版完全一致（0.337 / 1.40 / 1.42 / 1.55 / 0.2823 / 3.5426 / 690,150 等逐一核对无变化）。

**具体处理**：

| 类别 | 处理 |
|---|---|
| 破折号（em dash） | 稿件正文 37 处 → 1 处，投稿信 2 处 → 0 处；改为逗号、冒号或拆句 |
| 宣传性形容词 | 删除或替换 exaggerated 用词：`clinically entangled` → 直述关联；`coherent three-block` → `three blocks`；`exquisitely sensitive` → 直接说明依赖；`canonical` → `established`；`Crucially/`critically` 等强调词删除 |
| 加粗强调 | 正文内加粗（\*\*First/Second/Third\*\*、\*\*The pipeline is therefore...\*\* 等）改为自然行文；仅保留摘要标签、图注编号、表头等结构性加粗 |
| 斜体强调 | 删除用于"强调语气"的斜体（\*increased\*、\*same\*、\*quantity\*），保留物种/基因名与期刊名的规范斜体 |
| 三段式与同义轮换 | `not only...but`、`rule of three` 排比改写为直述；统一用词（BMD→OA 一律 `raise`，不再 increased/elevated/augmented 轮换） |
| 句尾 -ing 尾巴 | `..., reproducing/indicating/underscoring...` 类附着结构改为独立句 |
| 空泛结论 | 删除 `The pipeline is therefore directionally and quantitatively validated` 等口号式总结，改为陈述事实 |
| 弯引号 | 投稿信标题弯引号 → 直引号 |

**顺带修正的一处数据不一致**（非语言问题，需你知晓）：§4.3 原写 `β = +0.086, p = 0.402`，与 Table 3 及 §3.5 的 `β = +0.092, p = 0.395`（`mr_corrected_N.csv` 主口径）不一致，为旧版遗留值。已统一为 **+0.092 / 0.395**。方向、显著性与结论均不变。

## 8. 图片-数据一致性核验（v6 新增，2026-09-18）

全部 9 张图（Figure 1–7 + SupplFig S1–S2）经独立核验脚本 `Supplementary/verify_figures.py`
逐项比对底层统计文件，共 **120 项检查全部通过**：

- **Fig 1**：热图 169 个单元格的 r_g 逐位等于 LDSC 长表输出（`rg_results_long_13.csv`）；
  BH-FDR 星号用独立实现重算，与 LDSC 自带 q_fdr 列一致（最大差 <1e-9）。
- **Fig 2 / Fig 3**：每个点与误差棒 = 源 CSV 的估计值 ± 1.96×SE；图内硬编码的
  OR 1.40 (1.34–1.46) 与 1.42 (1.34–1.51) 由 β/SE 复算吻合。
  **本轮修正**：原硬编码坐标轴范围截断了部分 95% CI（Fig2a 6/10 条下缘、Fig2d 1/3、
  Fig3c 3/3 上缘被切）。已改为数据驱动轴范围并把实际范围落盘
  （`Supplementary/figure_axis_limits.json`），核验脚本据此确认全部误差棒完整可见。
- **Fig 4**：13 个 h² 与源 CSV 一致；三个易感性尺度 h²（0.0935/0.0726/0.0620）按
  Lee 2011 Eq 23 独立复算吻合。
  **本轮修正**：随包分发的 `Tables/oa_h2_liability.csv` 原为旧值（OA_selfrep 用已废弃
  N=350,000、K=0.25 → 0.0514，且缺 fracture 行），已按权威口径重写
  （OA_selfrep 0.0223/K=0.083 → 0.0726；fracture 0.0258/K=0.143 → 0.0620），与稿件
  Table 1 及图 4 完全一致。
- **Fig 6 / 7 / S1 / S2**：散点逐 SNP 来自重跑 harmonisation 的真实数据；从逐 SNP
  数据独立重算的 IVW/Egger/加权中位数与稿件 Table 3 逐位一致（最大 Δβ = 1.4e-15）；
  离群点数与 MR-PRESSO 结果一致；留一法 10 对均无符号翻转、无跨零。
  **本轮修正**：二分类性状轴标签精确为 "log-OR per allele"（原统称 "per allele"）。
- **Fig 5**：图上标注（OR 1.40–1.42、r_g +0.12~+0.25、OA→fracture ns）逐条复核。
  另：图 1 图注 "colour scale is symmetric about zero" 改为 "diverging and centred
  on zero"（实际数值范围 -0.65~+1.00 不对称，仅色相以 0 为中心）。

上述修正只影响**图形呈现与随包补充表**，未改动任何统计量、OR、p 值或结论。
