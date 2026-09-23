# Figures Package v17 — Nature/Science 顶刊规范配图

**稿件**: Bone Mineral Density and Osteoarthritis Share Genetic Architecture, but Evidence for an
Osteoarthritis-to-Fracture Causal Pathway Remains Limited

## 本包内容

| 文件 | 内容 |
|---|---|
| `Figure1`–`Figure10` | 10 张主图（每张含 PNG 600 dpi / 矢量 PDF / TIFF 300 dpi RGB-deflate） |
| `GraphicalAbstract` | 图文摘要 |
| `SupplFig1`–`SupplFig7` | 7 张补充图 |
| `plotted_data.json` | **每张图的每一个绘制数值**（机器可读，供独立复核） |
| `AUDIT_REPORT.txt` | 独立数值审计报告（**845 项比对，0 失败**（含 15 项图注数值断言核验）） |
| `Figures_Original_v14/` | 上一版全部原图（原样保留，未改动） |

## 相比上一版（Figures_Beautified）的变化

### 新增 4 张图（加强表达，全部基于已完成的真实统计输出）
- **Figure 10 — MRlap 样本重叠校正（§2.8）**：此前 v16 的 MRlap 结果只有文字与表，没有图。
  四面板：(a) 校正后 |β| 变化（中位 +20.5%，全部远离零）；(b) 跨性状 LDSC 截距 λ（|λ|≤0.036）；
  (c) 观测 vs 校正差异检验（9/10 显著）；(d) MRlap 独立 h² 与主 LDSC 一致性（r=0.999）。
- **SupplFig S5 — OA_knee h² 对代表 N 的敏感性**：Supplement S3 表的图形化（h² ∝ 1/N 确定性重标度）。
- **SupplFig S6 — fracture N 选择的证据链**：Table 1 §§ 脚注的图形化（4 个候选 N + 有效 N·K(1−K)）。
- **SupplFig S7 — 多效性/异质性诊断**：10 对 MR 的 I²、MR-Egger 截距检验、PRESSO 离群数全景。

### 全部 18 张图按顶刊规范重制
- **字体 Arial，最小 7 pt**（轴标签/刻度），面板标签小写粗体 **a, b, c**
- **Okabe–Ito 色盲安全配色**；蓝=保护、红=风险、灰=不显著（全文一致）
- 去除上/右边框、最小网格；双栏宽 183 mm 基准
- **PNG 600 dpi（线稿标准）+ 可编辑矢量 PDF（fonttype 42）+ TIFF 300 dpi RGB 无损压缩**
- 无缺失字形（Arial 全覆盖）

## 数据来源（无一编造）
所有数值直接读取 `code/output/genetics/` 与 `code/output/` 的分析输出：
`rg_matrix_13.csv`、`rg_results_long_13.csv`（=Table 2）、`rg_h2_trueN.csv`、
`N_sensitivity_decision.csv`、`oa_h2_liability.csv`、`AUDIT_mr_recomputed.csv`（=Table 3）、
`harmonised/*.csv`（逐 SNP 散点/漏斗）、`mr_power.csv`（=Table 4）、`overlap_statistics.csv`、
`mrlap_overlap_results.csv`、`gefos_benchmark.json`、`harmonised/_loo.csv`。

## 数值真实性核验（第 3 项要求）
`audit_figures_v17.py`（包外随代码提交）**不依赖绘图代码**，从源文件独立重算每个绘制值：
- Figure 1：169 格 r_g 与 BH-FDR q 值逐一比对
- Figure 3/6/7/S1：IVW 与 MR-Egger 从 harmonised 逐 SNP 数据重算，与 Table 3 完全一致；
  加权中位数取 AUDIT 管道已发表值
- Figure 8：功效曲线 z=2.80 处重算；10 对 power/MDE/K 与 Table 4 一致
- Figure 9：ρ_GE 与带符号 IVW z、10/10 符号一致性复核
- Figure 10：百分比变化、λ、p_difference、h² 一致性全部重算
- 稿件 Table 1–4 报告值与源数据交叉核对

**结果：830 项比对，0 失败（详见 AUDIT_REPORT.txt）。**

> 注：Figure 10 中 MRlap 的绝对效应量与 Table 3 不同（MRlap 输入使用代表 N），
> 图注与正文均已声明只解释"对内"的观测-校正对比，不与 Table 3 数值混用。
