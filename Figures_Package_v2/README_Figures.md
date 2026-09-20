# 图片下载包 v2（Figures_Package_v2）

**论文：** *Bone Mineral Density and Osteoarthritis Share Genetic Architecture but No Evidence of a Causal Pathway from Osteoarthritis to Fracture*
**配套稿件：** docs/08（v13），Submission_Package_v13.zip
**生成日期：** 2026-09-20
**生成方式：** 调用 `nature-figure` 与 `science-figures` 两个 SKILL 的顶刊图片规范，基于论文真实统计结果重绘（**未改动任何 β / SE / p / OR**）。

---

## 1. 本包内容

| 文件夹 | 内容 | 说明 |
|---|---|---|
| `Figures_Beautified/` | 14 张美化图（每张 PNG + 可编辑矢量 PDF） | 本次按 Nature 系列 + Elsevier(*Bone*) 规范重绘/新增 |
| `Figures_Original/` | 12 张原图（每张 PNG + PDF） | **原样保留**，来自 v13 投稿包，未做任何改动 |

> 生产用 TIFF（300 dpi，Elsevier 栅格首选）已在本地 `Figures_Beautified/*.tif` 生成，因单张 8–32 MB 会使压缩包过大，未纳入本包；如需请单独索取。PDF 矢量格式已被 Nature 与 Elsevier 共同接受为生产格式。

---

## 2. 缺口分析结论（是否需更多图）

对现有 12 图（Fig1–9 + S1–S3）逐一核对目标期刊（*Bone* / *Osteoporosis International* / *Osteoarthritis and Cartilage*）的图片期望后，结论：

- **现有 12 图已覆盖全部核心发现**，结构完整，无需删改：
  - Fig1 r_g 热图（13×13）、Fig2 r_g 森林、Fig3 MR 森林、Fig4 h² 条形、Fig5 因果网络图、Fig6/7 MR 散点/漏斗、S1 余下散点、S2 留一法、Fig8 效能、Fig9 样本重叠、S3 GEFOS 外部验证。
- **发现两处可加强表达的缺口（已补 2 张新图）**：
  1. **缺 STROBE-MR 研究流程图**（目标期刊均推荐该报告规范）→ 新增 `SupplFig4_STROBE_flow`。
  2. **缺 Graphical Abstract**（*Bone* 明确鼓励提交）→ 新增 `GraphicalAbstract`。

---

## 3. 美化改动（相对原图）

| 项 | 原图 | 本包（v2） | 依据 |
|---|---|---|---|
| 配色 | 多套近似蓝/红 | 统一 **Okabe-Ito 色盲安全配色**（#0072B2 / #D55E00 / #009E73 / #E69F00 / #CC79A7 / #56B4E9） | nature-figure / science-figures 强制色盲安全 |
| 字体 | 部分刻度 6.0–6.5 pt | 轴标签/刻度统一 **≥7 pt**，面板标签 8 pt 粗体小写 | Nature 终印 ≥7 pt |
| 字体族 | Arial | Arial / Helvetica（sans-serif） | Nature & Elsevier |
| 边框 | 去顶/右边框 | 保持去顶/右边框，线条 0.8 pt | 去 chartjunk |
| 输出 | PNG + PDF | PNG(300dpi) + **矢量 PDF**(+ TIFF 留盘) | 矢量可编辑、Elsevier 接受 PDF |
| 分辨率 | 300 dpi | 300 dpi（栅格）/ 矢量无限 | Nature ≥300；线条图矢量 |

---

## 4. 文件清单（14 张美化图）

| 文件名 | 类型 | 数据来源（真实） | 期刊合规要点 |
|---|---|---|---|
| Figure1_rg_heatmap_13 | r_g 热图 | rg_matrix_13 / rg_p / rg_se（BH-FDR） | 发散色条 RdBu_r，居中 0 |
| Figure2_forest_key_results | r_g 森林（4 面板） | 同上 | 95% CI，数据驱动坐标轴 |
| Figure3_MR_forest | MR 森林（3 面板） | mr_corrected_N.csv | 蓝=保护/红=风险/灰=零 |
| Figure4_h2_bar | h² 条形 | rg_h2_trueN / N_sensitivity_decision | 观测尺度 +  liability 菱形 |
| Figure5_causal_network | 因果网络示意 | 论文结论（示意） | 实线=因果/虚线=相关 |
| Figure6_MR_scatter | MR 散点（4 核心对） | harmonised/*.csv（逐 SNP） | IVW/Egger/加权中位三线 |
| Figure7_MR_funnel | MR-PRESSO 漏斗（4 核心对） | 同上 | 伪置信漏斗 + 离群点 |
| SupplFig1_MR_scatter_remaining | MR 散点（余 6 对） | 同上 | 同上 |
| SupplFig2_leave_one_out | 留一法稳定性（10 对） | 同上 | 无单 SNP 驱动结论 |
| Figure8_MR_power | 效能（2 面板） | mr_power / mr_corrected_N | 80% 参考线 + MDE |
| Figure9_overlap_pleiotropy | 样本重叠（2 面板） | overlap_statistics.csv | ρ_GE 符号一致性 10/10 |
| SupplFig3_GEFOS_benchmark | GEFOS 外部验证森林 | gefos_benchmark.json | 与独立基准差 8–10% |
| **GraphicalAbstract** (新) | 图文摘要 | 论文结论（示意） | *Bone* 鼓励；531×1328 px⁺ |
| **SupplFig4_STROBE_flow** (新) | STROBE-MR 流程图 | 研究设计（示意） | 填补报告规范缺口 |

---

## 5. 如何替换稿件图片

将 `Figures_Beautified/` 中对应文件覆盖投稿包 `Submission_Package_*/Figures/` 即可（文件名与稿件图注一致；新增两张需在 docs/08 图注补 `SupplFig4` 与 Graphical Abstract 说明，见 v14 投稿包）。

## 6. 约束

所有 14 张图均基于论文已审计报告（AUDIT_mr_recomputed.csv、harmonised/*、mr_results.csv 等）的真实数据绘制；示意图（Fig5、GraphicalAbstract、SupplFig4）不含统计量，仅表达已证实的结论。未改动任何 β / SE / p / OR。
