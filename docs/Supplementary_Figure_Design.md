# 补充图设计 · 对标 2 区以上顶刊（OAC / ARD / Nat Rheumatol / AJHG 级）

> 设计原则（依据 `nature-figure` skill 的 Nature 标准 + 本稿件已落地的 Okabe–Ito 配色）：
> - 字体 Arial ≥7 pt；面板标签 a/b/c 8 pt 粗体；去上/右边框；矢量 PDF + 600 dpi PNG + 300 dpi TIFF。
> - 色盲安全：蓝 `#0072B2` / 橙 `#E69F00` 为主对；红仅用于"风险升高"，不与绿同框。
> - 每张图只讲一个点；最后一面板做"对照/综合"。
> - **所有数值均来自一手 CSV，绝不编造**（见各设计的"数据来源"）。
>
> 以下 4 张图针对本稿当前"图已够、但综合表达不足"的短板（审计发现四）。

---

## 图 NEW-A · Causal × Correlational 三角验证象限图（首推，已出图）

**对应短板**：MR（因果）× LDSC r_g（相关）同向性目前只在文字里说，缺一张"三角验证"综合图。
**对标**：*Nature Genetics* / *AJHG* 常见的"因果效应 vs 遗传相关"对齐图。

- **面板**：单面板散点象限图。
  - X 轴 = LDSC r_g（exposure↔outcome，来自 `rg_results_long_13.csv`）；
  - Y 轴 = MR IVW β（N 校正尺度，来自 `AUDIT_mr_recomputed.csv`）；
  - 10 个点按 outcome 分组着色：fracture 控制（蓝）、knee-OA（橙）、reverse/BMD（灰）；
  - 四象限浅色阴影 + x=0 / y=0 参考线，标注"MR 与 r_g 同向"。
- **数据来源**：`code/output/genetics/rg_results_long_13.csv` + `code/output/genetics/AUDIT_mr_recomputed.csv`（已生成：`Figures_Proposed/FigNEW_A_triangulation.png`）。
- **图注要点**：n=10 MR 对；点 = IVW 估计；误差棒 = 95% CI；阴影 = MR 与 r_g 同号象限；说明 BMD→fracture 负向、BMD→OA 正向、反向对近零且 r_g 同向，支持"因果与相关一致"。

---

## 图 NEW-B · 共定位端点分裂条形图（高价值，待数据落盘）

**对应短板**：coloc "knee-OA 0/31 vs fracture 9/31" 是核心卖点，却只在文字 + 表，缺可视化。
**对标**：*Nature Genetics* coloc 结果常用"逐位点 PP.H4"条形/点阵。

- **面板**：横向条形图，y = 31 个 BMD 工具位点；每位点两根条 = PP.H4(knee-OA) vs PP.H4(fracture)；
  9 个 fracture 共定位位点高亮；附 Wilcoxon p=0.0043 标注。
- **数据来源（2026-11-02 复核）**：仓库中 coloc 产物有两套——
  ① **候选基因位点**级：`coloc_candidate_loci.csv`（BMD×OA_knee，含 `PP4`）与 `coloc_positive_control.csv`（BMD×fracture，含 `PP4`），按 `locus`（WNT16/SOX6/GDF5…）× `bmd_site` 组织；两者 PP.H4 已可对齐（如 WNT16：OA_knee≈0.002、fracture≈1.000），**可直接画"同一 BMD 位点 vs 两结局"的共定位分裂**（约 15 个命名候选位点）。
  ② **31 个基因组索引 SNP** 级屏幕：仅 `coloc_screen_all_loci.csv` / `supptable_S6_coloc.md` 提供了 **BMD×OA_knee** 的逐 SNP PP.H4（31 行，max 0.799 → 0/31）；而 `genetic_validation_coloc_contrast.csv` 已含 **BMD×fracture 的 31-SNP 逐位点 PP.H4**（`PP4_BMD_fracture`，9/31 ≥0.8）。
  → **2026-11-02 已导出**：`code/output/genetics/coloc_screen_fracture_all_loci.csv`（31 行，从 contrast 文件抽取骨折臂；provenance = 同一次 coloc 计算）。
     两版 NEW-B 均已生成：`FigNEW_B31_coloc_31snp`（31-SNP 分裂，0/31 vs 9/31，Wilcoxon exact 单侧 p=0.0043/双侧 0.0086）与 `FigNEW_Bcand_coloc_loci`（3 命名位点×5 部位候选版；含 GDF5·lumbar 反转 OA PP.H4=0.73 vs fracture 0.04 的诚实标注）。
- **图注要点**：n=31 仪器位点；条 = PP.H4；0/31 knee-OA vs 9/31 fracture 共定位；配对 Wilcoxon 单侧 p=0.0043；同一仪器集与骨而非关节共享因果变异。

---

## 图 NEW-C · SOS 机制假设示意图（面向 OAC 临床受众）

**对应短板**：SOS（Systemic Osteogenic Set-point）机制仅文字，临床读者（OAC 主力）缺一张可直觉理解的图。
**对标**：*Osteoarthritis and Cartilage* / *Ann Rheum Dis* 常用的"机制示意"（BioRender/Inkscape SVG）。

- **面板**：单张示意流程图。
  - 上游：多基因骨形成设定点（polygenic osteogenic set-point，可操作化为 areal BMD）；
  - 下游双分支：① 保护全身骨 → 骨折风险↓（正向对照）；② 关节软骨下重塑 → OA 风险↑；
  - 旁注临床对应：Hardcastle  triad（骨赘/附着点骨赘/高骨量，OR 3.00）与 Cui 双膦酸盐反例（RR 0.71，指向"形成侧"而非"密度"）。
- **数据来源**：无数值，纯概念图；按 `nature-figure` 建议用 SVG（BioRender/Inkscape）导出，嵌入 PDF。
- **图注要点**：SOS 为可检验假设而非已证机制；非负重 arms 结果（β=+0.075）将其与局部机械负荷区分；已通过 F6–F8。

---

## 图 NEW-D · 多方法稳健性三角热力图（综合收口）

**对应短板**：IVW/Egger/中位数/PRESSO/MRlap/MVMR/coloc 七条证据链分散在多处，缺一张"稳健性汇总"。
**对标**：*Nature* 方法学文章的"evidence triangulation"矩阵。

- **面板**：热力图，行 = 10 个 MR 对（或 3 类假设），列 = 7 种方法；单元格 = 方向一致性/显著性（绿=一致显著，灰=空，红=矛盾）；
  附底部"综合判读"一列。
- **数据来源**：`AUDIT_mr_recomputed.csv`（IVW/Egger/中位数/PRESSO/MRlap）、`mvmr_height_grid.csv` + `mvmr_bmi_grid.csv`（MVMR）、coloc 输出（coloc）。均为已算文件，汇总即可。
- **图注要点**：颜色 = 方向与显著性一致度；展示 BMD→OA 正向在 5 站点跨 5+ 方法一致；反向/ OA→fracture 空且弱功效；呼应正文"以估计量一致性而非固定效应 p 值下结论"。

---

## 落地优先级

| 图 | 价值 | 数据就绪度 | 建议 |
|---|---|---|---|
| NEW-A 三角验证 | 高 | ✅ 已生成 | 直接采用 |
| NEW-D 稳健性矩阵 | 高 | ✅ 汇总即可 | 下一轮生成 |
| NEW-B 共定位分裂 | 最高（卖点） | ⚠️ 需先导出逐位点 PP.H4 | 先补 `code/30*` 导出，再生成 |
| NEW-C SOS 示意 | 中（受众适配） | ✅ 概念图 | 用 BioRender 手绘 SVG |

> 所有新增图沿用 `Figures_v17_Nature` 的视觉规范（Okabe–Ito、Arial≥7pt、去上右边框），与正稿图保持同一视觉语言。
