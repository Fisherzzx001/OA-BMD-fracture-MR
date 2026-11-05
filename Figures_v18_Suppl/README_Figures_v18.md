# Figures_v18_Suppl — P1 补充图（S9–S12）

2026-10-16 依据 `docs/107_图件核验与补充图规划.md` 的 P1 清单绘制，视觉规范与
`Figures_v17_Nature/` 完全一致（Okabe–Ito、Arial ≥7pt、去上/右边框、PNG 600 dpi +
矢量 PDF + TIFF 300 dpi）。

| 文件 | 面板 | 生成脚本 | 一手数据源 |
|---|---|---|---|
| `SupplFig9_instrument_flow.*` | 五步流程 + 逐对流失表 | `code/112_fig_S9_instrument_flow.py` | `*_plink.clumped`、`strict_clump_summary.csv`、`harmonised/*.csv` |
| `SupplFig10_MVMR_decomposition.*` | (a) 配对森林 (b) 衰减%条形 | `code/110_fig_S10_mvmr_decomposition.py` | `mvmr_height_grid.csv` (spec=EUR\|n_h=100) |
| `SupplFig11_robustness_grid.*` | 10 对 × 7 估计器网格 | `code/111_fig_S11_robustness_grid.py` | `AUDIT_mr_recomputed.csv`、`harmonised/*.csv`、`mrlap_overlap_results.csv`、`mvmr_height_grid.csv` |
| `SupplFig12_literature_comparison.*` | (a) 已发表 (b) 本研究 | `code/114_fig_S12_literature_comparison.py` | `AUDIT_mr_recomputed.csv` + Hartley 2022 / Funck-Brentano 2018 / Zhou 2025 |

**核验**：`code/115_audit_figures_v18.py` — 187 项独立溯源核验，0 失败
（`code/output/genetics/AUDIT_115_figures_v18.json`）。每个绘制数值另落盘于
`plotted_data_v18.json` 供复核。

**勘误记录**：PLINK `.clumped` 文件尾部有 2 个空行；真实数据行 20–35/116
（此前 `wc -l` 误计为 22–37/118）。稿件 Methods 已同步更正
（`code/113_fix_instrument_counts.py`），knee-OA 流失披露为 47/116（41%）。

---

## 2026-11-02 增补 S18–S21（NEW-B / NEW-C / NEW-D 系列）

依据 `docs/Supplementary_Figure_Design.md` 的 NEW-B/C/D 规划，纳入正式投稿包
（docs/08 已补图注、S18–S21 列入 Supplementary Material 清单；门禁 code/148 仍
COMPLIANT：显示项 8/8，正文 3998/4000）。视觉规范同上（Okabe–Ito、Arial ≥7pt、
去上/右边框、PNG 300 dpi + 矢量 PDF + TIFF 300 dpi）。生成脚本 `code/263_build_proposed_figs.py`。

| 文件 | 内容 | 一手数据源 |
|---|---|---|
| `SupplFig18_coloc_31snp.*` | 31 仪器位点 PP.H4 分裂条形（BMD↔kneeOA 蓝 / BMD↔fracture 橙）；0/31 vs 9/31；Wilcoxon exact 单侧 0.0043 / 双侧 0.0086 | `genetic_validation_coloc_contrast.csv`（PP4_BMD_OA / PP4_BMD_fracture，31 行） |
| `SupplFig19_coloc_loci.*` | WNT16/SOX6/GDF5 × 5 BMD 部位 PP.H4；诚实标注 GDF5·lumbar 反转（OA 0.725 / fracture 0.036） | `coloc_candidate_loci.csv` + `coloc_positive_control.csv`（逐 locus×site PP4） |
| `SupplFig20_robustness_grid.*` | 10 MR 对 × 5 估计量（IVW/Egger/加权中位/留一法/MR-PRESSO）方向-显著性一致性热力图，A/B/C 三组 | `AUDIT_mr_recomputed.csv`（beta / egger_beta / wmed / loo_min-max / presso_Q_p / p） |
| `SupplFig21_sos_hypothesis.*` | SOS 机制示意图：设定点→全身骨/骨折↓ 与 软骨下重塑/膝OA↑；Hardcastle OR 3.00、Cui RR 0.71 旁注 | 稿件 §3/§4（β=−0.337、arms β=+0.075、ref [28][39][40]） |

> 注：S13–S17（r_g 森林、h² 条形、PRESSO 漏斗、MR 功效、ρ_GE 重叠）系第十二轮 P0-2 整改
> 由 `code/132_round12_display_items.py` 幂等重编号所得，源图见 `Figures_v17_Nature/`；本次未动其数据。

**数据闭合**：fracture 31-SNP 共定位屏幕此前缺独立 CSV，已自 `genetic_validation_coloc_contrast.csv`
抽取落盘 `code/output/genetics/coloc_screen_fracture_all_loci.csv`（非手填），与稿件 L79 的
9/31、Wilcoxon p=0.0043/0.0086 完全一致。
