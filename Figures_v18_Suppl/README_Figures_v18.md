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
