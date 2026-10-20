#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
14_figure3_mr_forest.py — 制作 Figure 3（MR 森林图）  [Item 4 高优先级]

读取 output/genetics/mr_results.csv，绘制三块面板：
  A. 阳性对照：BMD → 骨折（遗传高 BMD 应降低骨折）
  B. 主发现：   BMD → 膝 OA（遗传高 BMD 增加膝 OA 风险）
  C. 反向/无效：膝 OA → 骨折、膝 OA → BMD（应为空）
每点用 IVW(fixed) β 与 95% CI；竖线为 0（无效）。
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "output", "genetics")
CSV = os.path.join(OUT, "mr_results_corrected.csv")
PNG = os.path.join(OUT, "Figure3_MR_forest.png")

df = pd.read_csv(CSV)

def panel(rows, title, ax, color_pos="#c0392b", color_neg="#2874a6"):
    rows = rows.reset_index(drop=True)
    y = np.arange(len(rows))[::-1]
    labels, pts, lo, hi = [], [], [], []
    for _, r in rows.iterrows():
        b, se = r.ivw, r.ivw_se
        ci = 1.96 * se
        labels.append(f"{r.exposure} → {r.outcome}\n(n={int(r.n_snp)}, p={r.ivw_p:.1e})")
        pts.append(b); lo.append(b - ci); hi.append(b + ci)
    colors = [color_neg if b < 0 else color_pos for b in pts]
    ax.errorbar(pts, y, xerr=[np.array(pts)-np.array(lo), np.array(hi)-np.array(pts)],
                fmt="o", color="black", ecolor="gray", capsize=3, ms=5)
    for yi, b, c in zip(y, pts, colors):
        ax.scatter([b], [yi], color=c, s=42, zorder=5)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=7.5)
    ax.axvline(0, color="k", lw=1, ls="--")
    ax.set_title(title, fontsize=10, loc="left")
    ax.set_xlabel("IVW effect (log-OR per 1-SD genetically predicted BMD)", fontsize=8)
    return ax

fig, axes = plt.subplots(3, 1, figsize=(8.2, 9.0))
fig.subplots_adjust(hspace=0.35, left=0.34, right=0.96, top=0.94, bottom=0.07)

# A 阳性对照
pc = df[(df.outcome == "fracture") & (df.exposure.str.startswith("BMD"))]
panel(pc, "A. Positive control — genetically higher BMD protects against fracture",
      axes[0], color_neg="#2874a6", color_pos="#2874a6")

# B 主发现
main = df[(df.outcome == "OA_knee") & (df.exposure.str.startswith("BMD"))]
panel(main, "B. Main finding — genetically higher BMD increases knee OA risk",
      axes[1], color_neg="#c0392b", color_pos="#c0392b")

# C 反向/无效
rev = df[df.exposure == "OA_knee"]
panel(rev, "C. Reverse / null — knee OA does not raise BMD or fracture",
      axes[2], color_neg="#7f8c8d", color_pos="#7f8c8d")

axes[0].set_xlim(-1.6, 0.2)
axes[1].set_xlim(-0.05, 0.18)
axes[2].set_xlim(-0.1, 0.2)

# 注记
fig.text(0.02, 0.02,
         "Points: IVW (fixed-effects) estimate; bars: 95% CI. "
         "Blue = protective, Red = risk-increasing. Dashed line = null.\n"
         "Heterogeneity high (I² 51–85%); random-effects IVW in Table 3.",
         fontsize=7, color="#555")

fig.savefig(PNG, dpi=200)
plt.close(fig)
print(f"✅ Figure 3 -> {PNG}")
