#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
22_correct_figure3.py
=========================================================================
把 Figure 3 重画到 **N 校正后** 的效应量标度上（与稿件 Table 3 / §3.3 / §3.7 一致）。

标度规则（同 21_mr_corrected_N.py）：
  · BMD 为暴露           -> β, SE × 1/F   （F = sqrt(400000/31873) = 3.5426）
  · BMD 为结局           -> β, SE × F
  · OA_knee→fracture     -> × 1（两者都无需校正）
p 值与 n_snp 不变（β、SE 同因子缩放）。

输出：mr_results_corrected.csv；覆盖 Figure3_MR_forest.png
=========================================================================
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "output", "genetics")
CSV = os.path.join(OUT, "mr_results.csv")
CSV_OUT = os.path.join(OUT, "mr_results_corrected.csv")
PNG = os.path.join(OUT, "Figure3_MR_forest.png")

F = float(np.sqrt(400_000.0 / 31_873.0))     # 3.5426
INV = 1.0 / F                                 # 0.2823


def scale_factor(exp, out):
    exp_bmd, out_bmd = str(exp).startswith("BMD"), str(out).startswith("BMD")
    if exp_bmd and not out_bmd:
        return INV          # BMD 暴露被低估 -> MR β 高估 -> 缩小
    if out_bmd and not exp_bmd:
        return F            # BMD 结局被低估 -> MR β 低估 -> 放大
    return 1.0


def main():
    df = pd.read_csv(CSV)
    # 缩放所有 beta 类列
    beta_cols = [c for c in df.columns
                 if c.startswith("ivw") and c.endswith(("_re", "")) is False] \
                + [c for c in ("egger", "wmed") if c in df.columns]
    se_like = [c for c in df.columns if c.endswith("_se")]
    for i, r in df.iterrows():
        f = scale_factor(r.exposure, r.outcome)
        if f == 1.0:
            continue
        for c in df.columns:
            if c in ("exposure", "outcome", "note", "n_snp"):
                continue
        if c.endswith("_se") or c in ("ivw", "egger", "wmed", "ivw_radial"):
            df.at[i, c] = r[c] * f
    df.to_csv(CSV_OUT, index=False)
    print("✅ 校正后 MR 结果 ->", CSV_OUT)
    print(df[["exposure", "outcome", "n_snp", "ivw", "ivw_se", "ivw_p"]].to_string(index=False))

    # ---------------- 重画 Figure 3 ----------------
    def panel(rows, title, ax, color, xlabel):
        rows = rows.reset_index(drop=True)
        y = np.arange(len(rows))[::-1]
        labels, pts, lo, hi = [], [], [], []
        for _, r in rows.iterrows():
            b, se = r.ivw, r.ivw_se
            ci = 1.96 * se
            labels.append(f"{r.exposure} → {r.outcome}\n(n={int(r.n_snp)}, p={r.ivw_p:.1e})")
            pts.append(b); lo.append(b - ci); hi.append(b + ci)
        ax.errorbar(pts, y, xerr=[np.array(pts) - np.array(lo), np.array(hi) - np.array(pts)],
                    fmt="o", color="black", ecolor="gray", capsize=3, ms=5)
        for yi, b in zip(y, pts):
            ax.scatter([b], [yi], color=color, s=42, zorder=5)
        ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=7.5)
        ax.axvline(0, color="k", lw=1, ls="--")
        ax.set_title(title, fontsize=10, loc="left")
        ax.set_xlabel(xlabel, fontsize=8)
        # 动态坐标范围，确保 0 在内
        xmin, xmax = min(lo), max(hi)
        pad = 0.15 * (xmax - xmin if xmax > xmin else 1.0)
        ax.set_xlim(min(xmin - pad, -pad * 0.2), xmax + pad)

    fig, axes = plt.subplots(3, 1, figsize=(8.2, 9.0))
    fig.subplots_adjust(hspace=0.42, left=0.34, right=0.96, top=0.94, bottom=0.08)

    panel(df[(df.outcome == "fracture") & (df.exposure.str.startswith("BMD"))],
          "A. Positive control — genetically higher BMD protects against fracture",
          axes[0], "#2874a6",
          "IVW effect: log-OR of fracture per 1-SD higher genetically predicted BMD")
    panel(df[(df.outcome == "OA_knee") & (df.exposure.str.startswith("BMD"))],
          "B. Main finding — genetically higher BMD increases knee OA risk",
          axes[1], "#c0392b",
          "IVW effect on knee OA per 1-SD higher genetically predicted BMD")
    panel(df[df.exposure == "OA_knee"],
          "C. Reverse / null — knee OA does not raise BMD or fracture",
          axes[2], "#7f8c8d",
          "IVW effect per 1-SD genetically predicted knee OA liability")

    fig.text(0.02, 0.02,
             "Points: IVW (fixed-effects); bars: 95% CI. Blue = protective, "
             "Red = risk-increasing, Grey = null. Dashed line = null.\n"
             "Effect sizes are on the N-corrected scale (BMD N = 31,873; see §3.7). "
             "Heterogeneity high (I² 51–85%); random-effects IVW in Table 3.\n"
             "Panel A: OR per 1-SD lower BMD = 1.40 (total) / 1.42 (femoral neck), "
             "vs GEFOS benchmark 1.55.",
             fontsize=7, color="#555")

    fig.savefig(PNG, dpi=200)
    plt.close(fig)
    print("✅ Figure 3 (N-corrected) ->", PNG)


if __name__ == "__main__":
    main()
