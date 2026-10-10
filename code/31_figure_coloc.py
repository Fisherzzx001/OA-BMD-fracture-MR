#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
31_figure_coloc.py — Supplementary Figure S8：共定位结果与功效
=========================================================================
Panel a: 31 个 BMD 工具变量基因座的 PP.H4 vs PP.H3（0.8 阈值线）
Panel b: MR 预测的结局侧 |z| vs 实际观测 |z|（OA 灰 / fracture 红），附全基因组阈值

输出（300 dpi，三种格式）：
  Figures_v17_Nature/SupplFig8_coloc_power.{png,pdf,tif}
  submission_package_v16/Supplementary_Figures/SupplFig8_coloc_power.{png,pdf,tif}
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.abspath(__file__))
OUTG = os.path.join(ROOT, "output", "genetics")
FIGDIR = os.path.join(os.path.dirname(ROOT), "Figures_v17_Nature")
PKG = os.path.join(os.path.dirname(ROOT), "submission_package_v16", "Supplementary_Figures")

# Okabe–Ito（与 Figures_v17 其余图一致）
OI = {"blue": "#0072B2", "orange": "#E69F00", "vermillion": "#D55E00",
      "grey": "#999999", "black": "#000000", "sky": "#56B4E9"}

GW_Z = 5.45
MR_BETA_OA = 0.02312
MR_BETA_FX = -0.33713


def main():
    os.makedirs(FIGDIR, exist_ok=True)
    os.makedirs(PKG, exist_ok=True)

    scr = pd.read_csv(os.path.join(OUTG, "coloc_screen_all_loci.csv"))
    pw = pd.read_csv(os.path.join(OUTG, "coloc_power.csv"))
    cand = pd.read_csv(os.path.join(OUTG, "coloc_candidate_loci.csv"))

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.6))
    fig.patch.set_facecolor("white")

    # ---------------- Panel a ----------------
    ax = axes[0]
    ax.scatter(scr.PP3, scr.PP4, s=34, facecolors="none", edgecolors=OI["blue"],
               linewidths=1.2, zorder=3, label="BMD_total instrument loci (n=31)")
    w = scr[scr.index_snp == "rs3779381"]
    if len(w):
        ax.scatter(w.PP3, w.PP4, s=70, color=OI["vermillion"], zorder=4,
                   label="WNT16 (rs3779381)")
    # 三个点名基因座（其中 SOX6 / GDF5 不在 31 个 BMD 工具集中，用候选表补画）
    _lab_off = {"WNT16": (8, 4), "SOX6": (-8, -22), "GDF5": (24, 1)}
    for lab, rs in (("WNT16", "rs3779381"), ("SOX6", "rs7108738"), ("GDF5", "rs224333")):
        sub = cand[(cand.lead_rs == rs) & (cand.bmd_site == "total")]
        if len(sub):
            r = sub.iloc[0]
            ax.scatter(r.PP3, r.PP4, s=70, marker="D",
                       color=OI["vermillion"] if rs == "rs3779381" else OI["orange"],
                       zorder=5)
            ax.annotate(lab, (r.PP3, r.PP4), textcoords="offset points",
                        xytext=_lab_off[lab], fontsize=9, color=OI["black"])
    ax.axvline(0.8, color=OI["grey"], ls="--", lw=1)
    ax.axhline(0.8, color=OI["grey"], ls="--", lw=1)
    ax.text(0.82, 0.03, "PP.H3 = 0.8\n(distinct variants)", fontsize=8,
            color=OI["grey"], rotation=90, va="bottom")
    ax.text(0.03, 0.83, "PP.H4 = 0.8 (shared variant)", fontsize=8, color=OI["grey"])
    ax.set_xlabel("PP.H3  (both traits, distinct causal variants)")
    ax.set_ylabel("PP.H4  (shared causal variant)")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_title("a  Colocalisation of BMD with knee OA", loc="left", fontsize=11,
                 fontweight="bold")
    ax.legend(loc="center right", fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)

    # ---------------- Panel b ----------------
    ax = axes[1]
    m = scr.merge(pw[["snp", "z_oa_exp", "z_oa_obs", "z_fx_exp", "z_fx_obs"]],
                  left_on="index_snp", right_on="snp", how="left")
    ax.scatter(np.abs(m.z_oa_exp), np.abs(m.z_oa_obs), s=34, facecolors="none",
               edgecolors=OI["blue"], linewidths=1.2, label="knee OA (n=1.21M)")
    ax.scatter(np.abs(m.z_fx_exp), np.abs(m.z_fx_obs), s=34, color=OI["vermillion"],
               alpha=0.85, label="any fracture (n=265k, control)")
    lim = [0.01, 12]
    ax.plot(lim, lim, color=OI["black"], lw=0.8, ls=":")
    ax.axhline(GW_Z, color=OI["grey"], ls="--", lw=1)
    ax.axvline(GW_Z, color=OI["grey"], ls="--", lw=1)
    ax.text(0.015, GW_Z + 0.25, "|z| = 5.45  (p = 5×10⁻⁸)", fontsize=8, color=OI["grey"])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.08, 12)
    ax.set_ylim(0.08, 12)
    ax.set_xlabel("|z| expected in the outcome under the MR estimate")
    ax.set_ylabel("|z| observed in the outcome")
    ax.set_title("b  Power of locus-level colocalisation", loc="left", fontsize=11,
                 fontweight="bold")
    ax.legend(loc="lower right", fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)

    fig.suptitle("Supplementary Figure S8", fontsize=12, fontweight="bold", x=0.008,
                 ha="left", y=1.02)
    fig.tight_layout()

    for d in (FIGDIR, PKG):
        for ext in ("png", "pdf", "tif"):
            fig.savefig(os.path.join(d, "SupplFig8_coloc_power.%s" % ext), dpi=300,
                        facecolor="white", bbox_inches="tight")
    print("written to", FIGDIR, "and", PKG)


if __name__ == "__main__":
    main()
