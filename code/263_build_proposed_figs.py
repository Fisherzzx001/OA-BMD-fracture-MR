#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
263_build_proposed_figs.py — 生成 4 张补充图（NEW-B / NEW-C / NEW-D）+ 导出 fracture 31-SNP 屏幕
==========================================================================================
所有数值均回溯一手 CSV，绝不编造：
  · NEW-B-31       ← code/output/genetics/genetic_validation_coloc_contrast.csv (PP4_BMD_OA / PP4_BMD_fracture, 31 位仪器)
  · NEW-B-candidate← coloc_candidate_loci.csv (OA_knee) + coloc_positive_control.csv (fracture)，3 命名基因座×5 BMD 部位
  · NEW-D          ← AUDIT_mr_recomputed.csv (IVW/Egger/WM/LOO/PRESSO, 10 对)
  · NEW-C (SOS)    ← 概念图，无数值；标注数字均来自稿件已报告值（arms β=+0.075 p=1.7e-6; BMD→fracture β=−0.337;
                     Hardcastle triad OR 3.00; Cui RR 0.71; F6 τ²=0; F8 双性阳性）
  · 导出           ← coloc_screen_fracture_all_loci.csv（从 contrast 文件抽取骨折臂，附 provenance）

视觉规范：Nature / Okabe–Ito；Arial≥7pt；去上右边框；pdf.fonttype=42；双栏 7.2in。
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D
from scipy.stats import wilcoxon

ROOT = os.path.dirname(os.path.abspath(__file__))
OUTG = os.path.join(ROOT, "output", "genetics")
FIG = os.path.join(os.path.dirname(ROOT), "Figures_Proposed")
os.makedirs(FIG, exist_ok=True)

OI = {
    "blue": "#0072B2", "orange": "#E69F00", "vermillion": "#D55E00",
    "grey": "#999999", "black": "#000000", "sky": "#56B4E9",
    "lightblue": "#BDD7E9", "amber": "#EE854A", "green": "#009E73",
}

plt.rcParams.update({
    "font.family": "Arial", "font.size": 8, "axes.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "xtick.major.size": 3, "ytick.major.size": 3,
    "xtick.direction": "out", "ytick.direction": "out",
    "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05, "pdf.fonttype": 42, "ps.fonttype": 42,
})


# ============================================================ 导出 fracture 31-SNP 屏幕
def export_fracture_screen():
    src = pd.read_csv(os.path.join(OUTG, "genetic_validation_coloc_contrast.csv"))
    out = src[["index_snp", "chrom", "pos", "PP4_BMD_fracture", "PP3_BMD_fracture",
               "n_snps_fx", "top_snp_oa", "z_top_bmd_fx", "z_top_fx", "dir_fx",
               "delta_PP4"]].copy()
    out = out.rename(columns={
        "PP4_BMD_fracture": "PP4", "PP3_BMD_fracture": "PP3",
        "n_snps_fx": "n_snps", "z_top_bmd_fx": "z_bmd", "z_top_fx": "z_fracture",
    })
    out = out.sort_values("PP4", ascending=False).reset_index(drop=True)
    out.insert(0, "rank", np.arange(1, len(out) + 1))
    out_path = os.path.join(OUTG, "coloc_screen_fracture_all_loci.csv")
    out.to_csv(out_path, index=False)
    n_ge = int((out.PP4 >= 0.8).sum())
    print("[export] wrote %s  (n=31, PP.H4>=0.8: %d/31)" % (out_path, n_ge))
    return out


# ============================================================ NEW-B-31 : 31-SNP 分裂条形
def fig_new_b_31():
    c = pd.read_csv(os.path.join(OUTG, "genetic_validation_coloc_contrast.csv"))
    c = c.sort_values("PP4_BMD_fracture", ascending=True).reset_index(drop=True)
    y = np.arange(len(c))
    oa = c.PP4_BMD_OA.values
    fx = c.PP4_BMD_fracture.values
    fig, ax = plt.subplots(figsize=(7.2, 9.0))
    h = 0.38
    ax.barh(y - h / 2, oa, height=h, color=OI["blue"], label="PP.H4 (BMD ↔ knee OA)")
    ax.barh(y + h / 2, fx, height=h, color=OI["orange"], label="PP.H4 (BMD ↔ fracture)")
    ax.axvline(0.8, color=OI["grey"], ls="--", lw=1.0)
    ax.text(0.81, 0.6, "PP.H4 = 0.8\n(shared variant)", fontsize=7.5,
            color=OI["grey"], va="center")
    fx_ge = (fx >= 0.8)
    n_fx = int(fx_ge.sum()); n_oa = int((oa >= 0.8).sum())
    ax.set_yticks(y)
    labels = ["%s (%s:%d)" % (r.index_snp, r.chrom, r.pos) for r in c.itertuples()]
    ax.set_yticklabels(labels, fontsize=6.5)
    ax.set_xlim(0, 1.05)
    ax.set_ylim(-1.2, len(c) + 0.2)
    ax.set_xlabel("Posterior probability of shared causal variant (PP.H4)")
    ax.set_title("a  Endpoint-specific colocalisation at 31 BMD instrument loci",
                 loc="left", fontsize=11, fontweight="bold")
    ax.legend(loc="center right", fontsize=8, frameon=False)
    # 计数 + Wilcoxon 放左上空白区（顶部骨折高 PP.H4 条左侧 x<0.8 为空）
    ax.text(0.02, len(c) - 1.2, "knee OA: %d/31   |   fracture: %d/31   (≥ PP.H4 0.8)"
            % (n_oa, n_fx), fontsize=8.5, ha="left", color=OI["black"], fontweight="bold")
    w = wilcoxon(c.PP4_BMD_OA, c.PP4_BMD_fracture, alternative="two-sided", method="exact")
    wo = wilcoxon(c.PP4_BMD_OA, c.PP4_BMD_fracture, alternative="less", method="exact")
    ax.text(0.02, len(c) - 2.6,
            "paired Wilcoxon signed-rank (exact): one-sided p = %.4f, two-sided p = %.4f"
            % (wo.pvalue, w.pvalue), fontsize=7.5, ha="left", color=OI["black"],
            bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", pad=2.0))
    fig.subplots_adjust(left=0.30, right=0.97, top=0.95, bottom=0.10)
    for ext in ("pdf", "png"):
        ax.get_figure().savefig(os.path.join(FIG, "FigNEW_B31_coloc_31snp.%s" % ext),
                                dpi=300, facecolor="white")
    print("[NEW-B-31] kneeOA %d/31, fracture %d/31, Wilcoxon one-sided p=%.4f"
          % (n_oa, n_fx, wo.pvalue))


# ============================================================ NEW-B-candidate : 3 基因座×5 部位
def fig_new_b_cand():
    oa = pd.read_csv(os.path.join(OUTG, "coloc_candidate_loci.csv"))
    fx = pd.read_csv(os.path.join(OUTG, "coloc_positive_control.csv"))
    sites = ["total", "femur_neck", "lumbar", "arms", "leg"]
    loci = ["WNT16", "SOX6", "GDF5"]
    rows = []
    for loc in loci:
        for s in sites:
            o = oa[(oa.locus == loc) & (oa.bmd_site == s)].PP4.values
            f = fx[(fx.locus == loc) & (fx.bmd_site == s)].PP4.values
            rows.append((loc, s, o[0] if len(o) else np.nan, f[0] if len(f) else np.nan))
    df = pd.DataFrame(rows, columns=["locus", "site", "oa", "fx"])
    df["label"] = df.locus + " · " + df.site
    df = df.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(df))
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    h = 0.38
    ax.barh(y - h / 2, df.oa, height=h, color=OI["blue"], label="PP.H4 (BMD ↔ knee OA)")
    ax.barh(y + h / 2, df.fx, height=h, color=OI["orange"], label="PP.H4 (BMD ↔ fracture)")
    ax.axvline(0.8, color=OI["grey"], ls="--", lw=1.0)
    ax.text(0.81, 0.5, "PP.H4 = 0.8", fontsize=7.5, color=OI["grey"], va="center")
    # 分组分隔线
    bounds = [df[df.locus == l].index.max() for l in loci]
    for b in bounds[:-1]:
        ax.axhline(b + 0.5, color="white", lw=1.5)
    ax.set_yticks(y)
    ax.set_yticklabels(df.label, fontsize=7.0)
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("Posterior probability of shared causal variant (PP.H4)")
    ax.set_title("b  Candidate shared-effector loci (WNT16 · SOX6 · GDF5) × 5 BMD sites",
                 loc="left", fontsize=10.5, fontweight="bold")
    ax.legend(loc="lower right", fontsize=8, frameon=False)
    # GDF5 lumbar 反转诚实标注
    gli = df[(df.locus == "GDF5") & (df.site == "lumbar")].index
    if len(gli):
        gi = gli[0]
        ax.annotate("GDF5·lumbar: OA PP.H4=%.2f vs fracture %.2f\n(site-specific reversal)"
                    % (df.oa[gi], df.fx[gi]), (df.fx[gi], gi + h / 2),
                    xytext=(0.30, gi + 1.4), fontsize=6.8, color=OI["vermillion"],
                    arrowprops=dict(arrowstyle="->", color=OI["vermillion"], lw=0.7))
    fig.subplots_adjust(left=0.22, right=0.97, top=0.93, bottom=0.10)
    for ext in ("pdf", "png"):
        ax.get_figure().savefig(os.path.join(FIG, "FigNEW_Bcand_coloc_loci.%s" % ext),
                                dpi=300, facecolor="white")
    print("[NEW-B-cand] rows=%d loci=%s" % (len(df), loci))


# ============================================================ NEW-D : 多方法稳健性矩阵
def fig_new_d():
    m = pd.read_csv(os.path.join(OUTG, "AUDIT_mr_recomputed.csv"))

    def sgn(x):
        return 0 if x == 0 else (1 if x > 0 else -1)

    order = [
        "BMD_arms__OA_knee", "BMD_femur_neck__OA_knee", "BMD_leg__OA_knee",
        "BMD_lumbar__OA_knee", "BMD_total__OA_knee",
        "BMD_femur_neck__fracture", "BMD_total__fracture",
        "OA_knee__BMD_femur_neck", "OA_knee__BMD_total", "OA_knee__fracture",
    ]
    grp = {
        "BMD_arms__OA_knee": "A", "BMD_femur_neck__OA_knee": "A", "BMD_leg__OA_knee": "A",
        "BMD_lumbar__OA_knee": "A", "BMD_total__OA_knee": "A",
        "BMD_femur_neck__fracture": "B", "BMD_total__fracture": "B",
        "OA_knee__BMD_femur_neck": "C", "OA_knee__BMD_total": "C", "OA_knee__fracture": "C",
    }
    m = m.set_index("pair").loc[order].reset_index()

    methods = ["IVW", "MR-Egger", "Weighted\nMedian", "Leave-one-\nout", "MR-PRESSO\n(outlier)"]
    grid = np.zeros((len(m), 5), dtype=int)   # 0 blue,1 light,2 amber,3 verm
    txt = np.full((len(m), 5), "", dtype=object)

    for i, r in m.iterrows():
        ivw_dir = sgn(r.beta)
        # IVW
        if r.p < 0.05:
            grid[i, 0] = 0; txt[i, 0] = ("↑*" if ivw_dir > 0 else "↓*")
        else:
            grid[i, 0] = 1; txt[i, 0] = ("↑" if ivw_dir > 0 else "↓")
        # Egger
        ed = sgn(r.egger_beta)
        if ed == ivw_dir:
            grid[i, 1] = 0 if r.egger_p < 0.05 else 1
        else:
            grid[i, 1] = 3
        txt[i, 1] = (("↑*" if ed > 0 else "↓*") if r.egger_p < 0.05 else ("↑" if ed > 0 else "↓"))
        # Weighted median (direction only)
        wd = sgn(r.wmed)
        grid[i, 2] = 0 if wd == ivw_dir else 3
        txt[i, 2] = ("↑" if wd > 0 else "↓")
        # Leave-one-out stability
        stable = (sgn(r.loo_min) == sgn(r.loo_max)) and (r.loo_flips == 0)
        ld = sgn(r.loo_min)
        grid[i, 3] = 0 if (stable and ld == ivw_dir) else (2 if not stable else 1)
        txt[i, 3] = ("↑" if ld > 0 else "↓") + ("" if stable else "≈")
        # MR-PRESSO outlier test
        het = r.presso_Q_p < 0.05
        grid[i, 4] = 2 if het else 1
        txt[i, 4] = "het" if het else "clean"

    cmap = matplotlib.colors.ListedColormap(
        [OI["blue"], OI["lightblue"], OI["amber"], OI["vermillion"]])
    fig, ax = plt.subplots(figsize=(7.2, 6.0))
    ax.imshow(grid, cmap=cmap, vmin=0, vmax=3, aspect="auto")
    for i in range(len(m)):
        for j in range(5):
            ax.text(j, i, txt[i][j], ha="center", va="center", fontsize=7.5,
                    color="white" if grid[i, j] in (0, 3) else OI["black"],
                    fontweight="bold")
    # 组分隔
    ax.set_xticks(range(5))
    ax.set_xticklabels(methods, fontsize=8)
    short = {
        "BMD_arms__OA_knee": "BMD_arms → kneeOA",
        "BMD_femur_neck__OA_knee": "BMD_femur → kneeOA",
        "BMD_leg__OA_knee": "BMD_leg → kneeOA",
        "BMD_lumbar__OA_knee": "BMD_lumbar → kneeOA",
        "BMD_total__OA_knee": "BMD_total → kneeOA",
        "BMD_femur_neck__fracture": "BMD_femur → fracture",
        "BMD_total__fracture": "BMD_total → fracture",
        "OA_knee__BMD_femur_neck": "kneeOA → BMD_femur",
        "OA_knee__BMD_total": "kneeOA → BMD_total",
        "OA_knee__fracture": "kneeOA → fracture",
    }
    ax.set_yticks(range(len(m)))
    prefix = {p: g for g, ps in (("A", order[:5]), ("B", order[5:7]), ("C", order[7:])) for p in ps}
    ax.set_yticklabels(["%s  %s" % (prefix[p], short[p]) for p in m.pair], fontsize=7.5)
    # 组分隔横线
    for b in (4.5, 6.5):
        ax.axhline(b, color=OI["grey"], lw=0.8)
    # 图例
    leg = [Line2D([0], [0], marker="s", color="w", markerfacecolor=OI["blue"], markersize=9,
                  label="concordant & significant"),
           Line2D([0], [0], marker="s", color="w", markerfacecolor=OI["lightblue"], markersize=9,
                  label="concordant, not significant / direction-only"),
           Line2D([0], [0], marker="s", color="w", markerfacecolor=OI["amber"], markersize=9,
                  label="unstable / heterogeneity (PRESSO)"),
           Line2D([0], [0], marker="s", color="w", markerfacecolor=OI["vermillion"], markersize=9,
                  label="discordant direction")]
    ax.legend(handles=leg, loc="upper center", bbox_to_anchor=(0.5, -0.09), ncol=2,
              fontsize=7.5, frameon=False)
    ax.text(0.5, -0.20, "Groups: A = BMD → knee OA (5 sites)   B = BMD → fracture (2 sites)   "
            "C = reverse / OA → outcome (3 pairs).  Cell arrow = direction vs IVW; * = p<0.05.",
            transform=ax.transAxes, ha="center", fontsize=6.8, color=OI["grey"])
    ax.set_title("Evidence triangulation: 10 MR pairs × 5 estimators",
                 loc="left", fontsize=11, fontweight="bold")
    fig.subplots_adjust(left=0.22, right=0.97, top=0.93, bottom=0.24)
    for ext in ("pdf", "png"):
        ax.get_figure().savefig(os.path.join(FIG, "FigNEW_D_robustness.%s" % ext),
                                dpi=300, facecolor="white")
    print("[NEW-D] matrix built for %d pairs" % len(m))


# ============================================================ NEW-C : SOS 机制示意图
def fig_new_c():
    fig, ax = plt.subplots(figsize=(7.2, 6.0))
    ax.set_xlim(0, 10); ax.set_ylim(0, 11); ax.axis("off")

    def box(x, y, w, h, text, fc, ec=OI["black"], tc=OI["black"], fs=8.5, bold=False):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08",
                     fc=fc, ec=ec, lw=1.2))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
                color=tc, fontweight=("bold" if bold else "normal"))

    def arrow(x1, y1, x2, y2, color=OI["black"], ls="-", lw=1.5):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                     mutation_scale=12, color=color, lw=lw, linestyle=ls))

    # 上游设定点
    box(2.7, 8.9, 4.6, 1.3,
        "Polygenic osteogenic set-point\n(operationalised as areal BMD)",
        fc=OI["sky"], tc=OI["black"], fs=9, bold=True)
    arrow(4.0, 8.9, 2.1, 7.6, color=OI["blue"])
    arrow(6.0, 8.9, 7.9, 7.6, color=OI["orange"])

    # 左分支：全身骨 → 骨折↓
    box(0.4, 5.7, 3.3, 1.4, "Systemic bone strength ↑\n(positive control)",
        fc=OI["blue"], tc="white", fs=8.5)
    box(0.4, 3.4, 3.3, 1.4, "Fracture risk ↓\n(OR per SD lower BMD = 1.40;\nIVW beta = -0.337)",
        fc="#cfe2f3", tc=OI["black"], fs=8)
    arrow(2.05, 5.7, 2.05, 4.8, color=OI["blue"])
    arrow(2.05, 3.4, 2.05, 2.2, color=OI["blue"], ls="--")

    # 右分支：软骨下重塑 → 膝OA↑
    box(6.3, 5.7, 3.3, 1.4, "Subchondral remodelling ↑\n(formation-side)",
        fc=OI["orange"], tc="white", fs=8.5)
    box(6.3, 3.4, 3.3, 1.4, "Knee OA risk ↑\n(non-weight-bearing arms\nbeta = +0.075, p = 1.7×10^-6)",
        fc="#fde9c8", tc=OI["black"], fs=8)
    arrow(7.95, 5.7, 7.95, 4.8, color=OI["orange"])
    arrow(7.95, 3.4, 7.95, 2.2, color=OI["orange"], ls="--")

    # 旁注：Hardcastle + Cui（均指向膝OA/形成侧，避免交叉）
    box(0.2, 0.9, 4.4, 1.5,
        "Hardcastle triad (osteophytic /\nentheseal ossification, high BMD)\n→ OA OR 3.00 [1.96–4.58]",
        fc="#f2f2f2", tc=OI["black"], fs=7.3)
    box(5.4, 0.9, 4.4, 1.5,
        "Cui: bisphosphonate (antiresorptive)\n→ knee-replacement RR 0.71\nimplicates FORMATION side, not density",
        fc="#f2f2f2", tc=OI["black"], fs=7.3)
    arrow(4.6, 2.4, 6.3, 3.4, color=OI["grey"], ls=":", lw=1.0)
    arrow(7.95, 2.4, 7.95, 3.4, color=OI["grey"], ls=":", lw=1.0)

    ax.text(5.0, 10.7, "Systemic Osteogenic Set-point (SOS) hypothesis", fontsize=11,
            ha="center", fontweight="bold")
    ax.text(5.0, 7.2, "↑ systemic  /  ↓ local-mechanical", fontsize=7.5, ha="center",
            color=OI["grey"])
    ax.text(5.0, 0.35, "SOS = falsifiable hypothesis (F1–F10), not proven mechanism. "
            "Supported by F6 (τ²=0, cross-site homogeneity), F7 (locus exclusion), F8 (both sexes).",
            fontsize=6.8, ha="center", color=OI["grey"], va="bottom")
    fig.subplots_adjust(left=0.02, right=0.98, top=0.95, bottom=0.10)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, "FigNEW_C_sos.%s" % ext), dpi=300, facecolor="white")
    print("[NEW-C] SOS schematic built")


def emit_tif():
    """Add 300-dpi LZW TIFF versions for parity with the S9-S17 suppl set."""
    from PIL import Image
    base = {
        "FigNEW_B31_coloc_31snp": "FigNEW_B31_coloc_31snp",
        "FigNEW_Bcand_coloc_loci": "FigNEW_Bcand_coloc_loci",
        "FigNEW_D_robustness": "FigNEW_D_robustness",
        "FigNEW_C_sos": "FigNEW_C_sos",
    }
    for src in base:
        png = os.path.join(FIG, src + ".png")
        tif = os.path.join(FIG, src + ".tif")
        im = Image.open(png)
        im.save(tif, compression="tiff_lzw", resolution=(300, 300), dpi=(300, 300))
        print("[TIF] wrote", tif)


if __name__ == "__main__":
    export_fracture_screen()
    fig_new_b_31()
    fig_new_b_cand()
    fig_new_d()
    fig_new_c()
    emit_tif()
    print("ALL DONE ->", FIG)
