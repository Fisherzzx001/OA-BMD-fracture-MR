#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate Supplementary Figure S23: SOS hypothesis -> predictions (P1-P6)
-> falsifying tests (F1-F10) framework, colour-coded by outcome.

Outputs (reproducible):
  Figures_v18_Suppl/SupplFig23_SOS_framework.{png,pdf,svg}
  Figures/SupplFig23_SOS_framework.{png,pdf}
  submission_package/Figures/SupplFig23_SOS_framework.png

All numbers are taken verbatim from docs/08_投稿稿_OA_BMD_fracture.md
(§4, Table S11 / S26-S32).
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.lines import Line2D

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_DIRS = [
    os.path.join(ROOT, "Figures_v18_Suppl"),
    os.path.join(ROOT, "Figures"),
    os.path.join(ROOT, "submission_package", "Figures"),
]

# ---- outcome palette -------------------------------------------------------
GREEN  = "#2e7d32"; GREEN_F  = "#e7f1e7"
AMBER  = "#c77700"; AMBER_F  = "#fbf0d9"
GREY   = "#5f6b73"; GREY_F   = "#eceff1"
RED    = "#c62828"; RED_F    = "#f8e3e3"
INK    = "#1c1c1c"

# ---- hypothesis ------------------------------------------------------------
HYP_TITLE = "Systemic Osteogenic Set-point (SOS) Hypothesis"
HYP_BODY = [
    "One polygenic osteogenic set-point, operationalised as areal BMD,",
    "protects the whole bone (lower fracture risk) yet raises knee-OA",
    "risk at the joint via formation-side subchondral remodelling.",
]

# ---- predictions P1-P6 ----------------------------------------------------
# (id, title, [body lines], colour, face)
PRED = [
    ("P1", "Opposing endpoints", ["Fracture-protect +", "OA-raise (5/5 +)"], GREEN, GREEN_F),
    ("P2", "Systemic, not\nlocal-mechanical",
        ["Non-wt-bear arms", "b=+0.075*;", "ordering ns"], AMBER, AMBER_F),
    ("P3", "Diffuse polygenic", ["0/31 colocalise;", "direction-consistent,", "non-discriminating"], GREY, GREY_F),
    ("P4", "Size a separable\npartial component",
        ["67-91% non-size;", "volume arm", "untested"], AMBER, AMBER_F),
    ("P5", "Constitutional\nasymmetry", ["Reverse null;", "underpowered", "(0.14-0.26)"], AMBER, AMBER_F),
    ("P6", "Site-pathway\ndominance", ["Restored by F1;", "cross-site t2=0"], GREEN, GREEN_F),
]

# ---- falsifying tests F1-F10 ----------------------------------------------
# (id, name, tests_Px, [body], colour, face, target_pred_index or None)
# target_pred_index -> which P (0..5) the near-vertical connector attaches to
TOP = [
    ("F5", "Osteophyte-stratified outcome", "tests P1/P3",
        ["Site-strat MR;", "formation-side"], GREEN, GREEN_F, 0),
    ("F8", "Sex stratification", "tests P2",
        ["Positive in", "both sexes"], GREEN, GREEN_F, 1),
    ("F3", "Bone vs cartilage enrichment", "tests P3",
        ["Cartilage > bone", "(gene-set)"], GREEN, GREEN_F, 2),
    ("F2", "Mass-geometry decomposition", "tests P4",
        ["Knee-OA cond. +;", "volume arm", "untested"], AMBER, AMBER_F, 3),
    ("F1", "Hip-OA replication", "tests P6",
        ["b=+0.097;", "OR 1.102", "(1.068-1.136)"], GREEN, GREEN_F, 5),
]
BOT = [
    ("F9", "Site-gradient ordering", "tests P2/P4",
        ["arms-fem.neck", "+0.017; p=0.228", "non-discriminating"], GREY, GREY_F, 1),
    ("F4", "eQTL colocalisation", "tests P3",
        ["Fracture 9/31 vs", "OA 0/31", "(p=0.0043)"], GREEN, GREEN_F, 2),
    ("F7", "Joint-locus exclusion", "tests P3",
        ["GDF5 not a", "BMD instrument"], GREEN, GREEN_F, None),
    ("F6", "Cross-site homogeneity", "tests P6",
        ["Q=6.92, p=0.14;", "t2=0"], GREEN, GREEN_F, 5),
    ("F10", "Ancestry sensitivity", "external validity",
        ["Registered test", "unexecuted;", "proxy: S32"], RED, RED_F, None),
]

OUTCOME_LEGEND = [
    (GREEN, "Supported / Pass"),
    (AMBER, "Partial (arm untested)"),
    (GREY, "Null / non-discriminating"),
    (RED, "Not executed"),
]


def add_box(ax, x0, y0, w, h, title, body_lines, edge, face,
            title_fs=8.5, body_fs=6.6, title_weight="bold", tstep=0.020):
    ax.add_patch(FancyBboxPatch((x0, y0), w, h,
                 boxstyle="round,pad=0.004,rounding_size=0.010",
                 linewidth=1.2, edgecolor=edge, facecolor=face, zorder=3))
    tlines = title.split("\n")
    ty = y0 + h - 0.012
    for tl in tlines:
        ax.text(x0 + w / 2, ty, tl, ha="center", va="top",
                fontsize=title_fs, fontweight=title_weight, color=edge, zorder=4)
        ty -= tstep
    start_y = ty - 0.006
    step = 0.0205
    for i, line in enumerate(body_lines):
        ax.text(x0 + w / 2, start_y - i * step, line, ha="center", va="top",
                fontsize=body_fs, color=INK, zorder=4)


def connector(ax, x1, y1, x2, y2, color="#9aa3ab"):
    ax.add_line(Line2D([x1, x2], [y1, y2], color=color, lw=0.9,
                       zorder=1, alpha=0.85))
    # small arrowhead near the lower end
    ax.annotate("", xy=(x2, y2), xytext=(x2, y2 + 0.018),
                arrowprops=dict(arrowstyle="-|>", color=color, lw=0.9),
                zorder=1)


def main():
    fig, ax = plt.subplots(figsize=(15.0, 11.2))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # heading
    ax.text(0.5, 0.985, "Supplementary Figure S23.  SOS hypothesis \u2192 predictions (P1\u2013P6) "
            "\u2192 falsifying tests (F1\u2013F10)", ha="center", va="top",
            fontsize=12, fontweight="bold", color=INK)

    # ---- Hypothesis box ----
    hx0, hy0, hw, hh = 0.06, 0.845, 0.88, 0.115
    add_box(ax, hx0, hy0, hw, hh, HYP_TITLE, HYP_BODY, INK, "#dfe7ef",
            title_fs=11, body_fs=8.2)
    hcx = hx0 + hw / 2
    hbottom = hy0 + hh

    # ---- Predictions ----
    pw, ph = 0.135, 0.155
    p_cx = [0.1075, 0.2645, 0.4215, 0.5785, 0.7355, 0.8925]
    py0 = 0.625
    for cx, (pid, ptitle, pbody, edge, face) in zip(p_cx, PRED):
        x0 = cx - pw / 2
        title = f"{pid}.  {ptitle}"
        add_box(ax, x0, py0, pw, ph, title, pbody, edge, face)
        # connector H -> P
        connector(ax, cx, py0 + ph, cx, hbottom)

    # ---- Tests top row ----
    tw, th = 0.165, 0.145
    t_top_cx = [0.1225, 0.31125, 0.5, 0.68875, 0.8775]
    ttop_y0 = 0.405
    for cx, (fid, fname, tag, fbody, edge, face, tgt) in zip(t_top_cx, TOP):
        x0 = cx - tw / 2
        add_box(ax, x0, ttop_y0, tw, th, fid + "\n" + fname, [tag] + fbody,
                edge, face, title_fs=7.4, body_fs=6.3, tstep=0.017)
        if tgt is not None:
            connector(ax, p_cx[tgt], py0, cx, ttop_y0 + th)

    # ---- Tests bottom row ----
    tbot_y0 = 0.185
    for cx, (fid, fname, tag, fbody, edge, face, tgt) in zip(t_top_cx, BOT):
        x0 = cx - tw / 2
        add_box(ax, x0, tbot_y0, tw, th, fid + "\n" + fname, [tag] + fbody,
                edge, face, title_fs=7.4, body_fs=6.3, tstep=0.017)
        if tgt is not None:
            connector(ax, p_cx[tgt], py0, cx, tbot_y0 + th)

    # ---- Legend ----
    lx0 = 0.04
    ly0 = 0.075
    lh = 0.045
    step = 0.215
    for i, (col, lab) in enumerate(OUTCOME_LEGEND):
        x = lx0 + i * step
        ax.add_patch(FancyBboxPatch((x, ly0), 0.035, lh,
                     boxstyle="round,pad=0.002,rounding_size=0.006",
                     linewidth=1.0, edgecolor=col, facecolor=col, zorder=3))
        ax.text(x + 0.045, ly0 + lh / 2, lab, ha="left", va="center",
                fontsize=7.6, color=INK, zorder=4)

    # ---- Verdict ----
    # 2026-11-08 (R1 theory rectification): state the discriminating kernel explicitly
    # and replace "SOS survives 9/10" with the weaker, non-completeness claim.
    verdict_lines = [
        "Verdict: 9/10 executed tests do not refute SOS (F2 partial \u2014 volume arm untested; "
        "F9 non-discriminating; F10 registered test unexecuted, ancestry-composition proxy in S32).",
        "Kernel: one inherited variant set must fix opposite signs at whole-bone and joint endpoints; "
        "distinct variant sets or concordant signs refute it.",
        "Priority closing tests: F2 volumetric arm (true vBMD instruments) and F10 cross-ancestry "
        "replication; P4 is discriminating and F2 decisive, with the volume arm unexecuted.",
        "Supported, not demonstrated. * non-weight-bearing arms \u03b2 = +0.075, p = 1.7\u00d710\u207b\u2076 "
        "(systemic, not local-mechanical clause).",
    ]
    vy = 0.062
    for ln in verdict_lines:
        ax.text(0.04, vy, ln, ha="left", va="top", fontsize=7.4,
                color="#333333", zorder=4)
        vy -= 0.015

    # ---- save ----
    for d in OUT_DIRS:
        os.makedirs(d, exist_ok=True)
        base = os.path.join(d, "SupplFig23_SOS_framework")
        fig.savefig(base + ".png", dpi=300, bbox_inches="tight",
                   pad_inches=0.04)
        fig.savefig(base + ".pdf", dpi=300, bbox_inches="tight",
                   pad_inches=0.04)
        if d == OUT_DIRS[0]:
            try:
                fig.savefig(base + ".svg", bbox_inches="tight", pad_inches=0.04)
            except Exception as e:
                print("svg skip:", e)
    print("WROTE SupplFig23_SOS_framework to:", OUT_DIRS)


if __name__ == "__main__":
    main()
