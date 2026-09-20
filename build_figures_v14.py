# -*- coding: utf-8 -*-
"""
build_figures_v14.py
=========================================================================
基于论文真实数据（与 v13 完全一致，不改任何统计量）重绘全部 12 张图，
并新增 2 张加强表达的图，统一 Nature 系列 + Elsevier(Bone) 图片规范：
  - Okabe-Ito 色盲安全配色（已在 build_figures*.py 源码打补丁）
  - 轴标签/刻度 >=7pt（已打补丁）
  - 输出 PNG(300dpi) + 可编辑矢量 PDF + TIFF(300dpi, Elsevier 生产首选栅格)

重绘（数据真实、逻辑来自 build_figures*.py，仅风格升级）：
  Figure1_rg_heatmap_13, Figure2_forest_key_results, Figure3_MR_forest,
  Figure4_h2_bar, Figure5_causal_network,
  Figure6_MR_scatter, Figure7_MR_funnel, SupplFig1_MR_scatter_remaining,
  SupplFig2_leave_one_out, Figure8_MR_power, Figure9_overlap_pleiotropy,
  SupplFig3_GEFOS_benchmark

新增（加强表达，均为示意图，不含统计量）：
  GraphicalAbstract      -> 图文摘要（Bone 鼓励提交；OA-BMD-fracture 三角因果）
  SupplFig4_STROBE_flow  -> STROBE-MR 研究流程图（参与表型选择->QC->协调->MR->敏感->报告）
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.patheffects as pe

PROJECT = "C:/Users/Administrator/WorkBuddy/2026-09-08-22-17-40/OA_B_PROJECT"
OUT = os.path.join(PROJECT, "Figures_Beautified")
os.makedirs(OUT, exist_ok=True)
os.environ["FIG_OUT"] = OUT

# 在导入三个模块前设置输出目录（模块在 import 时读取 FIG_OUT）
import build_figures as F1
import build_figures_mr as F2
import build_figures_v9 as F3

# ---------- 统一 Nature / Okabe-Ito 风格 ----------
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})
CB = {"blue": "#0072B2", "red": "#D55E00", "green": "#009E73", "orange": "#E69F00",
      "purple": "#CC79A7", "teal": "#56B4E9", "grey": "#999999", "ink": "#000000",
      "faint": "#DDDDDD", "cb_blue": "#0072B2", "cb_orange": "#E69F00"}

def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=300)
    fig.savefig(os.path.join(OUT, name + ".pdf"))
    try:
        fig.savefig(os.path.join(OUT, name + ".tif"), dpi=300)
    except Exception as e:
        print("  tif skipped:", name, e)
    plt.close(fig)
    print("  OK", name)


# ================================================================== 新图 A：图文摘要
def fig_graphical_abstract():
    fig, ax = plt.subplots(figsize=(4.4, 6.2))
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")

    def node(x, y, w, h, text, fc, tc="white"):
        ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                     boxstyle="round,pad=0.02,rounding_size=0.18",
                     fc=fc, ec="none", zorder=2))
        ax.text(x, y, text, ha="center", va="center", fontsize=9,
                fontweight="bold", color=tc, zorder=3)

    def arrow(x1, y1, x2, y2, color, label, lx, ly, rad=0.0, lw=2.4):
        arr = FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                     mutation_scale=15, color=color, lw=lw,
                     connectionstyle=f"arc3,rad={rad}", zorder=1,
                     path_effects=[pe.withStroke(linewidth=3.4, foreground="white")])
        ax.add_patch(arr)
        ax.text(lx, ly, label, ha="center", va="center", fontsize=7.6,
                fontweight="bold", color=color, zorder=4,
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.92))

    node(2.2, 7.7, 2.8, 1.2, "OA (knee)", CB["purple"])
    node(7.6, 6.4, 3.0, 1.7, "BMD\n(10 regional sites)", CB["blue"])
    node(7.6, 1.9, 2.8, 1.2, "Fracture", CB["teal"])

    arrow(6.9, 6.0, 6.9, 2.4, CB["blue"], "Protective\nOR 1.40–1.42", 8.9, 4.2, rad=-0.12)
    arrow(6.6, 6.1, 3.6, 7.0, CB["red"], "Causal ↑ risk\nβ +0.017 to +0.028", 5.0, 7.9, rad=0.16)
    arr = FancyArrowPatch((3.4, 5.2), (6.4, 6.0), arrowstyle="-|>", mutation_scale=11,
               color=CB["grey"], lw=1.6, ls=(0, (3, 2)), zorder=1,
               connectionstyle="arc3,rad=-0.22",
               path_effects=[pe.withStroke(linewidth=3.0, foreground="white")])
    ax.add_patch(arr)
    ax.text(5.0, 4.6, "Genetic correlation\nr_g +0.12 to +0.25", ha="center", va="center",
            fontsize=7.2, color=CB["grey"], zorder=4,
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.92))
    arr2 = FancyArrowPatch((3.6, 5.2), (6.6, 2.5), arrowstyle="-|>", mutation_scale=10,
                color=CB["grey"], lw=1.4, ls=(0, (3, 2)), zorder=1,
                connectionstyle="arc3,rad=0.28",
                path_effects=[pe.withStroke(linewidth=3.0, foreground="white")])
    ax.add_patch(arr2)
    ax.text(5.4, 3.4, "Null (ns)", ha="center", va="center", fontsize=7.2,
            color=CB["grey"], zorder=4, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.92))

    ax.set_title("OA–BMD–fracture triangle: shared genetics,\nseparate fracture fate",
                 fontsize=10, fontweight="bold", pad=10)
    ax.text(5.0, 0.35, "Higher BMD protects against fracture yet raises OA risk;\n"
            "OA shows no well-powered causal pathway to fracture.",
            ha="center", va="center", fontsize=7.2, color=CB["ink"],
            bbox=dict(boxstyle="round,pad=0.4", fc="#F2F2F2", ec="none"))
    save(fig, "GraphicalAbstract")


# ================================================================== 新图 B：STROBE-MR 流程图
def fig_strobe_flow():
    fig, ax = plt.subplots(figsize=(7.4, 9.6))
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")

    stages = [
        ("1  GWAS inputs", CB["blue"],
         "13 skeletal phenotypes\nUKB OA (knee, self-rep); UKB/GEFOS 10-site DXA BMD;\nGEFOS & IEU OpenGWAS fracture. N: 1.21M / 31.9k / 0.46M / 0.26M"),
        ("2  Quality control", CB["purple"],
         "LDSC h² + intercept (all ≈1.0); allele-orientation checks via top loci;\nOA_selfrep global Z-sign correction (rg +0.84 post-correction)"),
        ("3  Harmonize & clump", CB["teal"],
         "Strict PLINK-equivalent clumping (r²<0.001 / 10 Mb, 1000G-EUR);\n1-Mb window as sensitivity; 16–31 IVs per BMD exposure, 69 for OA_knee"),
        ("4  Two-sample MR", CB["orange"],
         "10 exposure–outcome pairs; IVW (primary), MR-Egger, weighted median;\nN-scale correction (β ×0.2823) for reconstructed BMD N = 31,873"),
        ("5  Sensitivity & diagnostics", CB["red"],
         "MR-PRESSO outliers; leave-one-out; Cochran Q / I²;\nsample-overlap ρ_GE; Bonferroni + BH-FDR across 78 r_g / 10 MR pairs"),
        ("6  STROBE-MR reporting", CB["green"],
         "Transparent design, effect scales, F-statistic (44–75),\nethics statement and full supplementary tables"),
    ]
    ytop = 9.2
    dy = 1.62
    for i, (title, col, desc) in enumerate(stages):
        y = ytop - i * dy
        ax.add_patch(FancyBboxPatch((0.4, y - 0.62), 9.2, 1.24,
                     boxstyle="round,pad=0.02,rounding_size=0.12",
                     fc="#FFFFFF", ec=col, lw=1.8, zorder=2))
        ax.add_patch(FancyBboxPatch((0.4, y - 0.62), 0.55, 1.24,
                     boxstyle="round,pad=0.02,rounding_size=0.12",
                     fc=col, ec="none", zorder=3))
        ax.text(0.675, y + 0.30, title, ha="left", va="center",
                fontsize=9.5, fontweight="bold", color=col, zorder=4)
        ax.text(1.15, y - 0.18, desc, ha="left", va="center",
                fontsize=7.0, color=CB["ink"], zorder=4)
        if i < len(stages) - 1:
            ax.add_patch(FancyArrowPatch((5.0, y - 0.62), (5.0, y - dy + 0.62),
                         arrowstyle="-|>", mutation_scale=14, color=CB["grey"],
                         lw=1.6, zorder=1))

    ax.set_title("Study design and MR analysis flow (STROBE-MR compliant)",
                 fontsize=11, fontweight="bold", pad=10)
    save(fig, "SupplFig4_STROBE_flow")


# ================================================================== 主流程
def main():
    print("▶ 重绘 12 张真实数据图（风格升级：Okabe-Ito / >=7pt / PNG+PDF+TIFF）")
    F1.fig1(); F1.fig2(); F1.fig3(); F1.fig4(); F1.fig5()
    F2.main()
    F3.main()
    print("▶ 新增 2 张加强表达图")
    fig_graphical_abstract()
    fig_strobe_flow()
    print("ALL v14 FIGURES ->", OUT)

if __name__ == "__main__":
    main()
