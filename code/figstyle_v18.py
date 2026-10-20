# -*- coding: utf-8 -*-
"""
figstyle_v18.py — 补充图（S9–S12）统一顶刊版式
与 build_figures_v17.py 完全一致的视觉规范，避免在补充图里另起一套配色。
规范：Arial/Helvetica ≥7 pt；面板标签 a/b/c 小写粗体；去上/右边框；
      Okabe–Ito 色盲安全配色；PNG 600 dpi + 矢量 PDF + TIFF 300 dpi（Adobe-deflate）。
"""
import os
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJECT = "C:/Users/Administrator/WorkBuddy/2026-09-08-22-17-40/OA_B_PROJECT"
GEN = os.path.join(PROJECT, "code", "output", "genetics")
HARM = os.path.join(GEN, "harmonised")
OUT = os.path.join(PROJECT, "Figures_v18_Suppl")

try:
    plt.style.use(["science", "nature", "no-latex"])
except Exception:
    pass

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.labelsize": 8,
    "axes.titlesize": 9,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "legend.fontsize": 7,
    "legend.frameon": False,
    "figure.dpi": 300,
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.04,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

# Okabe–Ito 色盲安全配色（与 v17 完全一致）
C = {
    "blue": "#0072B2", "orange": "#E69F00", "vermillion": "#D55E00",
    "green": "#009E73", "sky": "#56B4E9", "purple": "#CC79A7",
    "yellow": "#F0E442", "grey": "#8C8C8C", "ink": "#000000",
    "faint": "#DCDCDC",
}

SAVED = []


def save(fig, name, tiff_dpi=300):
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=600)
    fig.savefig(os.path.join(OUT, name + ".pdf"))
    try:
        import io
        from PIL import Image
        buf = io.BytesIO()
        fig.savefig(buf, format="tiff", dpi=tiff_dpi)
        buf.seek(0)
        Image.open(buf).convert("RGB").save(
            os.path.join(OUT, name + ".tif"), compression="tiff_adobe_deflate")
    except Exception as e:
        print("   [tif skipped]", name, e)
    plt.close(fig)
    SAVED.append(name)
    print("   OK", name)


def panel_label(ax, lab, x=-0.14, y=1.06):
    ax.text(x, y, lab, transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="top", ha="right")


def dump_plotted(payload, fname="plotted_data_v18.json"):
    """把每张新图的每一个绘制数值落盘，供 104 系列同款溯源核验。"""
    os.makedirs(OUT, exist_ok=True)
    p = os.path.join(OUT, fname)
    if os.path.exists(p):
        with open(p, encoding="utf-8") as f:
            old = json.load(f)
        old.update(payload)
        payload = old
    with open(p, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=1, ensure_ascii=False)
    print("   OK", fname)
