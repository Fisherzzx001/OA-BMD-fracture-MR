# -*- coding: utf-8 -*-
"""Fig NEW-A: Causal (MR IVW beta) x Correlational (LDSC r_g) triangulation.
All numbers from real CSVs: AUDIT_mr_recomputed.csv + rg_results_long_13.csv. No invented values.
"""
import os, json
import pandas as pd, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GEN = os.path.join(ROOT, "code", "output", "genetics")
OUT = os.path.join(ROOT, "Figures_Proposed")
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "font.family": "Arial", "font.size": 8, "axes.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "xtick.direction": "out", "ytick.direction": "out",
    "figure.dpi": 300, "savefig.dpi": 600, "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05, "pdf.fonttype": 42,
})
CB_BLUE="#0072B2"; CB_ORANGE="#E69F00"; GREY="#8C8C8C"; RED="#D65F5F"

mr = pd.read_csv(os.path.join(GEN, "AUDIT_mr_recomputed.csv"))
rg = pd.read_csv(os.path.join(GEN, "rg_results_long_13.csv"))
def rg_of(a, b):
    r = rg[((rg.trait1==a)&(rg.trait2==b))|((rg.trait1==b)&(rg.trait2==a))]
    return None if r.empty else float(r.iloc[0].rg)

rows = []
for _, r in mr.iterrows():
    e, o = r["exposure"], r["outcome"]
    rgv = rg_of(e, o)
    if rgv is None: continue
    if str(e).startswith("OA"):
        grp = "OA->fracture" if o=="fracture" else "reverse OA->BMD"
    elif o=="fracture":
        grp = "fracture control"
    else:
        grp = "BMD->OA knee"
    rows.append((f"{e}->{o}", float(r["beta"]), float(r["se"]), rgv, grp))
df = pd.DataFrame(rows, columns=["pair","beta","se","rg","grp"])
print(df.to_string(index=False))

colmap = {"fracture control":CB_BLUE, "BMD->OA knee":CB_ORANGE,
          "reverse OA->BMD":GREY, "OA->fracture":RED}

fig, ax = plt.subplots(figsize=(7.2,5.0))
# quadrant shading for "MR and r_g same sign"
ax.add_patch(Rectangle((0,0),1,1, color="#F2F2F2", zorder=0))      # top-right (both +)
ax.add_patch(Rectangle((-1,-1),1,1, color="#F2F2F2", zorder=0))   # bottom-left (both -)
ax.axhline(0, color="k", lw=0.6); ax.axvline(0, color="k", lw=0.6)
# leader-line label positions (data coords) to avoid cluster overlap
LBL = {
 "BMD_lumbar->OA_knee": (0.34, 0.108, "left"),
 "BMD_total->OA_knee": (0.34, 0.085, "left"),
 "BMD_arms->OA_knee": (0.34, 0.062, "left"),
 "BMD_leg->OA_knee": (0.34, 0.039, "left"),
 "BMD_femur_neck->OA_knee": (0.34, 0.016, "left"),
 "OA_knee->BMD_total": (-0.16, 0.062, "right"),
 "OA_knee->BMD_femur_neck": (-0.16, 0.040, "right"),
 "OA_knee->fracture": (-0.16, 0.018, "right"),
 "BMD_total->fracture": (-0.28, -0.30, "left"),
 "BMD_femur_neck->fracture": (-0.28, -0.42, "left"),
}
for _, r in df.iterrows():
    c = colmap[r["grp"]]
    ax.errorbar(r["rg"], r["beta"], xerr=0, yerr=1.96*r["se"],
                fmt="o", color=c, ecolor=c, capsize=2, ms=6, elinewidth=1)
    lx, ly, ha = LBL[r["pair"]]
    ax.annotate(r["pair"].replace("->", "→"), (r["rg"], r["beta"]),
                xytext=(lx, ly), textcoords="data", ha=ha, va="center",
                fontsize=6.3, arrowprops=dict(arrowstyle="-", color="#999999", lw=0.5))
ax.set_xlabel("LDSC genetic correlation r_g  (exposure ↔ outcome)")
ax.set_ylabel("MR IVW β  (N-corrected per-SD scale)")
ax.set_xlim(-0.72, 0.58); ax.set_ylim(-0.45, 0.14)
ax.text(0.5, 0.965, "Shaded quadrants = MR and r_g same sign (causal parallel correlational)",
        transform=ax.transAxes, fontsize=6.8, color="#333333", va="top", ha="center")
from matplotlib.lines import Line2D
leg = [Line2D([0],[0], marker="o", color="w", markerfacecolor=v, markersize=7, label=k)
       for k,v in colmap.items()]
ax.legend(handles=leg, loc="upper left", fontsize=6.8, frameon=False)
ax.text(-0.15, 1.05, "a", transform=ax.transAxes, fontsize=8, fontweight="bold", va="top", ha="right")
fig.suptitle("Causal–correlational triangulation of the OA–BMD–fracture triangle", fontsize=9, y=0.995)
fig.savefig(os.path.join(OUT, "FigNEW_A_triangulation.png"))
fig.savefig(os.path.join(OUT, "FigNEW_A_triangulation.pdf"))
print("SAVED:", os.path.join(OUT, "FigNEW_A_triangulation.png"))
