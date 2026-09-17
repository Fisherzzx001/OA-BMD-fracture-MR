# -*- coding: utf-8 -*-
"""
Build improved submission figures following the nature-figure skill style spec:
  - Arial font, >=7pt, colorblind-safe palette, no top/right spines,
    lowercase-bold panel labels, 300 DPI, PNG + vector PDF.
Outputs (into submission_package/Figures/):
  Figure1_rg_heatmap_13.png/.pdf   (regenerated, with BH-FDR asterisks + block lines)
  Figure2_forest_key_results.png/.pdf  (regenerated, 4 panels, CIs)
  Figure3_MR_forest.png/.pdf       (regenerated, 3 panels, colour-coded)
  Figure5_causal_network.png/.pdf  (NEW conceptual diagram)
(Figure4_h2_bar is produced by the science-figures skill CLI in a separate step.)
"""
import os, csv, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe

GEN = "C:/Users/Administrator/WorkBuddy/2026-09-08-22-17-40/OA_B_PROJECT/code/output/genetics"
OUT = "C:/Users/Administrator/WorkBuddy/2026-09-08-22-17-40/OA_B_PROJECT/submission_package/Figures"

# ---------- Nature style ----------
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
CB = {  # colorblind-safe Nature palette
    "blue":   "#4878CF",
    "red":    "#D65F5F",
    "green":  "#6ACC65",
    "orange": "#EE854A",
    "purple": "#956CB4",
    "teal":   "#82C6E2",
    "grey":   "#8C8C8C",
    "cb_blue": "#0072B2",
    "cb_orange": "#E69F00",
    "ink": "#222222",
}

ORDER = ["OA_knee","OA_selfrep","fracture","BMD_total","BMD_femur_neck","BMD_femur_total",
         "BMD_lumbar","BMD_spine","BMD_arms","BMD_leg","BMD_pelvis","BMD_rib","BMD_trunk"]
PRETTY = {
    "OA_knee":"OA (knee)","OA_selfrep":"OA (self-rep)","fracture":"Fracture",
    "BMD_total":"BMD total","BMD_femur_neck":"BMD femoral neck","BMD_femur_total":"BMD femoral total",
    "BMD_lumbar":"BMD lumbar","BMD_spine":"BMD spine","BMD_arms":"BMD arms",
    "BMD_leg":"BMD legs","BMD_pelvis":"BMD pelvis","BMD_rib":"BMD rib","BMD_trunk":"BMD trunk",
}

def read_matrix(fname):
    path = os.path.join(GEN, fname)
    with open(path) as f:
        r = csv.reader(f)
        header = next(r)
        rows = [row for row in r if row]
    labels = [h for h in header if h]
    M = np.array([[float(x) if x != "" else np.nan for x in row[1:]] for row in rows])
    return labels, M

def panel_label(ax, txt, x=-0.02, y=1.05):
    """面板标签置于轴左上角、与标题同高，避开 y 轴刻度文字。"""
    ax.text(x, y, txt, transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="bottom", ha="right", color=CB["ink"])


def limits(los, his, include_zero=True, pad=0.06):
    """由数据本身推出 x 轴范围，保证每条误差棒都完整可见（避免坐标轴截断）。

    之前的做法是把范围写死（如 -0.75~0.35 / -0.30~0.15），实测会切掉多条
    95% CI 的端线：Fig2a 有 6/10 条下缘被截（最远 -0.921），Fig3c 有 3/3
    条上缘被截（最远 +0.354）。改为数据驱动后不再丢失任何信息。
    """
    lo, hi = float(min(los)), float(max(his))
    if include_zero:
        lo, hi = min(lo, 0.0), max(hi, 0.0)
    span = hi - lo
    if span <= 0:
        span = max(abs(hi), abs(lo), 1e-6)
    return lo - pad * span, hi + pad * span

def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=300)
    fig.savefig(os.path.join(OUT, name + ".pdf"))
    print("saved", name)

# ---------- BH-FDR q across the 78 unique off-diagonal pairs ----------
def bh_q(pairs_p):
    p = np.asarray(pairs_p, float)
    m = len(p)
    order = np.argsort(p)
    ps = p[order]
    q = ps * m / (np.arange(m) + 1)
    # enforce monotonicity from the end
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(m)
    out[order] = q
    return out

# ================= FIGURE 1 : r_g heatmap =================
def fig1():
    labels, M = read_matrix("rg_matrix_13.csv")
    _, P = read_matrix("rg_p_matrix_13.csv")
    n = len(labels)
    # unique off-diagonal p for BH
    pu, idx = [], []
    for i in range(n):
        for j in range(i+1, n):
            pu.append(P[i, j]); idx.append((i, j))
    q = bh_q(pu)
    Q = np.ones((n, n))
    for k, (i, j) in enumerate(idx):
        Q[i, j] = Q[j, i] = q[k]

    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    norm = TwoSlopeNorm(vmin=-0.65, vcenter=0.0, vmax=1.0)
    im = ax.imshow(M, cmap="RdBu_r", norm=norm, aspect="auto")
    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    ax.set_xticklabels([PRETTY[l] for l in labels], rotation=45, ha="right", fontsize=7)
    ax.set_yticklabels([PRETTY[l] for l in labels], fontsize=7)
    # cell text
    for i in range(n):
        for j in range(n):
            v = M[i, j]
            if i == j:
                ax.text(j, i, "1.00", ha="center", va="center", fontsize=6.2,
                        color="white", fontweight="bold")
                continue
            col = "white" if abs(v) > 0.55 or v < -0.35 else CB["ink"]
            star = ""
            if Q[i, j] < 0.001: star = "***"
            elif Q[i, j] < 0.01: star = "**"
            elif Q[i, j] < 0.05: star = "*"
            txt = f"{v:+.2f}{star}"
            ax.text(j, i, txt, ha="center", va="center", fontsize=6.0, color=col)
    # block separators: OA(0,1) | fracture(2) | BMD(3..12)
    for b in (2, 3):
        ax.axvline(b - 0.5, color="black", lw=1.1)
        ax.axhline(b - 0.5, color="black", lw=1.1)
    ax.set_xticks(np.arange(-0.5, n, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n, 1), minor=True)
    ax.tick_params(which="minor", length=0)
    # colorbar
    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label("Genetic correlation (r_g)", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    ax.set_title("Genetic correlation matrix (LDSC, 13 skeletal phenotypes)",
                 fontsize=9, fontweight="bold", pad=10)
    fig.text(0.5, 0.005, "Asterisks: BH-FDR q < 0.05 (*), < 0.01 (**), < 0.001 (***).  "
             "Black lines separate the three phenotype blocks (OA / fracture / BMD).",
             ha="center", fontsize=6.5, color=CB["grey"])
    save(fig, "Figure1_rg_heatmap_13")
    plt.close(fig)

# ================= FIGURE 2 : r_g forest =================
def fig2(LIMITS=None):
    labels, M = read_matrix("rg_matrix_13.csv")
    _, SE = read_matrix("rg_se_matrix_13.csv")
    _, P = read_matrix("rg_p_matrix_13.csv")
    n = len(labels)
    pu, idx = [], []
    for i in range(n):
        for j in range(i+1, n):
            pu.append(P[i, j]); idx.append((i, j))
    Q = np.ones((n, n)); q = bh_q(pu)
    for k, (i, j) in enumerate(idx):
        Q[i, j] = Q[j, i] = q[k]

    bmd = [l for l in labels if l.startswith("BMD_")]
    panels = [
        ("a", "Fracture vs 10 BMD sites (all negative)", "fracture", bmd, CB["cb_blue"]),
        ("b", "Knee OA vs 10 BMD sites (all positive)", "OA_knee", bmd, CB["cb_blue"]),
        ("c", "Self-reported OA vs 10 BMD sites (all positive)", "OA_selfrep", bmd, CB["cb_blue"]),
    ]
    fig = plt.figure(figsize=(7.2, 7.4))
    gs = fig.add_gridspec(4, 1, height_ratios=[1, 1, 1, 0.9], hspace=0.55)
    axes = [fig.add_subplot(gs[k]) for k in range(3)]
    axd = fig.add_subplot(gs[3])

    for ax, (lab, title, row, cols, color) in zip(axes, panels):
        ys = list(range(len(cols)))[::-1]
        los, his = [], []
        for y, c in zip(ys, cols):
            i, j = labels.index(row), labels.index(c)
            rg, se = M[i, j], SE[i, j]
            sig = Q[i, j] < 0.05
            col = color if sig else CB["grey"]
            ax.errorbar(rg, y, xerr=1.96*se, fmt="o", color=col, ecolor=col,
                        markersize=4, elinewidth=1.1, capsize=2.5, zorder=3)
            los.append(rg - 1.96*se); his.append(rg + 1.96*se)
        ax.axvline(0, color=CB["grey"], ls="--", lw=0.8, zorder=1)
        ax.set_yticks(ys); ax.set_yticklabels([PRETTY[c] for c in cols], fontsize=7)
        ax.set_xlabel("Genetic correlation (r_g)", fontsize=7.5)
        ax.set_xlim(*limits(los, his))
        if LIMITS is not None: LIMITS[f"fig2_{lab}"] = list(ax.get_xlim())
        panel_label(ax, lab)
        ax.set_title(title, fontsize=8, fontweight="bold", loc="left")

    # panel d: OA-OA + OA-fracture
    drows = [("OA_knee","OA_selfrep","OA (knee) vs OA (self-rep)"),
             ("OA_knee","fracture","OA (knee) vs Fracture"),
             ("OA_selfrep","fracture","OA (self-rep) vs Fracture")]
    ys = list(range(len(drows)))[::-1]
    los, his = [], []
    for y, (a, b, _) in zip(ys, drows):
        i, j = labels.index(a), labels.index(b)
        rg, se = M[i, j], SE[i, j]
        sig = Q[i, j] < 0.05
        col = CB["cb_blue"] if sig else CB["grey"]
        axd.errorbar(rg, y, xerr=1.96*se, fmt="o", color=col, ecolor=col,
                     markersize=4, elinewidth=1.1, capsize=2.5, zorder=3)
        los.append(rg - 1.96*se); his.append(rg + 1.96*se)
    axd.axvline(0, color=CB["grey"], ls="--", lw=0.8, zorder=1)
    axd.set_yticks(ys); axd.set_yticklabels([t for _,_,t in drows], fontsize=7)
    axd.set_xlabel("Genetic correlation (r_g)", fontsize=7.5)
    axd.set_xlim(*limits(los, his))
    if LIMITS is not None: LIMITS["fig2_d"] = list(axd.get_xlim())
    panel_label(axd, "d")
    axd.set_title("OA\u2013OA and OA\u2013fracture pairs", fontsize=8, fontweight="bold", loc="left")
    save(fig, "Figure2_forest_key_results")
    plt.close(fig)

# ================= FIGURE 3 : MR forest =================
def fig3(LIMITS=None):
    rows = []
    with open(os.path.join(GEN, "mr_corrected_N.csv")) as f:
        for row in csv.DictReader(f):
            rows.append(row)
    def get(exp, out):
        for r in rows:
            if r["exposure"] == exp and r["outcome"] == out:
                return r
    def num(x):
        try: return float(x)
        except: return np.nan

    panels = [
        ("a", "BMD \u2192 fracture (positive control)", [
            ("BMD total \u2192 fracture\nOR 1.40 (1.34\u20131.46)", "BMD_total","fracture", CB["cb_blue"]),
            ("BMD femoral neck \u2192 fracture\nOR 1.42 (1.34\u20131.51)", "BMD_femur_neck","fracture", CB["cb_blue"]),
        ]),
        ("b", "BMD \u2192 knee OA (causally risk-increasing)", [
            ("BMD total \u2192 OA knee", "BMD_total","OA_knee", CB["red"]),
            ("BMD femoral neck \u2192 OA knee", "BMD_femur_neck","OA_knee", CB["red"]),
            ("BMD lumbar \u2192 OA knee", "BMD_lumbar","OA_knee", CB["red"]),
            ("BMD arms \u2192 OA knee", "BMD_arms","OA_knee", CB["red"]),
            ("BMD legs \u2192 OA knee", "BMD_leg","OA_knee", CB["red"]),
        ]),
        ("c", "Reverse / null pairs (non-significant)", [
            ("OA knee \u2192 fracture", "OA_knee","fracture", CB["grey"]),
            ("OA knee \u2192 BMD total", "OA_knee","BMD_total", CB["grey"]),
            ("OA knee \u2192 BMD femoral neck", "OA_knee","BMD_femur_neck", CB["grey"]),
        ]),
    ]
    fig = plt.figure(figsize=(7.2, 6.6))
    gs = fig.add_gridspec(3, 1, height_ratios=[0.85, 1.30, 0.95], hspace=0.55)
    for k, (lab, title, items) in enumerate(panels):
        ax = fig.add_subplot(gs[k])
        ys = list(range(len(items)))[::-1]
        los, his = [], []
        for y, (ylab, e, o, col) in zip(ys, items):
            r = get(e, o)
            b = num(r["beta"]); se = num(r["se"])
            ax.errorbar(b, y, xerr=1.96*se, fmt="o", color=col, ecolor=col,
                        markersize=4.5, elinewidth=1.2, capsize=2.5, zorder=3)
            los.append(b - 1.96*se); his.append(b + 1.96*se)
        ax.axvline(0, color=CB["grey"], ls="--", lw=0.8, zorder=1)
        ax.set_yticks(ys); ax.set_yticklabels([t for t,_,_,_ in items], fontsize=7)
        ax.set_xlabel("MR effect (\u03b2 per 1 SD higher exposure)", fontsize=7.5)
        ax.set_xlim(*limits(los, his))
        if LIMITS is not None: LIMITS[f"fig3_{lab}"] = list(ax.get_xlim())
        panel_label(ax, lab)
        ax.set_title(title, fontsize=8, fontweight="bold", loc="left")
    fig.text(0.5, 0.005, "Points: IVW fixed-effects estimate with 95% CI.  "
             "Blue = protective, red = risk-increasing, grey = null.  "
             "Scale is N-corrected (BMD N = 31,873; \u00a73.7).",
             ha="center", fontsize=6.5, color=CB["grey"])
    save(fig, "Figure3_MR_forest")
    plt.close(fig)

# ================= FIGURE 4 : SNP heritability bar =================
def fig4(LIMITS=None):
    # load h2 from rg_h2_trueN (OA + BMD) and N_sensitivity_decision (fracture)
    h2, se = {}, {}
    with open(os.path.join(GEN, "rg_h2_trueN.csv")) as f:
        for row in csv.DictReader(f):
            h2[row["phenotype"]] = float(row["h2"])
            se[row["phenotype"]] = float(row["h2_se"])
    with open(os.path.join(GEN, "N_sensitivity_decision.csv")) as f:
        for row in csv.DictReader(f):
            if row["source"] == "264973_GEFOS发现集":
                h2["fracture"] = float(row["h2_obs"])
                se["fracture"] = float(row["h2_se"])
                break

    labels = ["OA_knee","OA_selfrep","fracture","BMD_total","BMD_femur_neck","BMD_femur_total",
              "BMD_lumbar","BMD_spine","BMD_arms","BMD_leg","BMD_pelvis","BMD_rib","BMD_trunk"]
    pretty = ["OA (knee)","OA (self-rep)","Fracture","BMD total","BMD femoral neck","BMD femoral total",
              "BMD lumbar","BMD spine","BMD arms","BMD legs","BMD pelvis","BMD rib","BMD trunk"]
    # group colours
    colors = [CB["purple"], CB["purple"], CB["teal"]] + [CB["cb_blue"]] * 10
    y = np.arange(len(labels))
    vals = [h2[l] for l in labels]
    err = [1.96 * se[l] for l in labels]
    # liability-scale h2 for binary traits
    liability = {"OA_knee":0.0935, "OA_selfrep":0.0726, "fracture":0.0620}
    liab_vals = [liability.get(l, np.nan) for l in labels]

    fig, ax = plt.subplots(figsize=(6.0, 5.6))
    bars = ax.barh(y, vals, color=colors, height=0.65, zorder=2, edgecolor="white", lw=0.4)
    ax.errorbar(vals, y, xerr=err, fmt="none", ecolor=CB["ink"], elinewidth=0.8,
                capsize=2.5, zorder=3)
    # diamond markers for liability h2
    ax.scatter([liability.get(l, np.nan) for l in labels], y, marker="D",
               color=CB["orange"], s=35, zorder=4, edgecolors="white", linewidths=0.6,
               label="Liability-scale h²")
    # annotate exact observed h2 values
    for yy, v, l in zip(y, vals, labels):
        ax.text(v + err[labels.index(l)] + 0.012, yy, f"{v:.3f}", va="center", ha="left",
                fontsize=6.5, color=CB["ink"])
    # add group legend manually, placed in empty space to the right
    from matplotlib.patches import Patch
    legend_handles = [
        Patch(facecolor=CB["purple"], edgecolor="white", label="OA definitions"),
        Patch(facecolor=CB["teal"], edgecolor="white", label="Fracture"),
        Patch(facecolor=CB["cb_blue"], edgecolor="white", label="BMD sites"),
        Line2D([0],[0], marker="D", color="w", markerfacecolor=CB["orange"],
               markersize=6, label="Liability-scale h²"),
    ]
    ax.legend(handles=legend_handles, loc="upper left", bbox_to_anchor=(1.02, 0.95),
              fontsize=7, frameon=True, framealpha=0.95, edgecolor="none")
    ax.set_yticks(y); ax.set_yticklabels(pretty, fontsize=7.5)
    ax.set_xlabel("SNP-based heritability (h², observed scale)", fontsize=8)
    ax.set_xlim(0, 0.52)
    if LIMITS is not None: LIMITS["fig4"] = list(ax.get_xlim())
    ax.set_title("SNP-based heritability of the 13 skeletal phenotypes", fontsize=9,
                 fontweight="bold", loc="left")
    ax.grid(axis="x", ls="-", lw=0.4, color="#E0E0E0", zorder=0)
    ax.invert_yaxis()
    fig.text(0.5, 0.005, "Bars: observed-scale h² with 95% CI; diamond markers: liability-scale h² for OA and fracture. "
             "BMD is a quantitative trait and is not liability-converted.",
             ha="center", fontsize=6.5, color=CB["grey"])
    save(fig, "Figure4_h2_bar")
    plt.close(fig)

# ================= FIGURE 5 : causal network =================
def fig5():
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")

    # nodes
    def node(x, y, w, h, text, fc, tc="white"):
        ax.add_patch(FancyBboxPatch((x-w/2, y-h/2), w, h,
                     boxstyle="round,pad=0.02,rounding_size=0.18",
                     fc=fc, ec="none", zorder=2))
        ax.text(x, y, text, ha="center", va="center", fontsize=8.5,
                fontweight="bold", color=tc, zorder=3)

    # OA block (left)
    node(2.0, 7.4, 2.4, 1.0, "OA (knee)", CB["purple"])
    node(2.0, 5.6, 2.4, 1.0, "OA (self-rep)", CB["purple"])
    # BMD cluster (right centre)
    node(7.6, 6.5, 2.8, 1.6, "BMD\n(10 regional sites)", CB["cb_blue"])
    # Fracture (bottom right)
    node(7.6, 1.8, 2.4, 1.0, "Fracture", CB["teal"])

    def arrow(x1, y1, x2, y2, color, label, lx, ly, rad=0.0):
        arr = FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                     mutation_scale=14, color=color, lw=2.2,
                     connectionstyle=f"arc3,rad={rad}", zorder=1,
                     path_effects=[pe.withStroke(linewidth=3.2, foreground="white")])
        ax.add_patch(arr)
        ax.text(lx, ly, label, ha="center", va="center", fontsize=7.2,
                fontweight="bold", color=color, zorder=4,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.9))

    # BMD -> Fracture (protective, blue)
    arrow(7.0, 6.0, 7.0, 2.3, CB["cb_blue"], "Protective\nOR 1.40\u20131.42", 8.7, 4.1, rad=-0.15)
    # BMD -> OA (causal risk, red)
    arrow(6.6, 6.1, 3.4, 6.9, CB["red"], "Causal \u2191 risk", 5.0, 7.7, rad=0.18)
    # OA <-> BMD genetic correlation (dashed grey)
    arr = FancyArrowPatch((3.3, 5.2), (6.4, 6.0), arrowstyle="-|>", mutation_scale=11,
               color=CB["grey"], lw=1.6, ls=(0, (3, 2)), zorder=1,
               connectionstyle="arc3,rad=-0.2",
               path_effects=[pe.withStroke(linewidth=3.0, foreground="white")])
    ax.add_patch(arr)
    ax.text(4.9, 4.7, "Genetic correlation\nr_g +0.12 to +0.25", ha="center", va="center",
            fontsize=7.0, color=CB["grey"], zorder=4,
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.9))
    # OA -> Fracture null (grey dashed, crossed)
    arr2 = FancyArrowPatch((3.4, 5.2), (6.6, 2.4), arrowstyle="-|>", mutation_scale=10,
                color=CB["grey"], lw=1.4, ls=(0, (3, 2)), zorder=1,
                connectionstyle="arc3,rad=0.25",
                path_effects=[pe.withStroke(linewidth=3.0, foreground="white")])
    ax.add_patch(arr2)
    ax.text(5.2, 3.4, "Null (ns)", ha="center", va="center", fontsize=7.0,
            color=CB["grey"], zorder=4, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.9))

    ax.set_title("Causal architecture of the OA\u2013BMD\u2013fracture triangle",
                 fontsize=9.5, fontweight="bold", pad=8)
    save(fig, "Figure5_causal_network")
    plt.close(fig)

if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    LIMITS = {}
    fig1()
    fig2(LIMITS)
    fig3(LIMITS)
    fig4(LIMITS)
    fig5()
    # 把每张图实际使用的坐标轴范围落盘，供 verify_figures.py 做端到端校验
    with open(os.path.join(GEN, "figure_axis_limits.json"), "w", encoding="utf-8") as fh:
        json.dump(LIMITS, fh, indent=2)
    print("axis limits ->", os.path.join(GEN, "figure_axis_limits.json"))
    print("ALL FIGURES DONE")
