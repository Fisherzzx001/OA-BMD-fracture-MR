# -*- coding: utf-8 -*-
"""
build_figures_v17.py
=============================================================================
OA–BMD–fracture 稿件 · 顶刊级配图重建 + 新增图
-----------------------------------------------------------------------------
规范依据：
  * nature-figure skill（Nature portfolio 视觉标准：Arial ≥7pt、无上/右框、
    面板小写粗体标签 a/b/c、去冗余、色盲安全、矢量 PDF、线稿 ≥600 dpi）
  * scienceplots（'science' + 'nature' 版式基线）
  * Okabe–Ito 色盲安全调色板

数据来源（全部为稿件真实统计输出，无任何编造数值）：
  code/output/genetics/rg_matrix_13.csv            Figure 1
  code/output/genetics/rg_results_long_13.csv      Figure 1(FDR)/2, Table 2
  code/output/genetics/rg_h2_trueN.csv             Figure 4
  code/output/genetics/N_sensitivity_decision.csv  Figure 4/S6
  code/output/genetics/oa_h2_liability.csv         Figure 4
  code/output/genetics/AUDIT_mr_recomputed.csv     Figure 3/6/7/S1/S2/S7（= Table 3）
  code/output/genetics/harmonised/*.csv            Figure 6/7/S1/S7（逐 SNP）
  code/output/genetics/mr_power.csv                Figure 8（= Table 4）
  code/output/genetics/overlap_statistics.csv      Figure 9
  code/output/mrlap_overlap_results.csv            Figure 10（新增，§2.8）
  code/output/genetics/gefos_benchmark.json        SupplFig S3
  Supplement_S3（h² scales as 1/N 量表）                   SupplFig S5（新增）

新增图（加强表达）：
  Figure 10        MRlap 样本重叠校正（§2.8，此前仅有文字与表，无图）
  SupplFig S5      OA_knee h² 对代表 N 的敏感性（Supplement S3 的图形化）
  SupplFig S6      fracture N 选择的证据链（Table 1 §§ 的图形化）
  SupplFig S7      多效性 / 异质性诊断（10 对 MR 的稳健性全景）

输出：Figures_v17_Nature/  （PNG 600 dpi + 矢量 PDF + TIFF 300 dpi）
      同时在包内保留 Figures_Original_v14/（原图，未改动）
=============================================================================
"""
import os
import json
import shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.colors import TwoSlopeNorm
import matplotlib.patheffects as pe
from scipy import stats

# ============================================================ 路径与输出
PROJECT = "C:/Users/Administrator/WorkBuddy/2026-09-08-22-17-40/OA_B_PROJECT"
GEN = os.path.join(PROJECT, "code", "output", "genetics")
HARM = os.path.join(GEN, "harmonised")
OUT = os.path.join(PROJECT, "Figures_v17_Nature")
ORIG = os.path.join(OUT, "Figures_Original_v14")
os.makedirs(OUT, exist_ok=True)

# ============================================================ 顶刊版式
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

# Okabe–Ito 色盲安全配色
C = {
    "blue": "#0072B2", "orange": "#E69F00", "vermillion": "#D55E00",
    "green": "#009E73", "sky": "#56B4E9", "purple": "#CC79A7",
    "yellow": "#F0E442", "grey": "#8C8C8C", "ink": "#000000",
    "faint": "#DCDCDC",
}

PHEN = {
    "OA_knee": "OA (knee)", "OA_selfrep": "OA (self-reported)", "fracture": "Any fracture",
    "BMD_total": "BMD total body", "BMD_femur_neck": "BMD femoral neck",
    "BMD_femur_total": "BMD femoral total", "BMD_lumbar": "BMD lumbar spine",
    "BMD_spine": "BMD spine", "BMD_arms": "BMD arms", "BMD_leg": "BMD legs",
    "BMD_pelvis": "BMD pelvis", "BMD_rib": "BMD ribs", "BMD_trunk": "BMD trunk",
}
ORDER = ["OA_knee", "OA_selfrep", "fracture", "BMD_total", "BMD_femur_neck",
         "BMD_femur_total", "BMD_lumbar", "BMD_spine", "BMD_arms", "BMD_leg",
         "BMD_pelvis", "BMD_rib", "BMD_trunk"]

PLOTTED = {}          # 记录每张图的每一个绘制数值，供独立审计
SAVED = []


def save(fig, name, tiff_dpi=300):
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=600)
    fig.savefig(os.path.join(OUT, name + ".pdf"))
    # TIFF：期刊生产格式 — 扁平化为 RGB、Adobe-deflate 压缩（体积小、无损）
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


# ============================================================ 数据加载
def load_all():
    d = {}
    d["rg"] = pd.read_csv(os.path.join(GEN, "rg_matrix_13.csv"), index_col=0)
    d["rg_long"] = pd.read_csv(os.path.join(GEN, "rg_results_long_13.csv"))
    d["h2"] = pd.read_csv(os.path.join(GEN, "rg_h2_trueN.csv"))
    d["liab"] = pd.read_csv(os.path.join(GEN, "oa_h2_liability.csv"))
    d["nsens"] = pd.read_csv(os.path.join(GEN, "N_sensitivity_decision.csv"))
    d["mr"] = pd.read_csv(os.path.join(GEN, "AUDIT_mr_recomputed.csv"))
    d["power"] = pd.read_csv(os.path.join(GEN, "mr_power.csv"))
    d["overlap"] = pd.read_csv(os.path.join(GEN, "overlap_statistics.csv"))
    d["mrlap"] = pd.read_csv(os.path.join(PROJECT, "code", "output", "mrlap_overlap_results.csv"))
    d["loo"] = pd.read_csv(os.path.join(HARM, "_loo.csv"))
    with open(os.path.join(GEN, "gefos_benchmark.json"), encoding="utf-8") as f:
        d["gefos"] = json.load(f)
    # q(FDR) 矩阵
    lg = d["rg_long"]
    q = pd.DataFrame(np.nan, index=d["rg"].index, columns=d["rg"].columns)
    for _, r in lg.iterrows():
        q.loc[r.trait1, r.trait2] = r.q_fdr
        q.loc[r.trait2, r.trait1] = r.q_fdr
    d["q"] = q
    return d


D = load_all()


def harm(e, o):
    return pd.read_csv(os.path.join(HARM, f"{e}__{o}.csv"))


def mr_row(e, o):
    m = D["mr"][(D["mr"].exposure == e) & (D["mr"].outcome == o)]
    return m.iloc[0]


def ivw_fit(x, y, se):
    w = 1.0 / se ** 2
    return float(np.sum(w * x * y) / np.sum(w * x * x))


def ivw_se(x, y, se):
    w = 1.0 / se ** 2
    return float(np.sqrt(1.0 / np.sum(w * x * x)))


def egger_fit(x, y, se):
    w = 1.0 / se ** 2
    X = np.column_stack([np.ones_like(x), x])
    b = np.linalg.solve(X.T @ (w[:, None] * X), X.T @ (w * y))
    return float(b[0]), float(b[1])


def plab(e, o):
    return f"{PHEN[e].replace('BMD ', '')} → {PHEN[o].replace(' (knee)', ' knee')}"


# ======================================================================
# Figure 1 — 13 表型遗传相关热图
# ======================================================================
def fig1():
    rg = D["rg"].loc[ORDER, ORDER].copy()
    q = D["q"].loc[ORDER, ORDER].copy()
    M = rg.values.astype(float)
    Q = q.values.astype(float)
    n = len(ORDER)
    show = M.copy()
    np.fill_diagonal(show, np.nan)

    fig = plt.figure(figsize=(7.2, 6.5))
    gs = fig.add_gridspec(1, 2, width_ratios=[100, 4], wspace=0.06)
    ax = fig.add_subplot(gs[0])
    cb = fig.add_subplot(gs[1])

    norm = TwoSlopeNorm(vcenter=0.0, vmin=-1.0, vmax=1.0)
    im = ax.imshow(show, cmap="RdBu_r", norm=norm, aspect="equal")

    for i in range(n):
        for j in range(n):
            if i == j:
                ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, fc="#F2F2F2", ec="white", lw=.5))
                ax.text(j, i, "—", ha="center", va="center", fontsize=6.5, color=C["grey"])
                continue
            v, qq = M[i, j], Q[i, j]
            star = ""
            if not np.isnan(qq):
                if qq < 0.001:
                    star = "***"
                elif qq < 0.01:
                    star = "**"
                elif qq < 0.05:
                    star = "*"
            txt = f"{v:.2f}{star}"
            ax.text(j, i, txt, ha="center", va="center", fontsize=5.4,
                    color="white" if abs(v) > 0.55 else C["ink"])

    for pos in (1.5, 2.5):
        ax.axhline(pos, color="black", lw=1.2)
        ax.axvline(pos, color="black", lw=1.2)
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels([PHEN[t] for t in ORDER], rotation=45, ha="right", fontsize=6.5)
    ax.set_yticklabels([PHEN[t] for t in ORDER], fontsize=6.5)
    ax.set_xticks(np.arange(-.5, n, 1), minor=True)
    ax.set_yticks(np.arange(-.5, n, 1), minor=True)
    ax.grid(which="minor", color="white", lw=0.6)
    ax.tick_params(which="minor", length=0)
    ax.set_title("LDSC genetic correlation (r_g), 13 skeletal phenotypes", fontsize=9, pad=8)

    cbar = fig.colorbar(im, cax=cb)
    cbar.set_label("r_g", fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    cbar.outline.set_linewidth(0.8)

    PLOTTED["Figure1"] = {
        "order": ORDER,
        "labels": [PHEN[t] for t in ORDER],
        "rg": [[None if np.isnan(show[i, j]) else float(show[i, j]) for j in range(n)] for i in range(n)],
        "q_fdr": [[None if np.isnan(Q[i, j]) else float(Q[i, j]) for j in range(n)] for i in range(n)],
        "source": "rg_matrix_13.csv + rg_results_long_13.csv(q_fdr)",
    }
    save(fig, "Figure1_rg_heatmap_13")


# ======================================================================
# Figure 2 — 关键遗传相关森林图（4 面板）
# ======================================================================
def _forest(ax, rows, title, xlab="r_g"):
    y = np.arange(len(rows))[::-1]
    for k, (lab, v, se, qq) in enumerate(rows):
        lo, hi = v - 1.96 * se, v + 1.96 * se
        sig = (not np.isnan(qq)) and qq < 0.05
        col = C["grey"] if not sig else (C["vermillion"] if v > 0 else C["blue"])
        ax.plot([lo, hi], [y[k], y[k]], color=col, lw=1.3, solid_capstyle="butt")
        ax.plot([lo, lo], [y[k] - .12, y[k] + .12], color=col, lw=1.0)
        ax.plot([hi, hi], [y[k] - .12, y[k] + .12], color=col, lw=1.0)
        ax.scatter([v], [y[k]], s=26, facecolor=col if sig else "white",
                   edgecolor=col, zorder=3, linewidths=1.0)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=6.8)
    ax.axvline(0, color=C["ink"], lw=0.8, ls=(0, (4, 2)))
    ax.set_xlabel(xlab, fontsize=8)
    ax.set_title(title, fontsize=8.5, pad=6)
    ax.set_ylim(-0.6, len(rows) - 0.4)
    return rows


def fig2():
    lg = D["rg_long"]
    BMD = [t for t in ORDER if t.startswith("BMD")]

    def pick(t1, t2):
        m = lg[((lg.trait1 == t1) & (lg.trait2 == t2)) | ((lg.trait1 == t2) & (lg.trait2 == t1))]
        return m.iloc[0] if len(m) else None

    def rows_for(a, others):
        out = []
        for b in others:
            r = pick(a, b)
            if r is None:
                continue
            out.append((PHEN[b].replace("BMD ", ""), float(r.rg), float(r.se), float(r.q_fdr)))
        return sorted(out, key=lambda t: t[1])

    ra = rows_for("fracture", BMD)
    rb = rows_for("OA_knee", BMD)
    rc = rows_for("OA_selfrep", BMD)
    rd = []
    for a, b, lab in [("OA_knee", "OA_selfrep", "OA (knee) – OA (self-rep)"),
                      ("OA_knee", "fracture", "OA (knee) – fracture"),
                      ("OA_selfrep", "fracture", "OA (self-rep) – fracture")]:
        r = pick(a, b)
        rd.append((lab, float(r.rg), float(r.se), float(r.q_fdr)))
    rd = sorted(rd, key=lambda t: t[1], reverse=True)

    fig = plt.figure(figsize=(7.2, 6.6))
    gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.34)
    axs = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]),
           fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]
    _forest(axs[0], ra, "Any fracture vs 10 BMD sites")
    _forest(axs[1], rb, "Knee OA vs 10 BMD sites")
    _forest(axs[2], rc, "Self-reported OA vs 10 BMD sites")
    _forest(axs[3], rd, "OA–OA and OA–fracture")
    for ax, lab in zip(axs, "abcd"):
        panel_label(ax, lab, x=-0.30, y=1.10)

    PLOTTED["Figure2"] = {
        "a_fracture_vs_BMD": [{"pair": r[0], "rg": r[1], "se": r[2], "q_fdr": r[3]} for r in ra],
        "b_kneeOA_vs_BMD": [{"pair": r[0], "rg": r[1], "se": r[2], "q_fdr": r[3]} for r in rb],
        "c_selfrepOA_vs_BMD": [{"pair": r[0], "rg": r[1], "se": r[2], "q_fdr": r[3]} for r in rc],
        "d_OA_OA_fracture": [{"pair": r[0], "rg": r[1], "se": r[2], "q_fdr": r[3]} for r in rd],
        "ci": "95% CI = r_g ± 1.96 × SE (block jackknife)",
        "source": "rg_results_long_13.csv",
    }
    save(fig, "Figure2_forest_key_results")


# ======================================================================
# Figure 3 — MR 森林图（N 校正效应尺度，= Table 3）
# ======================================================================
def fig3():
    mr = D["mr"]
    ctrl = mr[(mr.outcome == "fracture") & (mr.exposure.str.startswith("BMD"))]
    fwd = mr[(mr.outcome == "OA_knee")]
    null = mr[(mr.exposure == "OA_knee")]

    def rows(sub):
        out = []
        for _, r in sub.iterrows():
            out.append((plab(r.exposure, r.outcome), float(r.beta), float(r.se),
                        float(r.p), int(r.K)))
        return sorted(out, key=lambda t: t[1])

    ra, rb, rc = rows(ctrl), rows(fwd), rows(null)

    fig = plt.figure(figsize=(7.2, 6.0))
    hs = [max(1, len(ra)), max(1, len(rb)), max(1, len(rc))]
    gs = fig.add_gridspec(3, 1, height_ratios=hs, hspace=0.52)
    axs = [fig.add_subplot(gs[i]) for i in range(3)]

    titles = ["BMD → any fracture (positive control)",
              "BMD → knee OA (five exposure sites)",
              "Knee OA → fracture / BMD (non-significant, underpowered)"]
    for ax, dat, ti, lab in zip(axs, [ra, rb, rc], titles, ["a", "b", "c"]):
        y = np.arange(len(dat))[::-1]
        for k, (nm, b, se, p, K) in enumerate(dat):
            lo, hi = b - 1.96 * se, b + 1.96 * se
            if p >= 0.05:
                col = C["grey"]
            elif b < 0:
                col = C["blue"]
            else:
                col = C["vermillion"]
            ax.plot([lo, hi], [y[k], y[k]], color=col, lw=1.4, solid_capstyle="butt")
            ax.plot([lo, lo], [y[k] - .12, y[k] + .12], color=col, lw=1.1)
            ax.plot([hi, hi], [y[k] - .12, y[k] + .12], color=col, lw=1.1)
            ax.scatter([b], [y[k]], s=34, facecolor=col, edgecolor="white",
                       zorder=3, linewidths=0.6)
            ax.text(hi, y[k] + 0.30, f"{b:+.3f} ({b-1.96*se:+.3f}, {b+1.96*se:+.3f})  K={K}",
                    fontsize=6.2, color=C["grey"], va="bottom")
        ax.set_yticks(y)
        ax.set_yticklabels([d[0] for d in dat], fontsize=7)
        ax.axvline(0, color=C["ink"], lw=0.8, ls=(0, (4, 2)))
        ax.set_ylim(-0.7, len(dat) - 0.1)
        ax.set_title(ti, fontsize=8.5, pad=6)
        panel_label(ax, lab, x=-0.22, y=1.12)
    axs[2].set_xlabel("IVW β (N-corrected scale) ± 95% CI", fontsize=8)

    PLOTTED["Figure3"] = {
        "a_positive_control": [{"pair": d[0], "beta": d[1], "se": d[2], "p": d[3], "K": d[4]} for d in ra],
        "b_BMD_to_kneeOA": [{"pair": d[0], "beta": d[1], "se": d[2], "p": d[3], "K": d[4]} for d in rb],
        "c_null_pairs": [{"pair": d[0], "beta": d[1], "se": d[2], "p": d[3], "K": d[4]} for d in rc],
        "ci": "95% CI = β ± 1.96 × SE",
        "source": "AUDIT_mr_recomputed.csv (identical to manuscript Table 3)",
    }
    save(fig, "Figure3_MR_forest")


# ======================================================================
# Figure 4 — SNP 遗传力（观测尺度 + 责任尺度）
# ======================================================================
def fig4():
    h2 = D["h2"].set_index("phenotype")
    ns = D["nsens"]
    liab = D["liab"].set_index("trait")
    frac = ns[(ns.trait == "fracture") & (ns.source.str.contains("264973"))].iloc[0]

    rows = []
    for t in ORDER:
        if t == "fracture":
            rows.append((t, float(frac.h2_obs), float(frac.h2_se)))
        else:
            r = h2.loc[t]
            rows.append((t, float(r.h2), float(r.h2_se)))
    rows_sorted = sorted(rows, key=lambda r: r[1])
    y = np.arange(len(rows_sorted))
    vals = [r[1] for r in rows_sorted]
    errs = [1.96 * r[2] for r in rows_sorted]

    fig, ax = plt.subplots(figsize=(7.2, 3.9))
    cols = [C["purple"] if t in ("OA_knee", "OA_selfrep") else
            (C["sky"] if t == "fracture" else C["blue"]) for t, _, _ in rows_sorted]
    ax.barh(y, vals, xerr=errs, color=cols, edgecolor="white", height=0.68,
            error_kw=dict(ecolor=C["ink"], elinewidth=0.9, capsize=2.6, capthick=0.9))
    for k, (t, v, se) in enumerate(rows_sorted):
        if t in liab.index:
            hl = float(liab.loc[t, "h2_liab_primary"])
            ax.scatter([hl], [k], marker="D", s=42, color=C["vermillion"],
                       edgecolor="white", linewidths=0.6, zorder=4)
    ax.set_yticks(y)
    ax.set_yticklabels([PHEN[t] for t, _, _ in rows_sorted], fontsize=7)
    ax.set_xlabel("SNP heritability h² (observed scale)", fontsize=8)
    ax.set_xlim(0, 0.52)
    ax.xaxis.grid(True, color=C["faint"], lw=0.6)
    ax.set_axisbelow(True)

    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    leg = [Patch(facecolor=C["blue"], label="BMD sites (quantitative, DXA)"),
           Patch(facecolor=C["purple"], label="OA traits (binary)"),
           Patch(facecolor=C["sky"], label="Any fracture (binary)"),
           Line2D([], [], marker="D", color="w", markerfacecolor=C["vermillion"],
                  markeredgecolor="white", markersize=6,
                  label="h² on liability scale (Lee et al. 2011, Eq 23)"),
           Line2D([], [], color=C["ink"], lw=1.0, label="95% CI (h² ± 1.96 × SE)")]
    ax.legend(handles=leg, loc="lower right", fontsize=6.6)
    ax.set_title("SNP-based heritability of the 13 skeletal phenotypes", fontsize=9, pad=7)
    panel_label(ax, "", x=-0.10, y=1.06)

    PLOTTED["Figure4"] = {
        "traits": [PHEN[t] for t, _, _ in rows_sorted],
        "h2_obs": [float(v) for _, v, _ in rows_sorted],
        "se": [float(s) for _, _, s in rows_sorted],
        "h2_liability": {t: float(liab.loc[t, "h2_liab_primary"]) for t in liab.index},
        "source": "rg_h2_trueN.csv; fracture from N_sensitivity_decision.csv (N=264,973); "
                  "liability from oa_h2_liability.csv",
    }
    save(fig, "Figure4_h2_bar")


# ======================================================================
# Figure 5 — 概念性因果架构
# ======================================================================
def fig5():
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")

    def node(x, y, w, h, t, fc, fs=8.5):
        ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                     boxstyle="round,pad=0.03,rounding_size=0.30",
                     fc=fc, ec="none", zorder=2))
        ax.text(x, y, t, ha="center", va="center", fontsize=fs,
                fontweight="bold", color="white", zorder=3)

    def arrow(x1, y1, x2, y2, col, lab, lx, ly, rad=0.0, ls="-", lw=2.6):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                     mutation_scale=16, color=col, lw=lw, ls=ls, zorder=1,
                     connectionstyle=f"arc3,rad={rad}",
                     path_effects=[pe.withStroke(linewidth=3.6, foreground="white")]))
        ax.text(lx, ly, lab, ha="center", va="center", fontsize=6.8, color=col,
                zorder=5, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.30", fc="white", ec="none", alpha=0.95))

    node(5.0, 8.0, 3.6, 1.15, "BMD (10 regional DXA sites)", C["blue"])
    node(8.2, 2.0, 2.9, 1.15, "Any fracture", C["sky"], fs=8.5)
    node(1.8, 2.0, 2.9, 1.15, "Knee OA", C["purple"], fs=8.5)

    arrow(6.1, 7.4, 7.6, 2.7, C["blue"],
          "Protective\nβ = −0.337 (total), −0.351 (femoral neck)\nOR 1.40–1.42 per 1-SD lower BMD",
          7.55, 4.85, rad=-0.10)
    arrow(3.9, 7.4, 2.4, 2.7, C["vermillion"],
          "Increases risk\nβ = +0.058 to +0.095\n5 BMD sites, 5/5 significant",
          3.72, 4.55, rad=0.10)
    arrow(3.35, 2.25, 6.65, 2.25, C["grey"],
          "No well-powered evidence (power 0.14; excludes |β| ≥ 0.09)",
          5.0, 1.25, rad=0.0, ls=(0, (4, 2.2)), lw=1.8)
    arrow(2.30, 3.10, 3.75, 6.95, C["grey"],
          "Genetic correlation only\nr_g = +0.12 to +0.25 (20/20 q<0.05)",
          1.30, 5.60, rad=-0.42, ls=(0, (4, 2.2)), lw=1.8)

    ax.text(5.0, 9.35, "Shared genetic architecture, diverging causal fates",
            ha="center", va="center", fontsize=9.5, fontweight="bold")
    ax.text(5.0, 0.12, "Solid arrows: MR causal effects (survive MR-PRESSO, leave-one-out and\n"
            "MRlap sample-overlap correction, §2.8).   Dashed arrows: non-causal / uninformative.",
            ha="center", va="center", fontsize=6.2, color=C["ink"],
            bbox=dict(boxstyle="round,pad=0.35", fc="#F4F4F4", ec="none"))
    save(fig, "Figure5_causal_network")


# ======================================================================
# Figure 6 / SupplFig S1 — 逐 SNP 散点图
# ======================================================================
def _scatter(ax, e, o, lab):
    d = harm(e, o)
    r = mr_row(e, o)
    x = d.beta_e.values; y = d.beta_o.values
    sx = 1.96 * d.se_e.values; sy = 1.96 * d.se_o.values
    out = d.presso_outlier.values.astype(bool)

    ax.errorbar(x, y, xerr=sx, yerr=sy, fmt="none", ecolor="#C9C9C9",
                elinewidth=0.7, alpha=0.85, zorder=0)
    ax.scatter(x[~out], y[~out], s=20, facecolor=C["blue"], edgecolor="white",
               linewidths=0.5, zorder=3)
    if out.sum():
        ax.scatter(x[out], y[out], s=42, facecolor="none", edgecolor=C["vermillion"],
                   linewidths=1.2, zorder=4)

    b = ivw_fit(x, y, d.se_o.values)
    b0, b1 = egger_fit(x, y, d.se_o.values)
    wm = float(r.wmed)
    xs = np.linspace(x.min(), x.max(), 100)
    ax.plot(xs, b * xs, color="black", lw=1.6, zorder=5)
    ax.plot(xs, b0 + b1 * xs, color=C["orange"], lw=1.5, ls=(0, (5, 2.5)), zorder=5)
    ax.plot(xs, wm * xs, color=C["purple"], lw=1.5, ls=(0, (1.5, 2)), zorder=5)

    ax.axhline(0, color=C["faint"], lw=0.7, zorder=0)
    ax.axvline(0, color=C["faint"], lw=0.7, zorder=0)

    def pv(p):
        return "P < 1.0e-50" if p < 1e-50 else (f"P = {p:.2g}" if p >= 1e-4 else f"P = {p:.1e}")
    ax.text(0.03, 0.97, f"IVW  β = {b:+.3f}  ({pv(float(r.p))})\n"
                        f"MR-Egger  β = {b1:+.3f}\nWeighted median  β = {wm:+.3f}\n"
                        f"n = {int(r.K)} instruments",
            transform=ax.transAxes, ha="left", va="top", fontsize=6.2,
            bbox=dict(boxstyle="round,pad=0.30", fc="white", ec=C["faint"], lw=0.7))

    ax.set_xlabel(f"SNP effect on {PHEN[e]}", fontsize=7.5)
    ax.set_ylabel(f"SNP effect on {PHEN[o]}", fontsize=7.5)
    ax.set_title(lab, fontsize=8.2, pad=5)
    return dict(pair=f"{e}->{o}", n=int(len(d)), ivw=float(b), egger=float(b1),
                wmed=float(wm), p=float(r.p), K=int(r.K), n_outlier=int(out.sum()))


def fig6():
    pairs = [("BMD_total", "fracture"), ("BMD_femur_neck", "fracture"),
             ("BMD_total", "OA_knee"), ("OA_knee", "BMD_total")]
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 6.6))
    recs = []
    for k, (e, o) in enumerate(pairs):
        ax = axs[k // 2, k % 2]
        lab = f"({chr(97+k)})  {plab(e, o)}"
        r = _scatter(ax, e, o, lab)
        ax.set_title(lab[5:], fontsize=8.2, pad=5)
        panel_label(ax, chr(97 + k), x=-0.20, y=1.08)
        recs.append(r)
    from matplotlib.lines import Line2D
    leg = [Line2D([], [], color="black", lw=1.6, label="IVW (fixed effects)"),
           Line2D([], [], color=C["orange"], lw=1.5, ls=(0, (5, 2.5)), label="MR-Egger"),
           Line2D([], [], color=C["purple"], lw=1.5, ls=(0, (1.5, 2)), label="Weighted median"),
           Line2D([], [], marker="o", color="w", markerfacecolor=C["blue"], markersize=5, label="Instrument"),
           Line2D([], [], marker="o", color="w", markerfacecolor="none",
                  markeredgecolor=C["vermillion"], markersize=6, label="MR-PRESSO outlier")]
    fig.legend(handles=leg, loc="lower center", ncol=5, fontsize=6.8,
               frameon=False, bbox_to_anchor=(0.5, -0.012), columnspacing=1.4,
               handletextpad=0.5)
    fig.suptitle("Per-SNP Mendelian randomisation scatter plots", fontsize=9.5, y=0.995)
    PLOTTED["Figure6"] = {"panels": recs, "source": "harmonised/*.csv + AUDIT_mr_recomputed.csv"}
    save(fig, "Figure6_MR_scatter")


def suppl1():
    pairs = [("BMD_femur_neck", "OA_knee"), ("BMD_lumbar", "OA_knee"),
             ("BMD_arms", "OA_knee"), ("BMD_leg", "OA_knee"),
             ("OA_knee", "fracture"), ("OA_knee", "BMD_femur_neck")]
    fig, axs = plt.subplots(3, 2, figsize=(7.2, 9.4))
    recs = []
    for k, (e, o) in enumerate(pairs):
        ax = axs[k // 2, k % 2]
        lab = f"({chr(97+k)})  {plab(e, o)}"
        r = _scatter(ax, e, o, lab)
        ax.set_title(lab[5:], fontsize=8.2, pad=5)
        panel_label(ax, chr(97 + k), x=-0.20, y=1.08)
        recs.append(r)
    fig.suptitle("Per-SNP MR scatter plots — remaining six pairs", fontsize=9.5, y=0.995)
    PLOTTED["SupplFig1"] = {"panels": recs, "source": "harmonised/*.csv + AUDIT_mr_recomputed.csv"}
    save(fig, "SupplFig1_MR_scatter_remaining")


# ======================================================================
# Figure 7 — MR-PRESSO 漏斗图
# ======================================================================
def _funnel(ax, e, o, lab):
    d = harm(e, o)
    r = mr_row(e, o)
    th = d.ratio.values
    prec = 1.0 / d.se_ratio.values
    out = d.presso_outlier.values.astype(bool)
    b = float(r.beta)
    b0, b1 = egger_fit(d.beta_e.values, d.beta_o.values, d.se_o.values)

    ax.scatter(th[~out], prec[~out], s=20, facecolor=C["blue"], edgecolor="white",
               linewidths=0.5, zorder=3)
    if out.sum():
        ax.scatter(th[out], prec[out], s=42, facecolor="none", edgecolor=C["vermillion"],
                   linewidths=1.2, zorder=4)
    ax.axvline(b, color="black", lw=1.6, zorder=5)
    ax.axvline(b1, color=C["orange"], lw=1.5, ls=(0, (5, 2.5)), zorder=5)
    yy = np.linspace(prec.min() * 0.85, prec.max() * 1.05, 200)
    ax.plot(b - 1.96 / yy, yy, color=C["grey"], lw=1.0, ls=(0, (4, 2)), zorder=2)
    ax.plot(b + 1.96 / yy, yy, color=C["grey"], lw=1.0, ls=(0, (4, 2)), zorder=2)

    ax.set_xlabel("SNP-specific Wald ratio", fontsize=7.5)
    ax.set_ylabel("Precision (1 / SE)", fontsize=7.5)
    ax.set_title(lab, fontsize=8.2, pad=5)
    ax.text(0.03, 0.97, f"IVW β = {b:+.3f}\nMR-Egger β = {b1:+.3f}\n"
                        f"Cochran I² = {float(r.I2)*100:.0f}%",
            transform=ax.transAxes, ha="left", va="top", fontsize=6.2,
            bbox=dict(boxstyle="round,pad=0.30", fc="white", ec=C["faint"], lw=0.7))
    return dict(pair=f"{e}->{o}", n=int(len(d)), ivw=b, egger=b1,
                I2=float(r.I2), n_outlier=int(out.sum()))


def fig7():
    pairs = [("BMD_total", "fracture"), ("BMD_femur_neck", "fracture"),
             ("BMD_total", "OA_knee"), ("OA_knee", "BMD_total")]
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 6.6))
    recs = []
    for k, (e, o) in enumerate(pairs):
        ax = axs[k // 2, k % 2]
        lab = f"({chr(97+k)})  {plab(e, o)}"
        r = _funnel(ax, e, o, lab)
        ax.set_title(lab[5:], fontsize=8.2, pad=5)
        panel_label(ax, chr(97 + k), x=-0.20, y=1.08)
        recs.append(r)
    fig.suptitle("MR-PRESSO funnel plots (Wald ratio vs precision)", fontsize=9.5, y=0.995)
    PLOTTED["Figure7"] = {"panels": recs, "source": "harmonised/*.csv + AUDIT_mr_recomputed.csv"}
    save(fig, "Figure7_MR_funnel")


# ======================================================================
# Figure 8 — 统计功效
# ======================================================================
def fig8():
    pw = D["power"].copy()
    pw["pair"] = [plab(e, o) for e, o in zip(pw.exposure, pw.outcome)]
    pw = pw.sort_values("power_observed")

    fig = plt.figure(figsize=(7.2, 3.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.05, 1.0], wspace=0.44)
    axa, axb = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])

    z = np.linspace(0, 17, 500)
    power = stats.norm.sf(1.96 - z) + stats.norm.cdf(-1.96 - z)
    axa.plot(z, power, color="black", lw=1.6)
    for _, r in pw.iterrows():
        col = C["blue"] if r.verdict == "powered" else C["vermillion"]
        axa.scatter([r.abs_z], [r.power_observed], s=42,
                    facecolor=col if r.verdict == "powered" else "none",
                    edgecolor=col, linewidths=1.2, zorder=4)
    axa.axhline(0.8, color=C["grey"], lw=1.0, ls=(0, (5, 2.5)))
    axa.axvline(2.80, color=C["orange"], lw=1.1, ls=(0, (1.5, 2)))
    axa.text(0.55, 0.83, "80% power", fontsize=6.5, color=C["grey"])
    axa.text(2.95, 0.30, "|z| = 2.80\n(MDE)", fontsize=6.2, color=C["orange"])
    axa.set_xlabel("|β| / SE  (standardised effect size)", fontsize=8)
    axa.set_ylabel("Achieved power (two-sided α = 0.05)", fontsize=8)
    axa.set_ylim(0, 1.05)
    axa.set_title("Achieved power vs effect size", fontsize=8.5, pad=6)
    panel_label(axa, "a")

    y = np.arange(len(pw))
    cols = [C["blue"] if v >= 0.8 else C["vermillion"] for v in pw.power_observed]
    axb.barh(y, pw.power_observed.values, color=cols, edgecolor="white", height=0.68)
    for k, (_, r) in enumerate(pw.iterrows()):
        axb.text(min(r.power_observed + 0.02, 1.0), k,
                 f"{r.power_observed:.2f}  (K={int(r.K)})", va="center", fontsize=6.0,
                 color=C["ink"])
    axb.set_yticks(y)
    axb.set_yticklabels(pw.pair.values, fontsize=6.4)
    axb.axvline(0.8, color=C["orange"], lw=1.1, ls=(0, (1.5, 2)))
    axb.set_xlim(0, 1.28)
    axb.set_xlabel("Achieved power at the observed point estimate", fontsize=8)
    axb.set_title("Per-pair achieved power", fontsize=8.5, pad=6)
    panel_label(axb, "b", x=-0.22)

    PLOTTED["Figure8"] = {
        "panel_a_curve_power_at_z2.8": float(stats.norm.sf(1.96 - 2.80) + stats.norm.cdf(-1.96 - 2.80)),
        "pairs": [{"pair": r.pair, "K": int(r.K), "abs_z": float(r.abs_z),
                   "power": float(r.power_observed), "MDE": float(r.MDE_beta_80pct),
                   "verdict": r.verdict} for _, r in pw.iterrows()],
        "source": "mr_power.csv (identical to manuscript Table 4)",
    }
    save(fig, "Figure8_MR_power")


# ======================================================================
# Figure 9 — 样本重叠诊断
# ======================================================================
def fig9():
    ov = D["overlap"].copy()
    ov["exposure"] = ov.pair.str.split("__").str[0]
    ov["outcome"] = ov.pair.str.split("__").str[1].str.replace(".csv", "", regex=False)
    ov["pair_lab"] = [plab(e, o) for e, o in zip(ov.exposure, ov.outcome)]
    zs = {}
    for _, r in D["mr"].iterrows():
        zs[(r.exposure, r.outcome)] = float(r.z)
    ov["z"] = [zs[(e, o)] for e, o in zip(ov.exposure, ov.outcome)]
    ov["K"] = [int(mr_row(e, o).K) for e, o in zip(ov.exposure, ov.outcome)]
    ov = ov.sort_values("rho_GE_zcorr")

    fig = plt.figure(figsize=(7.2, 3.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1.0], wspace=0.48)
    axa, axb = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])

    y = np.arange(len(ov))
    cols = [C["vermillion"] if v < 0 else C["blue"] for v in ov.rho_GE_zcorr]
    axa.barh(y, ov.rho_GE_zcorr.values, color=cols, edgecolor="white", height=0.68)
    for k, (_, r) in enumerate(ov.iterrows()):
        if r.rho_GE_zcorr >= 0:
            axa.text(r.rho_GE_zcorr + 0.03, k, f"{r.rho_GE_zcorr:+.2f} (K={int(r.K)})",
                     va="center", ha="left", fontsize=6.0)
        else:
            axa.text(r.rho_GE_zcorr + 0.03, k, f"{r.rho_GE_zcorr:+.2f} (K={int(r.K)})",
                     va="center", ha="left", fontsize=6.0, color="white", fontweight="bold")
    axa.set_yticks(y)
    axa.set_yticklabels(ov.pair_lab.values, fontsize=6.4)
    axa.axvline(0, color=C["ink"], lw=0.8)
    axa.set_xlim(-0.95, 0.72)
    axa.set_xlabel("ρ_GE  (z-statistic correlation across instruments)", fontsize=7.5)
    axa.set_title("Sample-overlap proxy", fontsize=8.5, pad=6)
    panel_label(axa, "a", x=-0.42)

    axb.fill_between([-1, 0], 0, 20, color=C["faint"], alpha=0.35, zorder=0)
    axb.fill_between([0, 1], -20, 0, color=C["faint"], alpha=0.35, zorder=0)
    for _, r in ov.iterrows():
        col = C["blue"] if r.z > 0 else C["vermillion"]
        axb.scatter([r.rho_GE_zcorr], [r.z], s=42, facecolor=col, edgecolor="white",
                    linewidths=0.6, zorder=3)
    axb.axhline(0, color=C["ink"], lw=0.8)
    axb.axvline(0, color=C["ink"], lw=0.8)
    axb.set_xlim(-1, 1); axb.set_ylim(-17, 17)
    axb.set_xlabel("ρ_GE", fontsize=7.5)
    axb.set_ylabel("Signed IVW z-statistic", fontsize=7.5)
    n_conc = int(np.sum(np.sign(ov.rho_GE_zcorr) == np.sign(ov.z)))
    axb.set_title(f"Sign concordance: {n_conc}/10 pairs", fontsize=8.5, pad=6)
    panel_label(axb, "b", x=-0.22)

    PLOTTED["Figure9"] = {
        "pairs": [{"pair": r.pair_lab, "rho_GE": float(r.rho_GE_zcorr),
                   "rho_beta": float(r.rho_beta), "K": int(r.K),
                   "signed_ivw_z": float(r.z)} for _, r in ov.iterrows()],
        "n_sign_concordant": n_conc,
        "source": "overlap_statistics.csv + AUDIT_mr_recomputed.csv",
    }
    save(fig, "Figure9_overlap_pleiotropy")


# ======================================================================
# Figure 10（新增）— MRlap 样本重叠校正  §2.8
# ======================================================================
def fig10():
    m = D["mrlap"].copy()
    m["pair"] = [plab(e, o) for e, o in zip(m.exposure, m.outcome)]
    m["pct"] = (m.corrected_effect.abs() - m.observed_effect.abs()) / m.observed_effect.abs() * 100.0
    m = m.sort_values("pct")

    # ---- (d) 用：MRlap 自身 h² vs 主 LDSC h²（7 个唯一表型）
    h2_primary = {}
    for _, r in D["h2"].iterrows():
        h2_primary[r.phenotype] = float(r.h2)
    ns = D["nsens"]
    h2_primary["fracture"] = float(
        ns[(ns.trait == "fracture") & (ns.source.str.contains("264973"))].iloc[0].h2_obs)
    mr_h2 = {}
    for _, r in m.iterrows():
        mr_h2[r.exposure] = float(r.h2_exp)
        mr_h2[r.outcome] = float(r.h2_out)
    traits = sorted(mr_h2.keys())
    xs = [h2_primary[t] for t in traits]
    ys = [mr_h2[t] for t in traits]
    r_pear = float(np.corrcoef(xs, ys)[0, 1])

    fig = plt.figure(figsize=(7.2, 6.8))
    gs = fig.add_gridspec(2, 2, hspace=0.52, wspace=0.52)
    axa, axb = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    axc, axd = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])

    # (a) 相对变化
    y = np.arange(len(m))
    cols = [C["vermillion"] if v < 5 else C["blue"] for v in m.pct]
    axa.barh(y, m.pct.values, color=cols, edgecolor="white", height=0.68)
    for k, (_, r) in enumerate(m.iterrows()):
        if r.pct > 8:
            axa.text(r.pct - 0.6, k, f"{r.pct:.1f}%", va="center", ha="right",
                     fontsize=6.0, color="white", fontweight="bold")
        else:
            axa.text(r.pct + 0.5, k, f"{r.pct:.1f}%", va="center", ha="left", fontsize=6.0)
    med = float(np.median(m.pct))
    axa.axvline(med, color=C["grey"], lw=1.1, ls=(0, (5, 2.5)))
    axa.text(med + 0.5, 0.05, f"median {med:.1f}%", fontsize=6.2, color=C["grey"])
    axa.set_yticks(y)
    axa.set_yticklabels(m.pair.values, fontsize=6.2)
    axa.set_xlabel("Change in |β| after MRlap correction (%)", fontsize=7.5)
    axa.set_title("Change in |β| after correction", fontsize=8.5, pad=6)
    panel_label(axa, "a", x=-0.44)

    # (b) cross-trait LDSC intercept
    yb = np.arange(len(m))
    mb = m.sort_values("int_crosstrait")
    cols = [C["blue"] if v < 0 else C["vermillion"] for v in mb.int_crosstrait]
    axb.barh(yb, mb.int_crosstrait.values, color=cols, edgecolor="white", height=0.68)
    axb.set_yticks(yb)
    axb.set_yticklabels(mb.pair.values, fontsize=6.2)
    axb.axvline(0, color=C["ink"], lw=0.8)
    axb.set_xlabel("Cross-trait LDSC intercept (λ)", fontsize=7.5)
    axb.set_title("Estimated sample overlap (|λ| < 0.04)", fontsize=8.5, pad=6)
    panel_label(axb, "b", x=-0.44)

    # (c) MRlap difference test
    mc = m.sort_values("p_difference")
    yc = np.arange(len(mc))
    nlp = -np.log10(mc.p_difference.values.clip(1e-300))
    cols = [C["blue"] if v > -np.log10(0.05) else C["grey"] for v in nlp]
    axc.barh(yc, nlp, color=cols, edgecolor="white", height=0.68)
    axc.axvline(-np.log10(0.05), color=C["vermillion"], lw=1.1, ls=(0, (5, 2.5)))
    axc.text(-np.log10(0.05) + 0.35, len(mc) - 0.6, "α = 0.05", fontsize=6.2, color=C["vermillion"])
    axc.set_yticks(yc)
    axc.set_yticklabels(mc.pair.values, fontsize=6.2)
    axc.set_xlabel("−log10(p) for observed vs corrected difference", fontsize=7.5)
    n_sig = int(np.sum(mc.p_difference < 0.05))
    axc.set_title(f"Difference test: {n_sig}/10 significant", fontsize=8.5, pad=6)
    panel_label(axc, "c", x=-0.44)

    # (d) h² 一致性
    axd.scatter(xs, ys, s=46, facecolor=C["blue"], edgecolor="white", linewidths=0.6, zorder=3)
    lim = [0, max(max(xs), max(ys)) * 1.12]
    axd.plot(lim, lim, color=C["grey"], lw=1.0, ls=(0, (5, 2.5)), zorder=2)
    # arms 与 legs 的两组 h² 几乎重合，合并为一个标签避免文字重叠
    offs = {"BMD_total": (7, -2), "BMD_lumbar": (-62, -12), "BMD_femur_neck": (6, 4),
            "BMD_leg": None, "BMD_arms": (-62, -4), "OA_knee": (-50, 10),
            "fracture": (7, -10)}
    for t, x, yv in zip(traits, xs, ys):
        if t == "BMD_leg":
            continue
        lab = "arms / legs" if t == "BMD_arms" else PHEN[t].replace("BMD ", "")
        axd.annotate(lab, (x, yv), fontsize=5.8,
                     xytext=offs.get(t, (4, 4)), textcoords="offset points")
    axd.set_xlim(lim); axd.set_ylim(lim)
    axd.set_xlabel("h² — primary LDSC (Table 1)", fontsize=7.5)
    axd.set_ylabel("h² — MRlap (independent run)", fontsize=7.5)
    axd.set_title(f"Heritability concordance, r = {r_pear:.3f}", fontsize=8.5, pad=6)
    panel_label(axd, "d", x=-0.20)

    fig.suptitle("MRlap sample-overlap and winner's-curse correction (§2.8)", fontsize=9.5, y=0.995)

    PLOTTED["Figure10"] = {
        "pairs": [{"pair": r.pair, "observed": float(r.observed_effect),
                   "corrected": float(r.corrected_effect),
                   "pct_change": float(r.pct),
                   "p_difference": float(r.p_difference),
                   "int_crosstrait": float(r.int_crosstrait),
                   "rg": float(r.rg), "m_IVs": int(r.m_IVs),
                   "observed_p": float(r.observed_effect_p),
                   "corrected_p": float(r.corrected_effect_p)} for _, r in m.iterrows()],
        "median_pct_change": med,
        "n_difference_significant": n_sig,
        "h2_concordance": {"traits": traits,
                           "primary_LDSC": [float(v) for v in xs],
                           "MRlap": [float(v) for v in ys],
                           "pearson_r": r_pear},
        "source": "mrlap_overlap_results.csv + rg_h2_trueN.csv + N_sensitivity_decision.csv",
        "note": "MRlap input files use per-trait representative N, so its absolute β scale differs "
                "from Table 3; only the within-pair observed-vs-corrected contrast is interpreted.",
    }
    save(fig, "Figure10_MRlap_overlap_correction")


# ======================================================================
# SupplFig S2 — 留一法稳定性
# ======================================================================
def suppl2():
    loo = D["loo"]
    recs = []
    order = []
    for (e, o), g in loo.groupby(["exposure", "outcome"], sort=False):
        r = mr_row(e, o)
        recs.append((e, o, g, float(r.beta), float(r.se)))
    recs = sorted(recs, key=lambda t: t[3])

    n = len(recs)
    fig, axs = plt.subplots(2, 5, figsize=(7.2, 4.6))
    fig.subplots_adjust(hspace=0.62, wspace=0.30)
    out = []
    for k, (e, o, g, b, se) in enumerate(recs):
        ax = axs[k // 5, k % 5]
        y = np.arange(len(g))
        ax.fill_betweenx([-0.6, len(g) - 0.4], g.loo_beta.min(), g.loo_beta.max(),
                         color=C["sky"], alpha=0.35, zorder=1)
        ax.scatter(g.loo_beta.values, y, s=9, facecolor=C["sky"], edgecolor="none", zorder=2)
        ax.errorbar([b], [len(g) / 2], xerr=1.96 * se, fmt="o", color="black",
                    ecolor="black", elinewidth=1.3, capsize=3, markersize=5, zorder=4)
        ax.axvline(0, color=C["vermillion"], lw=0.9, ls=(0, (4, 2)), zorder=3)
        ax.set_yticks([])
        ax.set_ylim(-0.6, len(g) - 0.4)
        ax.tick_params(axis="x", labelsize=5.4)
        ax.set_title(plab(e, o), fontsize=5.6, pad=3)
        ax.text(0.03, 0.97, chr(97 + k), transform=ax.transAxes, fontsize=8,
                fontweight="bold", va="top", ha="left",
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
        out.append(dict(pair=f"{e}->{o}", full_beta=float(b), se=float(se),
                        loo_min=float(g.loo_beta.min()), loo_max=float(g.loo_beta.max()),
                        n_loo=int(len(g)), null_in_range=bool(g.loo_beta.min() <= 0 <= g.loo_beta.max())))
    fig.text(0.5, -0.015,
             "IVW estimate: black point = all instruments (± 95% CI); shaded band = min–max "
             "range over leave-one-out estimates; dashed red line = null.",
             ha="center", fontsize=6.6)
    fig.suptitle("Leave-one-out stability of the IVW estimate (all ten pairs)",
                 fontsize=9.5, y=1.03)
    PLOTTED["SupplFig2"] = {"panels": out, "source": "harmonised/_loo.csv + AUDIT_mr_recomputed.csv"}
    save(fig, "SupplFig2_leave_one_out")


# ======================================================================
# SupplFig S3 — GEFOS 外部验证
# ======================================================================
def suppl3():
    g = D["gefos"]
    labs = [l.replace("This study", "This study") for l in g["labels"]]
    y = np.arange(len(labs))[::-1]
    est = np.array(g["estimates"]); lo = np.array(g["ci_low"]); hi = np.array(g["ci_high"])
    fig, ax = plt.subplots(figsize=(4.6, 2.9))
    cols = [C["grey"]] + [C["blue"]] * (len(labs) - 1)
    ax.errorbar(est, y, xerr=[est - lo, hi - est], fmt="o", color=C["ink"],
                ecolor=C["ink"], elinewidth=1.0, capsize=3, markersize=5)
    for k, c in enumerate(cols):
        ax.scatter([est[k]], [y[k]], s=42, facecolor=c, edgecolor="white",
                   linewidths=0.6, zorder=3)
    for k, v in enumerate(est):
        ax.text(hi[k] + 0.015, y[k], f"{v:.2f}", va="center", fontsize=6.2)
    ax.set_yticks(y)
    ax.set_yticklabels(labs, fontsize=6.4)
    ax.axvline(1.0, color=C["vermillion"], lw=1.0, ls=(0, (4, 2)))
    ax.set_xlabel("Odds ratio for any fracture per 1-SD lower BMD", fontsize=7.5)
    ax.set_xlim(1.28, 1.72)
    ax.set_title("External validation against GEFOS 2018", fontsize=8.5, pad=6)
    PLOTTED["SupplFig3"] = dict(g)
    save(fig, "SupplFig3_GEFOS_benchmark")


# ======================================================================
# SupplFig S4 — STROBE-MR 流程
# ======================================================================
def suppl4():
    fig, ax = plt.subplots(figsize=(7.2, 9.2))
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
    stages = [
        ("1   GWAS inputs", C["blue"],
         "13 skeletal phenotypes — OA knee (N = 1,208,397), OA self-reported (462,933),\n"
         "any fracture (264,973), ten UKB DXA BMD sites (31,873)"),
        ("2   Data integrity & QC", C["purple"],
         "LDSC h² and intercept (all ≈ 1.0); allele-orientation checks via sentinel loci;\n"
         "OA_selfrep global Z-sign correction (yielding r_g = +0.84 with knee OA)"),
        ("3   Harmonisation & clumping", C["sky"],
         "Strict PLINK-equivalent clumping (r² < 0.001 over 10 Mb, 1000G-EUR reference);\n"
         "1-Mb window retained as a sensitivity analysis; F = 44–75"),
        ("4   Two-sample MR", C["orange"],
         "Ten exposure–outcome pairs; IVW (primary), MR-Egger, weighted median;\n"
         "BMD effect scale verified at N=31,873 (per-SNP SE check); no rescaling applied"),
        ("5   Sensitivity & diagnostics", C["vermillion"],
         "MR-PRESSO outliers, leave-one-out, Cochran Q / I², sample-overlap ρ_GE,\n"
         "MRlap overlap + winner's-curse correction (§2.8), Bonferroni and BH-FDR"),
        ("6   Reporting", C["green"],
         "STROBE-MR checklist; effect scales stated explicitly; power and MDE reported;\n"
         "code and summary statistics publicly archived"),
    ]
    ytop, dy = 9.25, 1.55
    for i, (ti, col, desc) in enumerate(stages):
        yy = ytop - i * dy
        ax.add_patch(FancyBboxPatch((0.35, yy - 0.58), 9.3, 1.16,
                     boxstyle="round,pad=0.02,rounding_size=0.14",
                     fc="white", ec=col, lw=1.6, zorder=2))
        ax.add_patch(FancyBboxPatch((0.35, yy - 0.58), 0.45, 1.16,
                     boxstyle="round,pad=0.02,rounding_size=0.10", fc=col, ec="none", zorder=3))
        ax.text(0.60, yy + 0.28, ti, ha="left", va="center", fontsize=9,
                fontweight="bold", color=col, zorder=4)
        ax.text(0.60, yy - 0.20, desc, ha="left", va="center", fontsize=6.8,
                color=C["ink"], zorder=4)
        if i < len(stages) - 1:
            ax.add_patch(FancyArrowPatch((5.0, yy - 0.58), (5.0, yy - dy + 0.58),
                         arrowstyle="-|>", mutation_scale=13, color=C["grey"], lw=1.5, zorder=1))
    ax.set_title("Study design and MR analysis flow (STROBE-MR compliant)",
                 fontsize=10, fontweight="bold", pad=8)
    save(fig, "SupplFig4_STROBE_flow")


# ======================================================================
# SupplFig S5（新增）— OA_knee h² 对代表 N 的敏感性
# ======================================================================
def suppl5():
    Nref, h2ref, seref = 1208397.0, 0.0354, 0.0007
    Ns = np.array([10167., 50000., 100000., 250000., 605000., 1208397., 1316504.])
    h2 = h2ref * Nref / Ns
    se = seref * Nref / Ns

    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    ax.plot(Ns, h2, marker="o", ms=5, lw=1.5, color=C["blue"])
    ax.fill_between(Ns, h2 - 1.96 * se, h2 + 1.96 * se, color=C["blue"], alpha=0.18, lw=0)
    ax.scatter([Nref], [h2ref], s=90, facecolor=C["vermillion"], edgecolor="white",
               linewidths=0.8, zorder=5)
    ax.annotate(f"Reported\nh² = {h2ref:.4f}\nN = {int(Nref):,}", (Nref, h2ref),
                xytext=(-64, 26), textcoords="offset points", fontsize=6.4,
                color=C["vermillion"], fontweight="bold",
                arrowprops=dict(arrowstyle="->", color=C["vermillion"], lw=1.0))
    ax.axhline(1.0, color=C["grey"], lw=0.9, ls=(0, (4, 2)))
    ax.text(Ns[0], 1.05, "h² = 1 (biologically impossible upper bound)",
            fontsize=5.8, color=C["grey"])
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Representative N used for OA_knee", fontsize=7.5)
    ax.set_ylabel("h² (observed scale)", fontsize=7.5)
    ax.set_title("OA_knee h² scales as 1 / N", fontsize=8.5, pad=6)
    ax.text(0.03, 0.05, "Deterministic rescaling of the reported estimate (h² scales as 1/N);\n"
                        "no new LDSC run. r_g and all MR directions are N-invariant (§5.2).",
            transform=ax.transAxes, fontsize=5.8, color=C["ink"], va="bottom",
            bbox=dict(boxstyle="round,pad=0.30", fc="#F4F4F4", ec="none"))

    PLOTTED["SupplFig5"] = {
        "N": [float(v) for v in Ns], "h2": [float(v) for v in h2], "se": [float(v) for v in se],
        "reported_N": Nref, "reported_h2": h2ref,
        "rule": "h2(N') = h2_reported × N_reported / N'",
        "source": "Manuscript Table 1 (OA_knee) + Supplement S3",
    }
    save(fig, "SupplFig5_OA_knee_N_sensitivity")


# ======================================================================
# SupplFig S6（新增）— fracture N 选择的证据链
# ======================================================================
def suppl6():
    ns = D["nsens"]
    f = ns[ns.trait == "fracture"].copy()
    f["label"] = f.source.str.split("_").str[0]
    f = f.sort_values("N")
    adopted = 264973

    fig, ax = plt.subplots(figsize=(5.2, 3.2))
    ax.plot(f.N.values, f.h2_obs.values, marker="o", ms=5, lw=1.5, color=C["grey"])
    for _, r in f.iterrows():
        col = C["vermillion"] if int(r.N) == adopted else C["blue"]
        ms = 9 if int(r.N) == adopted else 5
        ax.scatter([r.N], [r.h2_obs], s=ms * 9, facecolor=col, edgecolor="white",
                   linewidths=0.6, zorder=4)
        ax.annotate(f"{int(r.N):,}", (r.N, r.h2_obs), xytext=(0, 9),
                    textcoords="offset points", ha="center", fontsize=6.0, color=col)
    ax.axhline(0.0258, color=C["vermillion"], lw=0.9, ls=(0, (4, 2)))
    ax.text(f.N.max(), 0.0262, "adopted h² = 0.0258", fontsize=6.2,
            color=C["vermillion"], ha="right")
    ax.set_xlabel("Candidate N for the fracture GWAS", fontsize=7.5)
    ax.set_ylabel("h² (observed scale)", fontsize=7.5)
    ax.set_title("Fracture h² across candidate N", fontsize=8.5, pad=6)
    ax.text(0.02, 0.04,
            "Adopted N = 264,973 (GEFOS discovery: 37,857 cases / 227,116 controls).\n"
            "Measured effective N·K(1−K) = 31,605 (IQR 31,319–32,325) vs GEFOS\n"
            "theoretical 32,449 — 2.6% difference, closer than any UKB scenario.\n"
            "All r_g and MR conclusions are N-invariant.",
            transform=ax.transAxes, fontsize=5.8, va="bottom",
            bbox=dict(boxstyle="round,pad=0.30", fc="#F4F4F4", ec="none"))

    PLOTTED["SupplFig6"] = {
        "candidates": [{"N": int(r.N), "h2": float(r.h2_obs), "se": float(r.h2_se),
                        "source": r.source} for _, r in f.iterrows()],
        "adopted_N": adopted,
        "measured_effective_NK": 31605, "IQR": [31319, 32325],
        "GEFOS_theoretical": 32449,
        "source": "N_sensitivity_decision.csv",
    }
    save(fig, "SupplFig6_fracture_N_decision")


# ======================================================================
# SupplFig S7（新增）— 多效性 / 异质性诊断
# ======================================================================
def suppl7():
    SHORT = {"BMD_total": "total", "BMD_femur_neck": "fem. neck", "BMD_lumbar": "lumbar",
             "BMD_arms": "arms", "BMD_leg": "legs", "OA_knee": "knee OA", "fracture": "fracture"}

    def short(e, o):
        return f"{SHORT[e]} → {SHORT[o]}"

    mr = D["mr"].copy()
    mr["pair"] = [short(e, o) for e, o in zip(mr.exposure, mr.outcome)]
    mr["n_out"] = [int(harm(e, o).presso_outlier.sum()) for e, o in zip(mr.exposure, mr.outcome)]
    mr = mr.sort_values("I2")

    fig = plt.figure(figsize=(7.2, 4.2))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1.05, 0.85], wspace=0.62)
    axa, axb, axc = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2])

    y = np.arange(len(mr))
    axa.barh(y, mr.I2.values * 100, color=C["sky"], edgecolor="white", height=0.68)
    axa.axvline(50, color=C["grey"], lw=0.9, ls=(0, (4, 2)))
    axa.text(52, len(mr) - 0.6, "I² = 50%", fontsize=6.0, color=C["grey"],
             bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
    axa.set_yticks(y); axa.set_yticklabels(mr.pair.values, fontsize=6.2)
    axa.set_xlabel("Cochran I² (%)", fontsize=7.5)
    axa.set_title("Heterogeneity", fontsize=8.5, pad=6)
    panel_label(axa, "a", x=-0.46)

    mb = mr.sort_values("egger_int_p")
    yb = np.arange(len(mb))
    nlp = -np.log10(mb.egger_int_p.values.clip(1e-300))
    cols = [C["vermillion"] if v > -np.log10(0.05) else C["blue"] for v in nlp]
    axb.barh(yb, nlp, color=cols, edgecolor="white", height=0.68)
    axb.axvline(-np.log10(0.05), color=C["vermillion"], lw=1.0, ls=(0, (5, 2.5)))
    axb.set_xlim(0, max(nlp.max() * 1.18, 2.0))
    axb.text(-np.log10(0.05) + 0.2, len(mb) - 0.8, "α = 0.05", fontsize=6.0, color=C["vermillion"])
    axb.set_yticks(yb); axb.set_yticklabels(mb.pair.values, fontsize=6.2)
    axb.set_xlabel("−log10(MR-Egger intercept P)", fontsize=7.5)
    n_pleio = int(np.sum(mb.egger_int_p < 0.05))
    axb.set_title(f"Directional pleiotropy: {n_pleio}/10", fontsize=8.5, pad=6)
    panel_label(axb, "b", x=-0.44)

    mc = mr.sort_values("n_out")
    yc = np.arange(len(mc))
    axc.barh(yc, mc.n_out.values, color=C["orange"], edgecolor="white", height=0.68)
    for k, v in enumerate(mc.n_out.values):
        if v > 0:
            axc.text(v - 0.10, k, str(int(v)), va="center", ha="right", fontsize=6.0,
                     color="white", fontweight="bold")
        else:
            axc.text(0.08, k, "0", va="center", ha="left", fontsize=6.0)
    axc.set_yticks(yc); axc.set_yticklabels(mc.pair.values, fontsize=6.2)
    axc.set_xlabel("MR-PRESSO outliers", fontsize=7.5)
    axc.set_xlim(0, max(mc.n_out) + 0.6)
    axc.set_xticks(range(0, int(max(mc.n_out)) + 1))
    axc.set_title("Outliers flagged", fontsize=8.5, pad=6)
    panel_label(axc, "c", x=-0.26)

    fig.suptitle("Pleiotropy and heterogeneity diagnostics across the ten MR pairs",
                 fontsize=9.5, y=1.0)
    PLOTTED["SupplFig7"] = {
        "pairs": [{"pair": r.pair, "I2": float(r.I2), "Q_p": float(r.Q_p),
                   "egger_intercept_p": float(r.egger_int_p),
                   "n_presso_outlier": int(r.n_out), "K": int(r.K)}
                  for _, r in mr.iterrows()],
        "n_egger_intercept_significant": n_pleio,
        "source": "AUDIT_mr_recomputed.csv + harmonised/*.csv (presso_outlier)",
    }
    save(fig, "SupplFig7_pleiotropy_diagnostics")


# ======================================================================
# Graphical Abstract
# ======================================================================
def graphical_abstract():
    fig, ax = plt.subplots(figsize=(4.4, 6.2))
    ax.set_xlim(0, 10); ax.set_ylim(0, 10.4); ax.axis("off")

    def node(x, y, w, h, t, fc, fs=9):
        ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                     boxstyle="round,pad=0.03,rounding_size=0.22", fc=fc, ec="none", zorder=2))
        ax.text(x, y, t, ha="center", va="center", fontsize=fs, fontweight="bold",
                color="white", zorder=3)

    def arrow(x1, y1, x2, y2, col, lab, lx, ly, rad=0.0, ls="-", lw=2.4):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                     mutation_scale=15, color=col, lw=lw, ls=ls, zorder=1,
                     connectionstyle=f"arc3,rad={rad}",
                     path_effects=[pe.withStroke(linewidth=3.4, foreground="white")]))
        ax.text(lx, ly, lab, ha="center", va="center", fontsize=7.2, fontweight="bold",
                color=col, zorder=5,
                bbox=dict(boxstyle="round,pad=0.26", fc="white", ec="none", alpha=0.94))

    node(2.2, 7.9, 2.9, 1.25, "OA (knee)", C["purple"])
    node(7.6, 6.5, 3.1, 1.55, "BMD\n(10 regional sites)", C["blue"])
    node(7.6, 1.9, 2.9, 1.25, "Fracture", C["sky"])

    arrow(6.95, 6.0, 6.95, 2.5, C["blue"], "Protective\nOR 1.40–1.42", 8.95, 4.2, rad=-0.12)
    arrow(6.6, 6.2, 3.6, 7.0, C["vermillion"], "Causal ↑ risk\nβ +0.058 to +0.095", 4.95, 8.05, rad=0.16)
    ax.add_patch(FancyArrowPatch((3.4, 5.4), (6.4, 6.1), arrowstyle="-|>", mutation_scale=11,
                 color=C["grey"], lw=1.6, ls=(0, (3.5, 2.2)), zorder=1,
                 connectionstyle="arc3,rad=-0.22",
                 path_effects=[pe.withStroke(linewidth=3.0, foreground="white")]))
    ax.text(5.0, 4.75, "Genetic correlation\nr_g +0.12 to +0.25", ha="center", va="center",
            fontsize=7.0, color=C["grey"], zorder=4,
            bbox=dict(boxstyle="round,pad=0.26", fc="white", ec="none", alpha=0.94))
    ax.add_patch(FancyArrowPatch((3.6, 5.4), (6.6, 2.5), arrowstyle="-|>", mutation_scale=10,
                 color=C["grey"], lw=1.4, ls=(0, (3.5, 2.2)), zorder=1,
                 connectionstyle="arc3,rad=0.28",
                 path_effects=[pe.withStroke(linewidth=3.0, foreground="white")]))
    ax.text(5.5, 3.45, "No well-powered\ncausal pathway", ha="center", va="center",
            fontsize=7.0, color=C["grey"], fontweight="bold", zorder=4,
            bbox=dict(boxstyle="round,pad=0.26", fc="white", ec="none", alpha=0.94))

    ax.text(5.0, 9.85, "Shared genetic architecture,\ndiverging causal fates",
            ha="center", va="center", fontsize=10, fontweight="bold")
    ax.text(5.0, 0.35, "Higher BMD protects against fracture yet raises knee-OA risk;\n"
                       "OA shows no well-powered causal pathway to fracture.",
            ha="center", va="center", fontsize=6.8, color=C["ink"],
            bbox=dict(boxstyle="round,pad=0.38", fc="#F2F2F2", ec="none"))
    save(fig, "GraphicalAbstract")


# ====================================================================== 主流程
def main():
    print("▶ Building Nature/Science-standard figure suite →", OUT)
    fig1(); fig2(); fig3(); fig4(); fig5()
    fig6(); fig7(); fig8(); fig9(); fig10()
    suppl1(); suppl2(); suppl3(); suppl4(); suppl5(); suppl6(); suppl7()
    graphical_abstract()

    with open(os.path.join(OUT, "plotted_data.json"), "w", encoding="utf-8") as f:
        json.dump(PLOTTED, f, indent=1, ensure_ascii=False)
    print("   OK plotted_data.json")

    # 保留原图（只读复制，不改动 Figures_Beautified/）
    src = os.path.join(PROJECT, "Figures_Beautified")
    if os.path.isdir(src) and not os.path.isdir(ORIG):
        shutil.copytree(src, ORIG)
        print("   OK archived originals →", ORIG)

    print(f"\n✔ {len(SAVED)} figures written to {OUT}")


if __name__ == "__main__":
    main()
