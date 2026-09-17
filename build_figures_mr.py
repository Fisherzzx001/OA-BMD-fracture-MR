# -*- coding: utf-8 -*-
"""
build_figures_mr.py
=========================================================================
基于 **重跑 harmonisation 导出的逐 SNP 真实数据** 绘制 MR 诊断图。

数据来源：code/output/genetics/harmonised/<exposure>__<outcome>.csv
  （由 26_export_harmonised_snp.py 从原始 GWAS 汇总统计重建，
    已通过与 mr_corrected_N.csv 的逐位一致性校验，Δβ ~ 1e-17）

产出（Nature 风格：Arial、色盲安全配色、去顶/右边框、面板小写粗体标签、
      300 dpi PNG + 可编辑矢量 PDF）
  Figure6_MR_scatter.png/.pdf         2×2 散点图（4 个核心对）
  Figure7_MR_funnel.png/.pdf          2×2 MR-PRESSO 漏斗图（4 个核心对）
  SupplFig1_MR_scatter_all.png/.pdf   2×3 散点图（其余 6 对）
  SupplFig2_leave_one_out.png/.pdf    留一法稳定性（全部 10 对）
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import norm

HARM = "C:/Users/Administrator/WorkBuddy/2026-09-08-22-17-40/OA_B_PROJECT/code/output/genetics/harmonised"
OUT = "C:/Users/Administrator/WorkBuddy/2026-09-08-22-17-40/OA_B_PROJECT/submission_package/Figures"

# ---------- Nature style（与 build_figures.py 保持一致）----------
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
CB = {"blue": "#4878CF", "red": "#D65F5F", "orange": "#EE854A",
      "purple": "#956CB4", "grey": "#8C8C8C", "cb_blue": "#0072B2",
      "ink": "#222222"}

PRETTY = {
    "OA_knee": "knee OA", "OA_selfrep": "self-reported OA", "fracture": "fracture",
    "BMD_total": "total BMD", "BMD_femur_neck": "femoral-neck BMD",
    "BMD_femur_total": "femoral-total BMD", "BMD_lumbar": "lumbar BMD",
    "BMD_spine": "spine BMD", "BMD_arms": "arms BMD", "BMD_leg": "leg BMD",
    "BMD_pelvis": "pelvis BMD", "BMD_rib": "rib BMD", "BMD_trunk": "trunk BMD",
}

HEADLINE = [("BMD_total", "fracture"),
            ("BMD_femur_neck", "fracture"),
            ("BMD_total", "OA_knee"),
            ("OA_knee", "BMD_total")]
REST = [("BMD_femur_neck", "OA_knee"), ("BMD_lumbar", "OA_knee"),
        ("BMD_arms", "OA_knee"), ("BMD_leg", "OA_knee"),
        ("OA_knee", "fracture"), ("OA_knee", "BMD_femur_neck")]
ALLPAIRS = HEADLINE + REST


# ------------------------------------------------------------------ 统计
def ivw(be, bo, so):
    w = 1.0 / (so ** 2)
    return float(np.sum(w * be * bo) / np.sum(w * be ** 2))


def ivw_se(be, so):
    return float(1.0 / np.sqrt(np.sum(be ** 2 / so ** 2)))


def egger(be, bo, so):
    w = 1.0 / (so ** 2)
    mx = np.sum(w * be) / np.sum(w)
    my = np.sum(w * bo) / np.sum(w)
    b = np.sum(w * (be - mx) * (bo - my)) / np.sum(w * (be - mx) ** 2)
    return float(my - b * mx), float(b)


def wmed(be, bo, so):
    r = bo / be
    w = (be / so) ** 2
    o = np.argsort(r)
    r, w = r[o], w[o]
    cw = np.cumsum(w) / np.sum(w)
    i = min(int(np.searchsorted(cw, 0.5)), len(r) - 1)
    if i == 0:
        return float(r[0])
    if abs(cw[i - 1] - 0.5) < 1e-12:
        return float(0.5 * (r[i - 1] + r[i]))
    return float(r[i])


def load(exp, out):
    d = pd.read_csv(os.path.join(HARM, f"{exp}__{out}.csv"))
    d = d[np.isfinite(d.beta_e) & np.isfinite(d.beta_o)]
    d = d[(np.abs(d.beta_e) > 0) & (d.se_o > 0)]
    return d.reset_index(drop=True)


def stats(d):
    be, se_e = d.beta_e.values, d.se_e.values
    bo, so = d.beta_o.values, d.se_o.values
    b = ivw(be, bo, so)
    s = ivw_se(be, so)
    z = b / s
    p = float(2 * norm.sf(abs(z)))
    a, be2 = egger(be, bo, so)
    return dict(b=b, se=s, p=p, egger=be2, egger_int=a, wm=wmed(be, bo, so), K=len(d))


def pfmt(p):
    if p < 1e-300:
        return "<1e-300"
    if p < 1e-4:
        return f"{p:.1e}"
    return f"{p:.3f}"


def xlab(t):
    return (f"SNP effect on {PRETTY[t]} (SD per allele)"
            if t.startswith("BMD") else f"SNP effect on {PRETTY[t]} (per allele)")


def ylab(t):
    if t == "fracture":
        return "SNP effect on fracture (log-OR per allele)"
    if t.startswith("BMD"):
        return f"SNP effect on {PRETTY[t]} (SD per allele)"
    return f"SNP effect on {PRETTY[t]} (per allele)"


def panel_label(ax, txt, x=-0.14, y=1.06):
    ax.text(x, y, txt, transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="top", ha="left")


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(os.path.join(OUT, name + ".png"))
    fig.savefig(os.path.join(OUT, name + ".pdf"))
    plt.close(fig)
    print("  ✔", name)


# ------------------------------------------------------------------ 散点
def scatter(ax, exp, out, lab, title=None, show_ylab=True, show_xlab=True):
    d = load(exp, out)
    st = stats(d)
    be, se_e = d.beta_e.values, d.se_e.values
    bo, so = d.beta_o.values, d.se_o.values
    outl = d.presso_outlier.values.astype(bool)

    ax.errorbar(be, bo, xerr=1.96 * se_e, yerr=1.96 * so, fmt="none",
                ecolor="#C9C9C9", elinewidth=0.6, alpha=0.85, zorder=1)
    ax.scatter(be[~outl], bo[~outl], s=16, facecolor=CB["cb_blue"],
               edgecolors="white", linewidths=0.4, zorder=3)
    if outl.any():
        ax.scatter(be[outl], bo[outl], s=34, facecolor="none",
                   edgecolors=CB["red"], linewidths=1.1, zorder=4,
                   label="MR-PRESSO outlier")

    xs = np.linspace(min(0, be.min() * 1.1), max(0, be.max() * 1.1), 100)
    ax.plot(xs, st["b"] * xs, "-", color=CB["ink"], lw=1.4, zorder=5,
            label=f"IVW  β={st['b']:+.3f} (P={pfmt(st['p'])})")
    ax.plot(xs, st["egger_int"] + st["egger"] * xs, "--", color=CB["orange"],
            lw=1.1, zorder=5, label=f"MR-Egger  β={st['egger']:+.3f}")
    ax.plot(xs, st["wm"] * xs, ":", color=CB["purple"], lw=1.2, zorder=5,
            label=f"Weighted median  β={st['wm']:+.3f}")

    ax.axhline(0, color="#D8D8D8", lw=0.7, zorder=0)
    ax.axvline(0, color="#D8D8D8", lw=0.7, zorder=0)
    if title:
        ax.set_title(title, fontsize=8.5, fontweight="bold", loc="left", pad=3)
    if show_xlab:
        ax.set_xlabel(xlab(exp), fontsize=7)
    if show_ylab:
        ax.set_ylabel(ylab(out), fontsize=7)
    ax.tick_params(labelsize=6.5)
    ax.legend(fontsize=5.6, frameon=True, framealpha=0.92,
              edgecolor="#DDDDDD", loc="best", handlelength=1.4,
              borderpad=0.4, labelspacing=0.3)
    panel_label(ax, lab)
    return st


# ------------------------------------------------------------------ 漏斗
def funnel(ax, exp, out, lab, title=None, show_ylab=True, show_xlab=True):
    d = load(exp, out)
    st = stats(d)
    r = d.ratio.values
    prec = 1.0 / d.se_ratio.values
    outl = d.presso_outlier.values.astype(bool)

    # 稳健坐标范围（弱工具变量的比值可能极端，按分位数截断以免整图失真）
    lo, hi = np.nanpercentile(r, [2, 98])
    span = max(hi - lo, 4 * abs(st["b"]) + 0.05)
    c = 0.5 * (lo + hi)
    xlim = (min(c - 0.75 * span, st["b"] - 0.3 * span),
            max(c + 0.75 * span, st["b"] + 0.3 * span))

    ax.scatter(r[~outl], prec[~outl], s=16, facecolor=CB["cb_blue"],
               edgecolors="white", linewidths=0.4, zorder=3)
    if outl.any():
        ax.scatter(r[outl], prec[outl], s=38, facecolor="none",
                   edgecolors=CB["red"], linewidths=1.1, zorder=4,
                   label="MR-PRESSO outlier")

    yy = np.linspace(max(prec.min() * 0.85, 1e-6), prec.max() * 1.08, 300)
    ax.plot(st["b"] + 1.96 / yy, yy, ls="--", color="#9E9E9E", lw=0.8, zorder=1)
    ax.plot(st["b"] - 1.96 / yy, yy, ls="--", color="#9E9E9E", lw=0.8, zorder=1,
            label="95% pseudo-CI")
    ax.axvline(st["b"], color=CB["ink"], lw=1.3, zorder=2,
               label=f"IVW  β={st['b']:+.3f}")
    ax.axvline(st["egger"], color=CB["orange"], ls="--", lw=1.1, zorder=2,
               label=f"MR-Egger  β={st['egger']:+.3f}")

    ax.set_xlim(*xlim)
    ax.set_ylim(max(prec.min() * 0.85, 0), prec.max() * 1.08)
    if title:
        ax.set_title(title, fontsize=8.5, fontweight="bold", loc="left", pad=3)
    if show_xlab:
        ax.set_xlabel("SNP-specific causal estimate (Wald ratio)", fontsize=7)
    if show_ylab:
        ax.set_ylabel("1 / SE (ratio)", fontsize=7)
    ax.tick_params(labelsize=6.5)
    ax.legend(fontsize=5.6, frameon=True, framealpha=0.92,
              edgecolor="#DDDDDD", loc="upper right", handlelength=1.4,
              borderpad=0.4, labelspacing=0.3)
    panel_label(ax, lab)
    return st


# ------------------------------------------------------------------ 图 6 / 7
def fig_scatter(fname, pairs, shape, size, tag="Figure"):
    fig, axes = plt.subplots(shape[0], shape[1], figsize=size)
    axes = np.atleast_1d(axes).ravel()
    labs = "abcdefghij"
    for k, (exp, out) in enumerate(pairs):
        ax = axes[k]
        title = f"{PRETTY[exp]} → {PRETTY[out]}  (K={stats(load(exp, out))['K']})"
        scatter(ax, exp, out, labs[k], title=title,
                show_xlab=(k >= len(pairs) - shape[1]),
                show_ylab=(k % shape[1] == 0))
    for j in range(len(pairs), len(axes)):
        axes[j].axis("off")
    fig.tight_layout()
    save(fig, fname)


def fig_funnel(fname, pairs, shape, size):
    fig, axes = plt.subplots(shape[0], shape[1], figsize=size)
    axes = np.atleast_1d(axes).ravel()
    labs = "abcdefghij"
    for k, (exp, out) in enumerate(pairs):
        ax = axes[k]
        title = f"{PRETTY[exp]} → {PRETTY[out]}  (K={stats(load(exp, out))['K']})"
        funnel(ax, exp, out, labs[k], title=title,
               show_xlab=(k >= len(pairs) - shape[1]),
               show_ylab=(k % shape[1] == 0))
    for j in range(len(pairs), len(axes)):
        axes[j].axis("off")
    fig.tight_layout()
    save(fig, fname)


# ------------------------------------------------------------------ 留一法
def fig_loo(fname):
    rows = []
    for exp, out in ALLPAIRS:
        d = load(exp, out)
        st = stats(d)
        be, bo, so = d.beta_e.values, d.beta_o.values, d.se_o.values
        loo = []
        for i in range(len(d)):
            m = np.ones(len(d), bool)
            m[i] = False
            loo.append(ivw(be[m], bo[m], so[m]))
        loo = np.array(loo)
        rows.append(dict(exp=exp, out=out, b=st["b"], se=st["se"], p=st["p"],
                         loo=loo, flips=int(np.sum(np.sign(loo) != np.sign(st["b"])))))

    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    y = np.arange(len(rows))[::-1]
    for yy, r in zip(y, rows):
        lo, hi = r["loo"].min(), r["loo"].max()
        ax.plot([lo, hi], [yy, yy], "-", color=CB["blue"], lw=2.6,
                alpha=0.35, solid_capstyle="round", zorder=1)
        ax.scatter(r["loo"], np.full(len(r["loo"]), yy), s=5,
                   facecolor=CB["blue"], edgecolors="none", alpha=0.55, zorder=2)
        ax.errorbar(r["b"], yy, xerr=1.96 * r["se"], fmt="o",
                    color=CB["ink"], ecolor=CB["ink"], markersize=4.5,
                    elinewidth=1.3, capsize=3, zorder=4)
    ax.axvline(0, color=CB["grey"], ls="--", lw=0.9, zorder=0)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{PRETTY[r['exp']]} → {PRETTY[r['out']]}"
                        for r in rows], fontsize=6.5)
    ax.set_xlabel("MR effect (β, N-corrected scale)", fontsize=7.5)
    ax.tick_params(labelsize=6.5)
    ax.set_title("Leave-one-out stability of the IVW estimate",
                 fontsize=9, fontweight="bold", loc="left")

    handles = [
        Line2D([0], [0], marker="o", color=CB["ink"], markerfacecolor=CB["ink"],
               markersize=4.5, lw=1.3, label="IVW estimate (all SNPs), 95% CI"),
        Line2D([0], [0], color=CB["blue"], lw=2.6, alpha=0.5,
               label="Range of leave-one-out estimates"),
        Line2D([0], [0], marker="o", color=CB["blue"], lw=0, markersize=3,
               alpha=0.6, label="Individual leave-one-out estimate"),
    ]
    ax.legend(handles=handles, fontsize=6, loc="lower right", frameon=True,
              framealpha=0.95, edgecolor="#DDDDDD")
    fig.tight_layout()
    save(fig, fname)


def main():
    print("▶ 基于逐 SNP 真实数据绘制 MR 诊断图")
    fig_scatter("Figure6_MR_scatter", HEADLINE, (2, 2), (9.2, 7.6))
    fig_funnel("Figure7_MR_funnel", HEADLINE, (2, 2), (9.2, 7.6))
    fig_scatter("SupplFig1_MR_scatter_remaining", REST, (2, 3), (12.4, 7.4))
    fig_loo("SupplFig2_leave_one_out")
    print("ALL MR FIGURES DONE")


if __name__ == "__main__":
    main()
