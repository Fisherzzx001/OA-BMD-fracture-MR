# -*- coding: utf-8 -*-
"""
build_figures_v9.py
=========================================================================
在 v6/v8 已有的 Figure 1–7 / S1–S2 之外，**依据论文真实统计结果**补充三张图，
用于加强稿件在三处审稿人关切点上的表达：

  Figure8_MR_power.png/.pdf
      (a) MR 设计的检验效能函数 + 10 对观测点
      (b) 逐对实测效能条形图（80% 参考线）
      -> 支撑 §3.8 / Table 4；把 "无证据 ≠ 证据显示无效应" 可视化

  Figure9_overlap_pleiotropy.png/.pdf
      (a) 逐对 SNP–暴露 与 SNP–结局 z 统计量相关 ρ_GE（重叠/多效性量规）
      (b) ρ_GE 符号 vs IVW z 符号的一致性散点
      -> 支撑 §2.8（MR 样本重叠是预审报告的阻塞项之一）

  SupplFig3_GEFOS_benchmark.png/.pdf
      外部效度森林图：本研究 OR 与 GEFOS 2018 独立基准的对照
      -> 支撑 §3.7

风格：Nature 系列规范（Arial、≥7 pt、色盲安全配色、去顶/右边框、
小写粗体面板标签、300 dpi PNG + 可编辑矢量 PDF），与 build_figures.py /
build_figures_mr.py 完全一致。
"""
import os
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import norm

GEN = "C:/Users/Administrator/WorkBuddy/2026-09-08-22-17-40/OA_B_PROJECT/code/output/genetics"
OUT = os.environ.get(
    "FIG_OUT",
    "C:/Users/Administrator/WorkBuddy/2026-09-08-22-17-40/OA_B_PROJECT/submission_package_v2/Figures",
)

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
CB = {"blue": "#4878CF", "red": "#D65F5F", "orange": "#EE854A",
      "purple": "#956CB4", "grey": "#8C8C8C", "cb_blue": "#0072B2",
      "cb_orange": "#E69F00", "ink": "#222222", "faint": "#D8D8D8"}

PRETTY = {
    "OA_knee": "knee OA", "OA_selfrep": "self-reported OA", "fracture": "fracture",
    "BMD_total": "total BMD", "BMD_femur_neck": "femoral-neck BMD",
    "BMD_lumbar": "lumbar BMD", "BMD_arms": "arms BMD", "BMD_leg": "legs BMD",
}


def panel_label(ax, txt, x=-0.02, y=1.05):
    ax.text(x, y, txt, transform=ax.transAxes, fontsize=9,
            fontweight="bold", va="top", ha="right")


DUMP = {}


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(os.path.join(OUT, name + ".png"))
    fig.savefig(os.path.join(OUT, name + ".pdf"))
    plt.close(fig)
    print("  OK", name)


# ================================================================== 数据
def load_power():
    """|z| 直接由 mr_corrected_N.csv 的 β/SE 现算，避免读 mr_power.csv 的 2 位小数舍入值。"""
    d = pd.read_csv(os.path.join(GEN, "mr_power.csv"))
    src = pd.read_csv(os.path.join(GEN, "mr_corrected_N.csv"))
    exact = {(e, o): (b, s) for e, o, b, s in
             zip(src.exposure, src.outcome, src.beta, src.se)}
    d["abs_z"] = [abs(exact[(e, o)][0] / exact[(e, o)][1])
                  for e, o in zip(d.exposure, d.outcome)]
    d["power_curve"] = [float(norm.sf(1.96 - z) + norm.cdf(-1.96 - z))
                        for z in d.abs_z]
    d["pair"] = [f"{PRETTY.get(e, e)} → {PRETTY.get(o, o)}"
                 for e, o in zip(d.exposure, d.outcome)]
    return d


def load_overlap():
    o = pd.read_csv(os.path.join(GEN, "overlap_statistics.csv"))
    m = pd.read_csv(os.path.join(GEN, "mr_corrected_N.csv"))
    key = {(e, x): (b, s) for e, x, b, s in
           zip(m.exposure, m.outcome, m.beta, m.se)}
    rows = []
    for _, r in o.iterrows():
        stem = str(r["pair"]).replace(".csv", "")
        exp, out = stem.split("__")
        b, s = key[(exp, out)]
        rows.append(dict(pair=f"{PRETTY.get(exp, exp)} → {PRETTY.get(out, out)}",
                         rho=r.rho_GE_zcorr, K=int(r.n_iv),
                         z=b / s, beta=b))
    return pd.DataFrame(rows)


# ================================================================== 图 8
def fig_power():
    d = load_power()
    fig = plt.figure(figsize=(7.2, 3.9))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.05], wspace=0.42)

    # ---- (a) 效能函数 ----
    ax = fig.add_subplot(gs[0])
    az = np.linspace(0, 16.5, 600)
    pw = norm.sf(1.96 - az) + norm.cdf(-1.96 - az)
    ax.plot(az, pw, "-", color=CB["ink"], lw=1.4, zorder=3)
    ax.axhline(0.80, color=CB["grey"], ls="--", lw=0.9, zorder=1)
    ax.axvline(2.80, color=CB["cb_orange"], ls=":", lw=1.0, zorder=1)
    ax.text(2.80 + 0.25, 0.055, "80% power\n(|z| = 2.80)",
            fontsize=6, color=CB["cb_orange"], va="bottom")

    ok = d.verdict.eq("powered").values
    ax.scatter(d.abs_z[ok], d.power_curve[ok], s=34, facecolor=CB["cb_blue"],
               edgecolors="white", linewidths=0.5, zorder=5, label="Powered pairs (n=7)")
    ax.scatter(d.abs_z[~ok], d.power_curve[~ok], s=46, facecolor="none",
               edgecolors=CB["red"], linewidths=1.2, zorder=5,
               label="Underpowered nulls (n=3)")
    ax.set_xlim(0, 16.5)
    ax.set_ylim(0, 1.04)
    ax.set_xlabel("|β| / SE  (standardised effect size, |z|)", fontsize=7)
    ax.set_ylabel("Power at two-sided α = 0.05", fontsize=7)
    ax.set_xticks([0, 2.8, 5, 10, 15])
    ax.tick_params(labelsize=6.5)
    ax.legend(fontsize=6, loc="lower right", frameon=True, framealpha=0.95,
              edgecolor="#DDDDDD")
    panel_label(ax, "a")

    # ---- (b) 逐对实测效能 ----
    ax = fig.add_subplot(gs[1])
    dd = d.sort_values("power_observed", ascending=True).reset_index(drop=True)
    y = np.arange(len(dd))
    cols = [CB["cb_blue"] if v == "powered" else CB["red"]
            for v in dd.verdict]
    ax.barh(y, dd.power_observed, height=0.62, color=cols,
            edgecolor="white", linewidth=0.4, zorder=3)
    ax.axvline(0.80, color=CB["cb_orange"], ls=":", lw=1.0, zorder=4)
    for yy, v, k in zip(y, dd.power_observed, dd.K):
        ax.text(min(v + 0.015, 0.99), yy, f"{v:.2f}  (K={k})",
                va="center", ha="left" if v < 0.93 else "right",
                fontsize=5.8, color=CB["ink"],
                bbox=dict(facecolor="white", edgecolor="none", pad=0.6, alpha=0.75)
                if v >= 0.93 else None)
    ax.set_yticks(y)
    ax.set_yticklabels(dd.pair, fontsize=6)
    ax.set_xlim(0, 1.22)
    ax.set_ylim(-0.6, len(dd) - 0.4)
    ax.set_xlabel("Achieved power for the observed point estimate", fontsize=7)
    ax.tick_params(labelsize=6.5)
    ax.text(0.805, len(dd) - 0.35, "80%", fontsize=6, color=CB["cb_orange"])
    panel_label(ax, "b")

    # ---- 记录实际绘制值（供独立审计 92_audit_figures_v9.py 端到端复核）----
    DUMP["figure8"] = dict(
        panel_a=dict(
            curve_x=az.tolist(), curve_y=pw.tolist(),
            thresholds=dict(power80=0.80, z80=2.80),
            points=[dict(pair=p, abs_z=float(z), power=float(v), K=int(k),
                         verdict=vd)
                    for p, z, v, k, vd in zip(d.pair, d.abs_z, d.power_curve,
                                              d.K, d.verdict)],
            xlim=(0.0, 16.5), ylim=(0.0, 1.04)),
        panel_b=dict(
            order=dd.pair.tolist(),
            power=dd.power_observed.tolist(),
            K=dd.K.tolist(),
            verdict=dd.verdict.tolist(),
            xlim=(0.0, 1.22)),
    )

    fig.suptitle("Statistical power of the two-sample MR design",
                 fontsize=9, fontweight="bold", x=0.008, ha="left", y=1.0)
    save(fig, "Figure8_MR_power")


# ================================================================== 图 9
def fig_overlap():
    d = load_overlap().sort_values("rho").reset_index(drop=True)

    fig = plt.figure(figsize=(7.2, 4.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1.0], wspace=0.40)

    # ---- (a) ρ_GE 逐对 ----
    ax = fig.add_subplot(gs[0])
    y = np.arange(len(d))
    cols = [CB["cb_blue"] if v > 0 else CB["red"] for v in d.rho]
    ax.barh(y, d.rho, height=0.62, color=cols, edgecolor="white",
            linewidth=0.4, zorder=3)
    ax.axvline(0, color=CB["ink"], lw=0.9, zorder=2)
    for yy, v in zip(y, d.rho):
        off = 0.02 if v >= 0 else -0.02
        ax.text(v + off, yy, f"{v:+.2f}", va="center",
                ha="left" if v >= 0 else "right", fontsize=5.8)
    ax.set_yticks(y)
    ax.set_yticklabels(d.pair, fontsize=6)
    ax.set_xlim(-0.95, 0.72)
    ax.set_ylim(-0.6, len(d) - 0.4)
    ax.set_xlabel("ρ_GE: correlation of SNP–exposure with SNP–outcome z statistics",
                  fontsize=7)
    ax.tick_params(labelsize=6.5)
    handles = [Line2D([0], [0], color=CB["cb_blue"], lw=5, label="Positive (BMD→OA)"),
               Line2D([0], [0], color=CB["red"], lw=5, label="Negative (BMD→fracture)")]
    ax.legend(handles=handles, fontsize=6, loc="center right", frameon=True,
              framealpha=0.95, edgecolor="#DDDDDD")
    panel_label(ax, "a")

    # ---- (b) ρ_GE 符号 vs IVW z 符号 ----
    ax = fig.add_subplot(gs[1])
    ax.axhline(0, color=CB["faint"], lw=0.8, zorder=1)
    ax.axvline(0, color=CB["faint"], lw=0.8, zorder=1)
    # 一致象限（一、三象限）浅底
    xl, xh = -0.95, 0.72
    yl, yh = -16.5, 7.5
    ax.add_patch(plt.Rectangle((0, 0), xh, yh, facecolor=CB["cb_blue"],
                               alpha=0.06, zorder=0))
    ax.add_patch(plt.Rectangle((xl, yl), -xl, -yl, facecolor=CB["cb_blue"],
                               alpha=0.06, zorder=0))
    same = np.sign(d.rho) == np.sign(d.z)
    ax.scatter(d.rho[same], d.z[same], s=42, facecolor=CB["cb_blue"],
               edgecolors="white", linewidths=0.5, zorder=4,
               label="Sign-concordant (10/10)")
    if (~same).any():
        ax.scatter(d.rho[~same], d.z[~same], s=60, facecolor="none",
                   edgecolors=CB["red"], linewidths=1.3, zorder=5,
                   label="Sign-discordant")
    for _, r in d.iterrows():
        if abs(r.z) > 8 or abs(r.rho) > 0.6:
            ax.annotate(f"{r.z:+.1f}", (r.rho, r.z), fontsize=5.6,
                        xytext=(4, 4), textcoords="offset points")
    ax.set_xlim(xl, xh)
    ax.set_ylim(yl, yh)
    ax.set_xlabel("ρ_GE (overlap / pleiotropy gauge)", fontsize=7)
    ax.set_ylabel("IVW z statistic (signed)", fontsize=7)
    ax.tick_params(labelsize=6.5)
    ax.legend(fontsize=6, loc="lower left", frameon=True, framealpha=0.95,
              edgecolor="#DDDDDD")
    panel_label(ax, "b")

    # ---- 记录实际绘制值（供独立审计）----
    DUMP["figure9"] = dict(
        panel_a=dict(order=d.pair.tolist(), rho=d.rho.tolist(),
                     K=d.K.tolist(), xlim=(-0.95, 0.72)),
        panel_b=dict(rho=d.rho.tolist(), z=d.z.tolist(),
                     xlim=(xl, xh), ylim=(yl, yh),
                     n_sign_concordant=int(same.sum())),
    )

    fig.suptitle("Sample-overlap diagnostics for the ten MR pairs",
                 fontsize=9, fontweight="bold", x=0.008, ha="left", y=1.0)
    save(fig, "Figure9_overlap_pleiotropy")


# ================================================================== 补充图 S3
def fig_gefos():
    """外部效度森林图（布局取 science-figures 的 forest 规范，字体取 nature-figure 规范）"""
    g = json.load(open(os.path.join(GEN, "gefos_benchmark.json"), encoding="utf-8"))
    labs, est = g["labels"], np.array(g["estimates"])
    lo, hi = np.array(g["ci_low"]), np.array(g["ci_high"])

    fig, ax = plt.subplots(figsize=(7.2, 2.7))
    y = np.arange(len(labs))[::-1]          # 第一条在最上

    ax.axvspan(1.48, 1.63, facecolor=CB["grey"], alpha=0.13, zorder=0,
               label="GEFOS 95% CI band")
    ax.axvline(1.55, color=CB["grey"], ls=":", lw=1.0, zorder=1)

    for i, (yy, e, l, h) in enumerate(zip(y, est, lo, hi)):
        ours = i > 0
        c = CB["cb_blue"] if ours else CB["ink"]
        ax.plot([l, h], [yy, yy], "-", color=c, lw=1.6 if ours else 1.4,
                solid_capstyle="butt", zorder=2)
        ax.plot([l, l], [yy - 0.13, yy + 0.13], "-", color=c, lw=1.2, zorder=2)
        ax.plot([h, h], [yy - 0.13, yy + 0.13], "-", color=c, lw=1.2, zorder=2)
        ax.scatter(e, yy, s=78 if ours else 92,
                   facecolor=c if ours else "white",
                   edgecolors=c, linewidths=1.2,
                   marker="o" if ours else "s", zorder=4)
        ax.text(1.675, yy, f"{e:.2f} [{l:.2f}, {h:.2f}]", va="center",
                fontsize=6.5, color=CB["ink"])

    ax.set_yticks(y)
    ax.set_yticklabels(labs, fontsize=6.5)
    ax.set_xlim(1.28, 1.80)
    ax.set_ylim(-0.6, len(labs) - 0.4)
    ax.set_xlabel("Odds ratio for any fracture per 1 SD lower BMD", fontsize=7.5)
    ax.tick_params(labelsize=6.5)
    handles = [
        Line2D([0], [0], marker="s", color=CB["ink"], markerfacecolor="white",
               markersize=6, lw=1.2, label="External benchmark (GEFOS 2018)"),
        Line2D([0], [0], marker="o", color=CB["cb_blue"], markerfacecolor=CB["cb_blue"],
               markersize=6, lw=1.2, label="This study"),
    ]
    ax.legend(handles=handles, fontsize=6.2, loc="upper left", frameon=True,
              framealpha=0.95, edgecolor="#DDDDDD")
    fig.suptitle("External validation against the GEFOS benchmark",
                 fontsize=9, fontweight="bold", x=0.008, ha="left", y=1.0)

    DUMP["supplfig3"] = dict(labels=labs, estimates=est.tolist(),
                             ci_low=lo.tolist(), ci_high=hi.tolist(),
                             xlim=(1.28, 1.80), ref_band=(1.48, 1.63))
    save(fig, "SupplFig3_GEFOS_benchmark")


# ================================================================== 主流程
def main():
    print(">> 依据真实统计结果补充新图")
    fig_power()
    fig_overlap()
    fig_gefos()
    # 落盘实际绘制值，供 92_audit_figures_v9.py 做端到端独立复核
    with open(os.path.join(GEN, "figure_v9_plotted_data.json"), "w",
              encoding="utf-8") as fh:
        json.dump(DUMP, fh, ensure_ascii=False, indent=1)
    print("   plotted-data dump ->", os.path.join(GEN, "figure_v9_plotted_data.json"))
    print("NEW FIGURES DONE ->", OUT)


if __name__ == "__main__":
    main()
