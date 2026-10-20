# -*- coding: utf-8 -*-
"""
110 · Supplementary Figure S10 —— MVMR 分解：BMD→OA 效应中有多少来自骨骼尺寸（身高）

数据来源（一手，全部重算，不抄任何派生表）：
  code/output/genetics/mvmr_height_grid.csv
      spec == 'EUR|n_h=100'  （= Table S9 主配置，Table S10 网格中的主行）
      ref_uni_beta  : Table 3 的 univariable IVW β（主表基线）
      uni_beta/uni_se : 在 MVMR 工具集上重算的 univariable β（第二基线）
      mvmr_beta/mvmr_se : 条件化身高后的 MVMR-IVW β
      mvmr_h / mvmr_h_p : 身高的直接效应 θ_H

衰减定义（与 Table S10 脚注完全一致，脚本内闭式重算并断言）：
  Atten% vs Table 3     = 100 × (ref_uni_beta − mvmr_beta) / ref_uni_beta
  Atten% vs MVMR set    = 100 × (uni_beta     − mvmr_beta) / uni_beta

面板：(a) BMD→膝 OA 的 univariable / MVMR 配对森林图
      (b) BMD→骨折阳性对照的配对森林图
      (c) 衰减百分比条形图（条上标注 %；空心菱形 = MVMR 集基线；右列 = 独立于身高通路的份额）
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle_v18 import C, GEN, save, panel_label, dump_plotted, SAVED

SRC = os.path.join(GEN, "mvmr_height_grid.csv")

LABEL = {
    ("BMD_total", "OA_knee"): "BMD total body \u2192 knee OA",
    ("BMD_femur_neck", "OA_knee"): "Femoral-neck BMD \u2192 knee OA",
    ("BMD_lumbar", "OA_knee"): "Lumbar-spine BMD \u2192 knee OA",
    ("BMD_arms", "OA_knee"): "Arms BMD \u2192 knee OA",
    ("BMD_leg", "OA_knee"): "Legs BMD \u2192 knee OA",
    ("BMD_total", "fracture"): "BMD total body \u2192 fracture",
    ("BMD_femur_neck", "fracture"): "Femoral-neck BMD \u2192 fracture",
}
ORDER = [f"{e}|{o}" for e, o in [
    ("BMD_total", "OA_knee"), ("BMD_femur_neck", "OA_knee"),
    ("BMD_lumbar", "OA_knee"), ("BMD_arms", "OA_knee"),
    ("BMD_leg", "OA_knee"), ("BMD_total", "fracture"),
    ("BMD_femur_neck", "fracture")]]

# Table S10 已发表值（用于断言本脚本重算结果一致）
EXPECT_ATTEN_T3 = {"BMD_total|OA_knee": 19.1, "BMD_femur_neck|OA_knee": 14.6,
                   "BMD_lumbar|OA_knee": 27.6, "BMD_arms|OA_knee": 13.5,
                   "BMD_leg|OA_knee": 33.2, "BMD_total|fracture": -1.7,
                   "BMD_femur_neck|fracture": -3.0}
EXPECT_ATTEN_SET = {"BMD_total|OA_knee": 9.8, "BMD_femur_neck|OA_knee": 14.6,
                    "BMD_lumbar|OA_knee": 9.5, "BMD_arms|OA_knee": 11.0,
                    "BMD_leg|OA_knee": 21.2, "BMD_total|fracture": -2.8,
                    "BMD_femur_neck|fracture": -2.8}

Z = 1.959963985


def build_rows():
    g = pd.read_csv(SRC)
    s = g[g["spec"] == "EUR|n_h=100"].copy()
    assert len(s) == 7, f"expected 7 rows at EUR|n_h=100, got {len(s)}"
    rows = {}
    for _, r in s.iterrows():
        key = f"{r.exposure}|{r.outcome}"
        a_t3 = 100.0 * (r.ref_uni_beta - r.mvmr_beta) / r.ref_uni_beta
        a_set = 100.0 * (r.uni_beta - r.mvmr_beta) / r.uni_beta
        rows[key] = dict(
            exp=r.exposure, out=r.outcome,
            n_bmd=int(r.n_bmd), n_h_req=int(r.n_h),
            beta_t3=float(r.ref_uni_beta),
            beta_uni=float(r.uni_beta), se_uni=float(r.uni_se),
            beta_mv=float(r.mvmr_beta), se_mv=float(r.mvmr_se),
            atten_t3=float(a_t3), atten_set=float(a_set),
            theta_h=float(r.mvmr_h), theta_h_p=float(r.mvmr_h_p),
            fc_bmd=float(r.Fc_bmd), fc_h=float(r.Fc_h),
        )
        # 断言：与已发表 Table S10 / S9 一致（容差 0.05 个百分点）
        assert abs(a_t3 - EXPECT_ATTEN_T3[key]) < 0.06, \
            f"{key}: atten_vs_Table3 recomputed {a_t3:.3f} vs published {EXPECT_ATTEN_T3[key]}"
        assert abs(a_set - EXPECT_ATTEN_SET[key]) < 0.06, \
            f"{key}: atten_vs_MVMRset recomputed {a_set:.3f} vs published {EXPECT_ATTEN_SET[key]}"
    return rows


def paired_forest(ax, keys, rows, title, xlabel, xlim, xticks):
    """并排（univariable vs MVMR）森林图；同一工具集，位移仅来自条件化身高。"""
    n = len(keys)
    ys = np.arange(n)[::-1]
    for y, k in zip(ys, keys):
        r = rows[k]
        for dy, (b, se, col, mk) in [
                (0.19, (r["beta_uni"], r["se_uni"], C["blue"], "o")),
                (-0.19, (r["beta_mv"], r["se_mv"], C["orange"], "s"))]:
            lo, hi = b - Z * se, b + Z * se
            ax.plot([lo, hi], [y + dy, y + dy], lw=1.1, color=col, solid_capstyle="butt")
            ax.plot(b, y + dy, marker=mk, ms=4.6, color=col, mec="white", mew=0.5, zorder=3)
        # 位移箭头：univariable → MVMR
        ax.annotate("", xy=(r["beta_mv"], y - 0.19), xytext=(r["beta_uni"], y + 0.19),
                    arrowprops=dict(arrowstyle="-", lw=0.7, color=C["grey"], ls=(0, (2, 1.6))))
    ax.axvline(0, lw=0.8, color=C["ink"], ls="--", alpha=0.55)
    ax.set_yticks(ys)
    ax.set_yticklabels([f"{LABEL[(rows[k]['exp'], rows[k]['out'])]}\n({rows[k]['n_bmd']} IVs)"
                        for k in keys], fontsize=6.6)
    ax.set_ylim(-0.75, n - 0.25)
    ax.set_xlim(*xlim)
    ax.set_xticks(xticks)
    ax.set_xlabel(xlabel, fontsize=7.4)
    ax.set_title(title, fontsize=8.2, loc="left", pad=4)
    ax.tick_params(axis="y", length=0)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)


import matplotlib.transforms as mtransforms

def _right_outside(ax):
    """返回把 x 固定在轴外右侧、y 用数据坐标的 transform。"""
    return mtransforms.blended_transform_factory(ax.transAxes, ax.transData)


def fig_S10():
    rows = build_rows()
    fig = plt.figure(figsize=(7.20, 6.05))
    gs = GridSpec(2, 2, height_ratios=[1.0, 1.06], hspace=0.62, wspace=0.40,
                  left=0.215, right=0.975, top=0.925, bottom=0.085)

    oa = ORDER[:5]
    fx = ORDER[5:]

    # ---- (a) BMD → knee OA
    axa = fig.add_subplot(gs[0, 0])
    paired_forest(axa, oa, rows,
                  "Univariable vs height-conditioned MVMR",
                  "\u03b2 per 1 SD higher BMD\n(95% CI)",
                  (0.0, 0.12), [0, 0.02, 0.04, 0.06, 0.08, 0.10])
    panel_label(axa, "a")

    # ---- (b) BMD → fracture (positive control)
    axb = fig.add_subplot(gs[0, 1])
    paired_forest(axb, fx, rows,
                  "Positive control: BMD \u2192 fracture",
                  "\u03b2 per 1 SD higher BMD\n(95% CI)",
                  (-0.44, -0.26), [-0.44, -0.40, -0.36, -0.32, -0.28])
    panel_label(axb, "b")

    # 共享图例
    handles = [Line2D([], [], marker="o", ls="none", ms=4.6, color=C["blue"],
                      mec="white", mew=0.5, label="Univariable (MVMR instrument set)"),
               Line2D([], [], marker="s", ls="none", ms=4.6, color=C["orange"],
                      mec="white", mew=0.5, label="MVMR, conditioned on height"),
               Line2D([], [], color=C["grey"], lw=0.8, ls=(0, (2, 1.6)),
                      label="Shift attributable to conditioning"),
               Line2D([], [], marker="D", ms=4.0, mfc="white", mec=C["ink"], mew=0.8,
                      ls="none", label="Attenuation vs MVMR-set baseline (panel c)")]
    fig.legend(handles=handles, loc="upper center", ncol=4,
               bbox_to_anchor=(0.5, 0.998), fontsize=6.4, frameon=False,
               columnspacing=1.1, handletextpad=0.35)

    # ---- (c) 衰减百分比
    axc = fig.add_subplot(gs[1, :])
    tr_out = _right_outside(axc)
    n = len(ORDER)
    ys = np.arange(n)[::-1]
    for y, k in zip(ys, ORDER):
        r = rows[k]
        is_ctrl = r["out"] == "fracture"
        col = C["grey"] if is_ctrl else C["orange"]
        axc.barh(y, r["atten_t3"], height=0.56, color=col,
                 edgecolor="white", linewidth=0.6, zorder=2)
        axc.plot(r["atten_set"], y, marker="D", ms=4.0, mfc="white",
                 mec=C["ink"], mew=0.8, zorder=4)
        # 条上标注 %（负值条改放条上方，避免与空心菱形基线标记相撞）
        x = r["atten_t3"]
        if x >= 0:
            axc.text(x + 0.8, y, f"{x:+.1f}%", fontsize=6.6,
                     color=C["ink"], ha="left", va="center", fontweight="bold")
        else:
            axc.text(x, y - 0.46, f"{x:+.1f}%", fontsize=6.6,
                     color=C["ink"], ha="center", va="top", fontweight="bold")
        # 右列（轴外）：独立于身高通路的份额
        indep = 100.0 - x
        axc.text(1.015, y, f"{indep:.1f}%", transform=tr_out, fontsize=6.6,
                 color=C["ink"], ha="left", va="center")

    axc.axvline(0, lw=0.9, color=C["ink"])
    # OA 与对照之间的分隔
    axc.axhline(1.5, lw=0.7, color=C["grey"], ls=(0, (3, 2)), alpha=0.8)
    labels_c = []
    for k in ORDER:
        r = rows[k]
        lab = LABEL[(r["exp"], r["out"])]
        labels_c.append(lab + "  (positive control)" if r["out"] == "fracture" else lab)
    axc.set_yticks(ys)
    axc.set_yticklabels(labels_c, fontsize=6.6)
    for tick, k in zip(axc.get_yticklabels(), ORDER):
        if rows[k]["out"] == "fracture":
            tick.set_color(C["grey"])
            tick.set_fontstyle("italic")
    axc.set_ylim(-0.75, n - 0.25)
    axc.set_xlim(-6.2, 36.0)
    axc.set_xticks([-5, 0, 5, 10, 15, 20, 25, 30, 35])
    axc.set_xlabel("Attenuation of the BMD effect after conditioning on height (%)", fontsize=7.4)
    axc.tick_params(axis="y", length=0)
    for s in ("top", "right", "left"):
        axc.spines[s].set_visible(False)
    axc.text(1.015, n - 0.30 + 0.30, "Height-independent\nresidual",
             transform=tr_out, fontsize=6.2, ha="left", va="bottom", color=C["grey"])
    axc.set_title("Share of the BMD effect removed by conditioning on height",
                  fontsize=8.2, loc="left", pad=4)
    panel_label(axc, "c")

    save(fig, "SupplFig10_MVMR_decomposition")

    # ---- 溯源留痕
    dump_plotted({"SupplFig10_MVMR_decomposition": {
        "source": "code/output/genetics/mvmr_height_grid.csv :: spec=='EUR|n_h=100'",
        "atten_definition": {
            "vs_Table3": "100*(ref_uni_beta - mvmr_beta)/ref_uni_beta",
            "vs_MVMRset": "100*(uni_beta - mvmr_beta)/uni_beta",
        },
        "rows": rows,
        "CI": "beta +/- 1.959964*SE",
    }})


if __name__ == "__main__":
    print("\u25b6 Supplementary Figure S10 (file prefix SupplFig9) — MVMR decomposition")
    fig_S10()
    print(f"\n\u2714 {len(SAVED)} figure(s) written")
