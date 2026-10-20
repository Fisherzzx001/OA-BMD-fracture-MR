# -*- coding: utf-8 -*-
"""
114 · Supplementary Figure S12 —— 与既往文献的定量对照（BMD→OA）

设计要点（诚实性优先）：
  * 文献 OR 在 eBMD（足跟超声）尺度、本研究在 areal-BMD per-SD 重建尺度 ——
    绝对 OR 不可直接比较（Methods 已声明），因此拆成两个面板，各自带自己的横轴，
    并在图注强调"可比的是方向与髋/膝对比，不是绝对值"。
  * Zhou 2025 只报告 β=−2.17（P<0.01），未给 OR —— 只画方向标记，不伪造点估计。

数据来源：
  本研究者 5 对：code/output/genetics/AUDIT_mr_recomputed.csv（= Table 3，先断言再 exp 换算）
  Hartley 2022 IJE 51(4):1254-1267（1-sample MR，heel-ultrasound eBMD，per SD increase）：hip 1.28(1.05–1.57)、knee 1.40(1.20–1.63)
  Funck-Brentano 2019 Arthritis Rheumatol 71(10):1634-1641（2-sample MR，femoral-neck BMD，per SD increase）：all-OA OR 1.14(1.06–1.22)
  Zhou 2025 Front Genet 16:1597005：FN-BMD→hip OA β=−2.17（P<0.01），方向为保护
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle_v18 import C, GEN, save, panel_label, dump_plotted, SAVED

# ---- 文献值（已在正文 Discussion 引用，本轮独立核实）----
LIT = [  # (label, OR, lo, hi, joint, note)
    ("Hartley 2022 · hip OA (eBMD)", 1.28, 1.05, 1.57, "hip", ""),
    ("Hartley 2022 · knee OA (eBMD)", 1.40, 1.20, 1.63, "knee", ""),
    ("Funck-Brentano 2019 · all-OA (femoral neck BMD)", 1.14, 1.06, 1.22, "all", ""),
]
ZHOU_NOTE = ("Zhou 2025 \u00b7 femoral-neck BMD \u2192 hip OA: \u03b2=\u22122.17 (P<0.01), "
              "protective; OR not reported \u2014 direction only")

# ---- 本研究：Table 2 MR β（G6 尺度修正后，断言用）----
OURS = {  # key -> (beta, se)  —— log-OR 尺度，与 AUDIT_mr_recomputed 一致
    "BMD_total":      (0.078, 0.013),
    "BMD_lumbar":     (0.095, 0.016),
    "BMD_arms":       (0.075, 0.016),
    "BMD_leg":        (0.060, 0.014),
    "BMD_femur_neck": (0.058, 0.017),
}
OUR_LABEL = {
    "BMD_lumbar": "Lumbar-spine BMD \u2192 knee OA",
    "BMD_total": "BMD total \u2192 knee OA",
    "BMD_arms": "Arms BMD \u2192 knee OA",
    "BMD_leg": "Legs BMD \u2192 knee OA",
    "BMD_femur_neck": "Femoral-neck BMD \u2192 knee OA",
}
Z = 1.959963985
JOINT_C = {"hip": C["blue"], "knee": C["orange"], "all": C["purple"]}


def load_ours():
    aud = pd.read_csv(os.path.join(GEN, "AUDIT_mr_recomputed.csv")).set_index("pair")
    rows = {}
    for k, (b_ms, se_ms) in OURS.items():
        a = aud.loc[f"{k}__OA_knee"]
        assert abs(a["beta"] - b_ms) < 5e-4 and abs(a["se"] - se_ms) < 5e-4, \
            (k, a["beta"], b_ms)
        lo, hi = a["beta"] - Z * a["se"], a["beta"] + Z * a["se"]
        rows[k] = dict(beta=float(a["beta"]), se=float(a["se"]),
                       OR=float(np.exp(a["beta"])),
                       lo=float(np.exp(lo)), hi=float(np.exp(hi)))
    return rows


def fig_S12():
    ours = load_ours()
    fig, (axa, axb) = plt.subplots(
        2, 1, figsize=(7.20, 4.60),
        gridspec_kw=dict(height_ratios=[1.15, 1.0], hspace=0.62,
                         left=0.235, right=0.845, top=0.930, bottom=0.135))

    tr_a = mtransforms.blended_transform_factory(axa.transAxes, axa.transData)
    tr_b = mtransforms.blended_transform_factory(axb.transAxes, axb.transData)

    # ================= (a) 已发表估计 =================
    n = len(LIT) + 1          # + Zhou 方向行
    ys = np.arange(n)[::-1]
    for y, (lab, orv, lo, hi, joint, note) in zip(ys[:len(LIT)], LIT):
        col = JOINT_C[joint]
        axa.plot([lo, hi], [y, y], lw=1.3, color=col, solid_capstyle="butt")
        axa.plot(orv, y, marker="o", ms=5.0, color=col, mec="white", mew=0.6, zorder=3)
        axa.text(1.02, y, f"{orv:.2f} ({lo:.2f}\u2013{hi:.2f})", transform=tr_a,
                 fontsize=6.0, ha="left", va="center", color="#333333")
    # Zhou 方向行
    yz = ys[-1]
    axa.annotate("", xy=(0.862, yz), xytext=(0.98, yz),
                 arrowprops=dict(arrowstyle="-|>", lw=1.2, color=C["grey"]))
    axa.text(0.985, yz, ZHOU_NOTE, transform=tr_a, fontsize=5.8, ha="left",
             va="center", color=C["grey"], style="italic")
    axa.set_yticks(ys)
    axa.set_yticklabels([l for l, *_ in LIT] + ["Zhou 2025 \u00b7 FN-BMD \u2192 hip OA"],
                        fontsize=6.2)
    axa.axvline(1.0, lw=0.8, color=C["ink"], ls="--", alpha=0.55)
    axa.set_xscale("log")
    axa.set_xlim(0.84, 1.72)
    axa.set_xticks([0.9, 1.0, 1.2, 1.4, 1.6])
    axa.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:g}"))
    axa.set_xlabel("OR per 1 SD higher BMD (log scale)", fontsize=7.2)
    axa.set_ylim(-0.7, n - 0.3)
    axa.tick_params(axis="y", length=0)
    for s in ("top", "right", "left"):
        axa.spines[s].set_visible(False)
    axa.set_title("Published MR estimates (eBMD scale)", fontsize=8.2, loc="left", pad=4)
    panel_label(axa, "a")
    # Hartley 髋/膝 within-study 对照标注
    yk = ys[1]
    axa.annotate("", xy=(1.28, yk + 0.38), xytext=(1.40, yk + 0.38),
                 arrowprops=dict(arrowstyle="<->", lw=0.7, color=C["grey"]))
    axa.text(1.335, yk + 0.52, "hip 27% below knee (log-odds, same study)",
             fontsize=5.6, ha="center", va="bottom", color=C["grey"], style="italic")

    # ================= (b) 本研究 =================
    order = ["BMD_lumbar", "BMD_total", "BMD_arms", "BMD_leg", "BMD_femur_neck"]
    nb = len(order)
    ysb = np.arange(nb)[::-1]
    for y, k in zip(ysb, order):
        r = ours[k]
        axb.plot([r["lo"], r["hi"]], [y, y], lw=1.3, color=C["orange"],
                 solid_capstyle="butt")
        axb.plot(r["OR"], y, marker="o", ms=5.0, color=C["orange"], mec="white",
                 mew=0.6, zorder=3)
        axb.text(1.02, y, f'{r["OR"]:.3f} ({r["lo"]:.3f}\u2013{r["hi"]:.3f})',
                 transform=tr_b, fontsize=6.0, ha="left", va="center", color="#333333")
    axb.axvline(1.0, lw=0.8, color=C["ink"], ls="--", alpha=0.55)
    axb.set_yticks(ysb)
    axb.set_yticklabels([OUR_LABEL[k] for k in order], fontsize=6.2)
    axb.set_xscale("log")
    axb.set_xlim(0.995, 1.045)
    axb.set_xticks([0.995, 1.00, 1.01, 1.02, 1.03, 1.04])
    axb.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:g}"))
    axb.set_xlabel("OR per 1 SD higher areal BMD (log scale)", fontsize=7.2)
    axb.set_ylim(-0.7, nb - 0.3)
    axb.tick_params(axis="y", length=0)
    for s in ("top", "right", "left"):
        axb.spines[s].set_visible(False)
    axb.set_title("This study (5 BMD sites \u2192 knee OA)", fontsize=8.2,
                  loc="left", pad=4)
    panel_label(axb, "b")

    # 底注
    fig.text(0.010, 0.012,
             "ORs in (a) are per 1 SD higher BMD: Hartley et al. on the heel-ultrasound eBMD scale, Funck-Brentano et al. on the femoral-neck BMD scale. "
             "ORs in (b) are on this study's reconstructed "
             "areal-BMD per-SD scale.\nThe two scales are not directly comparable "
             "(Methods, Effect-scale reconstruction); direction and the hip-versus-knee "
             "contrast are the comparables. F1 criterion: hip effect positive AND \u2265 knee.",
             fontsize=5.8, ha="left", va="bottom", color="#444444", linespacing=1.5)

    save(fig, "SupplFig12_literature_comparison")
    dump_plotted({"SupplFig12_literature_comparison": {
        "sources": [
            "code/output/genetics/AUDIT_mr_recomputed.csv (this study; OR=exp(beta))",
            "Hartley 2022 IJE 51(4):1254-1267 (1SMR, eBMD per SD): hip 1.28 (1.05-1.57), knee 1.40 (1.20-1.63)",
            "Funck-Brentano 2019 Arthritis Rheumatol 71(10):1634-1641: femoral neck BMD->OA, all-OA OR 1.14 (1.06-1.22)",
            "Zhou 2025 Front Genet 16:1597005: FN-BMD->hip OA beta=-2.17, P<0.01 (no OR reported)",
        ],
        "ours": ours,
        "scale_caveat": "eBMD scale vs areal-BMD reconstructed per-SD scale; not directly comparable",
    }})


if __name__ == "__main__":
    print("\u25b6 Supplementary Figure S12 — literature comparison")
    fig_S12()
    print(f"\n\u2714 {len(SAVED)} figure(s) written")
