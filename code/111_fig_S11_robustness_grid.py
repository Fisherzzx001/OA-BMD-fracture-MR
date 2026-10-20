# -*- coding: utf-8 -*-
"""
111 · Supplementary Figure S11 —— 多估计器稳健性总览（10 对 MR × 7 种估计）

回答审稿最大风险点："I² 高达 68–86%，凭什么信？" —— 一图展示方向 10/10 稳定。

数据来源（全部一手，脚本内重算，不抄派生表）：
  code/output/genetics/AUDIT_mr_recomputed.csv   IVW fixed / IVW random / Egger /
                                                 weighted median（= Table 3，已经独立审计）
  code/output/genetics/harmonised/*.csv          PRESSO 离群剔除后 IVW（用 presso_outlier
                                                 列剔除后直接重算 IVW，Table 3 同尺度）
  code/output/mrlap_overlap_results.csv          MRlap 校正后 β/SE/p
  code/output/genetics/mvmr_height_grid.csv      MVMR（spec='EUR|n_h=100'，Table S9 主配置）

编码：底色 = 效应方向（暖 = 正 / 冷 = 负）；符号 = 显著性
      ● p<0.05   ○ p≥0.05   ◊ 仅方向（源输出无 SE，不作显著性声称）
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle_v18 import C, GEN, PROJECT, save, dump_plotted, SAVED

HARM = os.path.join(GEN, "harmonised")

# Table 3 印刷值（用于断言 AUDIT_mr_recomputed 就是 Table 3）
# 已同步 G6 尺度修正（β/SE 为 log-OR 尺度）与步近约定加权中位数（Bowden 2016）。
TABLE3 = {
    "BMD_total__fracture":      (-0.337, 0.022, 3.65e-51, -0.393, 0.558),
    "BMD_femur_neck__fracture": (-0.351, 0.030, 7.71e-31, -0.365, 0.018),
    "BMD_total__OA_knee":       (+0.078, 0.013, 3.86e-10, +0.065, 0.810),
    "BMD_lumbar__OA_knee":      (+0.095, 0.016, 7.70e-09, +0.115, 0.295),
    "BMD_arms__OA_knee":        (+0.075, 0.016, 1.71e-06, +0.026, 0.982),
    "BMD_leg__OA_knee":         (+0.060, 0.014, 2.29e-05, +0.026, 0.815),
    "BMD_femur_neck__OA_knee":  (+0.058, 0.017, 6.36e-04, +0.052, 0.875),
    "OA_knee__fracture":        (+0.027, 0.032, 3.95e-01, +0.019, 0.476),
    "OA_knee__BMD_total":       (+0.031, 0.032, 3.26e-01, +0.026, 0.376),
    "OA_knee__BMD_femur_neck":  (+0.042, 0.032, 1.90e-01, +0.023, 0.169),
}

ROW_ORDER = ["BMD_total__fracture", "BMD_femur_neck__fracture",
             "BMD_total__OA_knee", "BMD_lumbar__OA_knee", "BMD_arms__OA_knee",
             "BMD_leg__OA_knee", "BMD_femur_neck__OA_knee",
             "OA_knee__fracture", "OA_knee__BMD_total", "OA_knee__BMD_femur_neck"]

ROW_LABEL = {
    "BMD_total__fracture":      "BMD total \u2192 fracture  (control)",
    "BMD_femur_neck__fracture": "Femoral-neck BMD \u2192 fracture  (control)",
    "BMD_total__OA_knee":       "BMD total \u2192 knee OA",
    "BMD_lumbar__OA_knee":      "Lumbar-spine BMD \u2192 knee OA",
    "BMD_arms__OA_knee":        "Arms BMD \u2192 knee OA",
    "BMD_leg__OA_knee":         "Legs BMD \u2192 knee OA",
    "BMD_femur_neck__OA_knee":  "Femoral-neck BMD \u2192 knee OA",
    "OA_knee__fracture":        "Knee OA \u2192 fracture",
    "OA_knee__BMD_total":       "Knee OA \u2192 BMD total",
    "OA_knee__BMD_femur_neck":  "Knee OA \u2192 BMD femoral neck",
}

COLS = ["IVW (fixed)", "IVW (random)", "MR-Egger", "Weighted\nmedian",
        "PRESSO\noutlier-removed", "MRlap\ncorrected", "MVMR\n+ height"]

BG_POS = "#FCEBCB"   # 暖（正效应）
BG_NEG = "#D2E5F2"   # 冷（负效应）
BG_NA = "#F3F3F3"


def fmt_b(x):
    return f"{x:+.3f}"


def cell(sym, beta, sig_state):
    """sig_state: 'sig' | 'ns' | 'dir'；返回 dict，保留原始数值供一致性计数。"""
    return {"sym": sym, "beta": float(beta), "txt": fmt_b(beta), "st": sig_state}


def build_cells():
    aud = pd.read_csv(os.path.join(GEN, "AUDIT_mr_recomputed.csv"))
    aud = aud.set_index("pair")
    mrlap = pd.read_csv(os.path.join(PROJECT, "code", "output", "mrlap_overlap_results.csv"))
    mrlap["pair"] = mrlap["exposure"] + "__" + mrlap["outcome"]
    mrlap = mrlap.set_index("pair")
    g = pd.read_csv(os.path.join(GEN, "mvmr_height_grid.csv"))
    g = g[g["spec"] == "EUR|n_h=100"]
    g["pair"] = g["exposure"] + "__" + g["outcome"]
    g = g.set_index("pair")

    grid = {}
    details = {}
    for key in ROW_ORDER:
        a = aud.loc[key]
        e_beta, e_se, e_p, e_wm, e_ip = TABLE3[key]
        # 断言：AUDIT_mr_recomputed == Table 3
        assert abs(a["beta"] - e_beta) < 5e-4, (key, a["beta"], e_beta)
        assert abs(a["wmed"] - e_wm) < 5e-4, (key, a["wmed"], e_wm)
        assert abs(a["egger_int_p"] - e_ip) < 5e-4, (key, a["egger_int_p"], e_ip)

        # --- PRESSO 离群剔除后 IVW（harmonised 一手重算，与 Table 3 同尺度）---
        df = pd.read_csv(os.path.join(HARM, key + ".csv"))
        n_out = int(df["presso_outlier"].astype(float).sum())
        m = df[np.asarray(df["presso_outlier"].astype(float) == 0)]
        b = (m["beta_o"] / m["beta_e"]).to_numpy(float)
        se = np.abs(m["se_o"] / m["beta_e"]).to_numpy(float)
        w = 1.0 / se ** 2
        b_out = float(np.sum(w * b) / np.sum(w))
        se_out = float(np.sqrt(1.0 / np.sum(w)))
        p_out = float(2 * stats.norm.sf(abs(b_out / se_out)))
        details[key] = {"n_outliers": n_out, "n_after": int(len(m)),
                        "presso_out_beta": b_out, "presso_out_se": se_out,
                        "presso_out_p": p_out}

        # --- MRlap ---
        r = mrlap.loc[key]
        # 方向一致性断言：MRlap observed 与 IVW 同号
        assert np.sign(r["observed_effect"]) == np.sign(a["beta"]), key
        details[key].update({"mrlap_observed": float(r["observed_effect"]),
                             "mrlap_corrected": float(r["corrected_effect"]),
                             "mrlap_se": float(r["corrected_effect_se"]),
                             "mrlap_p": float(r["corrected_effect_p"])})

        # --- MVMR ---
        has_mvmr = key in g.index
        if has_mvmr:
            gg = g.loc[key]
            details[key].update({"mvmr_beta": float(gg["mvmr_beta"]),
                                 "mvmr_p": float(gg["mvmr_p"])})

        cells = [
            cell("\u25cf" if a["p"] < 0.05 else "\u25cb", a["beta"],
                 "sig" if a["p"] < 0.05 else "ns"),
            cell("\u25cf" if a["p_re"] < 0.05 else "\u25cb", a["beta"],
                 "sig" if a["p_re"] < 0.05 else "ns"),
            cell("\u25cf" if a["egger_p"] < 0.05 else "\u25cb", a["egger_beta"],
                 "sig" if a["egger_p"] < 0.05 else "ns"),
            cell("\u25ca", a["wmed"], "dir"),
            cell("\u25cf" if p_out < 0.05 else "\u25cb", b_out,
                 "sig" if p_out < 0.05 else "ns"),
            cell("\u25cf" if r["corrected_effect_p"] < 0.05 else "\u25cb",
                 r["corrected_effect"],
                 "sig" if r["corrected_effect_p"] < 0.05 else "ns"),
            cell("\u25cf" if has_mvmr and gg["mvmr_p"] < 0.05 else
                 ("\u25cb" if has_mvmr else ""),
                 gg["mvmr_beta"] if has_mvmr else np.nan,
                 "sig" if has_mvmr and gg["mvmr_p"] < 0.05 else
                 ("ns" if has_mvmr else "na")),
        ]
        grid[key] = cells

        # 方向一致性计数（MVMR 不适用时按 6 计）
        primary_sign = np.sign(a["beta"])
        n_est = sum(1 for c in cells if c["st"] != "na")
        n_conc = sum(1 for c in cells if c["st"] != "na"
                     and np.sign(c["beta"]) == primary_sign)
        grid[key] = (cells, f"{n_conc}/{n_est}")
    return grid, details


def fig_S11():
    grid, details = build_cells()
    n_r, n_c = len(ROW_ORDER), len(COLS)

    fig, ax = plt.subplots(figsize=(7.20, 4.55))
    gs_x0, gs_x1 = 0.262, 0.845
    x0, x1 = gs_x0, gs_x1
    y_top, y_bot = 0.865, 0.045
    cw = (x1 - x0) / n_c
    rh = (y_top - y_bot) / n_r

    group_bounds = [(0, 2), (2, 7), (7, 10)]

    # ---- 单元格
    for i, key in enumerate(ROW_ORDER):
        cells, conc = grid[key]
        y = y_top - (i + 0.5) * rh
        for j, c in enumerate(cells):
            sym, beta, st = c["sym"], c["beta"], c["st"]
            xc = x0 + (j + 0.5) * cw
            if st == "na":
                fc = BG_NA
            elif np.sign(beta) > 0:
                fc = BG_POS
            else:
                fc = BG_NEG
            ax.add_patch(Rectangle((x0 + j * cw, y - rh / 2), cw, rh,
                                   fc=fc, ec="white", lw=1.0, zorder=1))
            if st == "na":
                txt, col = "n/a", C["grey"]
            else:
                col = C["ink"] if st in ("sig", "dir") else "#6B6B6B"
                txt = f"{sym} {c['txt']}"
            ax.text(xc, y, txt, ha="center", va="center", fontsize=6.5,
                    color=col, zorder=3,
                    fontweight="bold" if st == "sig" else "normal")
        # 方向一致性列
        xc = x0 + n_c * cw + 0.008
        ax.add_patch(Rectangle((x0 + n_c * cw, y - rh / 2), 1.0 - 0.005 - (x0 + n_c * cw), rh,
                               fc="#FFFFFF", ec="white", lw=1.0, zorder=1))
        ax.text((x0 + n_c * cw + 0.995) / 2 + 0.0025, y, conc, ha="center",
                va="center", fontsize=6.6, fontweight="bold", color=C["ink"], zorder=3)

    # ---- 组分隔线（分组信息：control 行标签已自带 "(control)"，反向组见图注）
    for (a0, a1) in group_bounds:
        y_edge = y_top - a1 * rh
        if a1 < n_r:
            ax.plot([x0 - 0.09, 0.995], [y_edge, y_edge], lw=0.7,
                    color=C["grey"], ls=(0, (3, 2)), zorder=2)

    # ---- 行标签 / 列头
    for i, key in enumerate(ROW_ORDER):
        y = y_top - (i + 0.5) * rh
        ax.text(x0 - 0.012, y, ROW_LABEL[key], transform=ax.transAxes,
                fontsize=6.6, ha="right", va="center")
    for j, cname in enumerate(COLS):
        xc = x0 + (j + 0.5) * cw
        ax.text(xc, y_top + 0.018, cname, transform=ax.transAxes, fontsize=6.4,
                ha="center", va="bottom")
    ax.text((x0 + n_c * cw + 0.995) / 2 + 0.0025, y_top + 0.018,
            "Sign\nconcordance", transform=ax.transAxes, fontsize=6.4,
            ha="center", va="bottom")

    # ---- 图例说明
    ax.text(x0, y_bot - 0.028,
            "\u25cf p<0.05    \u25cb p\u22650.05    \u25ca direction only (no SE in source "
            "output; significance not assessed)    n/a estimator not applicable to the pair",
            transform=ax.transAxes, fontsize=5.9, ha="left", va="top", color="#444444")
    ax.text(x0, y_bot - 0.085,
            "Cell shading shows the sign of the estimate (warm = positive, cool = negative). "
            "IVW fixed/random, MR-Egger, weighted median: Table 3. "
            "PRESSO outlier-removed: IVW recomputed on harmonised SNP data after removing "
            "MR-PRESSO outliers (0\u20133 per pair). MRlap: overlap-corrected estimate. "
            "MVMR: additionally conditioned on height (n_h=100/EUR).",
            transform=ax.transAxes, fontsize=5.9, ha="left", va="top", color="#444444")

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    save(fig, "SupplFig11_robustness_grid")

    dump_plotted({"SupplFig11_robustness_grid": {
        "sources": [
            "code/output/genetics/AUDIT_mr_recomputed.csv (= Table 3)",
            "code/output/genetics/harmonised/*.csv (presso_outlier column; PRESSO-out IVW recomputed)",
            "code/output/mrlap_overlap_results.csv",
            "code/output/genetics/mvmr_height_grid.csv :: spec=='EUR|n_h=100'",
        ],
        "encoding": {"sig": "\u25cf p<0.05", "ns": "\u25cb p\u22650.05",
                     "dir": "\u25ca direction only", "na": "not applicable"},
        "details": details,
    }})


if __name__ == "__main__":
    print("\u25b6 Supplementary Figure S11 — multi-estimator robustness grid")
    fig_S11()
    print(f"\n\u2714 {len(SAVED)} figure(s) written")
