#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
12_final_figures.py — 合成 13×13 矩阵并出发表级图件

图 1：13×13 LDSC rg 热图（OA / fracture / BMD 三组分块）
图 2：森林图 —— OA×BMD（20 对）、fracture×BMD（10 对）、OA×fracture（2 对）
"""
import os
import sys
import importlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
m09 = importlib.import_module("09_rg_matrix")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "genetics")

GROUP_ORDER = ["OA_knee", "OA_selfrep", "fracture"] + [k for k in m09.PHENOS
                                                       if k.startswith("BMD_")]


def build_13x13():
    rg = pd.read_csv(os.path.join(OUT, "rg_matrix.csv"), index_col=0)
    se = pd.read_csv(os.path.join(OUT, "rg_se_matrix.csv"), index_col=0)
    pv = pd.read_csv(os.path.join(OUT, "rg_p_matrix.csv"), index_col=0)
    fr = pd.read_csv(os.path.join(OUT, "rg_fracture.csv"))

    names12 = list(rg.index)
    n13 = GROUP_ORDER
    R = pd.DataFrame(np.nan, index=n13, columns=n13).copy()
    S = pd.DataFrame(np.nan, index=n13, columns=n13).copy()
    P = pd.DataFrame(np.nan, index=n13, columns=n13).copy()
    for t in n13:
        R.loc[t, t] = 1.0

    R.loc[names12, names12] = rg.values
    S.loc[names12, names12] = se.values
    P.loc[names12, names12] = pv.values

    for _, r in fr.iterrows():
        t = r.trait
        R.loc["fracture", t] = R.loc[t, "fracture"] = r.rg
        S.loc["fracture", t] = S.loc[t, "fracture"] = r.se
        P.loc["fracture", t] = P.loc[t, "fracture"] = r.p
    R.to_csv(os.path.join(OUT, "rg_matrix_13.csv"))
    S.to_csv(os.path.join(OUT, "rg_se_matrix_13.csv"))
    P.to_csv(os.path.join(OUT, "rg_p_matrix_13.csv"))
    return R, S, P


# ------------------------------------------------------------------ 图 1
def heatmap(R, P, path):
    labels = [l.replace("BMD_", "").replace("OA_", "OA:")
              .replace("fracture", "FRACTURE") for l in R.index]
    n = len(labels)
    M = R.values.copy()
    off = M[~np.eye(n, dtype=bool)]
    vmax = float(np.nanmax(np.abs(off))) * 1.12

    fig, ax = plt.subplots(figsize=(11.5, 10))
    im = ax.imshow(M, cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(n)); ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_yticks(range(n)); ax.set_yticklabels(labels)

    # 分组框线
    for edge, lw in ((2, 2.2), (3, 2.2)):
        ax.axhline(edge - 0.5, color="black", lw=lw)
        ax.axvline(edge - 0.5, color="black", lw=lw)

    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            v, p = M[i, j], P.values[i, j]
            if not np.isfinite(v):
                continue
            mark = "" if (not np.isfinite(p) or p >= 0.05) else (
                "***" if p < 1e-3 else "**" if p < 0.01 else "*")
            ax.text(j, i, f"{v:.2f}{mark}", ha="center", va="center", fontsize=7.2,
                    color="white" if abs(v) > vmax * 0.55 else "black",
                    fontweight="bold" if mark else "normal")
    ax.set_title("LDSC genetic correlation, 13 phenotypes\n"
                 "(* q<0.05, ** q<0.01, *** q<0.001; BH-FDR across 78 pairs)",
                 fontsize=11.5)
    fig.colorbar(im, ax=ax, shrink=0.8, label="genetic correlation r$_g$")
    plt.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


# ------------------------------------------------------------------ 图 2
def forest(R, S, P, path):
    rows = []                       # (label, rg, se, p, group)
    for t in GROUP_ORDER[3:]:       # fracture × BMD
        rows.append((t.replace("BMD_", "BMD: "), R.loc["fracture", t],
                     S.loc["fracture", t], P.loc["fracture", t], "fracture×BMD"))
    for t in GROUP_ORDER[3:]:       # OA_knee × BMD
        rows.append((t.replace("BMD_", "BMD: "), R.loc["OA_knee", t],
                     S.loc["OA_knee", t], P.loc["OA_knee", t], "OA_knee×BMD"))
    for t in GROUP_ORDER[3:]:       # OA_selfrep × BMD
        rows.append((t.replace("BMD_", "BMD: "), R.loc["OA_selfrep", t],
                     S.loc["OA_selfrep", t], P.loc["OA_selfrep", t],
                     "OA_selfrep×BMD"))
    rows.append(("OA_selfrep", R.loc["OA_knee", "OA_selfrep"],
                 S.loc["OA_knee", "OA_selfrep"], P.loc["OA_knee", "OA_selfrep"],
                 "OA×OA"))
    rows.append(("any fracture", R.loc["OA_knee", "fracture"],
                 S.loc["OA_knee", "fracture"], P.loc["OA_knee", "fracture"],
                 "OA×fracture"))

    fig, ax = plt.subplots(figsize=(9.5, 11))
    colors = {"fracture×BMD": "#1f4e79", "OA_knee×BMD": "#b22222",
              "OA_selfrep×BMD": "#e07b39", "OA×OA": "#2e7d32",
              "OA×fracture": "#6a4c93"}
    y, ypos, ylab, dividers = [], 0, [], []
    seen = set()
    for lab, r, s, p, g in rows:
        if g not in seen:
            if seen:
                ypos -= 1.1
                dividers.append(ypos + 0.55)
            seen.add(g)
            ypos -= 1.2
            ax.text(0.02, ypos, g, fontsize=10.5, fontweight="bold",
                    transform=ax.get_yaxis_transform())
            ypos -= 1.0
        y.append(ypos)
        ylab.append(lab)
        lo, hi = r - 1.96 * s, r + 1.96 * s
        ax.plot([lo, hi], [ypos, ypos], color=colors[g], lw=1.6, alpha=0.85)
        ax.plot([lo, lo], [ypos - 0.18, ypos + 0.18], color=colors[g], lw=1.4)
        ax.plot([hi, hi], [ypos - 0.18, ypos + 0.18], color=colors[g], lw=1.4)
        ax.plot(r, ypos, "o", color=colors[g], ms=5.5, zorder=3)
        ax.text(1.005, ypos, f"{r:+.3f} ({s:.3f})", va="center", fontsize=8.2,
                transform=ax.get_yaxis_transform())
        ypos -= 1.0

    ax.axvline(0, color="black", lw=1.0, ls="--")
    ax.set_yticks(y); ax.set_yticklabels(ylab, fontsize=8.5)
    for d in dividers:
        ax.axhline(d, color="grey", lw=0.6, alpha=0.5)
    ax.set_xlabel("genetic correlation  r$_g$  (95% CI)")
    ax.set_title("LDSC genetic correlations with BMD sites, OA and fracture",
                 fontsize=11.5)
    ax.set_ylim(min(y) - 1.5, max(y) + 1.5)
    plt.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def main():
    os.makedirs(OUT, exist_ok=True)
    R, S, P = build_13x13()
    heatmap(R, P, os.path.join(OUT, "Figure1_rg_heatmap_13.png"))
    forest(R, S, P, os.path.join(OUT, "Figure2_forest_key_results.png"))
    print(f"✅ 图件已输出至 {OUT}")


if __name__ == "__main__":
    main()
