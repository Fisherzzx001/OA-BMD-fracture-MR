"""
04_make_outputs.py — 论文图件生成
================================================================================
Figure 1  中英 OA 共病网络并排图（节点大小∝患病率，边宽∝O/E，颜色=病因模块）
Figure 2  OA vs 匹配对照的平均度森林图（含置换检验 P）
Figure 3  模块内平均度：中国 vs 英国（含 bootstrap 95%CI）
================================================================================
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import DKEYS, DLABEL, DMODULE, MODULES, MODULE_COLOR, OUT_DIR

plt.rcParams.update({"font.size": 9, "figure.dpi": 200,
                     "savefig.bbox": "tight", "font.family": "DejaVu Sans"})


# ---------------------------------------------------------------- Figure 1
def figure1_networks(built, k):
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 6.2))
    order = [("charls", "CHARLS (China)"), ("ukb", "UK Biobank (UK)")]
    for ax, (arm, title) in zip(axes, order):
        b = built[arm]
        adj = b["adj_oa"]
        G = nx.from_numpy_array(adj)
        keep = b["keep_oa"]
        prev = b["prev_oa"]

        pos = nx.spring_layout(G, seed=7, k=0.9, iterations=200)
        nodes = [i for i in G.nodes() if keep[i]]
        sub = G.subgraph(nodes)

        sizes = [220 + 4200 * prev[i] for i in sub.nodes()]
        colors = [MODULE_COLOR[DMODULE[DKEYS[i]]] for i in sub.nodes()]

        st = b["st_oa"]
        weights = []
        for u, v in sub.edges():
            weights.append(min(4.5, max(0.6, 1.2 * (st["oe"][u, v] - 1) + 0.8)))

        nx.draw_networkx_edges(sub, pos, ax=ax, width=weights,
                               alpha=0.55, edge_color="#5D6D7E")
        nx.draw_networkx_nodes(sub, pos, ax=ax, nodelist=list(sub.nodes()),
                               node_size=sizes, node_color=colors,
                               edgecolors="white", linewidths=1.2)
        labels = {i: DLABEL[DKEYS[i]] for i in sub.nodes()}
        nx.draw_networkx_labels(sub, pos, labels, ax=ax, font_size=6.2,
                                font_weight="bold")
        ax.set_title(f"{title}\nOsteoarthritis arm · top-{k} excess co-occurrence edges",
                     fontsize=10)
        ax.axis("off")

    handles = [Line2D([0], [0], marker="o", color="w",
                      markerfacecolor=MODULE_COLOR[m], markersize=9,
                      label=m.replace("_", " ")) for m in MODULES]
    fig.legend(handles=handles, loc="lower center", ncol=6, frameon=False,
               fontsize=8, bbox_to_anchor=(0.5, -0.04))
    fig.suptitle("Figure 1. Harmonised comorbidity networks in OA", fontsize=11.5, y=1.0)
    p = os.path.join(OUT_DIR, "Figure1_network_CN_UK.png")
    fig.savefig(p); plt.close(fig)
    print("  →", p)


# ---------------------------------------------------------------- Figure 2
def figure2_forest():
    t = pd.read_csv(os.path.join(OUT_DIR, "table3_network_metrics.csv"))
    h = pd.read_csv(os.path.join(OUT_DIR, "table5_H1_permutation.csv"))
    fig, ax = plt.subplots(figsize=(6.6, 3.1))
    ypos, ytick = [], []
    for i, (_, r) in enumerate(t.iterrows()):
        y = len(t) - 1 - i
        ax.plot(r.average_degree, y, "o", ms=7,
                color="#C0392B" if r.group == "OA" else "#2471A3")
        ax.text(r.average_degree + 0.12, y, f"{r.average_degree:.2f}",
                va="center", fontsize=8)
        ypos.append(y)
        ytick.append(f"{r.cohort} — {r.group}")
    ax.set_yticks(ypos); ax.set_yticklabels(ytick, fontsize=8)
    ax.set_xlabel("Average degree of comorbidity network")
    ax.axvline(0, color="grey", lw=0.6)
    ax.set_title("Figure 2. Network complexity: OA vs matched controls", fontsize=10)
    ax.grid(axis="x", alpha=0.25)
    ps = "   ".join([f"{r.cohort}: P={r.p_permutation}" for _, r in h.iterrows()])
    ax.text(0.0, -0.30, ps, transform=ax.transAxes, fontsize=7.5, color="#555")
    p = os.path.join(OUT_DIR, "Figure2_average_degree_forest.png")
    fig.savefig(p); plt.close(fig)
    print("  →", p)


# ---------------------------------------------------------------- Figure 3
def figure3_modules():
    t = pd.read_csv(os.path.join(OUT_DIR, "table6_H2H4_cross_country.csv"))
    t = t[t.metric.str.startswith("Module")]
    labels = [m.replace("Module: ", "").replace("_", " ") for m in t.metric]
    delta = t.delta_CN_minus_UK.values
    lo = t.ci_low.values - delta
    hi = t.ci_high.values - delta
    colors = ["#C0392B" if d > 0 else "#2471A3" for d in delta]
    y = np.arange(len(t))
    fig, ax = plt.subplots(figsize=(6.8, 3.4))
    ax.barh(y, delta, color=colors, alpha=.85, height=.62)
    ax.errorbar(delta, y, xerr=[np.abs(lo), np.abs(hi)], fmt="none",
                ecolor="#333", elinewidth=1.1, capsize=3)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8.5)
    ax.axvline(0, color="black", lw=1)
    ax.set_xlabel("Δ module aggregation strength (mean O/E)  (CHARLS − UK Biobank)")
    ax.set_title("Figure 3. Cross-country difference by aetiological module", fontsize=10)
    ax.grid(axis="x", alpha=0.25)
    p = os.path.join(OUT_DIR, "Figure3_module_difference.png")
    fig.savefig(p); plt.close(fig)
    print("  →", p)


def main(built, k, built_k=None):
    print("▶ 出图")
    figure1_networks(built_k if built_k is not None else built, k)
    figure2_forest()
    figure3_modules()


if __name__ == "__main__":
    import importlib
    na = importlib.import_module("03_network_analysis")
    built, k, built_k = na.main()
    main(built, k, built_k)
