#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""335 · Candidate supplementary figures S33–S35 (revision assets).

Context (2026-11-11 full-text audit / figure-expansion plan)
------------------------------------------------------------
Three core 2026-11 results had NO graphical carrier:
  G1  the SOS discriminative kernel -- the SAME 31 instruments pointing in
      opposite directions at the joint (OA) vs whole-bone (fracture) endpoint;
  G6  r_g vs MR beta triangulation (genetic correlation is not causation);
  G5  cell-type enrichment, the cellular anchor of the formation-side SOS.

These are built here as REVISION-CANDIDATE supplementary figures, numbered
S33–S35 to avoid collision with the existing Figure S1–S24 and Table S25–S32.
They are NOT inserted into the submitted main-text Figure-legend list: OAC
caps displayed items at 8 (already 8/8), so adding main figures would breach
the format; the candidates are held for the revision / response-to-reviewers.

Hard rules
----------
* Every plotted value is read from an existing results file under
  code/output/genetics/. No statistic is hard-coded or invented.
* Each figure emits a plotted_data_SNN.json recording the exact numbers drawn,
  so the figure is independently auditable (same discipline as build_figures_v17).
* SOS sign kernel is drawn with the verified convention: a BMD-increasing
  allele raises OA (BMD vs OA SAME sign) and lowers fracture (BMD vs FX
  OPPOSITE sign). The count 21/31 is taken verbatim from three_trait_coloc.json.
* No change to any reported beta/SE/p/OR/CI; no change to the 9/10 count.

Outputs: Figures_v18_Suppl/SupplFig{33,34,35}_*.{png,pdf}
         Figures_v18_Suppl/plotted_data_S33-S35.json
"""
from __future__ import annotations

import json
import os
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    plt.style.use(["science", "nature", "no-latex"])
except Exception:
    pass

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8, "axes.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.labelsize": 8, "axes.titlesize": 9,
    "xtick.labelsize": 7, "ytick.labelsize": 7,
    "legend.fontsize": 7, "legend.frameon": False,
    "figure.dpi": 300, "savefig.dpi": 600,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.04,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
G = os.path.join(ROOT, "code", "output", "genetics")
OUT = os.path.join(ROOT, "Figures_v18_Suppl")
os.makedirs(OUT, exist_ok=True)

# Okabe–Ito colour-blind-safe palette
OI = {"blue": "#0072B2", "orange": "#E69F00", "green": "#009E73",
      "red": "#D55E00", "purple": "#CC79A7", "grey": "#999999",
      "sky": "#56B4E9", "black": "#000000"}

trace = {}


def loadj(name):
    return json.load(open(os.path.join(G, name), encoding="utf-8"))


def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=600)
    fig.savefig(os.path.join(OUT, name + ".pdf"))
    plt.close(fig)
    print(f"  wrote {name}.png / .pdf")


# ------------------------------------------------------------------ S33
def fig_s33():
    import csv
    rg = list(csv.DictReader(open(os.path.join(G, "rg_results_long_13.csv"), encoding="utf-8")))
    mr = list(csv.DictReader(open(os.path.join(G, "mr_results_corrected.csv"), encoding="utf-8")))

    def rg_row(a, b):
        for r in rg:
            if {r["trait1"], r["trait2"]} == {a, b}:
                return r
        raise KeyError((a, b))

    def mr_row(e, o):
        for r in mr:
            if r["exposure"] == e and r["outcome"] == o:
                return r
        raise KeyError((e, o))

    # Unicode subscripts / arrows are missing from the Arial metrics shipped
    # with matplotlib (render as tofu). Use mathtext for r_g and plain ASCII
    # for relation arrows so the vector PDF stays clean and portable.
    rels = [("BMD_total", "OA_knee", "BMD vs knee OA"),
            ("BMD_total", "fracture", "BMD vs any fracture"),
            ("OA_knee", "fracture", "knee OA vs fracture")]
    rg_pts, mr_pts = [], []
    for a, b, lab in rels:
        r = rg_row(a, b); rg_pts.append((lab, float(r["rg"]), float(r["se"]), float(r["p"])))
    mr_spec = [("BMD_total", "OA_knee", "BMD -> knee OA"),
               ("BMD_total", "fracture", "BMD -> fracture"),
               ("OA_knee", "fracture", "knee OA -> fracture")]
    for e, o, lab in mr_spec:
        r = mr_row(e, o); mr_pts.append((lab, float(r["ivw"]), float(r["ivw_se"]), float(r["ivw_p"])))

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.5))
    for ax, pts, title, color in [(axes[0], rg_pts, r"a  LDSC genetic correlation ($r_g$)", OI["blue"]),
                                  (axes[1], mr_pts, r"b  Two-sample MR (IVW $\beta$)", OI["orange"])]:
        y = np.arange(len(pts))[::-1]
        for yi, (lab, est, se, p) in zip(y, pts):
            ax.errorbar(est, yi, xerr=1.96 * se, fmt="o", color=color, capsize=3, ms=5)
            ptxt = f"p={p:.1e}"
            ax.annotate(ptxt, (est + 1.96 * se, yi), xytext=(4, -2),
                        textcoords="offset points", fontsize=6, va="center")
        ax.axvline(0, color=OI["grey"], lw=0.8, ls="--")
        ax.set_yticks(y); ax.set_yticklabels([p[0] for p in pts]); ax.set_title(title, loc="left", fontweight="bold")
        ax.set_xlabel("estimate (95% CI)")
    fig.suptitle("Genetic correlation is not causation", fontweight="bold", fontsize=10, x=0.02, ha="left")
    fig.text(0.5, 0.005,
             r"OA-fracture $r_g$ is significant while OA->fracture MR $\beta$ is null: "
             "shared polygenic background is not a causal effect.",
             ha="center", fontsize=6.5, color="#333333")
    fig.tight_layout(rect=(0, 0.05, 1, 0.95))
    save(fig, "SupplFig33_rg_vs_mr_triangulation")
    trace["S33"] = {"source": ["rg_results_long_13.csv", "mr_results_corrected.csv"],
                    "rg": [{"pair": l, "rg": e, "se": s, "p": p} for l, e, s, p in rg_pts],
                    "mr": [{"pair": l, "beta": e, "se": s, "p": p} for l, e, s, p in mr_pts]}


# ------------------------------------------------------------------ S34
def fig_s34():
    d = loadj("three_trait_coloc.json")
    loc = d["loci"]
    bmd = np.array([L["top_bmd_z"] for L in loc], float)
    oa = np.array([L["top_oa_z_at_bmd"] for L in loc], float)
    fx = np.array([L["top_fx_z_at_bmd"] for L in loc], float)
    sos = np.array([bool(L["sos_sign_pattern"]) for L in loc])
    h123 = np.array([L["H123"] for L in loc], float)

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))
    # (a) signed Z of the same instrument across the two endpoints
    ax = axes[0]
    order = np.argsort(bmd)
    xp = np.arange(len(loc))
    ax.axhline(0, color=OI["grey"], lw=0.8)
    ax.scatter(xp, oa[order], s=14, c=np.where(sos[order], OI["blue"], OI["grey"]),
               label="OA (joint)", zorder=3)
    ax.scatter(xp, fx[order], s=14, marker="^", c=np.where(sos[order], OI["red"], OI["grey"]),
               label="fracture (whole bone)", zorder=3)
    ax.plot(xp, bmd[order], color=OI["black"], lw=1.0, label="BMD (exposure)")
    ax.set_xlabel("31 BMD instrument loci (sorted by BMD Z)")
    ax.set_ylabel("lead-variant Z on BMD-increasing allele")
    ax.set_title("a  One variant set, opposite endpoint signs", loc="left", fontweight="bold")
    ax.legend(loc="upper left", fontsize=6)
    # (b) H123 posterior at loci meeting the 3-way threshold
    ax = axes[1]
    idx = np.where(h123 >= 0.8)[0]
    idx = idx[np.argsort(-h123[idx])]
    labs = [f"chr{L['chr']}:{L['bp']//1000}k" for L in [loc[i] for i in idx]]
    ax.barh(range(len(idx))[::-1], h123[idx], color=OI["green"])
    ax.axvline(0.8, color=OI["grey"], lw=0.8, ls="--")
    ax.set_yticks(range(len(idx))[::-1]); ax.set_yticklabels(labs)
    ax.set_xlim(0, 1); ax.set_xlabel("P(H123) shared across BMD + OA + fracture")
    ax.set_title(f"b  Three-way coloc: {len(idx)}/31 at H123≥0.8", loc="left", fontweight="bold")
    for yi, v in zip(range(len(idx))[::-1], h123[idx]):
        ax.annotate(f"{v:.2f}", (v, yi), xytext=(4, -2), textcoords="offset points", fontsize=6)
    fig.suptitle("SOS discriminative kernel: sign opposition at matched loci", fontweight="bold", fontsize=10, x=0.02, ha="left")
    fig.text(0.5, 0.005,
             f"SOS sign pattern {d['n_loci_SOS_sign_pattern']}/{d['n_loci']} genome-wide; "
             f"{len(idx)}/31 three-way-colocalised loci, all sign-concordant. Coloc does NOT replace the 0/31 two-locus result.",
             ha="center", fontsize=6.3, color="#333333")
    fig.tight_layout(rect=(0, 0.06, 1, 0.94))
    save(fig, "SupplFig34_SOS_sign_opposition")
    trace["S34"] = {"source": "three_trait_coloc.json", "n_loci": d["n_loci"],
                    "sos_sign_pattern": d["n_loci_SOS_sign_pattern"],
                    "n_three_way_H4_ge_0p8": d["n_three_way_H4_ge_0p8"],
                    "loci": [{"chr": L["chr"], "bp": L["bp"], "Z_BMD": L["top_bmd_z"],
                              "Z_OA": L["top_oa_z_at_bmd"], "Z_FX": L["top_fx_z_at_bmd"],
                              "sos": L["sos_sign_pattern"], "H123": L["H123"]} for L in loc]}


# ------------------------------------------------------------------ S35
def fig_s35():
    d = loadj("Chai2026_celltype_enrichment.json")
    rows = d["rows"]
    key = {k["cluster"]: k for k in d["key_skeletal_cells"]}
    fig, ax = plt.subplots(figsize=(5.4, 3.6))
    for r in rows:
        is_key = r["cluster"] in key
        sig = r["q"] < 0.05
        col = OI["red"] if (is_key and sig) else (OI["blue"] if sig else OI["grey"])
        ax.scatter(r["fold"], -np.log10(max(r["q"], 1e-300)), s=42 if is_key else 18,
                   color=col, marker="o" if is_key else "s", zorder=4 if is_key else 2,
                   edgecolor="black" if is_key else "none", linewidth=0.5 if is_key else 0)
    ax.axhline(-np.log10(0.05), color=OI["grey"], ls="--", lw=0.8)
    ax.axvline(1.0, color=OI["grey"], ls=":", lw=0.8)
    # Annotate only the two significant skeletal cells directly. The three
    # non-significant key cells (q~1.0, -log10 q~0) sit on top of each other;
    # listing them in a boxed note keeps the plot legible. Labelling all 34
    # clusters would pile ~30 q~1.0 points into an unreadable smear.
    for k in d["key_skeletal_cells"]:
        if k["q"] < 0.05:
            ax.annotate(f"{k['cell']}\n({k['fold']}x, q={k['q']:.1e})",
                        (k["fold"], -np.log10(max(k["q"], 1e-300))),
                        xytext=(6, 2), textcoords="offset points", fontsize=5.8)
    nonsig = [k for k in d["key_skeletal_cells"] if k["q"] >= 0.05]
    note = "Not enriched (q>=0.05): " + ", ".join(
        f"{k['cell']} {k['fold']}x" for k in nonsig)
    ax.text(0.02, 0.62, note, transform=ax.transAxes, fontsize=5.8,
            color="#444444", va="top",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=OI["grey"], alpha=0.9))
    ax.set_xlabel("enrichment fold (observed / expected)")
    ax.set_ylabel(r"$-\log_{10}(q)$")
    ax.set_title("Skeletal-cell enrichment of OA effector genes", loc="left", fontweight="bold")
    fig.text(0.5, 0.005,
             "Osteoblast and chondrocyte clusters are enriched (q<0.05); osteoclast / endothelial / SMC are not — "
             "a formation-side, not resorption-side, cellular anchor for SOS.",
             ha="center", fontsize=6.2, color="#333333")
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    save(fig, "SupplFig35_celltype_enrichment")
    trace["S35"] = {"source": "Chai2026_celltype_enrichment.json",
                    "n_clusters": d["n_clusters"], "n_clusters_q_lt_0p05": d["n_clusters_q_lt_0p05"],
                    "key_skeletal_cells": d["key_skeletal_cells"],
                    "rows": rows}


def main():
    print("=== building candidate supplementary figures S33–S35 ===")
    fig_s33(); fig_s34(); fig_s35()
    tp = os.path.join(OUT, "plotted_data_S33-S35.json")
    json.dump(trace, open(tp, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print(f"wrote {os.path.relpath(tp, ROOT)}")
    print("DONE")


if __name__ == "__main__":
    main()
