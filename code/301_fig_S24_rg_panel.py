# -*- coding: utf-8 -*-
"""
Supplementary Figure S24 - Site-stratified genetic-correlation panel.
Visualises the data in Supplementary Table S30 (S30a/BMD-OA r_g,
S30b/OA site-to-site r_g, S30c/SNP h^2). All numbers taken verbatim
from code/output/genetics/ldsc_rg_verified/verified_results.csv and the
manuscript Table S30.
Outputs to Figures_v18_Suppl/, Figures/, submission_package/Figures/.
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

try:
    import scienceplots  # noqa
    plt.style.use(["science", "nature"])
except Exception:
    pass

# LaTeX is not installed on this host; render math with matplotlib's
# built-in mathtext instead of usetex so the figure builds reproducibly.
plt.rcParams.update({
    "text.usetex": False,
    "mathtext.fontset": "dejavusans",
    "font.family": "sans-serif",
    "font.size": 9,
    "axes.labelsize": 9,
    "axes.titlesize": 10,
    "axes.linewidth": 0.8,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "legend.fontsize": 8.5,
    "figure.dpi": 300,
})

TEAL = "#2a7f9e"
TEAL_D = "#1d5e76"
ORANGE = "#d9822b"
GREY = "#9aa0a6"
INK = "#222222"

# ---- data (verified_results.csv, verified_N31986 + OA_site_vs_site_baseline) ----
# Panel A: BMD(total) x OA site r_g (S30a)
sites_a = ["Hip", "Knee", "Spine", "Hand"]
rg_a = [0.1309, 0.1523, 0.1957, 0.0713]
se_a = [0.0378, 0.0328, 0.0468, 0.0561]
p_a = [5.32e-4, 3.41e-6, 2.86e-5, 0.204]
sig_a = [True, True, True, False]

# Panel B: OA site-to-site r_g, knee reference (S30b)
pairs_b = ["Knee\u2013hip", "Knee\u2013spine", "Knee\u2013hand"]
rg_b = [0.5692, 0.6200, 0.5905]
se_b = [0.0191, 0.0233, 0.0392]

# Panel C: SNP h^2 (S30c)
traits_c = ["BMD\n(total)", "OA hip", "OA knee", "OA spine", "OA hand"]
h2_c = [0.3479, 0.0283, 0.0326, 0.0131, 0.0087]


def stars(p):
    if p < 1e-3:
        return "***"
    if p < 1e-2:
        return "**"
    if p < 0.05:
        return "*"
    return "n.s."


def build():
    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(11.2, 3.7), dpi=300)
    fig.subplots_adjust(left=0.06, right=0.985, top=0.80, bottom=0.16, wspace=0.32)

    # ---------- Panel A ----------
    x = range(len(sites_a))
    colors = [TEAL if s else GREY for s in sig_a]
    bars = axA.bar(x, rg_a, color=colors, width=0.62, zorder=3,
                   edgecolor="white", linewidth=0.6)
    for xi, (v, e) in enumerate(zip(rg_a, se_a)):
        axA.errorbar(xi, v, yerr=1.96 * e, fmt="none", ecolor=INK,
                     elinewidth=1.0, capsize=3.0, zorder=4)
    for xi, (v, s) in enumerate(zip(rg_a, sig_a)):
        lab = stars(p_a[xi]) if s else "n.s."
        axA.text(xi, v + 1.96 * se_a[xi] + 0.006, lab, ha="center",
                 va="bottom", fontsize=8.5, fontweight="bold", color=INK)
    axA.axhline(0, color=INK, linewidth=0.8, zorder=2)
    axA.set_xticks(list(x))
    axA.set_xticklabels(sites_a, fontsize=9)
    axA.set_ylim(0, 0.32)
    axA.set_ylabel(r"Genetic correlation $r_g$ (BMD $\times$ OA site)", fontsize=9)
    axA.set_title("A. BMD\u2013OA by joint site", fontsize=10, fontweight="bold",
                  loc="left", pad=8)
    axA.tick_params(labelsize=8.5)

    # ---------- Panel B ----------
    xb = range(len(pairs_b))
    axB.bar(xb, rg_b, color=TEAL_D, width=0.62, zorder=3,
            edgecolor="white", linewidth=0.6)
    for xi, (v, e) in enumerate(zip(rg_b, se_b)):
        axB.errorbar(xi, v, yerr=1.96 * e, fmt="none", ecolor=INK,
                     elinewidth=1.0, capsize=3.0, zorder=4)
        axB.text(xi, v + 1.96 * e + 0.008, f"{v:.3f}", ha="center",
                 va="bottom", fontsize=8.5, color=INK)
    axB.set_xticks(list(xb))
    axB.set_xticklabels(pairs_b, fontsize=9)
    axB.set_ylim(0, 0.72)
    axB.set_ylabel(r"Genetic correlation $r_g$ (OA site pairs)", fontsize=9)
    axB.set_title("B. OA site-to-site", fontsize=10, fontweight="bold",
                  loc="left", pad=8)
    axB.tick_params(labelsize=8.5)

    # ---------- Panel C ----------
    xc = range(len(traits_c))
    axC.bar(xc, h2_c, color=[TEAL] + [ORANGE] * 4, width=0.62, zorder=3,
            edgecolor="white", linewidth=0.6, log=True)
    for xi, v in enumerate(h2_c):
        axC.text(xi, v * 1.18, f"{v:.3f}", ha="center", va="bottom",
                 fontsize=8.0, color=INK)
    axC.set_xticks(list(xc))
    axC.set_xticklabels(traits_c, fontsize=8.5)
    axC.set_ylim(0.005, 0.6)
    axC.set_ylabel(r"SNP heritability $h^2$ (observed)", fontsize=9)
    axC.set_title("C. SNP heritability", fontsize=10, fontweight="bold",
                  loc="left", pad=8)
    axC.tick_params(labelsize=8.5)

    # ---------- figure caption / panel letters spacing ----------
    fig.suptitle("", fontsize=1)
    return fig


def save_all(fig):
    out_dirs = [
        "Figures_v18_Suppl",
        "Figures",
        "submission_package/Figures",
    ]
    name = "SupplFig24_rg_panel.png"
    for d in out_dirs:
        os.makedirs(d, exist_ok=True)
        p = os.path.join(d, name)
        fig.savefig(p, dpi=300, bbox_inches="tight")
        print("wrote", p)
    plt.close(fig)


if __name__ == "__main__":
    save_all(build())
    print("DONE")
