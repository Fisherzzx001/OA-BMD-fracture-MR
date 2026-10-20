# -*- coding: utf-8 -*-
"""
112 · Supplementary Figure S9 —— 工具变量筛选与流失流程（STROBE-MR 合规）

所有计数均从一手文件现场统计并断言，不沿用任何手抄数字：
  code/output/genetics/strict_clump_summary.csv   n_gws（全基因组显著 SNP 数）
  code/output/genetics/*_plink.clumped            严格 clump 后工具数（PLINK 基准版）
  code/output/genetics/harmonised/*.csv           与结局 harmonise 后进入分析的 K

勘误记录（2026-10-16）：PLINK .clumped 文件尾部含 2 个空行，此前用 `wc -l - 1`
误计为 22–37/118；真实数据行为 20–35/116（与 strict_clump_summary 一致），
稿件 Methods 已同步勘误（code/113_fix_instrument_counts.py）。
"""
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle_v18 import C, GEN, HARM, save, dump_plotted, SAVED

EXPOSURES = {
    "BMD_total": "BMD total", "BMD_femur_neck": "Femoral neck",
    "BMD_lumbar": "Lumbar spine", "BMD_arms": "Arms", "BMD_leg": "Legs",
    "OA_knee": "Knee OA",
}
PAIRS = [
    ("BMD_total", "knee OA", "BMD_total__OA_knee.csv"),
    ("BMD_total", "fracture", "BMD_total__fracture.csv"),
    ("BMD_femur_neck", "knee OA", "BMD_femur_neck__OA_knee.csv"),
    ("BMD_femur_neck", "fracture", "BMD_femur_neck__fracture.csv"),
    ("BMD_lumbar", "knee OA", "BMD_lumbar__OA_knee.csv"),
    ("BMD_arms", "knee OA", "BMD_arms__OA_knee.csv"),
    ("BMD_leg", "knee OA", "BMD_leg__OA_knee.csv"),
    ("OA_knee", "BMD total", "OA_knee__BMD_total.csv"),
    ("OA_knee", "femoral neck", "OA_knee__BMD_femur_neck.csv"),
    ("OA_knee", "fracture", "OA_knee__fracture.csv"),
]

TABLE3_K = {
    "BMD_total__OA_knee.csv": 31, "BMD_total__fracture.csv": 31,
    "BMD_femur_neck__OA_knee.csv": 25, "BMD_femur_neck__fracture.csv": 25,
    "BMD_lumbar__OA_knee.csv": 21, "BMD_arms__OA_knee.csv": 16,
    "BMD_leg__OA_knee.csv": 24, "OA_knee__BMD_total.csv": 69,
    "OA_knee__BMD_femur_neck.csv": 69, "OA_knee__fracture.csv": 69,
}


def count_rows(path):
    n = 0
    with open(path, encoding="utf-8") as f:
        for i, line in enumerate(f):
            if i > 0 and line.strip():
                n += 1
    return n


def build_stats():
    sc = pd.read_csv(os.path.join(GEN, "strict_clump_summary.csv")).set_index("exposure")
    stats = {}
    for exp in EXPOSURES:
        n_gws = int(sc.loc[exp, "n_gws"])
        n_clumped = count_rows(os.path.join(GEN, exp + "_plink.clumped"))
        # 与 python 端 strict 管道交叉验证（已知差异：BMD_lumbar 24 vs 23，
        # 1 个边界 LD SNP —— 见 docs/12_补充材料_PLINK_clump_benchmark.md）
        n_py = int(sc.loc[exp, "n_kept_strict"])
        assert abs(n_clumped - n_py) <= 1, \
            f"{exp}: clumped rows {n_clumped} vs strict_clump_summary {n_py}"
        if n_clumped != n_py:
            print(f"   [note] {exp}: PLINK clump {n_clumped} vs python strict {n_py} "
                  f"(documented 1-SNP boundary-LD difference)")
        stats[exp] = {"n_gws": n_gws, "n_clumped": n_clumped, "pairs": {}}
    for exp, out, fn in PAIRS:
        n = count_rows(os.path.join(HARM, fn))
        assert n == TABLE3_K[fn], f"{fn}: harmonised K={n} vs Table 3 {TABLE3_K[fn]}"
        stats[exp]["pairs"][f"{exp}->{out}"] = {"file": fn, "K": n}
    bmd_clumped = [stats[e]["n_clumped"] for e in EXPOSURES if e != "OA_knee"]
    assert min(bmd_clumped) == 20 and max(bmd_clumped) == 35, bmd_clumped
    assert stats["OA_knee"]["n_clumped"] == 116
    fwd = [p["K"] for e in EXPOSURES if e != "OA_knee" for p in stats[e]["pairs"].values()]
    assert min(fwd) == 16 and max(fwd) == 31, fwd
    return stats


def box(ax, x, y, w, h, text, fc="#FFFFFF", fs=6.6):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.004,rounding_size=0.008",
                                fc=fc, ec=C["ink"], lw=1.0, zorder=2))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
            fontsize=fs, zorder=3, linespacing=1.45)


def fig_S9():
    st = build_stats()
    fig, ax = plt.subplots(figsize=(7.20, 5.1))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(0.030, 0.972,
            "Instrument selection followed STROBE-MR; counts are tallied from the "
            "deposited PLINK clump files and harmonised instrument tables.",
            fontsize=6.0, ha="left", va="center", color="#444444", style="italic")

    # ============ 左列：五步流程 ============
    Lx, Lw = 0.030, 0.455
    hs = [0.135, 0.095, 0.115, 0.125, 0.105]
    gap = 0.045
    y = 0.915
    boxes = [
        ("Step 1 \u00b7 Exposure GWAS summary statistics\n"
         "5 regional BMD sites (Qian 2025; N=31,873 per site)\n"
         "+ radiographic knee OA (GO Consortium 2025; N\u22481.21M)", "#EFF5FB"),
        ("Step 2 \u00b7 Genome-wide significant SNPs ($p<5\\times10^{-8}$)\n"
         "170\u2013514 per BMD site;  1,980 knee OA", "#FFFFFF"),
        ("Step 3 \u00b7 Independent instruments after strict LD clumping\n"
         "r\u00b2<0.001 / 10 Mb (1000G EUR; PLINK-benchmarked)\n"
         "20\u201335 per BMD site;  116 knee OA", "#FFFFFF"),
        ("Step 4 \u00b7 Harmonised with each outcome GWAS\n"
         "drop palindromic / strand-ambiguous / absent-from-outcome\n"
         "16\u201331 (BMD\u2192knee OA);  31, 25 (BMD\u2192fracture);\n"
         "69 (knee OA as exposure; \u221241% from 116)", "#FDF3E3"),
        ("Step 5 \u00b7 10 MR pairs analysed\n"
         "IVW fixed/random, MR-Egger, weighted median,\n"
         "MR-PRESSO, leave-one-out; MVMR + height on 7 pairs", "#EDF5F0"),
    ]
    centers = []
    for (txt, fc), h in zip(boxes, hs):
        box(ax, Lx, y - h, Lw, h, txt, fc=fc)
        centers.append(y - h / 2)
        y -= h + gap
    for i in range(4):
        ax.add_patch(FancyArrowPatch((Lx + Lw / 2, centers[i] - hs[i] / 2 - 0.004),
                                     (Lx + Lw / 2, centers[i + 1] + hs[i + 1] / 2 + 0.004),
                                     arrowstyle="-|>", mutation_scale=9, lw=1.0,
                                     color=C["ink"], zorder=1))

    # ============ 右列：逐对流失表 ============
    Tx = 0.545
    ax.text(Tx, 0.945, "Instrument attrition by analysed pair", fontsize=7.4,
            fontweight="bold", ha="left", va="center")
    col_x = {"gws": Tx + 0.290, "clump": Tx + 0.335, "K": Tx + 0.385, "ret": Tx + 0.432}
    hy = 0.905
    ax.text(Tx, hy, "Exposure \u2192 outcome", fontsize=6.2, ha="left", va="center",
            color="#444444")
    ax.text(col_x["gws"], hy, "GWS", fontsize=6.2, ha="right", va="center",
            color="#444444")
    ax.text(col_x["clump"], hy, "Clumped", fontsize=6.2, ha="right", va="center",
            color="#444444")
    ax.text(col_x["K"], hy, "K", fontsize=6.2, ha="right", va="center",
            color="#444444")
    ax.text(col_x["ret"], hy, "Retained", fontsize=6.2, ha="right", va="center",
            color="#444444")

    row_y0, dy = 0.822, 0.0745
    rows = []
    for exp, out, fn in PAIRS:
        K = st[exp]["pairs"][f"{exp}->{out}"]["K"]
        rows.append((EXPOSURES[exp], out, st[exp]["n_gws"], st[exp]["n_clumped"], K))
    for i, (exp_l, out, ng, nc, K) in enumerate(rows):
        yr = row_y0 - i * dy
        if i in (2, 7):
            ax.plot([Tx - 0.008, Tx + 0.428], [yr + dy * 0.60, yr + dy * 0.60], lw=0.6,
                    color=C["grey"], ls=(0, (3, 2)))
        ax.text(Tx, yr, f"{exp_l} \u2192 {out}", fontsize=6.3, ha="left", va="center")
        ax.text(col_x["gws"], yr, f"{ng:,}", fontsize=6.3, ha="right", va="center")
        ax.text(col_x["clump"], yr, f"{nc}", fontsize=6.3, ha="right", va="center")
        ax.text(col_x["K"], yr, f"{K}", fontsize=6.3, ha="right", va="center",
                fontweight="bold")
        ret = 100.0 * K / nc
        col = C["vermillion"] if ret < 75 else "#3A3A3A"
        ax.text(col_x["ret"], yr, f"{ret:.0f}%", fontsize=6.3, ha="right",
                va="center", color=col)
    ax.text(Tx, row_y0 - 10 * dy - 0.006,
            "GWS = genome-wide significant SNPs ($p<5\\times10^{-8}$); Clumped = independent "
            "instruments after strict LD clumping; K = instruments entering MR after "
            "harmonisation; Retained = K / clumped.\nKnee OA as exposure loses 41% of "
            "instruments at harmonisation (116 \u2192 69); the attrition is reported in "
            "Methods and examined in the sensitivity analyses.",
            fontsize=5.9, ha="left", va="top", color="#444444")

    # 中缝分隔线
    ax.plot([0.516, 0.516], [0.06, 0.945], lw=0.6, color="#DDDDDD")

    save(fig, "SupplFig9_instrument_flow")
    dump_plotted({"SupplFig9_instrument_flow": {
        "sources": ["strict_clump_summary.csv (n_gws, n_kept_strict)",
                     "*_plink.clumped (strict clump counts; blank trailing lines skipped)",
                     "harmonised/*.csv (analysed K, asserted == Table 3)"],
        "erratum": "22-37/118 (wc -l artifact) corrected to 20-35/116 on 2026-10-16",
        "rows": [{"exposure": e, "outcome": o, "n_gws": ng, "n_clumped": nc, "K": K,
                  "retained_pct": round(100.0 * K / nc, 1)}
                 for (e, o, ng, nc, K) in rows],
    }})


if __name__ == "__main__":
    print("\u25b6 Supplementary Figure S9 — instrument selection flow")
    fig_S9()
    print(f"\n\u2714 {len(SAVED)} figure(s) written")
