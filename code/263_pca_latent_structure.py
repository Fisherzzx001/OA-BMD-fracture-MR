#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
263_pca_latent_structure.py  —  T1-A
=============================================================================
对 13 表型遗传相关矩阵 (rg_matrix_13.csv) 做谱分解 (相关矩阵 PCA),
把 SOS ("骨骼-关节对立 set-point") 从隐喻升级为数据衍生的潜结构:
  PC1 = BMD↔fracture 轴 (BMD 部位同正, fracture 负, OA 近 0)
  PC2 = OA 轴 (OA_knee/OA_selfrep 主导, BMD 近 0)
输出: 特征值/%方差, 载荷矩阵 (PC1/PC2/PC3), 补充图 S22 (双轴载荷散点)。
诚实标注: OA h² 低 (~0.035 liability), 含 OA 的 r_g SE 偏大 → 噪声, 本分析为探索性探针。
=============================================================================
"""
import os, csv, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

GEN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "genetics")
SRC = os.path.join(GEN, "rg_matrix_13.csv")
FIGDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "Figures_v18_Suppl")
os.makedirs(FIGDIR, exist_ok=True)


def load_mat(path):
    rows = [r for r in csv.reader(open(path)) if r and any(c.strip() for c in r)]
    hdr = [h.strip() for h in rows[0]]
    labels = [h for h in hdr if h]          # 去掉首列空行标签
    M = np.array([[float(x) for x in r[1:]] for r in rows[1:]], float)
    return labels, M


def main():
    labels, M = load_mat(SRC)
    n = len(labels)
    # 相关矩阵谱分解
    ev, evec = np.linalg.eigh(M)
    order = np.argsort(ev)[::-1]
    ev = ev[order]; evec = evec[:, order]
    # 定向: 使 PC1 上 BMD 簇均值载荷为正 (叙事: 高骨量 set-point)
    bmd = set(range(3, n))
    if np.mean([evec[j, 0] for j in bmd]) < 0:
        evec[:, 0] = -evec[:, 0]
    if evec[0, 1] < 0:   # PC2: OA_knee 为正
        evec[:, 1] = -evec[:, 1]
    tot = ev.sum()
    pct = 100.0 * ev / tot
    print("=== 特征值 / 累计方差 ===")
    cum = 0.0
    for i in range(n):
        cum += pct[i]
        print("  PC%d  λ=%.4f  %.2f%%  (累计 %.2f%%)" % (i + 1, ev[i], pct[i], cum))

    # 载荷矩阵 (标准化: 相关矩阵 PCA, 载荷=特征向量)
    load = evec  # 列=PC
    print("\n=== 载荷 (PC1, PC2, PC3) ===")
    header = "%-16s %8s %8s %8s" % ("phenotype", "PC1", "PC2", "PC3")
    print(header)
    for j, lab in enumerate(labels):
        print("%-16s %+8.3f %+8.3f %+8.3f" % (lab, load[j, 0], load[j, 1], load[j, 2]))

    # 分组
    oa = {0, 1}
    fx = {2}
    bmd = set(range(3, n))

    # 补充图 S22: PC1 (x) vs PC2 (y) 载荷散点
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    colors = {"OA": "#D55E00", "fracture": "#0072B2", "BMD": "#009E73"}
    for grp, idx in (("OA", oa), ("fracture", fx), ("BMD", bmd)):
        xs = [load[j, 0] for j in idx]
        ys = [load[j, 1] for j in idx]
        ax.scatter(xs, ys, s=110, color=colors[grp], edgecolor="black",
                   linewidth=0.6, zorder=3, label=grp)
    for j, lab in enumerate(labels):
        ax.annotate(lab.replace("BMD_", ""), (load[j, 0], load[j, 1]),
                    textcoords="offset points", xytext=(5, 4), fontsize=8)
    ax.axhline(0, color="#999999", lw=0.8, zorder=1)
    ax.axvline(0, color="#999999", lw=0.8, zorder=1)
    ax.set_xlabel("PC1 loading  (BMD \u2194 fracture axis)", fontsize=10)
    ax.set_ylabel("PC2 loading  (OA axis)", fontsize=10)
    ax.set_title("Supplementary Figure S22. Latent structure of the 13-phenotype "
                 "genetic-correlation matrix\nPCA of the 13\u00d713 LDSC r_g matrix", fontsize=10)
    ax.legend(loc="upper right", fontsize=9, framealpha=0.9)
    ax.grid(True, ls=":", alpha=0.4)
    fig.tight_layout()
    outfig = os.path.join(FIGDIR, "S22_dual_axis_latent_structure.png")
    fig.savefig(outfig, dpi=300)
    fig.savefig(outfig.replace(".png", ".pdf"), dpi=300)
    print("\nwritten", outfig)

    # 落盘载荷 CSV (供稿件引用)
    outcsv = os.path.join(GEN, "pca_loadings_13.csv")
    with open(outcsv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["phenotype", "PC1", "PC2", "PC3", "group"])
        for j, lab in enumerate(labels):
            grp = "OA" if j in oa else ("fracture" if j in fx else "BMD")
            w.writerow([lab, "%.4f" % load[j, 0], "%.4f" % load[j, 1],
                        "%.4f" % load[j, 2], grp])
    print("written", outcsv)

    # 关键读数回显 (供稿件直引)
    print("\n=== 关键读数 (稿件直引用) ===")
    print("PC1 方差解释 = %.1f%%" % pct[0])
    print("PC2 方差解释 = %.1f%%" % pct[1])
    print("PC1 上 BMD 簇均值载荷 = %.3f" % np.mean([load[j, 0] for j in bmd]))
    print("PC1 上 fracture 载荷 = %.3f" % load[2, 0])
    print("PC1 上 OA_knee 载荷 = %.3f" % load[0, 0])
    print("PC2 上 OA_knee 载荷 = %.3f / OA_selfrep 载荷 = %.3f" %
          (load[0, 1], load[1, 1]))
    print("PC2 上 BMD 簇均值载荷 = %.3f" % np.mean([load[j, 1] for j in bmd]))


if __name__ == "__main__":
    main()
