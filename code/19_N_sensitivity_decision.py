#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
19_N_sensitivity_decision.py
=========================================================================
投稿前整改 #1 收尾：为 **fracture** 与 **OA_selfrep** 建立「N 敏感性决策表」。

背景：LDSC 中 h² 与 N 不可同时识别（h² ∝ 1/N，斜率只给出 N·h²），
      因此 N 必须由**外部来源元数据**确定，本脚本只负责：
        给定候选 N -> 重算 h² / SE / 截距 -> 转 liability 尺度 -> 供定稿选择。

候选来源（均已联网核实）：
  · OA_selfrep —— IEU/OpenGWAS **ukb-b-14486**（MRC-IEU, Ben Elsworth, 2018）
        N = 462,933；cases = 38,472；controls = 424,461；SNPs = 9,851,867
        Build HG19/GRCh37，European。
        （两篇已发表文献 Table 1 交叉印证：Int J Gen Med 2021;14:6777
          doi:10.2147/IJGM.S331864；及 PMC8523904）
  · fracture   —— GEFOS/GENOMOS + 23andMe（Trajanoska et al., BMJ 2018;362:k3225,
        PMID 30158200, doi:10.1136/bmj.k3225）
        发现集 37,857 cases + 227,116 controls = 264,973
        复制集 up to 147,200 cases + 150,085 controls = 297,285
        合计（发现+复制）约 562,258
=========================================================================
"""
import os
import sys
import importlib
import numpy as np
import pandas as pd
from scipy.stats import norm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
m08 = importlib.import_module("08_genetic_correlation")
m09 = importlib.import_module("09_rg_matrix")
m10 = importlib.import_module("10_fracture_row")

OUT = os.path.join(m09.ROOT, "output", "genetics")

# ---- 候选 N ---------------------------------------------------------------
FRACTURE_N = {
    "400000_当前占位": 400_000,
    "264973_GEFOS发现集": 264_973,
    "562258_GEFOS合计": 562_258,
    "462933_UKB若同ukb-b": 462_933,
}
SELFREP_N = {
    "350000_当前文件值": 350_000,
    "462933_ukb-b-14486权威": 462_933,
}
# ukb-b-14486：cases / controls
SELFREP_CASES, SELFREP_CONTROLS = 38_472, 424_461


def liability_h2(h2_obs, K):
    """Lee et al. 2011 Eq.23（未 ascertainment 的总体样本，P≈K）。"""
    if K is None or not (0 < K < 1):
        return np.nan
    T = norm.isf(K)
    z = norm.pdf(T)
    return h2_obs * K * (1 - K) / (z ** 2)


def main():
    os.makedirs(OUT, exist_ok=True)
    ld = m09.load_ld_ref()
    M = m09.load_M()
    print("LD 参考 SNP=%d   M=%d" % (len(ld), M))

    # ---------------- fracture ----------------
    print("\n▶ fracture（对齐到 LD 参考）")
    bridge = m10.load_ld_bridge()
    pos2snp = {(int(r.CHR), int(r.BP)): r.SNP for r in bridge.itertuples()}
    fx = m10.load_fracture()
    fx["poskey"] = list(zip(fx.chromosome.astype(int),
                            fx.base_pair_location.astype(int)))
    fx = fx[fx.poskey.isin(pos2snp)].copy()
    fx["SNP"] = fx.poskey.map(pos2snp)
    snps = [s for s in fx.SNP if s in ld]
    fsub = fx.set_index("SNP").loc[snps]
    z_f = fsub.Z.values.astype(float)
    ell_f = np.array([ld[s] for s in fsub.index], dtype=float)
    print("  对齐 SNP = %d" % len(fsub))

    rows = []
    for label, N in FRACTURE_N.items():
        h2, se, inter, _ise = m08.compute_h2(z_f, ell_f, N, M, n_iter=3, n_blocks=100)
        print("  %-26s N=%7d  h2=%.4f  SE=%.4f  intercept=%.3f"
              % (label, N, h2, se, inter))
        rows.append(dict(trait="fracture", source=label, N=N,
                         h2_obs=h2, h2_se=se, intercept=inter))

    # ---------------- OA_selfrep ----------------
    print("\n▶ OA_selfrep")
    p = os.path.join(m09.FLAGSHIP_LDSC, "sumstats_oa_selfrep.sumstats.gz")
    d = pd.read_csv(p, sep="\t")
    d = d[d["SNP"].isin(ld)].copy()
    z_s = d["Z"].values.astype(float)
    ell_s = np.array([ld[s] for s in d["SNP"]], dtype=float)
    print("  对齐 SNP = %d" % len(d))

    K = SELFREP_CASES / (SELFREP_CASES + SELFREP_CONTROLS)
    for label, N in SELFREP_N.items():
        h2, se, inter, _ise = m08.compute_h2(z_s, ell_s, N, M, n_iter=3, n_blocks=100)
        h2l = liability_h2(h2, K if "462933" in label else None)
        print("  %-26s N=%7d  h2=%.4f  SE=%.4f  intercept=%.3f  h2_liab(K=%.4f)=%s"
              % (label, N, h2, se, inter, K,
                 ("%.4f" % h2l) if h2l == h2l else "NA"))
        rows.append(dict(trait="OA_selfrep", source=label, N=N,
                         h2_obs=h2, h2_se=se, intercept=inter,
                         K=K if "462933" in label else np.nan,
                         h2_liability=h2l if "462933" in label else np.nan))

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "N_sensitivity_decision.csv"), index=False)
    print("\n✅ -> output/genetics/N_sensitivity_decision.csv")
    print("\n提示：r_g 与 MR 的方向/显著性不随 N 改变（N 不变性），"
          "本表只影响 h² 绝对值与 liability 转换。")


if __name__ == "__main__":
    main()
