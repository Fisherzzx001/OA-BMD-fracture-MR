#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
15_liability_h2.py — 把 OA 的 SNP-h² 从观测尺度转换到易感性(liability)尺度  [Item 5]

公式（Lee et al. 2011, Am J Hum Genet 88:294-305, Eq 23；LDSC 采用）：
    h²_liab = h²_obs × [K(1-K)] / z² × [K(1-K)] / [P(1-P)]
其中 K = 人群患病率，P = 样本中病例比例，z = φ(T)，T = norm.isf(K)。
UKB 的 OA 为未分层(unascertained)人群样本，P≈K，简化为
    h²_liab = h²_obs × K(1-K) / z²。

注意：BMD 为连续性状，无需转换；仅 OA_knee / OA_selfrep（二分类）需要。
K（人群患病率）需确认：膝 OA 在 UKB 约 0.10–0.16，自评 OA(任一侧)约 0.20–0.30。
本脚本对一组 K 做敏感性展示，并采用主要 K 给出主值。
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "output", "genetics")
H2CSV = os.path.join(OUT, "rg_h2_trueN.csv")

# 主要 K（人群患病率）与敏感性网格
K_PRIMARY = {"OA_knee": 0.12, "OA_selfrep": 0.25}
K_GRID = [0.08, 0.10, 0.12, 0.15, 0.20, 0.25, 0.30]


def h2_liab(h2_obs, K, P=None):
    P = K if P is None else P
    T = norm.isf(K)
    z = norm.pdf(T)
    return h2_obs * (K * (1 - K)) / (z ** 2) * (K * (1 - K)) / (P * (1 - P))


def main():
    df = pd.read_csv(H2CSV)
    od = df.set_index("phenotype")["h2"].to_dict()
    rows = []
    for trait, K0 in K_PRIMARY.items():
        h2_obs = od.get(trait)
        if h2_obs is None:
            print(f"  ⚠ 缺少 {trait} 的 h²_obs（先跑 16_recompute）")
            continue
        main_val = h2_liab(h2_obs, K0)
        sens = {f"K={k}": round(h2_liab(h2_obs, k), 4) for k in K_GRID}
        print(f"\n{trait}: h²_obs={h2_obs:.4f}  ->  h²_liab(K={K0})={main_val:.4f}")
        print(f"   敏感性 {sens}")
        rows.append({"trait": trait, "h2_obs": h2_obs, "K_primary": K0,
                     "h2_liab_primary": main_val,
                     "h2_liab_min(K=0.08)": h2_liab(h2_obs, 0.08),
                     "h2_liab_max(K=0.30)": h2_liab(h2_obs, 0.30)})
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(OUT, "oa_h2_liability.csv"), index=False)
    print(f"\n✅ -> {os.path.join(OUT, 'oa_h2_liability.csv')}")


if __name__ == "__main__":
    main()
