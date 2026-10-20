#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""G6 探查 B：用真实膝 OA GWAS（GO Consortium 2025，Hatzikotoulas Nature）逐 SNP 估计 K(1-K)。

逻辑：
  munged 文件给出每 SNP 的 Z / MAF / N（N = 总样本量，线性重构用）
  真实文件给出每 SNP 的 beta_true / se_true（log-OR 尺度）
  则  se_true = 1/sqrt(N_eff * 2 f (1-f))， N_eff = N * K(1-K)
  =>  K(1-K) = N_eff / N
  若真实文件与 munged 同源，各 SNP 的 K(1-K) 应集中在同一个常数附近。
"""
import os
import numpy as np
import pandas as pd

LDSC = r"C:\Users\Administrator\WorkBuddy\2026-09-08-20-17-24\OA_FLAGSHIP\07_reference\ldsc"
TRUE = r"C:\Users\Administrator\WorkBuddy\2026-09-08-20-17-24\OA_FLAGSHIP\02_data\11_bmd_snp_cache\gcst_knee_union.csv.gz"

mun = pd.read_csv(os.path.join(LDSC, "sumstats_oa_knee.sumstats.gz"), sep=r"\s+")
tru = pd.read_csv(TRUE)
print("munged:", mun.shape, "true:", tru.shape)

tru[["CHR", "BP"]] = tru.key.str.split(":", expand=True)
tru["CHR"] = tru.CHR.astype(int)
tru["BP"] = tru.BP.astype(int)
tru = tru.rename(columns={"EA_o": "A1t", "OA_o": "A2t", "eaf_o": "MAFt"})

m = mun.set_index(["CHR", "BP"])
rows = []
for _, r in tru.iterrows():
    k = (r.CHR, r.BP)
    if k not in m.index:
        continue
    q = m.loc[k]
    if isinstance(q, pd.DataFrame):
        q = q.iloc[0]
    # Z 一致性检查（同源性）
    z_true = r.beta_o / r.se_o
    rows.append(dict(SNP=q.SNP, CHR=r.CHR, BP=r.BP,
                     Z_mun=q.Z, Z_true=z_true,
                     N_mun=q.N, MAF_mun=q.MAF, MAF_true=r.MAFt,
                     beta_true=r.beta_o, se_true=r.se_o,
                     A1m=q.A1, A2m=q.A2, A1t=r.A1t, A2t=r.A2t))
d = pd.DataFrame(rows)
print("matched:", len(d))
if len(d) == 0:
    raise SystemExit("no match")

d["z_ratio"] = d.Z_true / d.Z_mun
d["N_eff_true"] = 1.0 / (d.se_true ** 2 * 2 * d.MAF_true * (1 - d.MAF_true))
d["Kc"] = d.N_eff_true / d.N_mun
print("\n--- Z 一致性（应≈1，符号可能因等位基因方向差 -1）---")
print(d.z_ratio.abs().describe())
print("\n--- N_eff_true ---")
print(d.N_eff_true.describe())
print("\n--- K(1-K) 逐 SNP 估计 ---")
print(d.Kc.describe())
print("median %.6f  IQR %.6f-%.6f" % (d.Kc.median(), d.Kc.quantile(.25), d.Kc.quantile(.75)))
kc = d.Kc.median()
print("=> 隐含 K = %.4f (或 %.4f)" % ((1 - np.sqrt(1 - 4 * kc)) / 2,
                                      (1 + np.sqrt(1 - 4 * kc)) / 2))
print("=> 校正因子 1/sqrt(K(1-K)) = %.4f" % (1 / np.sqrt(kc)))

# 用 munged 的 Z 与真实尺度重算 beta/se（log-OR），与真实值比对
f = 1 / np.sqrt(kc)
d["beta_recon_from_Z"] = d.Z_mun * d.se_true
print("\n--- 用 |Z_mun| * se_true 比对真实 |beta|（应完全相等）---")
print((d.beta_recon_from_Z.abs() - d.beta_true.abs()).abs().describe())

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "genetics",
                   "_g6_k_estimation.csv")
d.to_csv(out, index=False)
print("\nwritten:", out)
