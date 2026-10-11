#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""16b_fracture_h2.py — 骨折 h² 重算（占位 N 敏感性）。h² 用 z²，无需等位基因翻转。"""
import os, sys, importlib, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
m10 = importlib.import_module("10_fracture_row")
m08 = importlib.import_module("08_genetic_correlation")
m09 = importlib.import_module("09_rg_matrix")

ld = m09.load_ld_ref(); M = m09.load_M()
bridge = m10.load_ld_bridge()
pos2snp = {(int(r.CHR), int(r.BP)): r.SNP for r in bridge.itertuples()}
fx = m10.load_fracture()
fx["poskey"] = list(zip(fx.chromosome.astype(int), fx.base_pair_location.astype(int)))
fx = fx[fx.poskey.isin(pos2snp)].copy()
fx["SNP"] = fx.poskey.map(pos2snp)
snps = [s for s in fx.SNP if s in ld]
fsub = fx.set_index("SNP").loc[snps]
z = fsub.Z.values.astype(float)
ellf = np.array([ld[s] for s in fsub.index], dtype=float)
print(f"骨折对齐 SNP={len(fsub):,}")
rows = []
for N in (400_000, 90_000, 32_000):
    r = m08.compute_h2(z, ellf, N, M, n_iter=3, n_blocks=100)
    print(f"  fracture N={N:>8,}  h2={r[0]:.4f} SE={r[1]:.4f} int={r[2]:.3f}")
    rows.append((N, r[0], r[1], r[2]))
pd.DataFrame(rows, columns=["N_assumed", "h2", "h2_se", "intercept"]).to_csv(
    os.path.join(m09.ROOT, "output", "genetics", "fracture_h2_sensitivity.csv"), index=False)
print("✅ -> output/genetics/fracture_h2_sensitivity.csv  (真实 N 待 Neale-lab manifest 确认)")
