#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
45_regen_selfrep_row.py  (v2 — 全局共同 SNP 集)
=========================================================================
2026-09-26 整改：OA_selfrep 文件经重新 munge（N 更正为 462,933，方向与独立
OA_knee GWAS 一致），故 09_rg_matrix.py 的 SIGN_FLIP 已置空。

本脚本**仅**重算受影响的 OA_selfrep 行/列，但为保证与矩阵其余 55 个单元格
（均在全局 12 表型交集上估计）一致，这里对所有 12 表型做统一 harmonise，
取全局共同 SNP 集，再：
  · 重算 OA_selfrep 的 h² / 截距（文件已变）；
  · 其余 11 个表型的 h² / 截距直接读 rg_h2.csv（它们未被重新 munge，不变）；
  · 对 11 个 (OA_selfrep × X) 对重算 gencov → rg / SE / P（均在全局集上）；
  · 对称回填 rg_matrix / rg_se_matrix / rg_p_matrix / rg_gencov_matrix / rg_h2。
=========================================================================
"""
import os, sys, importlib
import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)
m09 = importlib.import_module("09_rg_matrix")
OUT = m09.OUT_DIR
NB = 40

print("载入 LD 参考 / M ...")
ld = m09.load_ld_ref()
M = m09.load_M()

print("▶ 载入全部 12 表型并统一 harmonise（取全局共同 SNP 集）...")
dfs = {nm: m09.load_sumstats(m09.PHENOS[nm]) for nm in m09.PHENOS}
Zdf, ell, Ns, diag = m09.harmonise(dfs, ld)
print("  全局共同 SNP 数 = %d" % len(Zdf))

h2csv = pd.read_csv(os.path.join(OUT, "rg_h2.csv"), index_col=0)

print("▶ 重算 OA_selfrep h²（全局集）...")
zs = Zdf["OA_selfrep"].values
h_s = m09.compute_h2(zs, ell, Ns["OA_selfrep"], M, 3, NB)
print("  OA_selfrep h2 = %.4f (SE %.4f) intercept = %.3f" % (h_s[0], h_s[1], h_s[2]))

others = [p for p in m09.PHENOS if p != "OA_selfrep"]
rg_new, se_new, p_new, gc_new = {}, {}, {}, {}
print("▶ 重算 11 个 (OA_selfrep × X) 遗传相关（全局集）...")
for X in others:
    zx = Zdf[X].values
    h_x = (float(h2csv.loc[X, "h2"]), float(h2csv.loc[X, "intercept"]))
    gc = m09.compute_gencov(zs, zx, ell, Ns["OA_selfrep"], Ns[X], M,
                            h_s[0], h_x[0], h_s[2], h_x[1], 3, NB)
    gencov, gc_se = gc[0], gc[1]
    denom = np.sqrt(max(h_s[0], 1e-12) * max(h_x[0], 1e-12))
    rg = gencov / denom
    se = abs(gc_se / denom) if np.isfinite(gc_se) else np.nan
    pp = 2 * norm.sf(abs(rg / se)) if (np.isfinite(se) and se > 0) else np.nan
    rg_new[X], se_new[X], p_new[X], gc_new[X] = rg, se, pp, gencov
    print("  OA_selfrep ~ %-16s rg=%+.4f (SE %.4f) P=%.2e" % (X, rg, se, pp))

print("\n▶ 回填矩阵文件 ...")
rg_mat = pd.read_csv(os.path.join(OUT, "rg_matrix.csv"), index_col=0)
se_mat = pd.read_csv(os.path.join(OUT, "rg_se_matrix.csv"), index_col=0)
p_mat = pd.read_csv(os.path.join(OUT, "rg_p_matrix.csv"), index_col=0)
gc_mat = pd.read_csv(os.path.join(OUT, "rg_gencov_matrix.csv"), index_col=0)
for X in others:
    rg_mat.loc["OA_selfrep", X] = rg_mat.loc[X, "OA_selfrep"] = rg_new[X]
    se_mat.loc["OA_selfrep", X] = se_mat.loc[X, "OA_selfrep"] = se_new[X]
    p_mat.loc["OA_selfrep", X] = p_mat.loc[X, "OA_selfrep"] = p_new[X]
    gc_mat.loc["OA_selfrep", X] = gc_mat.loc[X, "OA_selfrep"] = gc_new[X]
h2csv.loc["OA_selfrep"] = [h_s[0], h_s[1], h_s[2], h_s[3],
                           float(np.nanmean(dfs["OA_selfrep"].loc[Zdf.index, "N"].values))]
rg_mat.to_csv(os.path.join(OUT, "rg_matrix.csv"))
se_mat.to_csv(os.path.join(OUT, "rg_se_matrix.csv"))
p_mat.to_csv(os.path.join(OUT, "rg_p_matrix.csv"))
gc_mat.to_csv(os.path.join(OUT, "rg_gencov_matrix.csv"))
h2csv.to_csv(os.path.join(OUT, "rg_h2.csv"))

print("\n==================== 关键结果 ====================")
print("rg(OA_knee, OA_selfrep) = %+.4f  (全文此前报 +0.838，来自旧 discordant 文件)" % rg_new["OA_knee"])
print("✅ 已回填 rg_matrix.csv / rg_se_matrix.csv / rg_p_matrix.csv / rg_gencov_matrix.csv / rg_h2.csv")
print("   其余 55 个非 OA_selfrep 单元格保持不变（与当前文件一致）。")
