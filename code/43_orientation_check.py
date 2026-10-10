#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
43_orientation_check.py
=========================================================================
决定 OA_selfrep 在当前（2026-09-26 re-munge）文件下的 *真实* 方向，并核对 live
管线（09_rg_matrix.py，SIGN_FLIP={"OA_selfrep":-1}）是否仍输出正确的正 rg。

两个 OA 表型（knee OA ↔ self-reported OA）遗传相关 **必为正**（文献预期
rg≈+0.7~+0.9）。因此：
  · 若 SIGN_FLIP=-1 下 LDSC rg ≈ -0.838  → 翻转在当前文件上已**错误**，必须删除；
  · 若 SIGN_FLIP={}  下 LDSC rg ≈ +0.838  → 当前文件已 concordant，翻转冗余。

同时报告 Pearson r(Z)（瞬时、足够判定符号），以及完整 LDSC rg（用于回填论文）。
=========================================================================
"""
import os, sys, importlib
import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)
m09 = importlib.import_module("09_rg_matrix")

print("载入 LD 参考 / M ...")
ld = m09.load_ld_ref()
M = m09.load_M()
dfs = {nm: m09.load_sumstats(m09.PHENOS[nm]) for nm in ["OA_knee", "OA_selfrep"]}
print("  OA_knee SNPs = %d, OA_selfrep SNPs = %d" % (len(dfs["OA_knee"]), len(dfs["OA_selfrep"])))


def rg_pair(flip, n_blocks=20):
    m09.SIGN_FLIP = flip
    Zdf, ell, Ns, diag = m09.harmonise(dfs, ld)
    zk = Zdf["OA_knee"].values
    zs = Zdf["OA_selfrep"].values
    pear = float(np.corrcoef(zk, zs)[0, 1])
    hk = m09.compute_h2(zk, ell, Ns["OA_knee"], M, 3, n_blocks)
    hs = m09.compute_h2(zs, ell, Ns["OA_selfrep"], M, 3, n_blocks)
    gc = m09.compute_gencov(zk, zs, ell, Ns["OA_knee"], Ns["OA_selfrep"], M,
                            hk[0], hs[0], hk[2], hs[2], 3, n_blocks)
    rg = gc[0] / np.sqrt(max(hk[0], 1e-12) * max(hs[0], 1e-12))
    return rg, pear, diag["n_final"]


print("\n▶ 计算（SIGN_FLIP=-1，当前 live 代码）...")
rg_cur, pear_cur, n_cur = rg_pair({"OA_selfrep": -1})
print("  LDSC rg(OA_knee, OA_selfrep) = %+.4f   Pearson r(Z) = %+.4f   n=%d"
      % (rg_cur, pear_cur, n_cur))

print("\n▶ 计算（SIGN_FLIP={}，删除翻转）...")
rg_nof, pear_nof, n_nof = rg_pair({})
print("  LDSC rg(OA_knee, OA_selfrep) = %+.4f   Pearson r(Z) = %+.4f   n=%d"
      % (rg_nof, pear_nof, n_nof))

print("\n==================== 结论 ====================")
truth = +0.838
if abs(rg_nof - truth) < abs(rg_cur - truth):
    print("  删除 SIGN_FLIP 后 rg 接近生物学真值 +0.838 → 当前文件已 concordant，翻转已冗余/错误。")
    print("  → 行动：将 09_rg_matrix.py 的 SIGN_FLIP 置空，并重生成 rg_matrix.csv。")
else:
    print("  SIGN_FLIP=-1 仍给出正确 rg → 当前文件仍为 discordant，"


          "翻转仍需保留（与 42 的 94%% 一致性结论矛盾，需排查）。")
