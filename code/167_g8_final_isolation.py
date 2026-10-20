# -*- coding: utf-8 -*-
"""
G8 追加审计 4（终极）：code/10 与 code/09  fracture 行差异的完全隔离
==================================================================
已确立事实：
  · code/10（原始文件桥接）可复现地给出 rg(fracture,OA_knee)=+0.148117, SE=0.149351
    （今日重跑与 9/16 备份逐位一致）；
  · code/09（munged 文件）给出 +0.118958, SE=0.020073；
  · test 165 显示：若两条路径都用 code/09 的 ld 字典 harmonise，所得 Zdf 完全一致
    （max|dZ|=8.9e-16, SNP 集对称差 0, ell 零差异）；
  · test 166 显示：compute_gencov 对参数顺序完全对称（两种顺序逐位相同）；
  · 13 个 h2 在两条路径间完全一致（1e-16）。

唯一未验证差异：code/10 用**桥接 ld 字典**（load_ld_bridge -> set_index("SNP")["L2"]），
code/09/test165/166 用 **load_ld_ref() 字典**。本脚本：
  1. 直接比对两个 ld 字典（SNP 集 + L2 值）；
  2. 用各自的原生 ld 字典构建两条路径的 Zdf/ell/Ns（完全复刻生产调用）；
  3. 逐层比对并分别计算 gencov（两种参数顺序），定位分歧精确位置。
"""
import os
import sys
import importlib
import numpy as np
import pandas as pd
from scipy.stats import norm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
m09 = importlib.import_module("09_rg_matrix")
m10 = importlib.import_module("10_fracture_row")

NB = 100

# ---------------- 1. 两个 ld 字典 ----------------
print("▶ [1] 比对两个 LD 字典")
ld09 = m09.load_ld_ref()
ld_df = m10.load_ld_bridge()
ld10 = ld_df.set_index("SNP")["L2"].to_dict()
print(f"  load_ld_ref : {len(ld09):,} SNP")
print(f"  桥接字典     : {len(ld10):,} SNP")
s9, s10 = set(ld09), set(ld10)
print(f"  集合差: |only09|={len(s9-s10):,}  |only10|={len(s10-s9):,}")
common_ld = s9 & s10
if common_ld:
    diffs = [(s, ld09[s], ld10[s]) for s in common_ld if ld09[s] != ld10[s]]
    print(f"  共同 SNP 上 L2 不一致: {len(diffs):,}")
    for s, a, b in diffs[:5]:
        print(f"    {s}: load_ld_ref={a}  bridge={b}")

# ---------------- 2. 两条路径的原生 Zdf ----------------
print("\n▶ [2] 构建两条路径的原生 Zdf（各自原生 ld 字典）")
dfs = {}
for nm, fn in m09.PHENOS.items():
    if nm == "fracture":
        continue
    dfs[nm] = m09.load_sumstats(fn)

fxA = m09.load_sumstats("sumstats_fracture.sumstats.gz")          # munged
dA = dict(dfs); dA["fracture"] = fxA
ZA, ellA, NsA, diagA = m09.harmonise(dA, ld09)                    # code/09 原生

fx = m10.load_fracture()
merged = fx.merge(ld_df[["CHR", "BP", "SNP", "L2"]],
                  left_on=["chromosome", "base_pair_location"],
                  right_on=["CHR", "BP"], how="inner")
merged = merged.drop_duplicates(subset=["SNP"], keep="first")
fxB = merged.set_index("SNP")[["A1", "A2", "Z", "effect_allele_frequency"]].rename(
    columns={"effect_allele_frequency": "MAF"})
fxB["N"] = float(m10.DEFAULT_N)
dB = dict(dfs); dB["fracture"] = fxB
ZB, ellB, NsB, diagB = m09.harmonise(dB, ld10)                    # code/10 原生（桥接 ld）
NsB["fracture"] = m10.DEFAULT_N

print(f"  A(munged)   n={diagA['n_final']:,}  B(bridge) n={diagB['n_final']:,}")
print(f"  索引集对称差: |A-B|={len(set(ZA.index)-set(ZB.index)):,} |B-A|={len(set(ZB.index)-set(ZA.index)):,}")
print(f"  索引顺序一致: {list(ZA.index) == list(ZB.index)}")

print("\n▶ [3] 逐层比对")
if list(ZA.index) == list(ZB.index):
    for c in ZA.columns:
        d = np.abs(ZA[c].values - ZB[c].values)
        d = d[np.isfinite(d)]
        print(f"  Z[{c:16s}] max|d|={d.max():.3e}" if len(d) else f"  Z[{c}] 无有限值")
    print(f"  ell max|d| = {np.abs(ellA - ellB).max():.3e}")
    for nm in NsA:
        if abs(NsA[nm] - NsB.get(nm, np.nan)) > 1e-6:
            print(f"  Ns[{nm}] A={NsA[nm]:,.4f} B={NsB[nm]:,.4f}  <-- 不一致")
    print("  Ns 其余一致" )

# ---------------- 4. h2 + gencov 分别计算 ----------------
print("\n▶ [4] 分别计算 h2 与 gencov（各自原生 Zdf）")
M = m09.load_M()
hA = {}
for nm in ("OA_knee", "fracture"):
    r = m09.compute_h2(ZA[nm].values, ellA, NsA[nm], M, 3, NB)
    hA[nm] = dict(h2=r[0], intercept=r[2])
hB = {}
for nm in ("OA_knee", "fracture"):
    r = m09.compute_h2(ZB[nm].values, ellB, NsB[nm], M, 3, NB)
    hB[nm] = dict(h2=r[0], intercept=r[2])
print(f"  A: h2_knee={hA['OA_knee']['h2']:.12f} h2_fx={hA['fracture']['h2']:.12f} "
      f"int_knee={hA['OA_knee']['intercept']:.9f} int_fx={hA['fracture']['intercept']:.9f}")
print(f"  B: h2_knee={hB['OA_knee']['h2']:.12f} h2_fx={hB['fracture']['h2']:.12f} "
      f"int_knee={hB['OA_knee']['intercept']:.9f} int_fx={hB['fracture']['intercept']:.9f}")

def rg_of(Z, ell, Ns, h, order):
    Zk, Zf = Z["OA_knee"].values, Z["fracture"].values
    Nk, Nf = Ns["OA_knee"], Ns["fracture"]
    hk, hf = h["OA_knee"]["h2"], h["fracture"]["h2"]
    ik, if_ = h["OA_knee"]["intercept"], h["fracture"]["intercept"]
    den = np.sqrt(max(hk, 1e-12) * max(hf, 1e-12))
    if order == "09":
        g = m09.compute_gencov(Zk, Zf, ell, Nk, Nf, M, h1=hk, h2=hf, i_h1=ik, i_h2=if_, n_iter=3, n_blocks=NB)
    else:
        g = m09.compute_gencov(Zf, Zk, ell, Nf, Nk, M, h1=hf, h2=hk, i_h1=if_, i_h2=ik, n_iter=3, n_blocks=NB)
    return g[0] / den, abs(g[1] / den), g[0], g[1]

print("\n▶ [5] gencov / rg 四象限")
for tag, Z, ell, Ns, h in (("A(munged)", ZA, ellA, NsA, hA), ("B(bridge)", ZB, ellB, NsB, hB)):
    for order in ("09", "10"):
        rg, se, gc, gse = rg_of(Z, ell, Ns, h, order)
        print(f"  {tag:10s} order={order}: gc={gc:.10f} gc_se={gse:.10f} -> rg={rg:+.6f} SE={se:.6f}")
print("\n  参考: code/09 生产 = +0.118958 (SE 0.020073) ; code/10 生产 = +0.148117 (SE 0.149351)")
