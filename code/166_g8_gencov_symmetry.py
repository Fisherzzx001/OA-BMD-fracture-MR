# -*- coding: utf-8 -*-
"""
G8 追加审计 3：compute_gencov 参数顺序对称性实证测试
==================================================
已证：两条 fracture 数据管线 harmonise 后 Zdf/ell/Ns 完全一致（max|dZ|=8.9e-16）；
code/08 自 2026-09-15 未改动；13 个 h2 在两轮运行间完全一致（1e-16）。
但 code/10（9/16，fracture 行 rg=0.148, SE=0.149）与 code/09 重跑（10/20，
rg=0.119, SE=0.020）输出不同。两轮的唯一调用差异是 compute_gencov 的参数顺序
（code/10: (fracture, trait)；code/09: (trait, fracture) 对 fracture 列序第 3 的表型）。

本脚本用当前数据对 (OA_knee, fracture) 按两种顺序各跑一次 compute_gencov
（n_iter=3, n_blocks=100，与两轮生产运行一致），检验：
  (1) 两种顺序是否给出相同结果（对称性）；
  (2) 各顺序复现 0.148 / 0.119 中的哪一个。
"""
import os
import sys
import importlib
import numpy as np
import pandas as pd
from scipy.stats import norm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
m09 = importlib.import_module("09_rg_matrix")

print("▶ 加载 LD 参考 + 13 表型（munged fracture）")
ld = m09.load_ld_ref()
M = m09.load_M()
print(f"  LD SNP={len(ld):,}  M={M:,.0f}")
dfs = {}
for nm, fn in m09.PHENOS.items():
    dfs[nm] = m09.load_sumstats(fn)
Zdf, ell, Ns, diag = m09.harmonise(dfs, ld)
print(f"  n_final={diag['n_final']:,}")

NB = 100
print(f"\n▶ compute_h2（n_blocks={NB}）")
h = {}
for nm in ("OA_knee", "fracture"):
    r = m09.compute_h2(Zdf[nm].values, ell, Ns[nm], M, 3, NB)
    h[nm] = dict(h2=r[0], se=r[1], intercept=r[2], intercept_se=r[3])
    print(f"  {nm:10s} h2={r[0]:.10f}  intercept={r[2]:.10f}")

Zk, Zf = Zdf["OA_knee"].values, Zdf["fracture"].values
Nk, Nf = Ns["OA_knee"], Ns["fracture"]
hk, hf = h["OA_knee"]["h2"], h["fracture"]["h2"]
ik, if_ = h["OA_knee"]["intercept"], h["fracture"]["intercept"]
den = np.sqrt(max(hk, 1e-12) * max(hf, 1e-12))

print("\n▶ 顺序 1：compute_gencov(OA_knee, fracture, Nk, Nf, h1=hk, h2=hf, i_h1=ik, i_h2=if)  [code/09 口径]")
g1 = m09.compute_gencov(Zk, Zf, ell, Nk, Nf, M,
                        h1=hk, h2=hf, i_h1=ik, i_h2=if_, n_iter=3, n_blocks=NB)
rg1 = g1[0] / den
se1 = abs(g1[1] / den)
print(f"  gencov={g1[0]:.10f}  gc_se={g1[1]:.10f}  -> rg={rg1:+.6f}  SE={se1:.6f}  "
      f"p={2*norm.sf(abs(rg1/se1)):.3e}")

print("\n▶ 顺序 2：compute_gencov(fracture, OA_knee, Nf, Nk, h1=hf, h2=hk, i_h1=if, i_h2=ik)  [code/10 口径]")
g2 = m09.compute_gencov(Zf, Zk, ell, Nf, Nk, M,
                        h1=hf, h2=hk, i_h1=if_, i_h2=ik, n_iter=3, n_blocks=NB)
rg2 = g2[0] / den
se2 = abs(g2[1] / den)
print(f"  gencov={g2[0]:.10f}  gc_se={g2[1]:.10f}  -> rg={rg2:+.6f}  SE={se2:.6f}  "
      f"p={2*norm.sf(abs(rg2/se2)):.3e}")

print("\n=== 裁决 ===")
print(f"顺序1 vs 顺序2:  rg 差={abs(rg1-rg2):.3e}  SE 差={abs(se1-se2):.3e}")
print(f"参考值:  code/09 重跑(10/20) rg=+0.118958 SE=0.020073 ; "
      f"code/10(9/16) rg=+0.148117 SE=0.149351")
