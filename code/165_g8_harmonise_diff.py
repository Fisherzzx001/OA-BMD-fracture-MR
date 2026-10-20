# -*- coding: utf-8 -*-
"""
G8 追加审计 2：两条 fracture 管线经 code/09.harmonise 后的 Zdf 逐格比对
====================================================================
审计 1（164_）已证：munged 与原始骨折文件彼此内部自洽（等位基因互换 <=> Z 反号，
45.58% SNP；无同 allele 反号）。若 harmonise 对齐正确，两条管线应产出同一 Zdf。
但 code/09 重跑（munged）与 code/10（原始桥接）的 fracture 行 rg 系统性不符
（max |delta| 0.155），而全部 13 个 h2 一致（1e-16）。

本脚本对同一套 12 表型，分别以 (A) munged 文件、(B) 原始文件 chr:pos 桥接
两种方式接入 fracture，跑同一个 harmonise，逐格比对：
  - SNP 集合（规模 + 对称差）
  - 对齐后 Z 矩阵（max |delta|）
  - ell / Ns / diag
从而定位分歧究竟在 SNP 集、Z 值、还是别处。
"""
import os
import sys
import importlib
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
m09 = importlib.import_module("09_rg_matrix")
m10 = importlib.import_module("10_fracture_row")

print("▶ 读 LD 参考（code/09 口径）")
ld = m09.load_ld_ref()
M = m09.load_M()
print(f"  LD SNP = {len(ld):,}  M = {M:,.0f}")

print("▶ 读 12 表型 munged")
dfs = {}
for nm, fn in m09.PHENOS.items():
    if nm == "fracture":
        continue
    dfs[nm] = m09.load_sumstats(fn)
print(f"  {len(dfs)} 个表型")

# ---------- (A) munged fracture ----------
fxA = m09.load_sumstats("sumstats_fracture.sumstats.gz")
print(f"  (A) munged fracture: {len(fxA):,} SNP, 重复 rsid {fxA.index.duplicated().sum()}")

# ---------- (B) 原始文件桥接（完全复刻 code/10） ----------
print("▶ (B) 原始文件桥接（复刻 code/10）")
ld_df = m10.load_ld_bridge()
fx = m10.load_fracture()
print(f"  LD 桥接参考 {len(ld_df):,} SNP; 原始骨折 {len(fx):,} 变异")
merged = fx.merge(ld_df[["CHR", "BP", "SNP", "L2"]],
                  left_on=["chromosome", "base_pair_location"],
                  right_on=["CHR", "BP"], how="inner")
print(f"  桥接成功 {len(merged):,}")
merged = merged.drop_duplicates(subset=["SNP"], keep="first")
fxB = merged.set_index("SNP")[["A1", "A2", "Z", "effect_allele_frequency"]].rename(
    columns={"effect_allele_frequency": "MAF"})
fxB["N"] = float(m10.DEFAULT_N)
print(f"  (B) 桥接 fracture: {len(fxB):,} SNP, 重复 rsid {fxB.index.duplicated().sum()}")

# ---------- 分别 harmonise ----------
print("\n▶ harmonise (A)")
dA = dict(dfs); dA["fracture"] = fxA
ZA, ellA, NsA, diagA = m09.harmonise(dA, ld)
print(f"  n_final={diagA['n_final']:,}  flipped={diagA['n_flipped'].get('fracture')}")

print("▶ harmonise (B)")
dB = dict(dfs); dB["fracture"] = fxB
ZB, ellB, NsB, diagB = m09.harmonise(dB, ld)
print(f"  n_final={diagB['n_final']:,}  flipped={diagB['n_flipped'].get('fracture')}")

# ---------- 比对 ----------
print("\n=== 比对 ===")
print(f"diag A: {diagA}")
print(f"diag B: {diagB}")
sa, sb = set(ZA.index), set(ZB.index)
print(f"SNP 集: |A|={len(sa):,} |B|={len(sb):,} |A&B|={len(sa&sb):,} "
      f"|A-B|={len(sa-sb):,} |B-A|={len(sb-sa):,}")
common = sorted(sa & sb)
za = ZA.loc[common, "fracture"].values
zb = ZB.loc[common, "fracture"].values
ok = np.isfinite(za) & np.isfinite(zb)
dz = za[ok] - zb[ok]
print(f"共同 SNP 上 fracture 对齐 Z: 双方均有限 {ok.sum():,}; "
      f"max|dZ|={np.abs(dz).max():.3e}; 中位|dZ|={np.median(np.abs(dz)):.3e}")
sign_dis = int((np.sign(za[ok]) != np.sign(zb[ok])).sum())
print(f"符号不一致 SNP: {sign_dis:,} / {ok.sum():,} = {100*sign_dis/max(ok.sum(),1):.3f}%")
# 12 表型对齐 Z 是否一致（应一致）
mx12 = 0.0
for nm in dfs:
    xa = ZA.loc[common, nm].values
    xb = ZB.loc[common, nm].values
    m = np.isfinite(xa) & np.isfinite(xb)
    if m.sum():
        mx12 = max(mx12, float(np.abs(xa[m] - xb[m]).max()))
print(f"12 表型对齐 Z max|dZ| = {mx12:.3e}")
print(f"ell max|d| = {float(np.abs(ellA - ellB).max()):.3e}" if len(ellA) == len(ellB)
      else f"ell 长度不同: {len(ellA):,} vs {len(ellB):,}")
print(f"Ns A: " + ", ".join(f"{k}={v:,.0f}" for k, v in list(NsA.items())[:4]))
print(f"Ns B: " + ", ".join(f"{k}={v:,.0f}" for k, v in list(NsB.items())[:4]))

# 若符号不一致，给出样例与量级
if sign_dis:
    idx = np.where((np.sign(za[ok]) != np.sign(zb[ok])))[0][:10]
    snp_ok = [s for s, o in zip(common, np.isfinite(za) & np.isfinite(zb)) if o]
    print("\n符号不一致样例（共同 SNP）:")
    for i in idx:
        s = snp_ok[i]
        print(f"  {s:16s} A={za[ok][i]:+.4f}  B={zb[ok][i]:+.4f}  "
              f"|A|={abs(za[ok][i]):.3f} |B|={abs(zb[ok][i]):.3f}")
    dis_abs = np.abs(za[ok][np.sign(za[ok]) != np.sign(zb[ok])])
    print(f"  不一致 SNP |Z| 中位数 = {np.median(dis_abs):.3f}, >2 的 "
          f"{int((dis_abs>2).sum()):,} 个")
