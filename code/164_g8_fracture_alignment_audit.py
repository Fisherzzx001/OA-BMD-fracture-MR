# -*- coding: utf-8 -*-
"""
G8 追加审计：munged 骨折 sumstats 与原始 fracture_anytype.tsv 逐 SNP 比对
====================================================================
背景：code/09（G8 修复后，munged 口径）重跑的 13x13 rg 矩阵中，fracture 行与
code/10（原始文件桥接口径）系统性不符（max |delta_rg| = 0.155），而两套流程的
全部 13 个 h2 完全一致（1e-16）-> h2 只依赖 z^2 对符号不敏感，gencov 依赖
z1*z2 对符号敏感 => 怀疑一部分 SNP 的 Z 符号在两套流程间不一致。

本脚本：按 (chr,pos) 直接比对两套文件的 A1/A2/Z（munged 文件自带 CHR/BP，
无需 LD 桥接），定位分歧 SNP 的确切模式（等位基因同/换/异 + Z 同号/反号），
并区分「原始文件同位置多行」与「munged 引入」两种可能机制。

用法：python 164_g8_fracture_alignment_audit.py
"""
import os
import numpy as np
import pandas as pd

MUNG = r"C:\Users\Administrator\WorkBuddy\2026-09-08-20-17-24\OA_FLAGSHIP\07_reference\ldsc\sumstats_fracture.sumstats.gz"
RAW = (r"C:\Users\Administrator\WorkBuddy\2026-09-09-14-17-27\Direction1_pipeline"
       r"\Direction1_pipeline\data\raw\fracture_anytype.tsv")

# ---------------------------------------------------------------- 1. munged
m = pd.read_csv(MUNG, sep=r"\s+",
                usecols=["CHR", "BP", "SNP", "A1", "A2", "Z"],
                dtype={"CHR": "int32", "BP": "int64", "SNP": str,
                       "A1": str, "A2": str, "Z": "float64"})
m["A1"] = m.A1.str.upper()
m["A2"] = m.A2.str.upper()
m["key"] = m.CHR.astype(np.int64) * 1_000_000_000 + m.BP.astype(np.int64)
print(f"munged: {len(m):,} SNP, 多字符等位基因 {(m.A1.str.len()>1).sum() + (m.A2.str.len()>1).sum()} 处")
m_single = m[(m.A1.str.len() == 1) & (m.A2.str.len() == 1)]
print(f"  其中单核苷酸等位基因: {len(m_single):,}")

mkey = set(m_single.key.values)
mpos = {k: (s, a1, a2, z) for k, s, a1, a2, z in
        zip(m_single.key, m_single.SNP, m_single.A1, m_single.A2, m_single.Z)}

# ---------------------------------------------------------------- 2. raw（分块，仅保留 munged 位置）
first_raw = {}      # key -> (ea, oa, Z)  文件顺序首行（code/10 keep='first' 口径）
n_raw_at = {}       # key -> 原始文件该位置的行数（单核苷酸行）
n_raw_rows = 0
n_raw_snp_rows = 0
for chunk in pd.read_csv(RAW, sep="\t", chunksize=1_000_000,
                         usecols=["chromosome", "base_pair_location",
                                  "effect_allele", "other_allele",
                                  "beta", "standard_error"],
                         dtype={"chromosome": "int32",
                                "base_pair_location": "int64",
                                "effect_allele": "category",
                                "other_allele": "category",
                                "beta": "float64",
                                "standard_error": "float64"}):
    n_raw_rows += len(chunk)
    c = chunk.dropna(subset=["beta", "standard_error"])
    c = c[(c.standard_error > 0) & (c.beta.notna())]
    ea = c.effect_allele.astype(str).str.upper()
    oa = c.other_allele.astype(str).str.upper()
    snp_row = (ea.str.len() == 1) & (oa.str.len() == 1)
    c = c[snp_row.values]
    ea, oa = ea[snp_row.values], oa[snp_row.values]
    n_raw_snp_rows += len(c)
    key = c.chromosome.astype(np.int64).values * 1_000_000_000 + c.base_pair_location.astype(np.int64).values
    z = (c.beta.values / c.standard_error.values)
    for k, e, o, zz in zip(key, ea.values, oa.values, z):
        if k in mkey:
            n_raw_at[k] = n_raw_at.get(k, 0) + 1
            if k not in first_raw:
                first_raw[k] = (e, o, float(zz))

print(f"raw: {n_raw_rows:,} 行, 单核苷酸行 {n_raw_snp_rows:,}; "
      f"落在 munged 位置上的行 {sum(n_raw_at.values()):,}, 涉及位置 {len(first_raw):,}")

# ---------------------------------------------------------------- 3. 比对
rows = []
no_raw = 0
for k, (snp, a1m, a2m, zm) in mpos.items():
    if k not in first_raw:
        no_raw += 1
        continue
    e, o, zr = first_raw[k]
    if (a1m, a2m) == (e, o):
        allele = "same"
    elif (a1m, a2m) == (o, e):
        allele = "swapped"
    else:
        allele = "different"
    dz = zm - zr
    if abs(dz) < 1e-9 * max(1.0, abs(zr)):
        zrel = "Z_same"
    elif abs(zm + zr) < 1e-9 * max(1.0, abs(zr)):
        zrel = "Z_flip"
    else:
        zrel = "Z_diff"
    rows.append((snp, allele, zrel, n_raw_at[k], abs(zm), abs(zr)))

df = pd.DataFrame(rows, columns=["SNP", "allele", "zrel", "n_raw_rows", "absZm", "absZr"])
print(f"\nmunged 单核苷酸 SNP {len(mpos):,}；原始文件无对应位置: {no_raw:,}；可比对: {len(df):,}")
print("\n=== 等位基因关系 x Z 关系 交叉表（仅计原始文件首行）===")
print(pd.crosstab(df.allele, df.zrel, margins=True))

print("\n=== 按原始文件同位置行数分层 ===")
print(pd.crosstab(df.n_raw_rows, df.zrel, margins=True))

bad = df[((df.allele == "same") & (df.zrel == "Z_flip")) |
         ((df.allele == "swapped") & (df.zrel == "Z_same")) |
         (df.zrel == "Z_diff")]
print(f"\n内部不一致 SNP（同 allele 反号 / 换 allele 同号 / Z 数值不符）: {len(bad):,}")
if len(bad):
    print("  |Z| 分布（munged）:", np.percentile(bad.absZm, [0, 25, 50, 75, 100]).round(3))
    print("  其中 |Z_munged|>2 的:", int((bad.absZm > 2).sum()))
    print("  样例:", bad.head(10).to_string(index=False))

# 分歧 SNP 占全体的比例与 rg 影响的量级评估
flip = df[df.zrel == "Z_flip"]
print(f"\nZ 反号 SNP 合计: {len(flip):,} / {len(df):,} = {100*len(flip)/len(df):.2f}%")
if len(flip):
    print("  反号 SNP 的 |Z| 中位数:", round(float(flip.absZm.median()), 3))
