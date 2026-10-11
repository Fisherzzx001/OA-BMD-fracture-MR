#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
42_selfrep_validation.py
=========================================================================
独立验证 OA_selfrep 全局 Z 符号翻转（投稿前整改 #8 的「基因组级独立验证」）。

背景
----
doc 11 用 **两个** 哨兵 OA 位点（GDF5 rs224333、FTO rs9930333）诊断 OA_selfrep
整体 Z 反号，并施加全局翻转 SIGN_FLIP={"OA_selfrep":-1}。但「仅 2 个哨兵」是全文
最脆弱的单点：翻转本身是用 OA_knee×OA_selfrep 这对的不一致推导出来的，再用这对
验证就是循环论证。

本脚本提供 **基因组级独立验证**：
  · 比较对象 = OA_knee（GO Consortium 2025，Nature，Ref 1）——一个**独立生成**的
    膝 OA GWAS，其方向即「已发表的 OA 风险方向」；
  · 在 OA_knee 全部基因组显著位点（|Z|>5.45，P<5e-8）上，计算 OA_selfrep 的 Z 符号
    与 OA_knee 的一致性，原始 vs 翻转后；
  · 一致性从 ~0%（原始）跃升至 ~100%（翻转后）即证明翻转正确，且验证不再依赖那
    两个哨兵 SNP。

输出
----
code/output/genetics/selfrep_validation.csv  +  stdout 汇总。
=========================================================================
"""
import os, gzip, csv
import numpy as np
import pandas as pd

LDSC = r"C:\Users\Administrator\WorkBuddy\2026-09-08-20-17-24\OA_FLAGSHIP\07_reference\ldsc"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "genetics")
os.makedirs(OUT, exist_ok=True)

SENTINELS = {"rs224333": "GDF5", "rs9930333": "FTO"}   # doc 11 诊断用哨兵


def load_z(path):
    """LDSC .sumstats.gz -> DataFrame[SNP, A1, Z]（仅读必要列）。"""
    cols = ["SNP", "A1", "Z"]
    d = pd.read_csv(path, sep=r"\s+", usecols=cols,
                    dtype={"SNP": "str", "A1": "str"})
    d["Z"] = d["Z"].astype(float)
    return d


print("载入 OA_knee / OA_selfrep 汇总统计 ...")
knee = load_z(os.path.join(LDSC, "sumstats_oa_knee.sumstats.gz"))
selfr = load_z(os.path.join(LDSC, "sumstats_oa_selfrep.sumstats.gz"))
print("  OA_knee 行数 = %d, OA_selfrep 行数 = %d" % (len(knee), len(selfr)))

# 对齐到 SNP；两文件均以 OA_knee 为基准表型做 A1 harmonisation（见 09_rg_matrix REF_PHENO）
m = knee.merge(selfr, on="SNP", suffixes=("_knee", "_self"), how="inner")
print("  交集 SNP = %d" % len(m))
a1_mismatch = int((m["A1_knee"] != m["A1_self"]).sum())
print("  A1 不一致数 = %d（应≈0，因均已对齐到 OA_knee）" % a1_mismatch)
if a1_mismatch:
    # 若 A1 不同，Z 符号不可直接比；翻转 A1 不一致者的 self Z 再比（保守）
    flipA1 = m["A1_knee"] != m["A1_self"]
    m.loc[flipA1, "Z_self"] = -m.loc[flipA1, "Z_self"]

# 基因组显著 OA_knee 位点
gws = m[m["Z_knee"].abs() > 5.45].copy()
print("\nOA_knee 基因组显著位点（|Z|>5.45）= %d" % len(gws))

# 原始一致性
raw_conc = int((np.sign(gws["Z_knee"]) == np.sign(gws["Z_self"])).sum())
# 翻转后一致性（SIGN_FLIP = -1）
flipped = np.sign(gws["Z_knee"]) == np.sign(-gws["Z_self"])
flip_conc = int(flipped.sum())

# 染色体级分解：从原文件载入 CHR
knee_chr = pd.read_csv(os.path.join(LDSC, "sumstats_oa_knee.sumstats.gz"),
                       sep=r"\s+", usecols=["SNP", "CHR"], dtype={"SNP": "str", "CHR": "str"})
gws = gws.merge(knee_chr, on="SNP", how="left")

rows = []
per_chr = gws.groupby("CHR")
for chr_, sub in per_chr:
    n = len(sub)
    rc = int((np.sign(sub["Z_knee"]) == np.sign(sub["Z_self"])).sum())
    fc = int((np.sign(sub["Z_knee"]) == np.sign(-sub["Z_self"])).sum())
    rows.append(dict(CHR=str(chr_), n_gws=n,
                     concord_raw=rc, concord_raw_pct=round(100*rc/n, 1),
                     concord_flipped=fc, concord_flipped_pct=round(100*fc/n, 1)))

summary = dict(
    n_knee_gws=len(gws),
    concord_raw=raw_conc, concord_raw_pct=round(100*raw_conc/len(gws), 2),
    concord_flipped=flip_conc, concord_flipped_pct=round(100*flip_conc/len(gws), 2),
    a1_mismatch=a1_mismatch,
)
print("\n==================== 结果 ====================")
print("OA_knee 基因组显著位点总数 : %d" % len(gws))
print("原始一致性  (OA_selfrep 未翻转) : %d/%d = %.2f%%"
      % (raw_conc, len(gws), summary["concord_raw_pct"]))
print("翻转后一致性(SIGN_FLIP=-1)     : %d/%d = %.2f%%"
      % (flip_conc, len(gws), summary["concord_flipped_pct"]))

# 哨兵位点明细
print("\n--- 哨兵位点明细（doc 11 诊断 SNP）---")
for snp, gene in SENTINELS.items():
    r = gws[gws["SNP"].str.contains(snp)]
    if len(r):
        r = r.iloc[0]
        print("  %-8s %-6s knee Z=%+.3f  self Z=%+.3f  raw_concord=%s  flipped_concord=%s"
              % (snp, gene, r["Z_knee"], r["Z_self"],
                  np.sign(r["Z_knee"]) == np.sign(r["Z_self"]),
                  np.sign(r["Z_knee"]) == np.sign(-r["Z_self"])))

print("\n--- 染色体级一致性（翻转后应为 ~100%%）---")
chr_tbl = pd.DataFrame(rows).sort_values("CHR")
print(chr_tbl.to_string(index=False))

# 写盘
with open(os.path.join(OUT, "selfrep_validation.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["metric", "value"])
    for k, v in summary.items():
        w.writerow([k, v])
    w.writerow([])
    w.writerow(["CHR", "n_gws", "concord_raw", "concord_raw_pct",
                "concord_flipped", "concord_flipped_pct"])
    for r in rows:
        w.writerow([r["CHR"], r["n_gws"], r["concord_raw"], r["concord_raw_pct"],
                    r["concord_flipped"], r["concord_flipped_pct"]])
print("\nSaved -> %s" % os.path.join(OUT, "selfrep_validation.csv"))
