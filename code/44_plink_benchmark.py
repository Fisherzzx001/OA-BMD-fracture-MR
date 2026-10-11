#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
44_plink_benchmark.py
=========================================================================
投稿前整改 #3 的「染色体级真实 PLINK 基准」。

自研纯 Python clumper (18_plink_clump_pure.py) 直接解析 1000G EUR .bed 复现
PLINK --clump。本脚本用**真实 PLINK 二进制**对同样的暴露 GWS SNP 集重做 clumping，
逐染色体比较保留的 lead SNP 集合，量化一致性。

步骤
----
1. 对每个 MR 暴露，取与 18 完全相同的 GWS 候选集（m13.load_all()，P<5e-8，
   列 CHR/BP/SNP/P），写出 PLINK --clump 输入；
2. 调用真实 PLINK：
     plink --bfile <EUR> --clump <in> --clump-p1 5e-8 --clump-p2 5e-8 \
           --clump-r2 0.001 --clump-kb 10000 --out <exp>_plink
3. 解析 PLINK .clumped 的 lead SNP 集，与 18 的 <exp>_clumped.clumped 保留集比较：
   · 总体 Jaccard / 命中率 / 召回率 / F1；
   · 逐染色体（chr1–22）的命中/召回与 Jaccard。
=========================================================================
"""
import os, sys, subprocess, importlib
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)
m13 = importlib.import_module("13_mr_analysis")
m09 = importlib.import_module("09_rg_matrix")

OUT = os.path.join(ROOT, "output", "genetics")
BFILE = r"H:/数据下载/1kg.v3/EUR"          # 不含扩展名
PLINK = r"H:/Direction1_LINCS/tools/plink.exe"
P_THRESH = 5e-8
R2 = 0.001
KB = 10000                                 # PLINK --clump-kb 单位 = kb → 10 Mb

# 与 manuscript / 18 一致的 MR 暴露
EXPOSURES = sorted({e for e, _, _ in m13.PAIRS})
print("MR 暴露（基准对象）：", ", ".join(EXPOSURES))


def write_plink_input(exp, dfs):
    d = dfs[exp]
    g = d[d["P"] < P_THRESH][["SNP", "CHR", "BP", "P"]].copy()
    g["SNP"] = g["SNP"].astype(str)
    g = g.sort_values("P")
    fp = os.path.join(OUT, f"{exp}_plink_clump_in.txt")
    g.to_csv(fp, sep="\t", index=False)
    return fp, len(g)


print("载入全部 GWAS（一次）...")
ALL = m13.load_all()


def run_plink(exp, inp):
    out = os.path.join(OUT, f"{exp}_plink")
    cmd = [PLINK, "--bfile", BFILE, "--clump", inp,
           "--clump-p1", str(P_THRESH), "--clump-p2", str(P_THRESH),
           "--clump-r2", str(R2), "--clump-kb", str(KB),
           "--out", out, "--silent"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    clumped = out + ".clumped"
    if not os.path.exists(clumped):
        print(f"  ⚠ {exp}: PLINK 未产生 .clumped；stderr 前 500 字：\n{r.stderr[:500]}")
        return set(), r.returncode
    c = pd.read_csv(clumped, sep=r"\s+")
    lead_col = "INDEX" if "INDEX" in c.columns else "SNP"
    leads = set(c[lead_col].astype(str))
    return leads, r.returncode


def read_self_clump(exp):
    fp = os.path.join(OUT, f"{exp}_clumped.clumped")
    c = pd.read_csv(fp, sep=r"\t")
    return set(c["SNP"].astype(str)), c


print("\n▶ 生成 PLINK clump 输入并运行真实 PLINK ...")
plink_leads = {}
self_sets = {}
for e in EXPOSURES:
    inp, n = write_plink_input(e, ALL)
    leads, rc = run_plink(e, inp)
    plink_leads[e] = leads
    self_sets[e], _ = read_self_clump(e)
    print(f"  {e:16s} GWS={n:4d}  PLINK保留={len(leads):4d}  自研保留={len(self_sets[e]):4d}")

print("\n==================== 染色体级一致性 ====================")
summary_rows = []
per_chr_rows = []
for e in EXPOSURES:
    A = plink_leads[e]          # PLINK lead 集
    B = self_sets[e]            # 自研 clump 保留集
    inter = A & B
    union = A | B
    precision = len(inter) / len(A) if A else float("nan")   # PLINK 视角召回自研
    recall = len(inter) / len(B) if B else float("nan")      # 自研被 PLINK 覆盖
    jac = len(inter) / len(union) if union else float("nan")
    f1 = (2 * precision * recall / (precision + recall)) if (precision + recall) else float("nan")
    summary_rows.append(dict(exposure=e, plink_kept=len(A), self_kept=len(B),
                             shared=len(inter), precision=round(precision, 4),
                             recall=round(recall, 4), jaccard=round(jac, 4),
                             f1=round(f1, 4)))
    # 逐染色体
    cA = pd.Series(list(A)).str  # 用 SNP→chr 映射
    # 建立 SNP->chr
    snp_chr = {}
    for src in (A, B):
        pass
    # 直接用 18 clumped 的 CHR 列 + PLINK clumped 的 CHR 列
    self_df = pd.read_csv(os.path.join(OUT, f"{e}_clumped.clumped"), sep=r"\t")
    plink_df = pd.read_csv(os.path.join(OUT, f"{e}_plink.clumped"), sep=r"\s+")
    self_chr = dict(zip(self_df["SNP"].astype(str), self_df["CHR"]))
    plink_chr = dict(zip(plink_df[("INDEX" if "INDEX" in plink_df.columns else "SNP")].astype(str),
                         plink_df["CHR"]))
    for chr_ in range(1, 23):
        a = {s for s in A if plink_chr.get(s) == chr_}
        b = {s for s in B if self_chr.get(s) == chr_}
        it = a & b
        un = a | b
        prec = len(it)/len(a) if a else float("nan")
        rec = len(it)/len(b) if b else float("nan")
        jc = len(it)/len(un) if un else float("nan")
        per_chr_rows.append(dict(exposure=e, CHR=chr_, plink_kept=len(a),
                                 self_kept=len(b), shared=len(it),
                                 precision=round(prec, 4) if a else "",
                                 recall=round(rec, 4) if b else "",
                                 jaccard=round(jc, 4) if un else ""))

sumdf = pd.DataFrame(summary_rows)
print(sumdf.to_string(index=False))
print("\n总体（跨暴露汇总）：")
allA = set().union(*plink_leads.values())
allB = set().union(*self_sets.values())
print("  PLINK 总保留=%d  自研总保留=%d  共享=%d" %
      (len(allA), len(allB), len(allA & allB)))
print("  总体 Jaccard = %.4f" % (len(allA & allB)/len(allA | allB)))

perdf = pd.DataFrame(per_chr_rows)
# 每个暴露的逐染色体摘要（前几行示例）
print("\n--- 逐染色体（示例：BMD_total 与 OA_knee）---")
print(perdf[perdf.exposure.isin(["BMD_total", "OA_knee"])].to_string(index=False))

sumdf.to_csv(os.path.join(OUT, "plink_benchmark_summary.csv"), index=False)
perdf.to_csv(os.path.join(OUT, "plink_benchmark_per_chr.csv"), index=False)
print("\n✅ 基准结果 -> plink_benchmark_summary.csv / plink_benchmark_per_chr.csv")
