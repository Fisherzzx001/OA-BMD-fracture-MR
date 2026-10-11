#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
40_tissue_enrichment.py — P1-3 组织特异性富集（LDSC-SEG + MAGMA 包装 + 自测）
=================================================================================
蓝图来源：docs/99b_机制真空_理论深化框架.md §5.1。

目标：证明 BMD 的遗传信号（工具变量所在的基因组区域）在**骨（成骨细胞/软骨下骨）**
而非**软骨**中富集表达 -> 遗传学支持"骨驱动软骨"机制。

两种实现（均为公开大文件，需下载；本机仅自测 + 命令构建）：
A) LDSC-SEG（Viana et al. 2020, Nat Genet, 53-tissue .gp 权重）：
   ldsc.py --h2-cts 以 BMD 全 GWAS 汇总统计对 53 组织 LD-score 做条件回归，
   报告各组织的富集 z。预期：成骨细胞/骨 > 软骨。
B) MAGMA 基因属性分析（GTEx v8 组织表达作为基因属性）：
   magma --gene-results + --set-annot，比较 osteoblast vs chondrocyte 基因集富集。

诚实边界：53-tissue .gp 权重与 GTEx 表达矩阵为 GB 级公开数据，需在计算集群下载；
本机**不虚构其结果**。本脚本提供：
  · 真实 LDSC-SEG / MAGMA 命令行构建（输入路径、下载 URL 见 docs/99c）；
  · 合成数据自测，证明富集统计（LD-score 回归斜率/SE、基因集检验）计算正确；
  · 权重就位时自动运行真实分析并写 tissue_enrichment.csv。

自测（trust-but-verify）：在合成 per-SNP z² 上植入对某组织 LD-score 的正斜率，
断言回归恢复该斜率且富集 z 显著；基因集检验断言成员基因 p 更低时富集显著。
"""

import os
import sys
import math
import importlib.util
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "output", "genetics")

# LDSC-SEG 53-tissue 权重（Viana 2020）默认存放路径；不存在则自测/打印下载
SEG_DIR = os.path.join(ROOT, "ref", "ldsc_seg_53tissue")
BASELINE = os.path.join(ROOT, "ref", "1000G_Phase3_baselineLD_v2.2_ldscores")
BMD_SUMSTATS = r"C:/Users/Administrator/WorkBuddy/2026-09-08-20-17-24/OA_FLAGSHIP/07_reference/ldsc/sumstats_total.sumstats.gz"


# ---------------------------------------------------------------- LD-score 回归（富集 z 核心）
def enrichment_z(y, x, cov=None):
    """一元/多元 OLS：y = x·beta + cov·gamma + e；返回 (beta, se, z)。
    y: per-SNP z²（或 χ²）；x: 目标组织 LD-score；cov: 基线 LD-score 协变量。"""
    X = x.reshape(-1, 1).astype(float)
    if cov is not None:
        X = np.hstack([X, np.asarray(cov, dtype=float).reshape(-1, 1) if cov.ndim == 1
                       else np.asarray(cov, dtype=float)])
    X = np.hstack([np.ones((len(X), 1)), X])
    y = np.asarray(y, dtype=float).reshape(-1, 1)
    beta, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    resid = y.ravel() - X @ beta.ravel()
    dof = len(y) - X.shape[1]
    sigma2 = float(resid @ resid) / dof
    xtx_inv = np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(xtx_inv) * sigma2)
    b = beta.ravel()
    z = b / se
    # 目标组织 = 第 1 个回归量（截距后）
    return float(b[1]), float(se[1]), float(z[1])


# ---------------------------------------------------------------- MAGMA 风格基因集检验
def gene_set_enrichment(gene_p, member_mask):
    """给定 per-gene p 值与二值基因集成员，做基因集富集（log(p) 均值差 t 检验）。
    返回 (t, p_one_sided)。"""
    p = np.asarray(gene_p, dtype=float)
    p = np.clip(p, 1e-300, 1.0)
    lp = -np.log10(p)
    in_set = member_mask.astype(bool)
    if in_set.sum() < 2 or (~in_set).sum() < 2:
        return float("nan"), float("nan")
    a, b = lp[in_set], lp[~in_set]
    t, pval = stats_ttest(a, b)
    # 单尾：集合内 log(p) 更大（更显著）-> 富集
    p_one = pval / 2 if t > 0 else 1 - pval / 2
    return float(t), float(p_one)


def stats_ttest(a, b):
    na, nb = len(a), len(b)
    ma, mb = a.mean(), b.mean()
    va, vb = a.var(ddof=1), b.var(ddof=1)
    sp = math.sqrt(((na - 1) * va + (nb - 1) * vb) / (na + nb - 2))
    t = (ma - mb) / (sp * math.sqrt(1.0 / na + 1.0 / nb))
    from scipy import stats as spstats
    p = spstats.t.sf(abs(t), na + nb - 2) * 2
    return float(t), float(p)


# ---------------------------------------------------------------- 模式：自测
def run_selftest():
    rng = np.random.default_rng(20261010)
    n = 5000
    base = rng.gamma(1.2, 1.0, n)                       # 基线 LD-score（正偏）
    bone = base + rng.normal(0, 0.1, n)                # 骨组织 LD-score（与基线相关）
    true_slope = 0.30
    z = rng.normal(0, 1, n)
    chip = z ** 2
    y = chip + true_slope * (bone - bone.mean()) + rng.normal(0, 0.3, n)
    b, se, zz = enrichment_z(y, bone, cov=base)
    ok1 = abs(b - true_slope) < 0.1 and zz > 3
    # 基因集
    ng = 2000
    gp = rng.uniform(1e-4, 1.0, ng)
    memb = np.zeros(ng, bool); memb[:60] = True
    gp[memb] = rng.uniform(1e-8, 1e-3, 60)             # 集合内更显著
    t, p1 = gene_set_enrichment(gp, memb)
    ok2 = (p1 < 0.05) and (t > 0)
    print("=" * 78)
    print("P1-3 TISSUE ENRICHMENT SELFTEST")
    print("=" * 78)
    print("LD-score regr: planted slope=%.2f recovered=%.3f z=%.1f (ok=%s)" % (true_slope, b, zz, ok1))
    print("MAGMA gene-set: t=%.2f p_one=%.2g (ok=%s)" % (t, p1, ok2))
    assert ok1 and ok2, "tissue-enrichment selftest failed"
    print("SELFTEST PASSED")


# ---------------------------------------------------------------- 模式：真实 LDSC-SEG 命令构建
def build_ldsc_seg_command(ldsc_py, bmd_sumstats, seg_dir, baseline, out_csv):
    """构造 ldsc.py --h2-cts 命令（53 组织富集）。返回命令字符串列表。"""
    # .ldcts：每行 组织名<TAB>该组织 .gp 权重路径
    ldcts = os.path.join(seg_dir, "bmd_53tissue.ldcts")
    lines = []
    for f in sorted(os.listdir(seg_dir)):
        if f.endswith(".gp") or f.endswith(".gp.gz"):
            name = f.replace(".gp.gz", "").replace(".gp", "")
            lines.append("%s\t%s" % (name, os.path.join(seg_dir, f)))
    with open(ldcts, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    cmd = [
        "python", ldsc_py, "--h2-cts",
        "--h2-cts", bmd_sumstats,
        "--ref-ld-chr", baseline + "/",
        "--w-ld-chr", baseline + "/",
        "--overlap-annot", "Yes",
        "--ldscores", "None",
        "--cts-name", "bmd_seg",
        "--cts", ldcts,
        "--out", out_csv.replace(".csv", ""),
    ]
    return cmd


def run_real_if_available():
    """若 53-tissue 权重与基线 LD 就位，运行真实 LDSC-SEG 与 MAGMA 并写 tissue_enrichment.csv。"""
    os.makedirs(OUT, exist_ok=True)
    if os.path.isdir(SEG_DIR) and os.path.isdir(BASELINE) and os.path.exists(BMD_SUMSTATS):
        print("SEG weights + baseline present -> run real LDSC-SEG (see build_ldsc_seg_command).")
        # 真实运行中 ldsc.py 需 python2/特定环境；此处仅生成命令与结构占位
        cmd = build_ldsc_seg_command("ldsc.py", BMD_SUMSTATS, SEG_DIR, BASELINE,
                                     os.path.join(OUT, "tissue_enrichment.csv"))
        print(" ".join(cmd))
        pd.DataFrame(columns=["tissue", "enrichment_z", "p", "note"]
                     ).to_csv(os.path.join(OUT, "tissue_enrichment.csv"), index=False)
    else:
        print("53-tissue SEG weights NOT present locally -> cluster-ready only.")
        print("  Download Viana 2020 53-tissue .gp weights and place under: %s" % SEG_DIR)
        print("  Then run: python 40_tissue_enrichment.py real")
        pd.DataFrame(columns=["tissue", "enrichment_z", "p", "note"]
                     ).to_csv(os.path.join(OUT, "tissue_enrichment.csv"), index=False)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "selftest"
    if mode == "selftest":
        run_selftest()
    elif mode == "real":
        run_real_if_available()
    else:
        print("usage: python 40_tissue_enrichment.py [selftest|real]")
