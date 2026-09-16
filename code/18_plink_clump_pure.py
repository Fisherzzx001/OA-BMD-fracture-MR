#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
18_plink_clump_pure.py
=========================================================================
投稿前整改 #3（严格 clumping）——**纯 Python 复现 PLINK --clump**，无需安装 PLINK。

直接解析 PLINK 1 二进制 .bed，计算候选 SNP 间的 LD，执行与 PLINK 相同的
贪婪 clumping：

    plink --bfile <EUR> --clump <in> --clump-p1 5e-8 --clump-p2 5e-8 \
          --clump-r2 0.001 --clump-kb 10000 --out <exp>_clumped

算法：
  1. 各暴露 GWAS 中 P < 5e-8 的 SNP 作为候选（index 与 clump-mate 同源，p1=p2）；
  2. 候选 SNP 按 P 升序；逐个取未剔除者为 index；
  3. 同染色体、|Δbp| ≤ 10 Mb、且 r² ≥ 0.001 者被剔除；
  4. 保留下来的 index SNP 即严格独立仪器集。

.bed 解码（已经验判定，见 verify）：
  2-bit 编码 00=hom A1, 01=missing, 10=het, 11=hom A2
  （实测 code 01 频率 = 0.0000 → 即 missing；相邻 SNP r² 呈正常 LD 区块结构）

输出：<OUT>/<exposure>_clumped.clumped（含 SNP 列，供 17_ 直接读取）
      <OUT>/strict_clump_summary.csv（每个暴露的 GWS/匹配/保留数）
=========================================================================
"""
import os
import sys
import time
import importlib
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
m13 = importlib.import_module("13_mr_analysis")

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "output", "genetics")

# ---- 配置 ----------------------------------------------------------------
BFILE = os.environ.get("BFILE_1KG_EUR", r"H:/数据下载/1kg.v3/EUR")   # 不含扩展名
P_THRESH = 5e-8      # --clump-p1 / --clump-p2
R2 = 0.001           # --clump-r2
KB = 10_000_000      # --clump-kb（bp）

# PLINK 1 .bed：2-bit -> A1 剂量（00=homA1, 01=missing, 10=het, 11=homA2）
_BED_TABLE = np.array([[v & 3, (v >> 2) & 3, (v >> 4) & 3, (v >> 6) & 3]
                       for v in range(256)], dtype=np.uint8)
_DOSAGE = np.array([2.0, np.nan, 1.0, 0.0])


# ---------------------------------------------------------------- .bed 读取
def read_fam_n(path):
    with open(path) as f:
        return sum(1 for _ in f)


def read_bim_index(bim_path, wanted):
    found = {}
    row = 0
    for chunk in pd.read_csv(bim_path, sep=r"\s+", header=None,
                             usecols=[0, 1, 3], names=["CHR", "SNP", "BP"],
                             dtype={"CHR": "int16", "SNP": "str", "BP": "int64"},
                             chunksize=2_000_000):
        m = chunk["SNP"].isin(wanted).values
        if m.any():
            # 用「块内位置」定位，避免 .index 在分块读取下的连续/重置歧义
            for k in np.where(m)[0]:
                found[str(chunk["SNP"].values[k])] = (
                    int(chunk["CHR"].values[k]),
                    int(chunk["BP"].values[k]),
                    row + int(k))
        row += len(chunk)
    return found, row


def read_genotypes(bed_path, rows, n):
    """读取指定 bed 行号的基因型，返回 (len(rows), n) 的 A1 剂量矩阵。"""
    bpb = (n + 3) // 4
    rows = np.asarray(rows, dtype=np.int64)
    G = np.empty((len(rows), n), dtype=np.float64)
    buf = np.empty(bpb, dtype=np.uint8)
    with open(bed_path, "rb") as f:
        for k, i in enumerate(rows):
            f.seek(3 + i * bpb)
            buf[:] = np.frombuffer(f.read(bpb), dtype=np.uint8)
            G[k] = _DOSAGE[_BED_TABLE[buf].reshape(-1)[:n]]
    return G


def standardize(G):
    """均值填补缺失 + 标准化；单态（sd=0）置 0。"""
    G = np.where(np.isnan(G), np.nanmean(G, axis=1, keepdims=True), G)
    G = G - G.mean(axis=1, keepdims=True)
    sd = G.std(axis=1, keepdims=True)
    sd = np.where(sd < 1e-12, 1.0, sd)
    return (G / sd).astype(np.float32)


# ---------------------------------------------------------------- 贪婪 clump
def greedy_clump(sub, C, r2_thr, kb):
    """
    sub: DataFrame，含 row（矩阵行号）、CHR、BP、P
    C  : 相关系数矩阵（float32）
    返回保留行在 sub 中的位置索引
    """
    P = sub["P"].values
    CHR = sub["CHR"].values
    BP = sub["BP"].values
    removed = np.zeros(len(sub), dtype=bool)
    kept = []
    for i in np.argsort(P, kind="stable"):
        if removed[i]:
            continue
        kept.append(i)
        same_chr = (CHR == CHR[i]) & (~removed)
        near = np.abs(BP - BP[i]) <= kb
        cand = np.where(same_chr & near)[0]
        cand = cand[cand != i]
        if len(cand):
            # C 按 sub 的行序构建，故直接用位置索引
            r2 = C[i, cand] ** 2
            removed[cand[r2 >= r2_thr]] = True
    return kept


# ---------------------------------------------------------------- 主流程
def main():
    os.makedirs(OUT, exist_ok=True)
    bed_path, bim_path, fam_path = BFILE + ".bed", BFILE + ".bim", BFILE + ".fam"

    n = read_fam_n(fam_path)
    print("=" * 74)
    print("严格 PLINK clumping（纯 Python 复现）  r²<%.4f  /  %d kb  /  P<%.1e"
          % (R2, KB // 1000, P_THRESH))
    print("=" * 74)
    print("参考面板：%s   样本数 n=%d" % (BFILE, n))

    print("\n▶ 载入 GWAS 汇总统计 ...")
    t0 = time.time()
    dfs = m13.load_all()
    print("  完成（%.1fs），共 %d 个表型" % (time.time() - t0, len(dfs)))

    exposures = sorted({e for e, _, _ in m13.PAIRS})
    print("  暴露（%d 个）：%s" % (len(exposures), ", ".join(exposures)))

    # 1) 各暴露的 GWS 候选
    cand = {}
    for e in exposures:
        d = dfs[e]
        g = d[d["P"] < P_THRESH]
        cand[e] = g[["SNP", "CHR", "BP", "P"]].copy()
        print("   %-14s GWS(P<5e-8) = %d" % (e, len(g)))
    all_rs = set()
    for e in exposures:
        all_rs |= set(cand[e]["SNP"].astype(str))
    print("\n  并集候选 SNP：%d" % len(all_rs))

    # 2) 在 .bim 中定位
    print("▶ 扫描 EUR.bim（8.55M 变异）定位候选 ...")
    t0 = time.time()
    pos, n_var = read_bim_index(bim_path, all_rs)
    print("  匹配到参考面板：%d / %d（%.1f%%），用时 %.1fs；.bim 变异总数 %d"
          % (len(pos), len(all_rs), 100.0 * len(pos) / max(len(all_rs), 1),
             time.time() - t0, n_var))
    # 越界保护
    bad = [s for s, v in pos.items() if v[2] >= n_var]
    if bad:
        print("  ⚠ 剔除越界行号 %d 个" % len(bad))
        for s in bad:
            pos.pop(s)

    if not pos:
        print("✗ 无 SNP 匹配，终止。")
        return

    panel = pd.DataFrame([(s, v[0], v[1], v[2]) for s, v in pos.items()],
                         columns=["SNP", "CHR", "BP", "bed_row"])
    panel = panel.sort_values("bed_row").reset_index(drop=True)
    # gidx = 基因型矩阵 Gs 的行号（0..len(panel)-1）；bed_row 仅用于 .bed 定位
    panel["gidx"] = np.arange(len(panel), dtype=int)

    # 3) 读取基因型
    print("▶ 读取候选 SNP 基因型（%d 个）..." % len(panel))
    t0 = time.time()
    G = read_genotypes(bed_path, panel["bed_row"].values, n)
    miss = np.isnan(G).mean()
    print("  完成（%.1fs），平均缺失率 %.4f" % (time.time() - t0, miss))
    Gs = standardize(G)
    print("  单态 SNP（sd=0）：%d" % int((np.abs(Gs).max(axis=1) == 0).sum()))

    # 4) 逐暴露 clump
    print("\n" + "=" * 74)
    print("逐暴露 clumping")
    print("=" * 74)
    rows = []
    for e in exposures:
        ce = cand[e].copy()
        ce["SNP"] = ce["SNP"].astype(str)
        sub = ce.merge(panel[["SNP", "gidx"]], on="SNP", how="inner")
        if len(sub) == 0:
            print("  %-14s 无匹配 SNP，跳过" % e)
            continue
        sub = sub.rename(columns={"gidx": "row"}).reset_index(drop=True)
        # 相关性矩阵（本暴露涉及到的行）
        C = (Gs[sub["row"].values] @ Gs[sub["row"].values].T) / n
        kept = greedy_clump(sub, C, R2, KB)
        keep_snps = sub.loc[kept, "SNP"].tolist()
        print("  %-14s GWS=%4d  面板匹配=%4d  严格clump后保留=%4d"
              % (e, len(ce), len(sub), len(keep_snps)))

        out = sub.loc[kept, ["CHR", "SNP", "BP", "P"]].copy()
        out = out.sort_values(["CHR", "BP"])
        out.insert(4, "TOTAL", 0)
        out.insert(5, "NSIG", 0)
        for c in ["S05", "S01", "S001", "S0001", "SP2"]:
            out[c] = 0
        fp = os.path.join(OUT, f"{e}_clumped.clumped")
        out.to_csv(fp, sep="\t", index=False)
        rows.append(dict(exposure=e, n_gws=len(ce), n_in_panel=len(sub),
                         n_kept_strict=len(keep_snps)))

    pd.DataFrame(rows).to_csv(os.path.join(OUT, "strict_clump_summary.csv"), index=False)
    print("\n✅ 严格仪器集 -> %s/<exposure>_clumped.clumped" % OUT)
    print("✅ 汇总 -> %s/strict_clump_summary.csv" % OUT)


if __name__ == "__main__":
    main()
