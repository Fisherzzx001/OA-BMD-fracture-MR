#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
47_mvmr_unify_recovery.py — M1 统一管线：按 chr:bp 找回身高 GWAS 中「rsID 缺失」的 BMD 工具
=============================================================================
46_ 的诊断结论
--------------
· 主 MR → MVMR 丢失的 12 个工具 **全部**是 NOT_IN_HEIGHT_GWAS（回文 0、比对失败 0）：
      BMD_total 31→27(丢4) | femur_neck 25→23(丢2) | lumbar 21→19(丢2)
      arms 16→15(丢1)     | leg 24→21(丢3)
· 用 harmonised 值在**完整**工具集上算 IVW，结果 == Table 3 已报告值（差 +0.01%）
      BMD_total: 31 个 → 0.02312（Table 3 = 0.02312）
      27 子集   → 0.01996（−13.69%）

假设
----
46_ 的 height loader 要求 `p[1].startswith("rs")`，会跳过所有 rsID 为空/非 rs 前缀的行。
Yengo 2022 EUR 文件含 CHR / BP 列，本项目自身的等位基因 harmonise 亦按 chr:bp 进行。
因此「rsID 找不到」≠「该位点真的不在文件中」。

本脚本
------
1. 重新加载身高 GWAS：**不**过滤 rsID，同时建立 rsid 索引与 (chr,bp) 索引；
2. 对 46_ 标记的 12 个 NOT_IN_HEIGHT_GWAS 位点，按 chr:bp 精确查回；
3. 报告回收率与每个位点的身高 β/SE/P；
4. 复核：补全后 BMD 工具数是否 == 主 MR 工具数（31/25/21/16/24）。

输出
----
  output/genetics/mvmr_unify_recovery.csv
=============================================================================
"""
import os, gzip, csv
import pandas as pd

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "output", "genetics")
MUNG = "C:/Users/Administrator/WorkBuddy/2026-09-08-20-17-24/OA_FLAGSHIP/07_reference/ldsc"
HG = os.path.join(ROOT, "..", "data", "height_gwas")
HEIGHT_EUR = "GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR.gz"

HN = {"BMD_total": "BMD_total__OA_knee.csv",
      "BMD_femur_neck": "BMD_femur_neck__OA_knee.csv",
      "BMD_lumbar": "BMD_lumbar__OA_knee.csv",
      "BMD_arms": "BMD_arms__OA_knee.csv",
      "BMD_leg": "BMD_leg__OA_knee.csv"}
REF_N = {"BMD_total": 31, "BMD_femur_neck": 25, "BMD_lumbar": 21,
         "BMD_arms": 16, "BMD_leg": 24}


def load_height_full(fn):
    """不按 rsID 过滤；同时建 rsid 与 (chr,bp) 两套索引。"""
    by_rs, by_pos, n_row = {}, {}, 0
    with gzip.open(os.path.join(HG, fn), "rt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        print("   header:", hdr[:11])
        for line in f:
            p = line.rstrip("\n").split("\t")
            n_row += 1
            if len(p) < 11:
                continue
            try:
                c, bp = int(p[2]), int(p[3])
                ea, oa = p[4].upper(), p[5].upper()
                fr, be, se, pv, n = (float(p[6]), float(p[7]),
                                     float(p[8]), float(p[9]), float(p[10]))
            except (ValueError, IndexError):
                continue
            rs = p[1]
            rec = (c, bp, ea, oa, fr, be, se, pv, n, rs)
            if rs.startswith("rs"):
                by_rs.setdefault(rs, rec)
            by_pos.setdefault((c, bp), rec)
    print("   行数 %d | rsid 索引 %d | (chr,bp) 索引 %d" % (n_row, len(by_rs), len(by_pos)))
    return by_rs, by_pos


def load_munge(fn):
    d = {}
    with gzip.open(os.path.join(MUNG, fn), "rt") as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 8:
                continue
            try:
                d[p[2]] = (int(p[0]), int(p[1]), p[3].upper(), p[4].upper(),
                           float(p[5]), float(p[6]), float(p[7]))
            except ValueError:
                continue
    return d


def main():
    print("=" * 92)
    print("重新加载身高 GWAS（不按 rsID 过滤）")
    print("=" * 92)
    by_rs, by_pos = load_height_full(HEIGHT_EUR)

    diag = pd.read_csv(os.path.join(OUT, "mvmr_unify_diag.csv"))
    miss = diag[diag["drop"] == "NOT_IN_HEIGHT_GWAS"]
    print("\n46_ 标记的 NOT_IN_HEIGHT_GWAS 位点数：%d" % len(miss))
    print("  按暴露：", miss.groupby("exposure").size().to_dict())

    rows = []
    print("\n" + "=" * 92)
    print("按 chr:bp 查回")
    print("=" * 92)
    print("  %-15s %-12s %6s %11s | %-10s %10s %9s %10s"
          % ("exposure", "SNP", "chr:bp", "rsID命中", "位置命中", "β_height", "SE", "P_height"))
    n_found = 0
    for r in miss.itertuples():
        c, bp = int(r.CHR), int(r.BP)
        hit_rs = r.SNP in by_rs
        hit_pos = by_pos.get((c, bp))
        rec = dict(exposure=r.exposure, SNP=r.SNP, CHR=c, BP=bp,
                   found_by_rsid=hit_rs, found_by_pos=hit_pos is not None)
        if hit_pos:
            n_found += 1
            ch, bph, eah, oah, frh, bh, sh, ph, nh, rs_alt = hit_pos
            rec.update(pos_rsid=rs_alt, beta_h=bh, se_h=sh, p_h=ph, eaf_h=frh,
                       ea_h=eah, oa_h=oah)
            print("  %-15s %-12s %d:%d %11s | %-10s %10.5f %9.5f %10.2e"
                  % (r.exposure, r.SNP, c, bp, hit_rs, "YES", bh, sh, ph))
        else:
            rec.update(pos_rsid=None, beta_h=None, se_h=None, p_h=None)
            print("  %-15s %-12s %d:%d %11s | %-10s %10s %9s %10s"
                  % (r.exposure, r.SNP, c, bp, hit_rs, "NO", "-", "-", "-"))
        rows.append(rec)

    print("\n回收：%d / %d" % (n_found, len(miss)))

    print("\n" + "=" * 92)
    print("补全后的 BMD 工具数 vs 主 MR（Table 3）")
    print("=" * 92)
    print("  %-15s %8s %10s %10s %8s" % ("exposure", "主MR", "rsID命中", "+位置命中", "差额"))
    for e, fn in HN.items():
        h_rows = list(csv.DictReader(open(os.path.join(OUT, "harmonised", fn), encoding="utf-8")))
        n_all = len(h_rows)
        n_rs = sum(1 for r in h_rows if r["SNP"] in by_rs)
        n_tot = 0
        for r in h_rows:
            if r["SNP"] in by_rs:
                n_tot += 1
            elif (int(r["CHR"]), int(r["BP"])) in by_pos:
                n_tot += 1
        print("  %-15s %8d %10d %10d %8d" % (e, REF_N[e], n_rs, n_tot, n_tot - REF_N[e]))

    pd.DataFrame(rows).to_csv(os.path.join(OUT, "mvmr_unify_recovery.csv"),
                              index=False, encoding="utf-8")
    print("\nwritten mvmr_unify_recovery.csv rows=%d" % len(rows))


if __name__ == "__main__":
    main()
