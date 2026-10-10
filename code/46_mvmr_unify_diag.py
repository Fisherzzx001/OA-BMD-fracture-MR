#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
46_mvmr_unify_diag.py — M1 根因诊断：为何 MVMR 管线的 BMD 工具数 < 主 MR 工具数
=============================================================================
背景（复审报告 102 · M1）
------------------------
同一对 BMD→OA_knee，稿件内部并存三套互不相等的「单变量基线」：
  · Table 3 主 MR β        （31/25/21/16/24 个严格 clump 工具）
  · S9 的 uni_beta         （在 K = BMD + 身高 个工具上做的 IVW，随 n_h 漂移）
  · MVMR BMD-only β        （27/23/19/15/21 个工具）

本脚本回答三个问题，不改动任何产出：
  Q1 每个暴露丢了几个工具？分别因 (a) 不在身高 GWAS (b) 回文链歧义 (c) 比对失败 而丢？
  Q2 `uni_beta` 的定义到底是什么？为何随 n_h 变化？
  Q3 若把回文 SNP 按「等位基因频率可判定」恢复（与主 MR 一致），能回收几个？

输出
----
  output/genetics/mvmr_unify_diag.csv   逐 SNP 丢因
  stdout 汇总表
=============================================================================
"""
import os, gzip, math, csv
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "output", "genetics")
MUNG = "C:/Users/Administrator/WorkBuddy/2026-09-08-20-17-24/OA_FLAGSHIP/07_reference/ldsc"
HG = os.path.join(ROOT, "..", "data", "height_gwas")
HEIGHT_EUR = "GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR.gz"
GW_P, CLUMP_KB = 5e-8, 1_000_000
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}

PAIRS = [
    ("BMD_total", "OA_knee", "sumstats_total.sumstats.gz", "sumstats_oa_knee.sumstats.gz"),
    ("BMD_femur_neck", "OA_knee", "sumstats_femur_neck.sumstats.gz", "sumstats_oa_knee.sumstats.gz"),
    ("BMD_lumbar", "OA_knee", "sumstats_lumbar.sumstats.gz", "sumstats_oa_knee.sumstats.gz"),
    ("BMD_arms", "OA_knee", "sumstats_arms.sumstats.gz", "sumstats_oa_knee.sumstats.gz"),
    ("BMD_leg", "OA_knee", "sumstats_leg.sumstats.gz", "sumstats_oa_knee.sumstats.gz"),
]
HN = {"BMD_total": "BMD_total__OA_knee.csv",
      "BMD_femur_neck": "BMD_femur_neck__OA_knee.csv",
      "BMD_lumbar": "BMD_lumbar__OA_knee.csv",
      "BMD_arms": "BMD_arms__OA_knee.csv",
      "BMD_leg": "BMD_leg__OA_knee.csv"}
# 主 MR（Table 3）已报告值，用于校验
REF = {"BMD_total": (0.02312, 0.00369, 31),
       "BMD_femur_neck": (0.01709, 0.00500, 25),
       "BMD_lumbar": (0.02791, 0.00483, 21),
       "BMD_arms": (0.02218, 0.00464, 16),
       "BMD_leg": (0.01772, 0.00418, 24)}


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
    print("   [munge] %-28s %d SNPs" % (fn, len(d)))
    return d


def load_height(fn):
    d = {}
    with gzip.open(os.path.join(HG, fn), "rt") as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 11 or not p[1].startswith("rs"):
                continue
            try:
                d[p[1]] = (int(p[2]), int(p[3]), p[4].upper(), p[5].upper(),
                           float(p[6]), float(p[7]), float(p[8]), float(p[9]), float(p[10]))
            except ValueError:
                continue
    print("   [height] %-46s %d SNPs" % (fn, len(d)))
    return d


def ivw(b_e, b_o, s_o):
    w = 1.0 / s_o ** 2
    den = float(np.sum(w * b_e ** 2))
    return float(np.sum(w * b_e * b_o) / den), float(math.sqrt(1.0 / den))


def align(a1r, a2r, a1, a2):
    if a1 == a1r and a2 == a2r:
        return 1.0
    if a1 == a2r and a2 == a1r:
        return -1.0
    if COMP.get(a1, "") == a1r and COMP.get(a2, "") == a2r:
        return 1.0
    if COMP.get(a1, "") == a2r and COMP.get(a2, "") == a1r:
        return -1.0
    return None


def align_by_freq(a1e, a2e, ea, oa, maf_e, maf_o, tol=0.10):
    """回文 SNP：若两侧 MAF 明显偏离 0.5（或一侧接近 0/1 而另一侧接近 1/0），
    按频率判定方向（TwoSampleMR 的标准做法）。返回 (g, 判定依据)。"""
    # 情形 1：两者同侧高频 → 同向
    if abs(maf_e - maf_o) < tol and (abs(maf_e - 0.5) > tol):
        return 1.0, "same-side (maf %.3f vs %.3f, |Δ|=%.3f)" % (maf_e, maf_o, abs(maf_e - maf_o))
    # 情形 2：互补（一个低频一个高频，且和为 ~1）→ 反向
    if abs((maf_e + maf_o) - 1.0) < tol and abs(maf_e - 0.5) > tol:
        return -1.0, "complementary (maf %.3f + %.3f = %.3f)" % (maf_e, maf_o, maf_e + maf_o)
    return None, "undecidable (maf %.3f vs %.3f)" % (maf_e, maf_o)


def main():
    dh = load_height(HEIGHT_EUR)
    rows = []
    print()
    print("=" * 96)
    print("Q1/Q3  逐暴露工具丢失诊断")
    print("=" * 96)
    for (e, o, me, mo) in PAIRS:
        pc = os.path.join(OUT, "harmonised", HN[e])
        h_rows = list(csv.DictReader(open(pc, encoding="utf-8")))
        de, do = load_munge(me), load_munge(mo)
        ref_b, ref_s, ref_n = REF[e]
        assert len(h_rows) == ref_n, "%s: harmonised %d vs REF %d" % (e, len(h_rows), ref_n)

        n_miss_h = n_amb = n_amb_ok = n_amb_fail = n_align = 0
        for r in h_rows:
            s = r["SNP"]
            rec = dict(exposure=e, SNP=s, CHR=r["CHR"], BP=r["BP"])
            if s not in de:
                rec.update(drop="NO_MUNGE_EXP"); rows.append(rec); continue
            if s not in do:
                rec.update(drop="NO_MUNGE_OUT"); rows.append(rec); continue
            ce, bpe, a1e, a2e, mafe, ne, ze = de[s]
            co, bpo, a1o, a2o, mafo, no_, zo = do[s]
            in_h = s in dh
            rec["in_height_gwas"] = in_h
            rec["palindromic"] = (COMP.get(a1e, None) == a2e)
            if not in_h:
                n_miss_h += 1
                rec.update(drop="NOT_IN_HEIGHT_GWAS")
                rows.append(rec)
                continue
            ch, bph, eah, oah, frh, bh, sh, ph, nh = dh[s]
            go = align(a1e, a2e, a1o, a2o)
            if go is None:
                n_align += 1
                rec.update(drop="ALIGN_FAIL_OUTCOME")
                rows.append(rec)
                continue
            if COMP.get(a1e, None) == a2e:
                n_amb += 1
                g2, why = align_by_freq(a1e, a2e, eah, oah, mafe, frh)
                rec["freq_rescue"] = why
                if g2 is None:
                    n_amb_fail += 1
                    rec.update(drop="AMBIGUOUS_UNDROPPABLE")
                    rows.append(rec)
                    continue
                n_amb_ok += 1
                rec.update(drop="AMBIGUOUS_RESCUABLE", rescue_g=g2)
                rows.append(rec)
                continue
            gh = align(a1e, a2e, eah, oah)
            if gh is None:
                n_align += 1
                rec.update(drop="ALIGN_FAIL_HEIGHT")
                rows.append(rec)
                continue
            rec.update(drop="KEEP")
            rows.append(rec)

        kept = sum(1 for r in rows if r["exposure"] == e and r["drop"] == "KEEP")
        rescuable = sum(1 for r in rows if r["exposure"] == e and r["drop"] == "AMBIGUOUS_RESCUABLE")
        print("  %-15s 主MR=%2d  MVMR现用=%2d  丢失=%2d  "
              "[不在身高GWAS %d | 回文可救 %d | 回文不可救 %d | 比对失败 %d]  → 恢复回文后可用 %d"
              % (e, ref_n, kept, ref_n - kept, n_miss_h, n_amb_ok, n_amb_fail, n_align,
                 kept + rescuable))

    pd.DataFrame(rows).to_csv(os.path.join(OUT, "mvmr_unify_diag.csv"),
                              index=False, encoding="utf-8")
    print("\nwritten mvmr_unify_diag.csv rows=%d" % len(rows))

    print()
    print("=" * 96)
    print("Q2  uni_beta 的定义：用 harmonised 值在不同工具子集上重算 IVW")
    print("=" * 96)
    print("  %-15s %-28s %4s %12s %12s" % ("exposure", "instrument set", "K", "beta", "se"))
    for (e, o, me, mo) in PAIRS:
        pc = os.path.join(OUT, "harmonised", HN[e])
        h_rows = list(csv.DictReader(open(pc, encoding="utf-8")))
        be = np.array([float(r["beta_e"]) for r in h_rows])
        bo = np.array([float(r["beta_o"]) for r in h_rows])
        so = np.array([float(r["se_o"]) for r in h_rows])
        rs = [r["SNP"] for r in h_rows]
        ref_b, ref_s, ref_n = REF[e]
        b1, s1 = ivw(be, bo, so)
        sub = [i for i, r in enumerate(rows) if r["exposure"] == e
               and r["drop"] in ("KEEP", "AMBIGUOUS_RESCUABLE")]
        b2, s2 = ivw(be[sub], bo[sub], so[sub])
        print("  %-15s %-28s %4d %12.5f %12.5f" % (e, "全 31/25/21/16/24 (=Table3)", len(be), b1, s1))
        print("  %-15s %-28s %4d %12.5f %12.5f" % ("", "MVMR 恢复回文后的子集", len(sub), b2, s2))
        print("  %-15s %-28s %4s %12.5f %12.5f" % ("", "Table 3 已报告", "-", ref_b, ref_s))
        print("  %-15s  差异 vs Table3: 全集 %+.2f%%   子集 %+.2f%%"
              % ("", 100 * (b1 - ref_b) / ref_b, 100 * (b2 - ref_b) / ref_b))
        print()


if __name__ == "__main__":
    main()
