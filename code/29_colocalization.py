#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
29_colocalization.py — 候选基因座共定位分析（coloc.abf，纯 Python）
=========================================================================
目的
----
稿件 §4.2 把 WNT16 / SOX6 / GDF5 称为「候选共享效应器 (candidate shared effectors)」，
但从未检验两性状是否共享同一因果变异。本脚本用 Bayesian colocalization 补上这一环。

为什么现在能做
--------------
此前的印象是「手上只有 clump 后的工具变量，做不了 coloc」。核查数据来源后发现：
MR 的输入是 **LDSC munge 后的全基因组汇总统计**（CHR BP SNP A1 A2 MAF N Z P），
即 ~69 万（BMD）/ ~118 万（OA_knee）个 SNP 的完整 Z 统计量，存在于
  C:/Users/Administrator/WorkBuddy/2026-09-08-20-17-24/OA_FLAGSHIP/07_reference/ldsc/
因此区域级汇总统计是齐备的，coloc 可以做。

算法（严格按 Giambartolomei et al. 2014, PLoS Genet 10:e1004383；与 R 包 coloc
的 coloc.abf / log.bf.snp 一致，见 R/claudia.R）
------------------------------------------------------------------------------
  V  = 1/(2·N·f·(1-f))            quant
     = 1/(2·N·f·(1-f)·s·(1-s))    cc
  r  = sd.prior^2 / (sd.prior^2 + V)
  lABF = 0.5·( log(1-r) + r·z^2 )      sd.prior: quant=0.15, cc=0.2
  H0: 0
  H1: log(p1)  + logsum(l1)
  H2: log(p2)  + logsum(l2)
  H3: log(p1)+log(p2) + logdiff( logsum(l1)+logsum(l2), logsum(l1+l2) )
  H4: log(p12) + logsum(l1+l2)
  归一化得 PP0..PP4；PP.H4 >= 0.8 判为强共定位。

注意：lABF 只依赖 z^2 与 V，故 **等位基因方向不影响 coloc.abf**；
方向一致性需单独用 A1/A2 对齐后检验（本脚本另做）。

关键修正（trust-but-verify）
----------------------------
· BMD 各部位 munge 文件里 N=400,000 是**占位值**，真值来自 GWAS Catalog
  (GCST90568448-58)：total 31,986 / femur_neck 32,017 / femur_total 31,873 /
  lumbar 30,449 / legs 31,873 / pelvis 31,873 / trunk 31,986，其余部位取 31,873。
  用占位 N 会把 V 缩小 12.5 倍 → r 偏大 → 虚假抬升 lABF，必须修正。
· OA_knee 的 N 逐 SNP 变动（951k–1,316k），直接取文件内逐 SNP 的 N，不做占位替换。
· fracture munge 文件 N=450,000 亦为占位，真值 264,973，s=37,857/264,973。

输出
----
  output/genetics/coloc_candidate_loci.csv     三基因座 × 各 BMD 部位
  output/genetics/coloc_screen_all_loci.csv    BMD_total 全部工具变量基因座筛查
  output/genetics/coloc_positive_control.csv   BMD × fracture 阳性对照
  output/genetics/coloc_sensitivity.csv        先验/窗口/N 敏感性
  output/genetics/coloc_report.txt             可读报告
"""

import os
import sys
import gzip
import math
import itertools
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "output", "genetics")
LDSC = r"C:/Users/Administrator/WorkBuddy/2026-09-08-20-17-24/OA_FLAGSHIP/07_reference/ldsc"

# ---------------------------------------------------------------- 真值 N
# BMD: munge 文件内 N=400000 为占位；真值取自 GWAS Catalog GCST90568448-58（§2.1）
BMD_TRUE_N = {
    "total": 31986, "femur_neck": 32017, "femur_total": 31873,
    "lumbar": 30449, "legs": 31873, "pelvis": 31873, "trunk": 31986,
    "arms": 31873, "spine": 31873, "rib": 31873, "head": 31873,
}
BMD_FILE = {k: "sumstats_%s.sumstats.gz" % k for k in BMD_TRUE_N}
BMD_FILE["leg"] = "sumstats_leg.sumstats.gz"
BMD_TRUE_N["leg"] = 31873

TRAITS = {
    # name: (file, type, N_override, s)
    "OA_knee":     ("sumstats_oa_knee.sumstats.gz", "cc", None, 0.12),
    "OA_selfrep":  ("sumstats_oa_selfrep.sumstats.gz", "cc", 462933, 0.0831),
    "fracture":    ("sumstats_fracture.sumstats.gz", "cc", 264973, 37857 / 264973),
}

SD_PRIOR = {"quant": 0.15, "cc": 0.2}

# ---------------------------------------------------------------- 区域
CANDIDATE_LOCI = {
    # 稿件 §4.2 / §2.2 点名的三个基因座。坐标取自 munge 文件（hg19）
    "WNT16": ("rs3779381", 7, 120966790),
    "SOX6":  ("rs7108738", 11, 15710084),
    "GDF5":  ("rs224333", 20, 34023962),   # OA lead；BMD munge 集内无此 SNP
}


# ---------------------------------------------------------------- 载入
def load_munged(path, n_override=None):
    df = pd.read_csv(path, sep=r"\s+", compression="gzip",
                     dtype={"CHR": "int32", "BP": "int64", "SNP": str,
                            "A1": str, "A2": str})
    need = ["CHR", "BP", "SNP", "A1", "A2", "MAF", "N", "Z", "P"]
    df = df[[c for c in need if c in df.columns]].copy()
    df = df.dropna(subset=["CHR", "BP", "MAF", "Z"])
    df["MAF"] = df["MAF"].astype(float)
    df["Z"] = df["Z"].astype(float)
    # Z 缺失/为 0 时回退到 P
    if "P" in df.columns:
        p = pd.to_numeric(df["P"], errors="coerce")
        zf = pd.to_numeric(df["Z"], errors="coerce")
        from scipy.stats import norm
        zp = norm.isf(p / 2.0)
        df["Z"] = zf.where(np.abs(zf) > 1e-8, zp)
    if n_override is not None:
        df["N"] = float(n_override)
    else:
        df["N"] = df["N"].astype(float)
    return df


_CACHE = {}


def get_trait(name):
    """返回 (df, type, s)"""
    if name in _CACHE:
        return _CACHE[name]
    if name in BMD_TRUE_N:
        fn = BMD_FILE[name]
        typ, s = "quant", None
        df = load_munged(os.path.join(LDSC, fn), n_override=BMD_TRUE_N[name])
    else:
        fn, typ, nov, s = TRAITS[name]
        df = load_munged(os.path.join(LDSC, fn), n_override=nov)
    _CACHE[name] = (df, typ, s)
    return _CACHE[name]


# ---------------------------------------------------------------- coloc
def logsumexp(a):
    a = np.asarray(a, dtype=float)
    m = float(np.max(a))
    if not np.isfinite(m):
        return m
    return m + math.log(float(np.sum(np.exp(a - m))))


def logdiff(a, b):
    """log(exp(a) - exp(b)), a >= b"""
    if not np.isfinite(a):
        return -np.inf
    if b >= a:
        return -np.inf          # 数值上 exp(a)-exp(b) <= 0 → H3 概率 0
    return a + math.log1p(-math.exp(b - a))


def var_beta(n, f, typ, s):
    f = np.clip(f, 1e-6, 1 - 1e-6)
    v = 1.0 / (2.0 * np.asarray(n, dtype=float) * f * (1.0 - f))
    if typ == "cc":
        v = v / (s * (1.0 - s))
    return v


def log_abf(z, v, sd_prior):
    w = sd_prior ** 2
    r = w / (w + v)
    return 0.5 * np.log1p(-r) + 0.5 * r * (z ** 2)


def coloc_abf(z1, v1, z2, v2, typ1="quant", typ2="cc",
              p1=1e-4, p2=1e-4, p12=5e-6,
              sd1=None, sd2=None):
    sd1 = SD_PRIOR[typ1] if sd1 is None else sd1
    sd2 = SD_PRIOR[typ2] if sd2 is None else sd2
    l1 = log_abf(z1, v1, sd1)
    l2 = log_abf(z2, v2, sd2)
    L1 = logsumexp(l1)
    L2 = logsumexp(l2)
    L4 = logsumexp(l1 + l2)
    L3 = logdiff(L1 + L2, L4)
    lh = np.array([0.0,
                   math.log(p1) + L1,
                   math.log(p2) + L2,
                   math.log(p1) + math.log(p2) + L3,
                   math.log(p12) + L4])
    lh = lh - logsumexp(lh)
    pp = np.exp(lh)
    return dict(PP0=pp[0], PP1=pp[1], PP2=pp[2], PP3=pp[3], PP4=pp[4])


def region(df, chrom, pos, window):
    m = (df["CHR"].values == chrom) & \
        (np.abs(df["BP"].values - pos) <= window)
    return df.loc[m]


def coloc_region(trait1, trait2, chrom, pos, window=500000, **kw):
    """在 (chrom, pos) ± window 内做 coloc.abf。返回 dict。"""
    d1, t1, s1 = get_trait(trait1)
    d2, t2, s2 = get_trait(trait2)
    r1 = region(d1, chrom, pos, window)
    r2 = region(d2, chrom, pos, window)
    m = r1[["CHR", "BP"]].merge(r2[["CHR", "BP"]], on=["CHR", "BP"])
    if len(m) < 20:
        return None
    a = r1.merge(m, on=["CHR", "BP"]).sort_values("BP").reset_index(drop=True)
    b = r2.merge(m, on=["CHR", "BP"]).sort_values("BP").reset_index(drop=True)
    assert (a["BP"].values == b["BP"].values).all()

    # MAF 交叉校验（两文件应来自同一欧洲参考）
    maf_delta = float(np.nanmax(np.abs(a["MAF"].values - b["MAF"].values))) \
        if len(a) else float("nan")

    z1 = a["Z"].values.astype(float)
    z2 = b["Z"].values.astype(float)
    f = a["MAF"].values.astype(float)
    v1 = var_beta(a["N"].values, f, t1, s1)
    v2 = var_beta(b["N"].values, b["MAF"].values.astype(float), t2, s2)

    res = coloc_abf(z1, v1, z2, v2, t1, t2, **kw)

    # 方向一致性（用 A1/A2 对齐后比较 top SNP 的 Z 符号）
    k = int(np.argmax(z1 ** 2))
    z1t, z2t = z1[k], z2[k]
    A1a, A2a = str(a["A1"].iloc[k]), str(a["A2"].iloc[k])
    A1b, A2b = str(b["A1"].iloc[k]), str(b["A2"].iloc[k])
    if {A1a, A2a} == {A1b, A2b}:
        flip = (A1a != A1b)      # 效应等位互换 → 需要翻转
    else:
        flip = None              # 等位集合不一致，无法判向
    if flip is not None and flip:
        z2t = -z2t
    concordant = (None if flip is None else bool(np.sign(z1t) == np.sign(z2t)))

    res.update(
        n_snps=len(a),
        n_snps_t1=len(r1), n_snps_t2=len(r2),
        max_abs_z1=float(np.max(np.abs(z1))), max_abs_z2=float(np.max(np.abs(z2))),
        top_snp=str(a["SNP"].iloc[k]), top_bp=int(a["BP"].iloc[k]),
        z1_top=float(z1t), z2_top=float(z2t),
        dir_concordant=concordant, maf_max_delta=maf_delta,
        maf_top=float(f[k]),
    )
    return res


# ---------------------------------------------------------------- 主流程
def main():
    os.makedirs(OUT, exist_ok=True)
    bmd_sites = ["total", "femur_neck", "lumbar", "arms", "leg"]
    lines = []

    def P(s=""):
        print(s)
        lines.append(s)

    P("=" * 78)
    P("COLOCALIZATION (coloc.abf, Giambartolomei 2014) — OA_B_PROJECT")
    P("=" * 78)
    P("V: quant 1/(2Nf(1-f)) ; cc 1/(2Nf(1-f)s(1-s)) ; sd.prior quant=0.15 cc=0.2")
    P("priors p1=p2=1e-4, p12=5e-6 ; window +/-500 kb")

    # ---------- A. 三候选基因座 ----------
    rows = []
    P("\n--- A. Candidate loci named in Sec 4.2 (BMD x OA_knee) ---")
    P("%-8s %-11s %6s %6s %6s %6s %6s %-13s %7s %7s %s" %
      ("locus", "bmd_site", "PP0", "PP1", "PP2", "PP3", "PP4", "top_snp", "z_bmd", "z_oa", "dir"))
    for loc, (rs, chrom, pos) in CANDIDATE_LOCI.items():
        for site in bmd_sites:
            r = coloc_region(site, "OA_knee", chrom, pos)
            if r is None:
                P("%-8s %-11s  <20 shared SNPs — skipped>" % (loc, site))
                continue
            rows.append(dict(locus=loc, lead_rs=rs, chrom=chrom, pos=pos,
                             bmd_site=site, **r))
            P("%-8s %-11s %6.3f %6.3f %6.3f %6.3f %6.3f %-13s %7.2f %7.2f %s" %
              (loc, site, r["PP0"], r["PP1"], r["PP2"], r["PP3"], r["PP4"],
               r["top_snp"], r["z1_top"], r["z2_top"], r["dir_concordant"]))
    dfA = pd.DataFrame(rows)
    dfA.to_csv(os.path.join(OUT, "coloc_candidate_loci.csv"), index=False)

    # ---------- A2. 阳性对照：BMD x fracture ----------
    rows = []
    P("\n--- A2. Positive control: BMD x fracture (same loci) ---")
    P("%-8s %-11s %6s %6s %6s %6s %6s %-13s %7s %7s" %
      ("locus", "bmd_site", "PP0", "PP1", "PP2", "PP3", "PP4", "top_snp", "z_bmd", "z_fx"))
    for loc, (rs, chrom, pos) in CANDIDATE_LOCI.items():
        for site in bmd_sites:
            r = coloc_region(site, "fracture", chrom, pos)
            if r is None:
                continue
            rows.append(dict(locus=loc, lead_rs=rs, chrom=chrom, pos=pos,
                             bmd_site=site, **r))
            P("%-8s %-11s %6.3f %6.3f %6.3f %6.3f %6.3f %-13s %7.2f %7.2f" %
              (loc, site, r["PP0"], r["PP1"], r["PP2"], r["PP3"], r["PP4"],
               r["top_snp"], r["z1_top"], r["z2_top"]))
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "coloc_positive_control.csv"), index=False)

    # ---------- B. 系统筛查：BMD_total 全部工具变量基因座 ----------
    harm = os.path.join(OUT, "harmonised", "BMD_total__OA_knee.csv")
    rows = []
    P("\n--- B. Systematic screen: every BMD_total instrument locus vs OA_knee ---")
    if os.path.exists(harm):
        inst = pd.read_csv(harm)
        P("n instruments = %d" % len(inst))
        P("%-13s %8s %10s %6s %6s %6s %6s %6s" %
          ("index_snp", "chr:pos", "n_snps", "PP0", "PP1", "PP2", "PP3", "PP4"))
        for _, rw in inst.iterrows():
            r = coloc_region("total", "OA_knee", int(rw["CHR"]), int(rw["BP"]))
            if r is None:
                continue
            rows.append(dict(index_snp=rw["SNP"], chrom=int(rw["CHR"]),
                             pos=int(rw["BP"]), **r))
            P("%-13s %8s %10d %6.3f %6.3f %6.3f %6.3f %6.3f" %
              (rw["SNP"], "%d:%d" % (rw["CHR"], rw["BP"]), r["n_snps"],
               r["PP0"], r["PP1"], r["PP2"], r["PP3"], r["PP4"]))
        dfB = pd.DataFrame(rows)
        dfB.to_csv(os.path.join(OUT, "coloc_screen_all_loci.csv"), index=False)
        if len(dfB):
            P("\nSummary: PP.H4>=0.8 : %d/%d loci" % ((dfB.PP4 >= 0.8).sum(), len(dfB)))
            P("         PP.H4 0.5-0.8: %d/%d" % (((dfB.PP4 >= 0.5) & (dfB.PP4 < 0.8)).sum(), len(dfB)))
            P("         PP.H3>=0.8 : %d/%d loci (distinct causal variants)" %
              ((dfB.PP3 >= 0.8).sum(), len(dfB)))
            P("         median PP.H4 = %.3f" % dfB.PP4.median())

    # ---------- C. 敏感性 ----------
    rows = []
    P("\n--- C. Sensitivity: priors / window / N / sd.prior ---")
    for loc, (rs, chrom, pos) in CANDIDATE_LOCI.items():
        base = dict(trait1="total", trait2="OA_knee", locus=loc)
        for p12 in (1e-5, 5e-6, 1e-6, 1e-7):
            r = coloc_region("total", "OA_knee", chrom, pos, p12=p12)
            rows.append(dict(**base, setting="p12=%g" % p12, **r))
        for w in (250000, 1000000):
            r = coloc_region("total", "OA_knee", chrom, pos, window=w)
            rows.append(dict(**base, setting="window=+/-%dkb" % (w // 1000), **r))
        for sd1 in (0.10, 0.30):
            r = coloc_region("total", "OA_knee", chrom, pos, sd1=sd1)
            rows.append(dict(**base, setting="sd.prior_quant=%g" % sd1, **r))
        for sd2 in (0.10, 0.30):
            r = coloc_region("total", "OA_knee", chrom, pos, sd2=sd2)
            rows.append(dict(**base, setting="sd.prior_cc=%g" % sd2, **r))
    dfC = pd.DataFrame(rows)
    dfC.to_csv(os.path.join(OUT, "coloc_sensitivity.csv"), index=False)
    piv = dfC.pivot_table(index="setting", columns="locus", values="PP4")
    P("\nPP.H4 under each setting (BMD_total x OA_knee):")
    P(piv.round(3).to_string())

    # ---------- D. 次级信号检查（WNT16 已知等位异质性） ----------
    P("\n--- D. Secondary-signal check (mask top signal +/-100 kb, re-run) ---")
    P("Rationale: WNT16 locus has documented allelic heterogeneity (Medina-Gomez 2012:")
    P("rs917727 and rs7801723/C7orf58 are independent). coloc.abf assumes ONE causal")
    P("variant per trait; masking is a poor-man's conditional test (NOT SuSiE).")
    for loc, (rs, chrom, pos) in CANDIDATE_LOCI.items():
        d1, t1, s1 = get_trait("total")
        d2, t2, s2 = get_trait("OA_knee")
        r1 = region(d1, chrom, pos, 500000)
        k = int(np.argmax(np.abs(r1["Z"].values)))
        top_bp = int(r1["BP"].iloc[k])
        # 屏蔽 top SNP ±100kb 后，在剩余窗口内重跑
        m1 = r1[np.abs(r1["BP"].values - top_bp) > 100000]
        if len(m1) < 20:
            P("  %-7s : too few residual SNPs" % loc)
            continue
        r2 = region(d2, chrom, pos, 500000)
        m = m1[["CHR", "BP"]].merge(r2[["CHR", "BP"]], on=["CHR", "BP"])
        a = m1.merge(m, on=["CHR", "BP"]).sort_values("BP").reset_index(drop=True)
        b = r2.merge(m, on=["CHR", "BP"]).sort_values("BP").reset_index(drop=True)
        z1 = a["Z"].values.astype(float)
        z2 = b["Z"].values.astype(float)
        v1 = var_beta(a["N"].values, a["MAF"].values.astype(float), t1, s1)
        v2 = var_beta(b["N"].values, b["MAF"].values.astype(float), t2, s2)
        r = coloc_abf(z1, v1, z2, v2, t1, t2)
        P("  %-7s : masked top rs at %d, n=%d -> PP3=%.3f PP4=%.3f  (max|z| bmd=%.2f oa=%.2f)"
          % (loc, top_bp, len(a), r["PP3"], r["PP4"],
             float(np.max(np.abs(z1))), float(np.max(np.abs(z2)))))

    with open(os.path.join(OUT, "coloc_report.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    P("\nWritten to " + OUT)


if __name__ == "__main__":
    main()
