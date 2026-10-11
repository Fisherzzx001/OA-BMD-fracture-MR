#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
41_eqtl_coloc.py — P1-3 遗传验证：组织富集 + eQTL coloc（骨驱动机制的遗传旁证）
===================================================================================
蓝图来源：docs/99b_机制真空_理论深化框架.md §5（P1-3，最高价值新增分析）。

目标：遗传学证明 BMD->OA 的工具变量"在骨（成骨细胞/软骨下骨）而非软骨中富集表达"，
从而把"骨驱动软骨"从文字假说转为有遗传证据支持的假说。

为什么能落地（trust-but-verify）
------------------------------
· 真实 BMD / OA / fracture 全基因组 munge 汇总统计已齐备，位于
  C:/Users/Administrator/WorkBuddy/2026-09-08-20-17-24/OA_FLAGSHIP/07_reference/ldsc/
  （code/29_colocalization.py 已用同一批数据完成 S8 共定位，方法学经过审计）。
· 31 个 BMD_total 工具位点已 harmonise 在
  code/output/genetics/harmonised/BMD_total__OA_knee.csv（与 BMD_total__fracture.csv 同位点）。
· 本脚本**复用 code/29 的 coloc_abf 引擎**（相同先验 p1=p2=1e-4, p12=5e-6；相同 BMD 真值 N 修正），
  保证 S8 与 P1-3 的共定位口径完全一致。

三种运行模式
------------
1) local   — 纯本地数据，立即可跑（无需任何外部下载）。
   对 31 个 BMD_total 工具位点，分别跑
      PP4(BMD_total, OA_knee)   与   PP4(BMD_total, fracture)
   的 coloc.abf。fracture 是"纯骨"结局；若 BMD 工具位点在骨折（骨）上的共定位
   系统性高于在 OA（关节）上的共定位，则遗传学支持"这些工具变量是骨侧（bone-side）
   位点"，与"骨驱动软骨"机制自洽。这是本回合交付的**真实新补充分析**。
   输出：output/genetics/genetic_validation_coloc_contrast.csv

2) eqtl    — 骨/软骨 eQTL 就位时运行（数据来自 eQTLGen / GTEx，需下载，见 docs/99c）。
   对每个位点跑 coloc.abf，比较
      PP4(bone_eQTL, OA_knee)  vs  PP4(chondro_eQTL, OA_knee)
   若骨组织 eQTL 的 PP4 系统性更高 -> 因果中介发生在骨而非软骨内源。
   输出：output/genetics/eqtl_coloc_bone.csv / eqtl_coloc_cartilage.csv

3) selftest— 合成数据验证"组织判别力"：在骨 eQTL 中植入共享因果变异、软骨 eQTL 中
   不植入，断言 PP4(bone) > PP4(chondro)。证明引擎能区分组织语境，而非恒输出同一值。

注意（诚实边界）
----------------
local 模式给出的是"骨侧信号"的遗传学旁证，不等于组织内 eQTL 中介的直接证据；
eQTL/LDSC-SEG 的精细定位需多 GB 级公开权重，在 docs/99c 中以"代码已验证、集群就绪"
方式给出，本机不虚构其结果。
"""

import os
import sys
import gzip
import math
import importlib.util
import numpy as np
import pandas as pd
from scipy import stats as spstats

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "output", "genetics")
HARM = os.path.join(OUT, "harmonised", "BMD_total__OA_knee.csv")

# ---------------------------------------------------------------- 载入 code/29 引擎
_spec = importlib.util.spec_from_file_location("c29", os.path.join(ROOT, "29_colocalization.py"))
C29 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(C29)

WINDOW = 500_000


# ---------------------------------------------------------------- eQTL 载入（eqtl 模式）
def load_eqtl(path, n_study, sd_prior=0.15):
    """读组织 eQTL 汇总：CHR BP SNP A1 A2 MAF N Z（空格/制表分隔，可 .gz）。"""
    df = pd.read_csv(path, sep=r"\s+", compression="infer" if path.endswith(".gz") else None,
                     dtype={"CHR": "int32", "BP": "int64", "SNP": str, "A1": str, "A2": str})
    df = df.dropna(subset=["CHR", "BP", "MAF", "Z", "N"]).copy()
    for c in ("MAF", "Z", "N"):
        df[c] = df[c].astype(float)
    df["N"] = float(n_study)
    return df


def coloc_pair(trait1_zv, trait2_zv):
    """trait?_zv = (z-array, v-array, type, s)；返回 dict(PP0..PP4)。"""
    z1, v1, t1, s1 = trait1_zv
    z2, v2, t2, s2 = trait2_zv
    return C29.coloc_abf(z1, v1, z2, v2, t1, t2)


# ---------------------------------------------------------------- 模式 1：local cross-trait contrast
def run_local():
    os.makedirs(OUT, exist_ok=True)
    inst = pd.read_csv(HARM)
    rows = []
    for _, rw in inst.iterrows():
        chr_, pos = int(rw["CHR"]), int(rw["BP"])
        r_oa = C29.coloc_region("total", "OA_knee", chr_, pos, window=WINDOW)
        r_fx = C29.coloc_region("total", "fracture", chr_, pos, window=WINDOW)
        if r_oa is None or r_fx is None:
            continue
        rows.append(dict(
            index_snp=rw["SNP"], chrom=chr_, pos=pos,
            PP4_BMD_OA=r_oa["PP4"], PP4_BMD_fracture=r_fx["PP4"],
            PP3_BMD_OA=r_oa["PP3"], PP3_BMD_fracture=r_fx["PP3"],
            n_snps_oa=r_oa["n_snps"], n_snps_fx=r_fx["n_snps"],
            top_snp_oa=r_oa["top_snp"], z_top_bmd_oa=r_oa["z1_top"],
            z_top_oa=r_oa["z2_top"], z_top_bmd_fx=r_fx["z1_top"],
            z_top_fx=r_fx["z2_top"],
            dir_oa=r_oa.get("dir_concordant"), dir_fx=r_fx.get("dir_concordant"),
        ))
    d = pd.DataFrame(rows)
    d["delta_PP4"] = d["PP4_BMD_fracture"] - d["PP4_BMD_OA"]
    d = d.sort_values("delta_PP4", ascending=False).reset_index(drop=True)
    d.to_csv(os.path.join(OUT, "genetic_validation_coloc_contrast.csv"), index=False)

    # ---- 聚合检验：BMD 工具位点在"纯骨"结局(fracture)上的共定位是否系统性高于 OA ----
    pp_oa = d["PP4_BMD_OA"].values
    pp_fx = d["PP4_BMD_fracture"].values
    # 单样本（配对）Wilcoxon：H1: median(pp_fx - pp_oa) > 0
    w = spstats.wilcoxon(pp_fx, pp_oa, alternative="greater")
    n = len(d)
    frac_fx_higher = float((pp_fx > pp_oa).mean())
    median_oa = float(np.median(pp_oa))
    median_fx = float(np.median(pp_fx))
    n_oa_h4 = int((pp_oa >= 0.8).sum())
    n_fx_h4 = int((pp_fx >= 0.8).sum())

    print("=" * 78)
    print("P1-3 LOCAL MODE — cross-trait coloc contrast at 31 BMD_total instruments")
    print("=" * 78)
    print("loci compared          : %d" % n)
    print("median PP.H4  BMD->OA   : %.3f" % median_oa)
    print("median PP.H4  BMD->fx   : %.3f" % median_fx)
    print("fraction fx > oa        : %.2f" % frac_fx_higher)
    print("PP.H4>=0.8  BMD->OA     : %d/%d" % (n_oa_h4, n))
    print("PP.H4>=0.8  BMD->fx     : %d/%d" % (n_fx_h4, n))
    print("paired Wilcoxon (fx>oa) : W=%d p=%.3g" % (w.statistic, w.pvalue))
    print("verdict                 : " + (
        "BMD instruments colocalize MORE with fracture (bone) than OA (joint) -> bone-side signal"
        if (w.pvalue < 0.05 and frac_fx_higher > 0.5) else
        "no systematic bone-side coloc signal at nominal threshold -> see docs/99c for interpretation"))
    print("\nTop 5 loci by delta_PP4 (fracture - OA):")
    print(d.head(5)[["index_snp", "PP4_BMD_OA", "PP4_BMD_fracture", "delta_PP4"]].to_string(index=False))
    print("\nWritten: output/genetics/genetic_validation_coloc_contrast.csv")
    return dict(n=n, median_oa=median_oa, median_fx=median_fx,
                frac_fx_higher=frac_fx_higher, w_p=w.pvalue,
                n_oa_h4=n_oa_h4, n_fx_h4=n_fx_h4)


# ---------------------------------------------------------------- 模式 2：eQTL coloc（需数据）
def run_eqtl(bone_path, chondro_path, bone_n, chondro_n):
    os.makedirs(OUT, exist_ok=True)
    inst = pd.read_csv(HARM)
    bone = load_eqtl(bone_path, bone_n)
    chondro = load_eqtl(chondro_path, chondro_n)
    bmd, t_bmd, s_bmd = C29.get_trait("total")
    oa, t_oa, s_oa = C29.get_trait("OA_knee")

    rows_b, rows_c = [], []
    for _, rw in inst.iterrows():
        chr_, pos = int(rw["CHR"]), int(rw["BP"])
        rb = C29.region(bone, chr_, pos, WINDOW)
        rc = C29.region(chondro, chr_, pos, WINDOW)
        r_oa = C29.region(oa, chr_, pos, WINDOW)
        r_bmd = C29.region(bmd, chr_, pos, WINDOW)
        # 交集对齐
        for label, re_A, re_B, tA, sA in [
            ("bone", rb, r_oa, "quant", None),
            ("chondro", rc, r_oa, "quant", None),
        ]:
            pass
        # bone-eQTL vs OA
        if len(rb) and len(r_oa):
            m = rb[["CHR", "BP"]].merge(r_oa[["CHR", "BP"]], on=["CHR", "BP"])
            a = rb.merge(m, on=["CHR", "BP"]).sort_values("BP")
            b = r_oa.merge(m, on=["CHR", "BP"]).sort_values("BP")
            if len(a) >= 20:
                v1 = C29.var_beta(a["N"].values, a["MAF"].values.astype(float), "quant", None)
                v2 = C29.var_beta(b["N"].values, b["MAF"].values.astype(float), t_oa, s_oa)
                res = coloc_pair((a["Z"].values.astype(float), v1, "quant", None),
                                 (b["Z"].values.astype(float), v2, t_oa, s_oa))
                rows_b.append(dict(index_snp=rw["SNP"], chrom=chr_, pos=pos, **res))
        # chondro-eQTL vs OA
        if len(rc) and len(r_oa):
            m = rc[["CHR", "BP"]].merge(r_oa[["CHR", "BP"]], on=["CHR", "BP"])
            a = rc.merge(m, on=["CHR", "BP"]).sort_values("BP")
            b = r_oa.merge(m, on=["CHR", "BP"]).sort_values("BP")
            if len(a) >= 20:
                v1 = C29.var_beta(a["N"].values, a["MAF"].values.astype(float), "quant", None)
                v2 = C29.var_beta(b["N"].values, b["MAF"].values.astype(float), t_oa, s_oa)
                res = coloc_pair((a["Z"].values.astype(float), v1, "quant", None),
                                 (b["Z"].values.astype(float), v2, t_oa, s_oa))
                rows_c.append(dict(index_snp=rw["SNP"], chrom=chr_, pos=pos, **res))
        # BMD vs bone-eQTL（骨内一致性）
        if len(rb) and len(r_bmd):
            m = rb[["CHR", "BP"]].merge(r_bmd[["CHR", "BP"]], on=["CHR", "BP"])
            a = rb.merge(m, on=["CHR", "BP"]).sort_values("BP")
            b = r_bmd.merge(m, on=["CHR", "BP"]).sort_values("BP")
            if len(a) >= 20:
                v1 = C29.var_beta(a["N"].values, a["MAF"].values.astype(float), "quant", None)
                v2 = C29.var_beta(b["N"].values, b["MAF"].values.astype(float), t_bmd, s_bmd)
                res = coloc_pair((a["Z"].values.astype(float), v1, "quant", None),
                                 (b["Z"].values.astype(float), v2, t_bmd, s_bmd))
                rows_b[-1].update(PP4_BMD_bone_eqtl=res["PP4"]) if rows_b else None
    pd.DataFrame(rows_b).to_csv(os.path.join(OUT, "eqtl_coloc_bone.csv"), index=False)
    pd.DataFrame(rows_c).to_csv(os.path.join(OUT, "eqtl_coloc_cartilage.csv"), index=False)
    print("Written eqtl_coloc_bone.csv / eqtl_coloc_cartilage.csv")


# ---------------------------------------------------------------- 模式 3：selftest
def run_selftest():
    """合成区域：在骨 eQTL 中植入与 BMD 共享的因果变异，软骨 eQTL 中不植入。
    断言 PP4(bone) > PP4(chondro) 且 PP4(bone) 显著 -> 证明引擎能区分组织语境。"""
    rng = np.random.default_rng(20261010)
    n = 400
    f = rng.uniform(0.05, 0.5, n)
    N = 30000.0
    z_bmd = rng.normal(0, 1, n)                       # BMD 背景
    idx = 200
    z_bmd[idx] = 5.5                                  # BMD 共享因果变异
    # 骨 eQTL：与 BMD 共享该变异
    z_bone = rng.normal(0, 1, n); z_bone[idx] = 5.3
    # 软骨 eQTL：独立噪声，无共享
    z_chondro = rng.normal(0, 1, n)
    v = 1.0 / (2 * N * f * (1 - f))
    pp_bone = C29.coloc_abf(z_bmd, v, z_bone, v, "quant", "quant")["PP4"]
    pp_chondro = C29.coloc_abf(z_bmd, v, z_chondro, v, "quant", "quant")["PP4"]
    ok = (pp_bone > pp_chondro) and (pp_bone > 0.5)
    print("=" * 78)
    print("P1-3 SELFTEST — tissue discrimination (synthetic planted shared variant in bone)")
    print("=" * 78)
    print("PP4(BMD, bone_eQTL)    = %.3f" % pp_bone)
    print("PP4(BMD, chondro_eQTL) = %.3f" % pp_chondro)
    print("assert PP4(bone)>PP4(chondro) & PP4(bone)>0.5 : %s" % ok)
    assert ok, "selftest failed: engine cannot discriminate tissue context"
    print("SELFTEST PASSED")
    return ok


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "local"
    if mode == "local":
        run_local()
    elif mode == "selftest":
        run_selftest()
    elif mode == "eqtl":
        # 用法: python 41_eqtl_coloc.py eqtl <bone.gz> <chondro.gz> <bone_N> <chondro_N>
        run_eqtl(sys.argv[2], sys.argv[3], float(sys.argv[4]), float(sys.argv[5]))
    else:
        print("usage: python 41_eqtl_coloc.py [local|selftest|eqtl ...]")
