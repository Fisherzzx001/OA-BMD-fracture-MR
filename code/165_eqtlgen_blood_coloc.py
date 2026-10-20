#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
165_eqtlgen_blood_coloc.py — F4-eqtl（血液语境）：eQTLGen whole-blood cis-eQTL × OA/BMD 共定位
=============================================================================================
动机（docs/118 F4 / docs/121 §5）
------------------------------
F4 字面设计是"骨组织 eQTL↔OA vs 软骨 eQTL↔OA 的 PP.H4 比较"，但**骨/软骨组织 eQTL
汇总统计无公开资源**（eQTLGen 主体为 whole blood；GTEx 无骨/软骨组织）→ 字面对比
不可执行。本脚本执行**血液语境 coloc**（诚实降级，结论按此口径表述）：

  对 31 个 BMD_total harmonised 仪器位点（±500 kb，与 code/41 同窗口），按基因分别做
      PP.H4(blood_eQTL, OA_knee)   与   PP.H4(blood_eQTL, BMD_total)
  若某位点的 OA 信号与某基因的血液 cis-eQTL 共享因果变异（PP.H4≥0.8），说明该仪器
  位点在血液中受调控且与 OA 关联共享因果变异 —— 是"骨驱动软骨"链条上**可验证的一环**，
  但不区分骨/软骨组织（组织特异性未解）。

方法学（复用经审计的 code/29 引擎，保证与 S8/P1-3 同口径）
--------------------------------------------------------
· coloc.abf：Giambartolomei 2014；先验 p1=p2=1e-4, p12=5e-6；quant sd=0.15, cc sd=0.2
  （OA_knee s=0.12 取自 code/29 TRAITS）。
· 等位基因对齐：eQTL {A1=AssessedAllele, Z 参照 A1} vs GWAS {A1, A2}；集合相等且
  A1 一致→保留；A1 互换→Z 取负；集合不等→剔除（code/29 的 merge 只按 CHR/BP，
  必须先做等位基因清洗，否则静默错向）。
· v_e = var_beta(NrSamples, eQTLGen MAF, "quant")；v_g = var_beta(GWAS N, GWAS MAF, type)。
· 每基因 ≥20 个等位一致 SNP 才跑；跨全部 (位点×基因×对比) 检验做 BH-FDR。

输入
----
data/gwas_stratified/eQTLGen_blood_loci.tsv.gz（code/166 生成：
  CHR BP SNP A1 A2 MAF N Z GENE GENE_SYMBOL）

输出
----
output/genetics/eqtlgen_blood_coloc_locus_gene.csv   全表（位点×基因×对比）
output/genetics/eqtlgen_blood_coloc_locus_summary.csv 每 locus 最优基因汇总
output/genetics/eqtlgen_blood_coloc_report.txt        可读报告

用法
----
python code/165_eqtlgen_blood_coloc.py
"""
import os
import sys
import importlib.util
import numpy as np
import pandas as pd
from scipy.stats import norm, wilcoxon

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "output", "genetics")
HARM = os.path.join(OUT, "harmonised", "BMD_total__OA_knee.csv")
EQTL = os.path.join(ROOT, "..", "data", "gwas_stratified", "eQTLGen_blood_loci.tsv.gz")
WINDOW = 500_000
MIN_SNP = 20

_spec = importlib.util.spec_from_file_location("c29", os.path.join(ROOT, "29_colocalization.py"))
C29 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(C29)


def bh_fdr(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(q, 1.0)
    return out


def align_and_coloc(e_g, g_g, typ_g, s_g):
    """e_g: eQTL 区域 (CHR BP A1 A2 Z MAF N)；g_g: GWAS 区域。等位基因对齐后 coloc.abf。"""
    m = e_g[["CHR", "BP", "A1", "A2", "Z", "MAF", "N"]].merge(
        g_g[["CHR", "BP", "A1", "A2", "Z", "MAF", "N"]],
        on=["CHR", "BP"], suffixes=("_e", "_g"))
    if m.empty:
        return None
    same = (m["A1_e"] == m["A1_g"]) & (m["A2_e"] == m["A2_g"])
    swap = (m["A1_e"] == m["A2_g"]) & (m["A2_e"] == m["A1_g"])
    m = m[same | swap].copy()
    if len(m) < MIN_SNP:
        return None
    z_e = np.where(m["A1_e"].values == m["A1_g"].values, m["Z_e"].values, -m["Z_e"].values)
    v_e = C29.var_beta(m["N_e"].values, m["MAF_e"].values, "quant", None)
    v_g = C29.var_beta(m["N_g"].values, m["MAF_g"].values, typ_g, s_g)
    res = C29.coloc_abf(z_e, v_e, m["Z_g"].values.astype(float), v_g, "quant", typ_g)
    return dict(n_snp=int(len(m)), n_allele_drop=int((~(same | swap)).sum()), **res)


def selftest():
    """合成数据验证等位基因对齐 + coloc 集成：植入共享因果变异，断言 PP4 高且剔除数正确。"""
    rng = np.random.default_rng(20261017)
    n = 40
    ch, pos0 = 7, 120966790
    bp = pos0 + np.arange(n) * 1000
    maf = rng.uniform(0.1, 0.45, n)
    z_g = rng.normal(0, 1, n); z_g[10] = 6.0      # GWAS 共享因果变异
    z_e = rng.normal(0, 1, n); z_e[10] = 5.5      # eQTL 共享因果变异（同方向）
    a1 = np.where(rng.random(n) < 0.5, "A", "G")
    a2 = np.where(a1 == "A", "G", "A")
    e = pd.DataFrame(dict(CHR=ch, BP=bp, A1=a1.copy(), A2=a2.copy(), Z=z_e.copy(),
                          MAF=maf, N=np.full(n, 30000.0)))
    g = pd.DataFrame(dict(CHR=ch, BP=bp, A1=a1.copy(), A2=a2.copy(), Z=z_g.copy(),
                          MAF=maf, N=np.full(n, 32000.0)))
    # 3 个 SNP：eQTL 等位基因互换且 Z 同步取负（同一物理效应）→ 对齐后应复原
    for i in (30, 31, 32):
        e.loc[i, "A1"], e.loc[i, "A2"] = a2[i], a1[i]
        e.loc[i, "Z"] = -z_e[i]
    # 2 个 SNP：等位基因不一致 → 应剔除
    for i in (35, 36):
        e.loc[i, "A1"], e.loc[i, "A2"] = "C", "T"
    r = align_and_coloc(e, g, "quant", None)
    assert r is not None, "selftest: align_and_coloc returned None"
    assert r["n_snp"] == n - 2, "selftest: n_snp=%d expected %d" % (r["n_snp"], n - 2)
    assert r["n_allele_drop"] == 2, "selftest: allele drop=%d expected 2" % r["n_allele_drop"]
    assert r["PP4"] > 0.8, "selftest: PP4=%.3f expected >0.8" % r["PP4"]
    print("=" * 78)
    print("165 SELFTEST — allele alignment + coloc integration (synthetic shared variant)")
    print("  n_snp=%d (dropped %d discordant)  PP4=%.3f  PP3=%.3f"
          % (r["n_snp"], r["n_allele_drop"], r["PP4"], r["PP3"]))
    print("SELFTEST PASSED — swap-aware alignment restores shared-signal PP4")
    return True


def main():
    if not os.path.exists(EQTL):
        sys.exit("FATAL: %s not found — run code/166_eqtlgen_loci_extract.py first" % EQTL)
    os.makedirs(OUT, exist_ok=True)
    eqtl = pd.read_csv(EQTL, sep="\t",
                       dtype={"CHR": "int32", "BP": "int64", "SNP": str,
                              "A1": str, "A2": str, "GENE": str, "GENE_SYMBOL": str})
    inst = pd.read_csv(HARM, usecols=["SNP", "CHR", "BP"])
    oa, t_oa, s_oa = C29.get_trait("OA_knee")
    bmd, t_bmd, s_bmd = C29.get_trait("total")
    print("eQTL 子集: %d 行 / %d SNP / %d 基因；仪器 %d 个；OA_knee %d SNP；BMD_total %d SNP"
          % (len(eqtl), eqtl.SNP.nunique(), eqtl.GENE.nunique(), len(inst), len(oa), len(bmd)))

    rows = []
    for _, rw in inst.iterrows():
        ch, pos = int(rw["CHR"]), int(rw["BP"])
        e_win = C29.region(eqtl, ch, pos, WINDOW)
        o_win = C29.region(oa, ch, pos, WINDOW)
        b_win = C29.region(bmd, ch, pos, WINDOW)
        genes = sorted(e_win["GENE"].unique())
        if not genes:
            print("  %-14s chr%d:%d — 窗口内无 eQTL 行" % (rw["SNP"], ch, pos))
            continue
        for g in genes:
            e_g = e_win[e_win["GENE"] == g]
            sym = e_g["GENE_SYMBOL"].iloc[0]
            for label, g_win, typ_g, s_g in [("OA_knee", o_win, t_oa, s_oa),
                                             ("BMD_total", b_win, t_bmd, s_bmd)]:
                r = align_and_coloc(e_g, g_win, typ_g, s_g)
                if r is None:
                    rows.append(dict(index_snp=rw["SNP"], chrom=ch, pos=pos, gene=g,
                                     gene_symbol=sym, contrast=label, n_snp=len(e_g),
                                     note="skip(<20 aligned SNPs or no overlap)"))
                    continue
                rows.append(dict(index_snp=rw["SNP"], chrom=ch, pos=pos, gene=g,
                                 gene_symbol=sym, contrast=label, **r))
        print("  %-14s chr%d:%d — %d genes" % (rw["SNP"], ch, pos, len(genes)))

    d = pd.DataFrame(rows)
    ok = d[d["PP4"].notna()].copy() if "PP4" in d.columns else d.iloc[0:0].copy()
    if len(ok):
        # 以 1-PP4 作"无共定位"p 值代理做 BH 校正（小 FDR = 共定位证据强）
        ok["FDR"] = bh_fdr(1.0 - ok["PP4"].values)
        d.loc[ok.index, "FDR"] = ok["FDR"]
    d.to_csv(os.path.join(OUT, "eqtlgen_blood_coloc_locus_gene.csv"), index=False)

    # ---- 每 locus 汇总：按 PP4(OA_knee) 取最优基因 ----
    oa_rows = ok[ok["contrast"] == "OA_knee"] if len(ok) else ok
    summ = []
    for (snp, ch, pos), grp in oa_rows.groupby(["index_snp", "chrom", "pos"], sort=False):
        if grp.empty:
            continue
        best = grp.sort_values(["PP4", "n_snp"], ascending=False).iloc[0]
        bmd_same = ok[(ok["contrast"] == "BMD_total") & (ok["index_snp"] == snp)
                      & (ok["gene"] == best["gene"])]
        summ.append(dict(index_snp=snp, chrom=ch, pos=pos,
                         best_gene=best["gene"], best_gene_symbol=best["gene_symbol"],
                         n_snp=int(best["n_snp"]),
                         PP4_blood_eqtl_OA=float(best["PP4"]),
                         PP4_blood_eqtl_BMD=float(bmd_same["PP4"].iloc[0]) if len(bmd_same) else np.nan,
                         n_genes_tested=int(grp["gene"].nunique())))
    s = pd.DataFrame(summ).sort_values("PP4_blood_eqtl_OA", ascending=False)
    s.to_csv(os.path.join(OUT, "eqtlgen_blood_coloc_locus_summary.csv"), index=False)

    # ---- 报告 ----
    lines = []
    lines.append("=" * 78)
    lines.append("F4-eqtl 血液语境 coloc（eQTLGen whole blood × OA_knee / BMD_total）")
    lines.append("=" * 78)
    lines.append("位点数                 : %d" % s["index_snp"].nunique())
    lines.append("(位点×基因×对比) 检验数 : %d" % len(ok))
    lines.append("median PP.H4 (eQTL×OA) : %.3f" % s["PP4_blood_eqtl_OA"].median())
    lines.append("PP.H4>=0.8 (eQTL×OA)   : %d/%d" % (int((s["PP4_blood_eqtl_OA"] >= 0.8).sum()), len(s)))
    lines.append("median PP.H4 (eQTL×BMD): %.3f" % s["PP4_blood_eqtl_BMD"].median())
    lines.append("PP.H4>=0.8 (eQTL×BMD)  : %d/%d" % (int((s["PP4_blood_eqtl_BMD"] >= 0.8).sum()), len(s)))
    if len(s) >= 5:
        w = wilcoxon(s["PP4_blood_eqtl_OA"].values - 0.5, alternative="greater")
        lines.append("Wilcoxon PP4>0.5       : W=%d p=%.3g" % (w.statistic, w.pvalue))
    lines.append("")
    lines.append("诚实边界：eQTLGen 为 whole blood；骨/软骨组织 eQTL 无公开汇总统计，")
    lines.append("字面 F4（骨 eQTL vs 软骨 eQTL）不可执行；本结果只说明仪器位点在血液中")
    lines.append("是否受调控并与 OA 共享因果变异，不区分骨/软骨组织。")
    lines.append("")
    lines.append("Top 10 loci by PP.H4(blood eQTL, OA_knee):")
    lines.append(s.head(10)[["index_snp", "best_gene_symbol", "n_snp",
                             "PP4_blood_eqtl_OA", "PP4_blood_eqtl_BMD"]].to_string(index=False))
    rep = "\n".join(lines)
    print("\n" + rep)
    with open(os.path.join(OUT, "eqtlgen_blood_coloc_report.txt"), "w") as f:
        f.write(rep + "\n")
    print("\nwritten: eqtlgen_blood_coloc_locus_gene.csv / _locus_summary.csv / _report.txt")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "selftest":
        selftest()
    else:
        main()
