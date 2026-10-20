#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
166_eqtlgen_loci_extract.py — F4 数据制备：eQTLGen 全量 cis-eQTL → 31 个 BMD 工具位点窗口子集
=============================================================================================
目的（docs/118 F4 / docs/121 §5）
------------------------------
eQTLGen 全量 cis 文件（2019-12-11，3.88 GB，~38M 行）只保留落在 31 个 BMD_total
harmonised 仪器 ±500 kb 窗口内的行（与 code/41 WINDOW=500_000 同口径），并连接
eQTLGen 自有等位基因频率文件取 MAF，输出 coloc 可直接消费的位点子集。

输入
----
1) data/gwas_stratified/eQTLGen_cis_full.tsv.gz
   列（README_cis 核实）：Pvalue SNP SNPChr SNPPos AssessedAllele OtherAllele Zscore
   Gene GeneSymbol GeneChr GenePos NrCohorts NrSamples FDR BonferroniP
   —— "all cis-eQTL results"：每个 SNP-gene 对（<1Mb，≥2 cohorts 检测过）一行。
2) data/gwas_stratified/eQTLGen_SNP_AF.txt.gz
   列（README_allele_frequency 核实）：SNP hg19_chr hg19_pos AlleleA AlleleB
   allA_total allAB_total allB_total AlleleB_all（AlleleB = assessed allele 的频率）
3) code/output/genetics/harmonised/BMD_total__OA_knee.csv（31 个仪器：SNP CHR BP）

输出
----
data/gwas_stratified/eQTLGen_blood_loci.tsv.gz
   列：CHR BP SNP A1 A2 MAF N Z GENE GENE_SYMBOL
   A1 = AssessedAllele（Z 的方向参照 A1）；MAF = min(AlleleB_all, 1-AlleleB_all)；
   N = NrSamples（该 SNP-gene 对在 eQTLGen 的样本数）。
   等位基因与 GWAS 的翻转/剔除在 coloc 封装（code/165）中按分析逐次进行。

诚实边界
--------
· eQTLGen 主体为 whole blood；骨/软骨组织 eQTL 无公开汇总统计，字面"骨 vs 软骨"
  F4 不可执行 —— 本子集只支撑血液语境 coloc（code/165）。
· 不做任何筛选/显著性过滤：保留窗口内全部检测过的 SNP-gene 对（coloc 需要完整
  区域统计量，按显著性筛选会向上偏倚 PP.H4）。
"""
import os
import sys
import gzip
import time
import pandas as pd

ROOT = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(ROOT, "..", "data", "gwas_stratified")
HARM = os.path.join(ROOT, "output", "genetics", "harmonised", "BMD_total__OA_knee.csv")
CIS = os.path.join(DATA, "eQTLGen_cis_full.tsv.gz")
AF = os.path.join(DATA, "eQTLGen_SNP_AF.txt.gz")
OUT = os.path.join(DATA, "eQTLGen_blood_loci.tsv.gz")
WINDOW = 500_000


def load_windows():
    inst = pd.read_csv(HARM, usecols=["SNP", "CHR", "BP"])
    win = {}
    for _, r in inst.iterrows():
        win.setdefault(int(r["CHR"]), []).append(int(r["BP"]))
    for c in win:
        win[c] = sorted(set(win[c]))
    return inst, win


def pass_cis(win):
    """流式扫描全量 cis 文件，保留窗口内行。"""
    rows = []
    n_seen = 0
    t0 = time.time()
    with gzip.open(CIS, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for col in ("Pvalue", "SNP", "SNPChr", "SNPPos", "AssessedAllele",
                    "OtherAllele", "Zscore", "Gene", "GeneSymbol", "NrSamples"):
            if col not in idx:
                sys.exit("FATAL: cis file missing column %s (header=%s)" % (col, header))
        i_ch, i_pos = idx["SNPChr"], idx["SNPPos"]
        i_snp, i_a1, i_a2 = idx["SNP"], idx["AssessedAllele"], idx["OtherAllele"]
        i_z, i_n = idx["Zscore"], idx["NrSamples"]
        i_gene, i_sym, i_p = idx["Gene"], idx["GeneSymbol"], idx["Pvalue"]
        for line in f:
            n_seen += 1
            if n_seen % 5_000_000 == 0:
                print("  cis scan: %d rows (%.0fs, in-window %d)" % (n_seen, time.time() - t0, len(rows)))
            p = line.rstrip("\n").split("\t")
            if len(p) < len(header):
                continue
            try:
                ch = int(p[i_ch]); pos = int(p[i_pos])
            except ValueError:
                continue
            bpl = win.get(ch)
            if not bpl:
                continue
            hit = False
            for b in bpl:
                if -WINDOW <= pos - b <= WINDOW:
                    hit = True
                    break
            if not hit:
                continue
            try:
                z = float(p[i_z]); n_s = float(p[i_n]); pv = float(p[i_p])
            except ValueError:
                continue
            if n_s <= 0 or z != z:
                continue
            rows.append(dict(CHR=ch, BP=pos, SNP=p[i_snp],
                             A1=p[i_a1].upper(), A2=p[i_a2].upper(),
                             N=n_s, Z=z, GENE=p[i_gene], GENE_SYMBOL=p[i_sym], P=pv))
    print("  cis scan done: %d rows total, %d in-window (%.0fs)"
          % (n_seen, len(rows), time.time() - t0))
    return pd.DataFrame(rows)


def pass_af(need_snps):
    """流式扫描 AF 文件，按 SNP rsID 连接窗口内 SNP 的频率。"""
    af = {}
    n_seen = 0
    with gzip.open(AF, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for col in ("SNP", "AlleleB_all"):
            if col not in idx:
                sys.exit("FATAL: AF file missing column %s (header=%s)" % (col, header))
        i_snp, i_afb = idx["SNP"], idx["AlleleB_all"]
        for line in f:
            n_seen += 1
            p = line.rstrip("\n").split("\t")
            if len(p) < len(header):
                continue
            s = p[i_snp]
            if s in need_snps and s not in af:
                try:
                    af[s] = float(p[i_afb])
                except ValueError:
                    pass
    print("  AF scan done: %d rows, %d window SNPs matched" % (n_seen, len(af)))
    return af


def selftest():
    """合成数据验证窗口过滤 + AF 连接 + MAF 计算（trust-but-verify，不触真实大文件）。"""
    import tempfile
    inst = pd.read_csv(HARM)
    r0 = inst.iloc[0]; ch, bp = int(r0["CHR"]), int(r0["BP"])
    tmp = tempfile.mkdtemp(prefix="eqtlgen166_")
    cis = os.path.join(tmp, "cis.tsv.gz")
    with gzip.open(cis, "wt") as f:
        f.write("Pvalue\tSNP\tSNPChr\tSNPPos\tAssessedAllele\tOtherAllele\tZscore\t"
                "Gene\tGeneSymbol\tGeneChr\tGenePos\tNrCohorts\tNrSamples\tFDR\tBonferroniP\n")
        f.write("1e-5\trsTEST1\t%d\t%d\tA\tG\t3.5\tENSGTEST1\tTESTG1\t%d\t%d\t30\t30000\t0.01\t0.02\n"
                % (ch, bp, ch, bp))                                   # in-window, 有 AF
        f.write("1e-4\trsTEST2\t%d\t%d\tC\tT\t-2.5\tENSGTEST1\tTESTG1\t%d\t%d\t30\t29900\t0.02\t0.03\n"
                % (ch, bp + 100, ch, bp))                             # in-window, AF=0.9 → MAF 0.1
        f.write("1e-3\trsTESTOUT\t%d\t%d\tG\tA\t1.5\tENSGTEST2\tTESTG2\t%d\t%d\t30\t29800\t0.03\t0.04\n"
                % (ch, bp + 10 ** 7, ch, bp + 10 ** 7))               # 窗外 → 剔除
        f.write("1e-3\trsTESTBAD\t%d\t%d\tT\tC\tNA\tENSGTEST2\tTESTG2\t%d\t%d\t30\t29800\t0.03\t0.04\n"
                % (ch, bp + 50, ch, bp))                              # Z=NA → 剔除
    aff = os.path.join(tmp, "af.tsv.gz")
    with gzip.open(aff, "wt") as f:
        f.write("SNP\thg19_chr\thg19_pos\tAlleleA\tAlleleB\tallA_total\tallAB_total\t"
                "allB_total\tAlleleB_all\n")
        f.write("rsTEST1\t%d\t%d\tG\tA\t100\t200\t300\t0.25\n" % (ch, bp))
        f.write("rsTEST2\t%d\t%d\tT\tC\t100\t200\t300\t0.9\n" % (ch, bp + 100))
        # rsTESTOUT / rsTESTBAD 不在 AF → 连接后剔除
    global CIS, AF, OUT
    cis0, af0, out0 = CIS, AF, OUT
    tmp_out = os.path.join(tmp, "out.tsv.gz")
    CIS, AF, OUT = cis, aff, tmp_out
    try:
        main()
    finally:
        CIS, AF, OUT = cis0, af0, out0
    d = pd.read_csv(tmp_out, sep="\t")
    assert len(d) == 2, "selftest: expected 2 rows, got %d" % len(d)
    assert set(d["SNP"]) == {"rsTEST1", "rsTEST2"}, "selftest: SNP set wrong: %s" % set(d["SNP"])
    m1 = float(d.loc[d.SNP == "rsTEST1", "MAF"].iloc[0])
    m2 = float(d.loc[d.SNP == "rsTEST2", "MAF"].iloc[0])
    assert abs(m1 - 0.25) < 1e-9 and abs(m2 - 0.1) < 1e-9, "selftest: MAF wrong %s %s" % (m1, m2)
    assert abs(float(d.loc[d.SNP == "rsTEST1", "Z"].iloc[0]) - 3.5) < 1e-6, "selftest: Z wrong"
    print("SELFTEST PASSED — window filter / AF join / MAF=min(AF,1-AF) all correct")
    return True


def main():
    for p in (CIS, AF, HARM):
        if not os.path.exists(p):
            sys.exit("FATAL: missing input %s (download first; see docs/121 §5)" % p)
    inst, win = load_windows()
    print("instruments: %d on %d chromosomes; window ±%d bp" % (len(inst), len(win), WINDOW))

    d = pass_cis(win)
    if d.empty:
        sys.exit("FATAL: no in-window eQTL rows — check instrument coordinates/window")
    af = pass_af(set(d["SNP"]))
    d["AF_B"] = d["SNP"].map(af)
    n0 = len(d)
    d = d[d["AF_B"].notna()].copy()
    d["MAF"] = d["AF_B"].clip(0, 1)
    d["MAF"] = d["MAF"].where(d["MAF"] <= 0.5, 1.0 - d["MAF"])
    d = d[(d["MAF"] > 0) & (d["MAF"] < 0.5)]
    d = d.sort_values(["CHR", "BP", "GENE"]).reset_index(drop=True)
    d[["CHR", "BP", "SNP", "A1", "A2", "MAF", "N", "Z", "GENE", "GENE_SYMBOL"]] \
        .to_csv(OUT, sep="\t", index=False, compression="gzip", float_format="%.6g")

    print("\n=== extract QC ===")
    print("in-window SNP-gene rows (pre-AF) : %d" % n0)
    print("after AF join / MAF filter       : %d" % len(d))
    print("unique SNPs                      : %d" % d["SNP"].nunique())
    print("unique genes                     : %d" % d["GENE"].nunique())
    print("chromosomes covered              : %s" % sorted(d["CHR"].unique()))
    print("per-chromosome rows:")
    print(d.groupby("CHR").size().to_string())
    print("\nwritten: %s" % OUT)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "selftest":
        selftest()
    else:
        main()
