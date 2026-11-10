#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
329_chai_celltype_enrichment.py  —  T2-B: cell-type enrichment of eBMD effector genes
======================================================================================
把本文 31 BMD 工具之外、Chai 2026 (448K eBMD) 的 1,236 个 GWS lead 的
mappedGenes 作为**独立效应基因集**，检验其是否富集于 Chai 单细胞图谱的
34 个骨/骨髓细胞簇，重点 5 类 SOS 机制相关细胞：
  27. Chondrocytes / 28. Osteoblasts / 30. Osteoclasts /
  31. Endothelial cells / 32. Vascular SMCs

数据（全部来自 Chai 2026 Nat Genet 补充表，MOESM4）：
  Table 2a (sheet '2a')  34 细胞簇 marker 基因程序 (CELL LINEAGE / CLUSTER / GENE)
  Table 3a (sheet '3a')  GO 骨骼基因：HGNC(人) <-> MGI(鼠)
  Table 3b (sheet '3b')  MGI 异常骨骼基因：HGNC(人) <-> MGI(鼠)
效应基因：code/328 输出的 Chai2026_leads_full.tsv 的 gene 列 (mappedGenes)。

方法：超几何检验（N=Table2a 全基因背景，K=某簇基因数，n=效应基因中在背景内的数，
      k=交集）；BH-FDR 跨 34 簇校正。人鼠同源：先用 3a/3b 建 HGNC->MGI 映射，
      未映射者以大写直配（同源基因符号多保守）。
诚实：这是**效应基因集层面的细胞富集**（不是 31 工具的效应基因，因为 areal-BMD
      lead 的效应基因未系统注释）；它检验「骨骼 eBMD 遗传信号定位于哪些细胞」，
      为 SOS 的细胞级机制对应物提供独立、公开数据的支持，不改动 docs/08。
输出：code/output/genetics/Chai2026_celltype_enrichment.json / .tsv
"""
import os, csv, json, math, collections
import openpyxl
from scipy.stats import hypergeom

ROOT=os.path.dirname(os.path.abspath(__file__))
GEN=os.path.join(ROOT,"output","genetics")
LEADS=os.path.join(GEN,"Chai2026_leads_full.tsv")
XLSX=os.path.join(GEN,"Chai2026_MOESM4.xlsx")
OUTJ=os.path.join(GEN,"Chai2026_celltype_enrichment.json")
OUTT=os.path.join(GEN,"Chai2026_celltype_enrichment.tsv")
KEY_CLUSTERS={"27. Chondrocytes":"chondrocyte","28. Osteoblasts":"osteoblast",
              "30. Osteoclasts":"osteoclast","31. Endothelial cells":"endothelial",
              "32. Vascular SMCs":"vascular_SMC"}

def load_effector_genes():
    genes=set()
    with open(LEADS) as f:
        for r in csv.DictReader(f,delimiter="\t"):
            for g in (r.get("gene") or "").split(","):
                g=g.strip().upper()
                if g: genes.add(g)
    return genes

def load_homology(wb):
    """HGNC(human, upper) -> set(MGI mouse symbols)."""
    h2m=collections.defaultdict(set)
    for sh in ("3a","3b"):
        ws=wb[sh]
        hdr=None
        for row in ws.iter_rows(values_only=True):
            vals=list(row)
            if hdr is None:
                # find header row containing HGNC Symbol / MGI Symbol
                low=[str(v).lower() if v is not None else "" for v in vals]
                if "hgnc symbol" in low and "mgi symbol" in low:
                    hdr={k:i for i,k in enumerate(low)}
                continue
            try:
                mgi=str(vals[hdr["mgi symbol"]]).strip()
                hgnc=str(vals[hdr["hgnc symbol"]]).strip()
            except (KeyError,IndexError,TypeError):
                continue
            if hgnc and hgnc!="None" and mgi and mgi!="None":
                h2m[hgnc.upper()].add(mgi)
    return h2m

def load_clusters(wb):
    cl=collections.defaultdict(set); lineage={}
    for row in wb["2a"].iter_rows(min_row=2,values_only=True):
        if not row[0]: continue
        clu=str(row[1]).strip() if len(row)>1 and row[1] else None
        gene=str(row[2]).strip() if len(row)>2 and row[2] else None
        if clu and gene:
            cl[clu].add(gene); lineage[clu]=str(row[0]).strip()
    return cl,lineage

def bh(pvals):
    m=len(pvals); order=sorted(range(m),key=lambda i:pvals[i])
    q=[0.0]*m; prev=1.0
    for rank,i in enumerate(reversed(order),1):
        val=pvals[i]*m/(m-rank+1); prev=min(prev,val); q[i]=prev
    return q

def main():
    eff_h=load_effector_genes()
    wb=openpyxl.load_workbook(XLSX,read_only=True,data_only=True)
    h2m=load_homology(wb)
    cl,lineage=load_clusters(wb)
    wb.close()
    # map effector genes (human) to mouse symbols for comparison with mouse scRNA clusters
    eff_m=set()
    for g in eff_h:
        if g in h2m: eff_m|=h2m[g]
        else: eff_m.add(g)          # unconverted -> uppercase direct match (mostly conserved)
    universe=set().union(*cl.values()) if cl else set()
    N=len(universe); tested=eff_m & universe; n=len(tested)
    print(f"[329] effector genes (human)={len(eff_h)}; mapped to mouse set={len(eff_m)}; "
          f"in scRNA universe={n}/{N}",flush=True)
    rows=[]
    for clu,genes in cl.items():
        K=len(genes & universe); k=len(tested & genes)
        if K==0 or N==0: continue
        p=float(hypergeom.sf(k-1,N,K,n)) if k>0 else 1.0
        exp=K*n/N
        rows.append(dict(cluster=clu,lineage=lineage.get(clu,""),K=K,k=k,expected=round(exp,2),
                         fold=round(k/exp,2) if exp>0 else None,p=p))
    qs=bh([r["p"] for r in rows])
    for r,q in zip(rows,qs): r["q"]=float(q)
    rows.sort(key=lambda r:r["p"])
    # write TSV
    with open(OUTT,"w",newline="") as f:
        w=csv.writer(f,delimiter="\t"); w.writerow(["cluster","lineage","K","k","expected","fold","p","q"])
        for r in rows: w.writerow([r["cluster"],r["lineage"],r["K"],r["k"],r["expected"],r["fold"],f"{r['p']:.3g}",f"{r['q']:.3g}"])
    key=[]
    for clu,label in KEY_CLUSTERS.items():
        hit=next((r for r in rows if r["cluster"]==clu),None)
        if hit: key.append(dict(cell=label,**hit))
    sig34=sum(1 for r in rows if r["q"]<0.05)
    out=dict(design="T2-B: Chai 2026 448K-eBMD effector genes vs 34-cell-type scRNA programs",
             effector_source="GWAS Catalog mappedGenes of 1,236 Chai GCST90726625 leads",
             cluster_source="Chai 2026 Nat Genet Suppl Table 2a (mouse endosteal scRNA)",
             homology="Chai Suppl Table 3a/3b HGNC->MGI; unconverted uppercase direct match",
             universe_genes=N,effector_in_universe=n,n_clusters=len(rows),
             n_clusters_q_lt_0p05=sig34,key_skeletal_cells=key,
             note="effector-gene-set enrichment (NOT the 31 areal-BMD instruments); "
                  "independent public-data cell-level anchor for SOS mechanism",
             rows=rows)
    json.dump(out,open(OUTJ,"w",encoding="utf-8"),ensure_ascii=False,indent=2)
    print(f"\n[329] {sig34}/{len(rows)} clusters q<0.05",flush=True)
    print("[329] key skeletal cells:",flush=True)
    for h in key:
        print(f"   {h['cell']:14s} k={h['k']}/{h['K']} fold={h['fold']} p={h['p']:.2g} q={h['q']:.2g}",flush=True)
    print(f"[write] {OUTJ}\n[write] {OUTT}")

if __name__=="__main__":
    main()
