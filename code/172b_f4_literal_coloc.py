#!/usr/bin/env python3
"""F4-literal colocalisation driver (reuses code/29 coloc_abf engine).

For a given tissue eQTL file (cartilage now; bone-proxy later), at each of the
31 BMD_total instrument loci, colocalise the GWAS (OA_knee / BMD_total / fracture)
against the tissue eQTL with proper allele alignment, and report PP.H4 per locus.

Allele alignment: eQTL Z is signed for its ALT (effect) allele; GWAS Z is signed
for its A1. We align the eQTL Z to the GWAS A1 orientation before coloc_abf
(coloc_abf assumes both Z vectors share the same allele at each SNP).

Usage:
  python code/172b_f4_literal_coloc.py <eqtl_tsv> <tissue_label> [bone_eqtl_tsv] [bone_label]
Outputs: G:/OA_B_PROJECT_cache/t2/f4_literal_<tissue>_<gwas>.csv  and a summary .txt
"""
import os, sys, json
import numpy as np
import pandas as pd
from scipy import stats as spstats
import importlib.util

ROOT = os.path.dirname(os.path.abspath(__file__))
C29spec = importlib.util.spec_from_file_location("c29", os.path.join(ROOT, "29_colocalization.py"))
C29 = importlib.util.module_from_spec(C29spec); C29spec.loader.exec_module(C29)

HARM = os.path.join(ROOT, "output", "genetics", "harmonised", "BMD_total__OA_knee.csv")
CACHE = "G:/OA_B_PROJECT_cache/t2"
WINDOW = 500_000
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
GWAS_NAMES = [("OA_knee", "OA_knee"), ("total", "BMD_total"), ("fracture", "fracture")]


def align_z(a1g, a2g, a1e, a2e, ze):
    """Return eQTL Z aligned to GWAS A1 orientation, or None if ambiguous."""
    a1g, a2g, a1e, a2e = str(a1g), str(a2g), str(a1e), str(a2e)
    if {a1g, a2g} == {a1e, a2e}:
        return ze if a1g == a1e else -ze
    ca1, ca2 = COMP.get(a1e, ""), COMP.get(a2e, "")
    if {a1g, a2g} == {ca1, ca2}:
        return ze if a1g == ca2 else -ze
    return None


def load_eqtl(path):
    df = pd.read_csv(path, sep="\t")
    df["CHR"] = df["CHR"].astype(int); df["BP"] = df["BP"].astype(int)
    return df


def run_tissue(eqtl_df, tissue_label, out_prefix):
    inst = pd.read_csv(HARM)
    all_rows = {g: [] for g, _ in GWAS_NAMES}
    for _, rw in inst.iterrows():
        chrom, pos, snp = int(rw["CHR"]), int(rw["BP"]), rw["SNP"]
        # eQTL region for this locus
        m = (eqtl_df["CHR"].values == chrom) & (np.abs(eqtl_df["BP"].values - pos) <= WINDOW)
        re = eqtl_df.loc[m]
        if len(re) < 20:
            continue
        for gname, glabel in GWAS_NAMES:
            dg, tg, sg = C29.get_trait(gname)
            rg = C29.region(dg, chrom, pos, WINDOW)
            if len(rg) < 20:
                continue
            mg = rg[["CHR", "BP"]].merge(re[["CHR", "BP"]], on=["CHR", "BP"])
            if len(mg) < 20:
                continue
            a = rg.merge(mg, on=["CHR", "BP"]).sort_values("BP").reset_index(drop=True)
            b = re.merge(mg, on=["CHR", "BP"]).sort_values("BP").reset_index(drop=True)
            # align eQTL Z (b) to GWAS A1 (a)
            zg = a["Z"].values.astype(float)
            ze_raw = b["Z"].values.astype(float)
            ze_al = np.array([
                align_z(a["A1"].iloc[i], a["A2"].iloc[i], b["A1"].iloc[i], b["A2"].iloc[i], ze_raw[i])
                for i in range(len(a))
            ], dtype=float)
            keep = ~np.isnan(ze_al)
            if keep.sum() < 20:
                continue
            zg = zg[keep]; ze = ze_al[keep]
            af = a["MAF"].values.astype(float)[keep]
            bf = b["MAF"].values.astype(float)[keep]
            ng = a["N"].values.astype(float)[keep]
            ne = b["N"].values.astype(float)[keep]
            vg = C29.var_beta(ng, af, tg, sg)
            ve = C29.var_beta(ne, bf, "quant", None)
            res = C29.coloc_abf(zg, vg, ze, ve, tg, "quant")
            # direction concordance of top SNP
            k = int(np.argmax(zg ** 2))
            concordant = (np.sign(zg[k]) == np.sign(ze[k]))
            all_rows[gname].append(dict(
                index_snp=snp, chrom=chrom, pos=pos, n_snps=len(a),
                PP0=res["PP0"], PP1=res["PP1"], PP2=res["PP2"],
                PP3=res["PP3"], PP4=res["PP4"],
                z_gwas_top=float(zg[k]), z_eqtl_top=float(ze[k]),
                dir_concordant=bool(concordant)))
    # write per-GWAS
    summary = {}
    for gname, glabel in GWAS_NAMES:
        rows = all_rows[gname]
        if not rows:
            summary[glabel] = dict(n_loci=0); continue
        d = pd.DataFrame(rows)
        d.to_csv(f"{out_prefix}_{glabel}.csv", index=False)
        pp = d["PP4"].values
        summary[glabel] = dict(
            n_loci=len(d),
            median_PP4=float(np.median(pp)),
            mean_PP4=float(np.mean(pp)),
            n_PP4_ge_0p8=int((pp >= 0.8).sum()),
            frac_ge_0p8=float((pp >= 0.8).mean()),
        )
    with open(f"{out_prefix}_summary.json", "w") as f:
        json.dump({"tissue": tissue_label, "gwas": summary}, f, indent=1)
    # print
    print("=" * 70)
    print(f"F4-literal — tissue eQTL: {tissue_label}")
    print("=" * 70)
    for glabel, s in summary.items():
        if s["n_loci"] == 0:
            print(f"  {glabel:10s}: no loci with >=20 aligned SNPs")
        else:
            print(f"  {glabel:10s}: n={s['n_loci']:2d} medianPP4={s['median_PP4']:.3f} "
                  f"meanPP4={s['mean_PP4']:.3f} PP4>=0.8: {s['n_PP4_ge_0p8']}/{s['n_loci']}")
    return summary


if __name__ == "__main__":
    eqtl_path = sys.argv[1]
    tissue = sys.argv[2]
    out_prefix = os.path.join(CACHE, f"f4_literal_{tissue}")
    eqtl = load_eqtl(eqtl_path)
    print(f"Loaded eQTL file: {eqtl_path}  n_snp={len(eqtl)}")
    run_tissue(eqtl, tissue, out_prefix)
