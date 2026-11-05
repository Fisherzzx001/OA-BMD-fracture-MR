#!/usr/bin/env python3
"""Extract skeletal-muscle (bone-proxy) eQTL signal at the 31 BMD_total instrument loci
from the Wilson 2025 skeletal-muscle eQTL meta-analysis (hg19), producing the same
loci TSV format consumed by code/172b (LOCI_IDX CHR BP SNP A1 A2 MAF N Z).

MAF is sourced from the cached eQTLGen European allele-frequency file
(G:/OA_B_PROJECT_cache/eQTLGen_SNP_AF.txt.gz) matched on rsID, so the muscle and
cartilage arms use a consistent, real MAF source. SNPs absent from eQTLGen get a
documented default MAF=0.25.

Usage: python code/173_muscle_eqtl_loci.py
"""
import os, sys, gzip, csv
import pandas as pd
import numpy as np

ROOT = "C:/Users/Administrator/WorkBuddy/2026-09-08-22-17-40/OA_B_PROJECT"
HARM = os.path.join(ROOT, "code", "output", "genetics", "harmonised", "BMD_total__OA_knee.csv")
MSDIR = "G:/OA_B_PROJECT_cache/t2/MarginalSumStats_hg19"
AF = "G:/OA_B_PROJECT_cache/eQTLGen_SNP_AF.txt.gz"
OUT = "G:/OA_B_PROJECT_cache/t2/muscle_eqtl_loci.tsv"
WINDOW = 500_000
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
MUSCLE_N = 1002   # Wilson et al. AJHG 2025 skeletal-muscle eQTL meta-analysis (GTEx+FUSION), N=1,002 individuals
DEFAULT_MAF = 0.25


def load_maf():
    """rsID -> minor allele frequency, from eQTLGen AF."""
    d = {}
    with gzip.open(AF, "rt") as f:
        h = f.readline().rstrip("\n").split("\t")
        i_snp = h.index("SNP"); i_a = h.index("AlleleA"); i_b = h.index("AlleleB"); i_ball = h.index("AlleleB_all")
        for line in f:
            c = line.rstrip("\n").split("\t")
            rs = c[i_snp]
            aA, aB = c[i_a], c[i_b]
            try:
                fB = float(c[i_ball])
            except ValueError:
                continue   # NA / missing AF -> skip; caller falls back to DEFAULT_MAF
            fA = 1.0 - fB
            # store both allele freqs so caller can pick minor regardless of orientation
            d[rs] = (aA, aB, fA, fB)
    return d


def maf_for(rsid, ea, nea, af):
    """Return minor allele frequency for the (ea, nea) pair given eQTLGen AF."""
    if rsid in af:
        aA, aB, fA, fB = af[rsid]
        if ea == aB:
            f_ea = fB
        elif ea == aA:
            f_ea = fA
        else:
            return DEFAULT_MAF
        return min(f_ea, 1.0 - f_ea)
    return DEFAULT_MAF


def main():
    inst = pd.read_csv(HARM)
    by_chrom = {}
    for i, rw in inst.iterrows():
        c = int(rw["CHR"]); p = int(rw["BP"])
        by_chrom.setdefault(c, []).append((i, p - WINDOW, p + WINDOW))
    af = load_maf()
    print(f"MAF lookup loaded: {len(af):,} rsIDs", flush=True)

    best = {}   # (loci_idx, snpkey) -> (abs_z, z, ea, nea, maf, rsid)
    n_in = 0
    for chrom, wins in by_chrom.items():
        fp = os.path.join(MSDIR, f"Muscle_eQTL_MetaAnalysis_1e5_fullMarginalSumStats_chr{chrom}_hg19.txt")
        if not os.path.exists(fp):
            print(f"  missing {fp}", flush=True); continue
        with open(fp) as f:
            hdr = f.readline().rstrip("\n").split("\t")
            idx = {name: k for k, name in enumerate(hdr)}
            i_pos = idx["pos"]; i_ne = idx["NEA"]; i_ea = idx["EA"]; i_b = idx["beta"]; i_se = idx["se"]; i_rs = idx["rsID"]; i_sh = idx["SNP_hg19"]
            n = 0
            for line in f:
                n += 1
                if n % 5000000 == 0:
                    print(f"  chr{chrom}: scanned {n:,} rows, in-locus kept {n_in:,}, entries {len(best):,}", flush=True)
                c = line.rstrip("\n").split("\t")
                try:
                    bp = int(c[i_pos])
                except ValueError:
                    continue
                for (li, lo, hi) in wins:
                    if lo <= bp <= hi:
                        try:
                            beta = float(c[i_b]); se = float(c[i_se])
                        except ValueError:
                            break
                        if se <= 0:
                            break
                        z = beta / se
                        rsid = c[i_rs].strip()
                        snpkey = rsid if rsid and rsid != "." else c[i_sh]
                        ea = c[i_ea]; nea = c[i_ne]
                        maf = maf_for(rsid, ea, nea, af)
                        key = (li, snpkey)
                        cur = best.get(key)
                        az = abs(z)
                        if cur is None or az > cur[0]:
                            best[key] = (az, z, ea, nea, maf, rsid, bp)
                            n_in += 1 if cur is None else 0
                        break
    # write
    rows = []
    for (li, snpkey), v in best.items():
        az, z, ea, nea, maf, rsid, bp = v
        chrom = int(inst.iloc[li]["CHR"])
        rows.append((li, chrom, int(bp), snpkey, ea, nea, round(maf, 4), MUSCLE_N, round(z, 4)))
    rows.sort(key=lambda r: (r[0], r[1]))
    with open(OUT, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["LOCI_IDX", "CHR", "BP", "SNP", "A1", "A2", "MAF", "N", "Z"])
        for r in rows:
            w.writerow(r)
    print(f"\nWrote {len(rows):,} muscle eQTL locus-SNP rows to {OUT} (in-locus kept {n_in:,})", flush=True)
    # sanity: per-locus counts
    from collections import Counter
    cnt = Counter(r[0] for r in rows)
    print("loci with >=1 SNP:", sum(1 for v in cnt.values() if v >= 1), "of", len(by_chrom))


if __name__ == "__main__":
    main()
