#!/usr/bin/env python3
"""F4-literal (cartilage arm): extract FunGen low-grade cartilage eQTL signal at the
31 BMD_total instrument loci, building a per-SNP cartilage eQTL summary
(CHR BP SNP A1 A2 MAF N Z) for downstream colocalisation.

For each SNP (genotype_id = chr:pos, hg19) we keep its STRONGEST cartilage
eQTL association (min nominal p across genes) — the cartilage regulatory
signal most likely to colocalise with a GWAS hit. Z = slope/slope_se (Wald,
signed for ALT = effect allele per FunGen header).

Output: G:/OA_B_PROJECT_cache/t2/cartilage_eqtl_loci.tsv  (N=200 documentned)
"""
import os, sys, subprocess, gzip

ZIP = "G:/OA_B_PROJECT_cache/t2/FunGen_LowGradeCartilage.zip"
INNER = "FunGen_eQTL_LowGradeCartilage.FastQTL_perm_nom_info.ForMSK-KP_16Jan2021.txt.gz"
INNER_GZ = "G:/OA_B_PROJECT_cache/t2/FunGen_inner.txt.gz"  # extracted to disk (zcat unavailable on Win)
HARM = "C:/Users/Administrator/WorkBuddy/2026-09-08-22-17-40/OA_B_PROJECT/code/output/genetics/harmonised/BMD_total__OA_knee.csv"
OUT = "G:/OA_B_PROJECT_cache/t2/cartilage_eqtl_loci.tsv"
WINDOW = 500_000
CART_N = 200  # FunGen cartilage donor count (documented; approximate)

def parse_info(s):
    for kv in s.split(";"):
        if kv.startswith("MAF="):
            try:
                return float(kv[4:])
            except ValueError:
                return float("nan")
    return float("nan")

def main():
    # load 31 loci
    import csv
    loci = []  # (idx, chrom, lo, hi, pos)
    with open(HARM) as f:
        r = csv.DictReader(f)
        for i, row in enumerate(r):
            ch = int(row["CHR"]); pos = int(row["BP"])
            loci.append((i, ch, pos - WINDOW, pos + WINDOW, pos))
    by_chrom = {}
    for (i, ch, lo, hi, pos) in loci:
        by_chrom.setdefault(ch, []).append((i, lo, hi))
    print(f"Loaded {len(loci)} loci across chroms {sorted(by_chrom)}", flush=True)

    # best[(idx, genotype_id)] = (pval, slope, slope_se, ref, alt, maf)
    best = {}
    n_in = 0
    with gzip.open(INNER_GZ, "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        idx = {name: k for k, name in enumerate(hdr)}
        i_gid = idx["genotype_id"]; i_pv = idx["pval"]; i_sl = idx["slope"]
        i_se = idx["slope_se"]; i_ref = idx["REF"]; i_alt = idx["ALT"]; i_info = idx["INFO"]
        n = 0
        for line in fh:
            n += 1
            if n % 5000000 == 0:
                print(f"  scanned {n:,} rows, in-locus kept {n_in:,}, best entries {len(best):,}", flush=True)
            c = line.rstrip("\n").split("\t")
            gid = c[i_gid]            # e.g. 1:100000012
            try:
                cc, pp = gid.split(":")
                chrom = int(cc); bp = int(pp)
            except ValueError:
                continue
            wins = by_chrom.get(chrom)
            if not wins:
                continue
            inside = False
            for (li, lo, hi) in wins:
                if lo <= bp <= hi:
                    inside = True; lidx = li
                    break
            if not inside:
                continue
            try:
                pv = float(c[i_pv]); sl = float(c[i_sl]); se = float(c[i_se])
            except ValueError:
                continue
            if se <= 0 or pv <= 0:
                continue
            maf = parse_info(c[i_info])
            if not (0 < maf < 1):
                continue
            key = (lidx, gid)
            cur = best.get(key)
            if cur is None or pv < cur[0]:
                best[key] = (pv, sl, se, c[i_ref], c[i_alt], maf)
                n_in += 1 if cur is None else 0
    print(f"Scanned {n:,} FunGen rows. Writing {len(best):,} in-locus SNP-gene-best records...", flush=True)
    # sort by (locus idx, bp)
    items = sorted(best.items(), key=lambda kv: (kv[0][0], int(kv[0][1].split(":")[1])))
    with open(OUT, "w") as f:
        f.write("LOCI_IDX\tCHR\tBP\tSNP\tA1\tA2\tMAF\tN\tZ\n")
        for (lidx, gid), (pv, sl, se, ref, alt, maf) in items:
            cc, pp = gid.split(":")
            z = sl / se
            f.write(f"{lidx}\t{cc}\t{pp}\t{gid}\t{ref}\t{alt}\t{maf:.5f}\t{CART_N}\t{z:.4f}\n")
    print(f"WROTE {OUT}")

if __name__ == "__main__":
    main()
