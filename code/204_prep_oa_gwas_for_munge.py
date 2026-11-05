#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
204_prep_oa_gwas_for_munge.py
=============================================================================
Convert Hatzikotoulas 2025 site-stratified OA GWAS (GWAS Catalog harmonised
TSV) into munge_sumstats.py-ready input for downstream cross-trait LDSC rg.

WHY THIS EXISTS
---------------
The OA GWAS files use POSITION-ENCODED variant IDs (e.g. "9_15867726_G_GC"),
not rsIDs. LDSC matches sumstats to LD scores by rsID, so every retained
variant must be mapped back to the rsID used by the 1000G phase-3 EUR panel
(G:/REF/1000G_EUR/EUR.bim). Coordinates are GRCh37 in BOTH files (verified:
35.1% of knee-GWAS (CHR,BP) pairs hit EUR.bim; EUR.bim carries 8,550,156
common SNPs vs ~24M GWAS variants, so the remaining ~65% are simply rarer
variants absent from the 1000G reference -- expected, not a build mismatch).

MUNGE INPUT WE EMIT  (tab-separated, gzip)
    SNP  A1  A2  P  OR  N
  - SNP : rsID from EUR.bim, matched on (CHR, BP) with allele-compatibility
  - A1  : GWAS effect_allele (upper-cased)
  - A2  : GWAS other_allele  (upper-cased)
  - P   : p_value  (munge derives Z = sqrt(chi2.isf(P,1)); sign from OR)
  - OR  : odds_ratio as the SIGNED_SUMSTAT (null value 1; OR<1 flips Z sign)
  - N   : NCASES + NCONTROLS, PER-ROW (this GWAS is a multi-cohort meta --
          N varies SNP by SNP; a single constant would be wrong)

NOTES / HONEST BOUNDARIES
-------------------------
* standard_error is "#NA" in these files -- harmless: munge does not use SE.
* Palindromic SNPs (A/T, C/G) are dropped: strand direction is ambiguous and
  could corrupt the sign of Z (and hence rg).
* Only biallelic ACGT SNPs are kept (indels / multiallelic dropped).
* Allele compatibility: the GWAS (effect, other) pair must be consistent with
  the EUR.bim (A1, A2) pair up to a swap of the two alleles (EUR.bim defines
  A1 as the minor allele, GWAS as the effect allele, so ~half are swapped --
  this is expected). Strand-complement cases are also accepted.
* --merge-alleles is NOT used: we perform the rsID/allele harmonisation here
  explicitly, so munge's HapMap3 whitelist is unnecessary. The regression
  SNP set is then determined by the intersection with the self-computed EUR
  LD-score reference at rg time.

USAGE
-----
  python3 code/204_prep_oa_gwas_for_munge.py --check       # verify inputs
  python3 code/204_prep_oa_gwas_for_munge.py --all         # all 4 OA sites
  python3 code/204_prep_oa_gwas_for_munge.py --sites knee spine
  python3 code/204_prep_oa_gwas_for_munge.py --sites knee --limit 200000
                                                        # smoke test subset

OUTPUT
------
  code/output/genetics/ldsc_munge_input/OA_<site>.munge_input.tsv.gz
  code/output/genetics/ldsc_munge_input/OA_<site>.prep_stats.json
"""
import argparse
import gzip
import json
import os
import sys
import time
from collections import Counter

# Python needs drive-letter paths on this Windows box (POSIX /g/ does NOT work
# inside Python, only inside the Git Bash shell itself).
HERE = "C:/Users/Administrator/WorkBuddy/2026-09-08-22-17-40/OA_B_PROJECT"
CACHE = "G:/OA_B_PROJECT_cache/l7"
BIM = "G:/REF/1000G_EUR/EUR.bim"
OUTDIR = os.path.join(HERE, "code", "output", "genetics", "ldsc_munge_input")

# Hatzikotoulas 2025, GWAS Catalog harmonised (GRCh37)
SITES = {
    "hip":   "GCST90566798.tsv.gz",
    "knee":  "GCST90566800.tsv.gz",
    "spine": "GCST90566801.tsv.gz",
    "hand":  "GCST90566803.tsv.gz",
}

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
VALID = set("ACGT")
PALINDROMIC = {frozenset(("A", "T")), frozenset(("C", "G"))}

GWAS_COLS = {
    "chromosome": 0, "base_pair_location": 1, "effect_allele": 2,
    "other_allele": 3, "odds_ratio": 4, "effect_allele_frequency": 6,
    "p_value": 7, "NCASES": 11, "NCONTROLS": 12,
}


def load_bim_map(path, verbose=True):
    """(CHR,BP) -> (rsID, A1, A2). Keep only biallelic ACGT SNPs.

    EUR.bim is plink format: CHR, SNP, CM(genetic dist, all 0), BP, A1, A2.
    Duplicate (CHR,BP) positions (multiallelic) are skipped to avoid
    ambiguous rsID assignment.
    """
    seen = {}
    dup = 0
    n = 0
    with open(path, "r") as fh:
        for line in fh:
            q = line.split()
            if len(q) < 6:
                continue
            a1, a2 = q[4].upper(), q[5].upper()
            if a1 not in VALID or a2 not in VALID or len(a1) != 1 or len(a2) != 1:
                continue
            try:
                key = (int(q[0]), int(q[3]))
            except ValueError:
                continue
            n += 1
            if key in seen:
                dup += 1
                seen[key] = None          # mark ambiguous, drop later
                continue
            seen[key] = (q[1], a1, a2)
    clean = {k: v for k, v in seen.items() if v is not None}
    if verbose:
        print(f"[bim] {n:,} biallelic SNPs -> {len(clean):,} unique (CHR,BP) "
              f"(dropped {dup:,} duplicate positions)", flush=True)
    return clean


def alleles_compatible(ea, oa, ref_a1, ref_a2):
    """Accept exact, swap, strand-flip, or swap+strand-flip.

    Reject palindromic (A/T, C/G) pairs: strand ambiguity would make the
    direction of the effect uninterpretable.
    """
    if frozenset((ea, oa)) in PALINDROMIC:
        return False
    if frozenset((ref_a1, ref_a2)) in PALINDROMIC:
        return False
    if ea == ref_a1 and oa == ref_a2:
        return True
    if ea == ref_a2 and oa == ref_a1:
        return True
    cea, coa = COMPLEMENT.get(ea, "?"), COMPLEMENT.get(oa, "?")
    if cea == ref_a1 and coa == ref_a2:
        return True
    if cea == ref_a2 and coa == ref_a1:
        return True
    return False


def prep_one(site, fname, bim, limit=None):
    src = os.path.join(CACHE, fname)
    dst = os.path.join(OUTDIR, f"OA_{site}.munge_input.tsv.gz")
    stats_path = os.path.join(OUTDIR, f"OA_{site}.prep_stats.json")
    stats = Counter()
    t0 = time.time()

    with gzip.open(src, "rt") as fin, gzip.open(dst, "wt") as fout:
        header = fin.readline().rstrip("\n").split("\t")
        idx = {c: header.index(c) for c in GWAS_COLS}
        fout.write("SNP\tA1\tA2\tP\tOR\tN\n")
        for line in fin:
            stats["rows_total"] += 1
            if limit and stats["rows_total"] > limit:
                break
            p = line.rstrip("\n").split("\t")
            if len(p) < 13:
                stats["drop_malformed"] += 1
                continue
            ea = p[idx["effect_allele"]].upper()
            oa = p[idx["other_allele"]].upper()
            if ea not in VALID or oa not in VALID or len(ea) != 1 or len(oa) != 1:
                stats["drop_non_biallelic"] += 1
                continue
            if frozenset((ea, oa)) in PALINDROMIC:
                stats["drop_palindromic"] += 1
                continue
            try:
                chrom = int(p[idx["chromosome"]])
                bp = int(p[idx["base_pair_location"]])
            except ValueError:
                stats["drop_bad_coord"] += 1
                continue
            ref = bim.get((chrom, bp))
            if ref is None:
                stats["drop_no_rsid"] += 1
                continue
            rsid, r1, r2 = ref
            if not alleles_compatible(ea, oa, r1, r2):
                stats["drop_allele_mismatch"] += 1
                continue
            try:
                pv = float(p[idx["p_value"]])
                OR = float(p[idx["odds_ratio"]])
                nc = float(p[idx["NCASES"]])
                nk = float(p[idx["NCONTROLS"]])
            except ValueError:
                stats["drop_bad_number"] += 1
                continue
            if not (0.0 < pv <= 1.0) or not (OR > 0) or OR == float("inf"):
                stats["drop_bad_p_or"] += 1
                continue
            if nc <= 0 or nk <= 0:
                stats["drop_bad_N"] += 1
                continue
            N = nc + nk
            fout.write(f"{rsid}\t{ea}\t{oa}\t{pv:.10g}\t{OR:.10g}\t{N:.0f}\n")
            stats["kept"] += 1

    stats["elapsed_sec"] = round(time.time() - t0, 1)
    stats["output"] = dst
    with open(stats_path, "w") as jf:
        json.dump(dict(stats), jf, indent=2, sort_keys=True)
    print(f"[{site}] rows={stats['rows_total']:,} kept={stats['kept']:,} "
          f"({100.0*stats['kept']/max(stats['rows_total'],1):.1f}%) "
          f"no_rsid={stats['drop_no_rsid']:,} "
          f"allele_mismatch={stats['drop_allele_mismatch']:,} "
          f"palindromic={stats['drop_palindromic']:,} "
          f"non_biallelic={stats['drop_non_biallelic']:,} "
          f"[{stats['elapsed_sec']}s]", flush=True)
    return dst


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sites", nargs="+", choices=list(SITES), default=None)
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--limit", type=int, default=None,
                    help="Only process first N rows per file (smoke test)")
    ap.add_argument("--check", action="store_true",
                    help="Verify inputs exist and report sizes, do not convert")
    args = ap.parse_args()

    os.makedirs(OUTDIR, exist_ok=True)
    sites = list(SITES) if args.all else (args.sites or ["knee"])

    print("=== 204 · OA GWAS -> munge_sumstats input ===", flush=True)
    for s in sites:
        f = os.path.join(CACHE, SITES[s])
        ok = os.path.exists(f)
        size = os.path.getsize(f) if ok else 0
        print(f"  {s:6s} {'OK ' if ok else 'MISS'} {size/1e6:8.1f} MB  {f}", flush=True)
    ok_bim = os.path.exists(BIM)
    print(f"  EUR.bim {'OK ' if ok_bim else 'MISS'} "
          f"{os.path.getsize(BIM)/1e6 if ok_bim else 0:.1f} MB  {BIM}", flush=True)
    if args.check:
        print("[check] done (no conversion performed)", flush=True)
        return 0
    if not ok_bim or any(not os.path.exists(os.path.join(CACHE, SITES[s])) for s in sites):
        print("ERROR: missing input(s)", file=sys.stderr)
        return 2

    print("[bim] loading rsID map ...", flush=True)
    bim = load_bim_map(BIM)
    for s in sites:
        prep_one(s, SITES[s], bim, limit=args.limit)
    print(f"[done] outputs in {OUTDIR}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
