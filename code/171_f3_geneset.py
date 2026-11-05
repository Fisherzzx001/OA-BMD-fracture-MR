#!/usr/bin/env python3
"""F3 — bone-versus-cartilage tissue enrichment (competitive gene-set test).

Self-contained MAGMA-style competitive test:
  1. Parse UCSC hg19 refGene -> per-gene genomic span.
  2. For each GWAS (OA_knee, BMD_total, fracture), assign every SNP to genes
     whose extended span (gene body +/- WINDOW) contains the SNP; gene-level
     statistic = min p among assigned SNPs (VEGAS-like).
  3. Competitive gene-set test per GO term: compare gene-level -log10(p) of
     set genes vs background (non-set) genes by Mann-Whitney U (one-sided:
     set MORE significant). Report mean ranks / mean -log10p and MWW p.
  4. Bone-vs-cartilage contrast: per trait, compare enrichment of bone-dev
     gene sets vs cartilage-dev gene sets.

Outputs:
  /g/OA_B_PROJECT_cache/t2/f3_geneset_results.csv
  /g/OA_B_PROJECT_cache/t2/f3_geneset_summary.txt
"""
import gzip, json, os, sys, math
import numpy as np
from scipy import stats

CACHE = "G:/OA_B_PROJECT_cache/t2"
REFGENE = os.path.join(CACHE, "refGene_hg19.txt.gz")
GOJSON = os.path.join(CACHE, "go_genesets.json")
LDSC = r"C:\Users\Administrator\WorkBuddy\2026-09-08-20-17-24\OA_FLAGSHIP\07_reference\ldsc"
WINDOW = 50000  # 50 kb up/downstream of gene body

TRAITS = {
    "OA_knee":      os.path.join(LDSC, "sumstats_oa_knee.sumstats.gz"),
    "BMD_total":    os.path.join(LDSC, "sumstats_total.sumstats.gz"),
    "fracture":     os.path.join(LDSC, "sumstats_fracture.sumstats.gz"),
}

# Bone-side and cartilage-side GO sets (from go_genesets.json labels)
BONE_SETS = ["bone_development", "ossification", "skeletal_system_development"]
CART_SETS = ["cartilage_development", "chondrocyte_differentiation",
             "extracellular_matrix_organization", "skeletal_joint_morphogenesis"]

# ----------------------------------------------------------------------------
def parse_refgene(path):
    """return dict gene_symbol -> (chrom, min_start, max_end) aggregated over tx."""
    genes = {}
    with gzip.open(path, "rt") if path.endswith(".gz") else open(path) as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 13:
                continue
            chrom = c[2]
            if not chrom.startswith("chr") or chrom in ("chrM", "chrUn"):
                continue
            try:
                start = int(c[4]); end = int(c[5])
            except ValueError:
                continue
            sym = c[12]
            if not sym or sym == "":
                continue
            prev = genes.get(sym)
            if prev is None:
                genes[sym] = [chrom, start, end]
            else:
                if start < prev[1]:
                    prev[1] = start
                if end > prev[2]:
                    prev[2] = end
    return genes

def load_gwas(path):
    """yield (chrom(int), bp, p) per SNP. chrom 'chr1'->1."""
    snps = {}  # chrom -> list of (bp, p)
    with gzip.open(path, "rt") if path.endswith(".gz") else open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ic = hdr.index("CHR"); ib = hdr.index("BP"); ip = hdr.index("P")
        for line in fh:
            c = line.rstrip("\n").split("\t")
            try:
                ch = c[ic]
                if ch.startswith("chr"):
                    ch = ch[3:]
                ch = int(ch)
                bp = int(float(c[ib]))
                p = float(c[ip])
            except (ValueError, IndexError):
                continue
            if p <= 0 or p > 1 or math.isnan(p):
                continue
            snps.setdefault(ch, []).append((bp, p))
    for ch in snps:
        snps[ch].sort()
    return snps

def gene_minp_per_trait(genes, gwas):
    """For each gene, min p among SNPs within [start-W, end+W]. Return dict gene-> -log10p."""
    out = {}
    for sym, (chrom, gstart, gend) in genes.items():
        if isinstance(chrom, str):
            chrom = chrom.replace("chr", "")
        try:
            ch = int(chrom)
        except ValueError:
            continue
        lo = gstart - WINDOW; hi = gend + WINDOW
        arr = gwas.get(ch)
        if not arr:
            continue
        # binary search
        import bisect
        i0 = bisect.bisect_left(arr, (lo, 0))
        i1 = bisect.bisect_right(arr, (hi, 2.0))
        if i1 <= i0:
            continue
        minp = min(p for (_, p) in arr[i0:i1])
        if minp <= 0:
            minp = 1e-300
        out[sym] = -math.log10(minp)
    return out

def mww_one_sided(metrics_set, metrics_back):
    """Mann-Whitney one-sided: set MORE significant than background.
    Returns (U, p_one_sided, mean_set, mean_back, n_set, n_back)."""
    if len(metrics_set) < 3 or len(metrics_back) < 3:
        return (None, None, float(np.mean(metrics_set)) if metrics_set else float('nan'),
                float(np.mean(metrics_back)) if metrics_back else float('nan'),
                len(metrics_set), len(metrics_back))
    # scipy: alternative='greater' tests set > back
    U, p = stats.mannwhitneyu(metrics_set, metrics_back, alternative="greater")
    return (U, p, float(np.mean(metrics_set)), float(np.mean(metrics_back)),
            len(metrics_set), len(metrics_back))

# ----------------------------------------------------------------------------
def main():
    print("Parsing refGene ...", flush=True)
    genes = parse_refgene(REFGENE)
    print(f"  {len(genes)} genes with coordinates", flush=True)

    print("Loading GO sets ...", flush=True)
    go = json.load(open(GOJSON))
    # map label -> set of symbols
    label2sym = {v["label"]: set(v["genes"]) for v in go.values()}

    results = []  # rows
    print("Processing traits ...", flush=True)
    for trait, path in TRAITS.items():
        print(f"  [{trait}] loading GWAS ...", flush=True)
        gwas = load_gwas(path)
        nsnps = sum(len(v) for v in gwas.values())
        print(f"    {nsnps} SNPs", flush=True)
        genep = gene_minp_per_trait(genes, gwas)
        print(f"    {len(genep)} genes mapped to >=1 SNP", flush=True)
        # background metric array (all mapped genes)
        back = list(genep.values())
        for label, syms in label2sym.items():
            in_syms = [s for s in syms if s in genep]
            if not in_syms:
                continue
            set_metrics = [genep[s] for s in in_syms]
            back_metrics = [genep[s] for s in genep if s not in syms]
            U, p, mset, mback, nset, nback = mww_one_sided(set_metrics, back_metrics)
            p_str = ("%.2e" % p) if (p is not None and isinstance(p, (int, float))) else "NA"
            results.append({
                "trait": trait, "set": label, "tissue": (
                    "bone" if label in BONE_SETS else
                    "cartilage" if label in CART_SETS else "other"),
                "n_genes_in_set": nset, "n_genes_background": nback,
                "mean_neglog10p_set": (round(float(mset), 4) if mset == mset else ""),
                "mean_neglog10p_background": (round(float(mback), 4) if mback == mback else ""),
                "enrichment_delta": (round(float(mset) - float(mback), 4) if (mset == mset and mback == mback) else ""),
                "mww_U": (round(U, 1) if U is not None else ""),
                "mww_p_one_sided": p_str,
            })
            print(f"    {label:38s} n={nset:3d} delta={(mset-mback) if (mset==mset and mback==mback) else float('nan'):+.3f} p={p_str}", flush=True)

    # write CSV
    import csv
    csv_path = os.path.join(CACHE, "f3_geneset_results.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(results[0].keys()))
        w.writeheader(); w.writerows(results)

    # summary
    summ = os.path.join(CACHE, "f3_geneset_summary.txt")
    with open(summ, "w") as f:
        f.write("F3 — bone-versus-cartilage competitive gene-set enrichment\n")
        f.write(f"refGene hg19, window +/-{WINDOW//1000}kb, gene-level stat = min p in window\n")
        f.write(f"Traits: OA_knee, BMD_total, fracture\n\n")
        for trait in TRAITS:
            f.write(f"=== {trait} ===\n")
            for r in results:
                if r["trait"] == trait:
                    f.write(f"  {r['set']:38s} [{r['tissue']:9s}] n={r['n_genes_in_set']:3d} "
                            f"mean-set={r['mean_neglog10p_set']:.3f} mean-bg={r['mean_neglog10p_background']:.3f} "
                            f"delta={r['enrichment_delta']:+.3f} MWW p={r['mww_p_one_sided']}\n")
            f.write("\n")
    print(f"WROTE {csv_path}\nWROTE {summ}")
    # also print summary
    print(open(summ).read())

if __name__ == "__main__":
    main()
