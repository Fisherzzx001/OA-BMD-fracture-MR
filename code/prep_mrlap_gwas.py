#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Prepare full GWAS summary stats for MRlap (n-mounier/MRlap).

All input GWAS (BMD sites, fracture, OA_knee) are aligned to 1000G-phase3-hg19
(verified: allele-identity + nearest-position recovery gives median offset 0 vs
EUR.bim). rsIDs are therefore recovered by ALLELE-IDENTITY + NEAREST POSITION,
not by exact (chr,pos) lookup (the GWAS positions are slightly shifted vs the
1000G bim for many SNPs, and the leading rows are indels that are absent from
the bim).

Steps per trait:
  1. Load 1000G EUR.bim (hg19) -> per-chromosome sorted (pos, A1, A2, rsid).
  2. For each GWAS SNP (autosome, biallelic):
       - recover rsid by matching {effect,other} alleles (both strands) to the
         nearest EUR.bim SNP within +/-W bp; assign that rsid.
       - skip if no allele-matching SNP found (indels / non-1000G variants).
  3. Write clean TSV: rsid, chr, pos, a1(effect), a2(other), beta|or, se, N.
  4. Deduplicate by rsid (keep row with smallest p).

Output: code/ref/mrlap_gwas/<trait>.tsv
"""
import gzip, os, sys, time, math, re, bisect

EUR_BIM   = "H:/数据下载/1kg.v3/EUR.bim"
OUTDIR    = "C:/Users/Administrator/WorkBuddy/2026-09-08-22-17-40/OA_B_PROJECT/code/ref/mrlap_gwas"
SAMPLE    = "--sample" in sys.argv          # process only first SAMPLE_N SNPs per file
SAMPLE_N  = 80000
WIN       = 200                             # +/- bp window for nearest-position rsID recovery

TRAITS = {
    # name          : (path, type, N_or_Ncol, hg38?)
    # All traits are hg19 (no liftover needed).
    # N values are the manuscript Table-1 values (docs/08_投稿稿_OA_BMD_fracture.md):
    #   BMD sites 31,873 ; fracture 264,973 ; OA_knee 1,208,397.
    # OA_knee: the per-variant NCASES/NCONTROLS column in GCST90566800 is internally
    # inconsistent (ranges 605-1,316,504 with no mode), so the study-level N is used.
    "BMD_arms":       ("H:/UBK/Arms_BMD.tsv",          "cont", 31873, False),
    "BMD_femur_neck": ("H:/UBK/Femur_neck_BMD.tsv",    "cont", 31873, False),
    "BMD_leg":        ("H:/UBK/Leg_BMD.tsv",           "cont", 31873, False),
    "BMD_lumbar":     ("H:/UBK/LS_BMD.tsv",            "cont", 31873, False),
    "BMD_total":      ("H:/UBK/Total_BMD.tsv",         "cont", 31873, False),
    "fracture":       ("H:/UBK/fracture_any-type.tsv", "cont", 264973, False),
    "OA_knee":        ("H:/UBK/oa_gwas_catalog/GCST90566800.tsv.gz", "cc", 1208397, False),
}

def log(m): print(m, flush=True)

# ---------------------------------------------------------------------------
# 1. EUR.bim -> per-chromosome sorted (pos, A1, A2, rsid)  (hg19)
# ---------------------------------------------------------------------------
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}

def load_eur():
    log("Loading EUR.bim (hg19) for allele+position rsID recovery ...")
    chr_pos, chr_a1, chr_a2, chr_rs = {}, {}, {}, {}
    with open(EUR_BIM) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6: continue
            try:
                c = int(p[0]); pos = int(p[3]); rsid = p[1]
                a1 = p[4].upper(); a2 = p[5].upper()
            except: continue
            if not rsid.startswith("rs"): continue
            if a1 in ("A", "C", "G", "T") and a2 in ("A", "C", "G", "T"):
                chr_pos.setdefault(c, []).append(pos)
                chr_a1.setdefault(c, []).append(a1)
                chr_a2.setdefault(c, []).append(a2)
                chr_rs.setdefault(c, []).append(rsid)
    for c in list(chr_pos.keys()):
        idx = sorted(range(len(chr_pos[c])), key=lambda i: chr_pos[c][i])
        chr_pos[c] = [chr_pos[c][i] for i in idx]
        chr_a1[c]  = [chr_a1[c][i]  for i in idx]
        chr_a2[c]  = [chr_a2[c][i]  for i in idx]
        chr_rs[c]  = [chr_rs[c][i]  for i in idx]
    log(f"  EUR SNPs: {sum(len(v) for v in chr_pos.values()):,}")
    return chr_pos, chr_a1, chr_a2, chr_rs

def eur_rsid(c, pos, eff, oth, CP, CA1, CA2, CRS, W=WIN):
    """Recover rsID by allele-identity + nearest position within +/-W bp.
    Tries both strands. Returns rsid or None (picks closest allele-matching SNP)."""
    if c not in CP: return None
    P, A1, A2, RS = CP[c], CA1[c], CA2[c], CRS[c]
    lo = bisect.bisect_left(P, pos - W)
    hi = bisect.bisect_right(P, pos + W)
    ce = COMP.get(eff, eff); co = COMP.get(oth, oth)
    best = None; best_off = None
    for i in range(lo, hi):
        a = A1[i]; b = A2[i]
        if (a == eff and b == oth) or (a == oth and b == eff) or \
           (a == ce and b == co) or (a == co and b == ce):
            off = P[i] - pos
            if off < 0: off = -off
            if best is None or off < best_off:
                best = RS[i]; best_off = off
    return best

# ---------------------------------------------------------------------------
# 2. process one trait
# ---------------------------------------------------------------------------
def process(name, path, gtype, Nval, hg38, eur, C):
    log(f"\n=== {name} ({gtype}, hg38={hg38}) ===")
    CP, CA1, CA2, CRS = eur
    is_gz = path.endswith(".gz")
    fh = gzip.open(path, "rt") if is_gz else open(path)
    header = fh.readline().rstrip("\n").split("\t")
    ci = header.index("chromosome"); pi = header.index("base_pair_location")
    ai = header.index("effect_allele"); oi = header.index("other_allele")
    bi = header.index("odds_ratio") if gtype == "cc" else header.index("beta")
    si = header.index("standard_error")
    ei = header.index("effect_allele_frequency"); pvi = header.index("p_value")
    ncases_i = header.index("NCASES") if "NCASES" in header else None
    ncontrols_i = header.index("NCONTROLS") if "NCONTROLS" in header else None

    out = {}   # rsid -> best row (by p)
    total = matched = lifted = written = 0
    ns = 0     # SNP rows processed (sample mode stops after SAMPLE_N SNPs)
    t0 = time.time()
    for i, line in enumerate(fh):
        parts = line.rstrip("\n").split("\t")
        if len(parts) <= bi: continue
        total += 1
        try:
            c = int(parts[ci]); pos = int(parts[pi])
        except: continue
        eff = parts[ai].upper(); oth = parts[oi].upper()
        if eff in ("", "NA") or oth in ("", "NA"): continue
        # skip indels / non-biallelic SNPs (absent from 1000G bim)
        if eff not in ("A", "C", "G", "T") or oth not in ("A", "C", "G", "T"):
            continue
        if SAMPLE and ns >= SAMPLE_N: break
        ns += 1
        # beta / or
        beta_raw = parts[bi]
        se_raw = parts[si]
        try:
            beta = float(beta_raw)
        except:
            continue
        # se: may be missing (#NA) -> compute from CI for cc
        se = None
        try:
            se = float(se_raw)
        except:
            se = None
        p_str = parts[pvi]
        try:
            pval = float(p_str)
        except:
            pval = 1.0
        if gtype == "cc":
            if (se is None or se != se or se <= 0):
                try:
                    cu = float(parts[header.index("ci_upper")]); cl = float(parts[header.index("ci_lower")])
                    if cu > 0 and cl > 0:
                        se = (math.log(cu) - math.log(cl)) / (2 * 1.96)
                except:
                    se = None
            if se is None or se <= 0 or beta <= 0:
                continue
            # MRlap's tidy_inputGWAS has a bug in its "or" branch: it calls
            #   dplyr::mutate(Z = log(.data$beta)/.data$se, beta=NULL, se=NULL)
            # WITHOUT piping the data in, which throws
            # "Can't subset `.data` outside of a data mask context".
            # Workaround: pre-convert OR -> log(OR) here and label the column "beta",
            # so MRlap takes the (correct) continuous branch.
            effcol = "beta"; effval = math.log(beta)
            N = Nval          # study-level N; per-variant N column is unreliable
        else:
            if se is None or se <= 0:
                continue
            effcol = "beta"; effval = beta
            N = Nval
        # recover rsid (allele-identity + nearest position within +/-WIN bp)
        hg19_c, hg19_pos = c, pos
        rsid = eur_rsid(hg19_c, hg19_pos, eff, oth, CP, CA1, CA2, CRS)
        if rsid is None:
            continue
        matched += 1
        row = (rsid, hg19_c, hg19_pos, eff, oth, effcol, effval, se, N)
        if rsid not in out or pval < out[rsid][-1]:
            out[rsid] = row + (pval,)
    fh.close()
    # write
    os.makedirs(OUTDIR, exist_ok=True)
    outpath = os.path.join(OUTDIR, name + ".tsv")
    eff_header = "beta"   # always: case-control traits are pre-converted to log-OR
    with open(outpath, "w") as o:
        o.write("rsid\tchr\tpos\ta1\ta2\t%s\tse\tN\n" % eff_header)
        n = 0
        for rsid, hg19_c, hg19_pos, eff, oth, effcol, effval, se, N, pval in out.values():
            o.write(f"{rsid}\t{hg19_c}\t{hg19_pos}\t{eff}\t{oth}\t{effval}\t{se}\t{N}\n")
            n += 1
    written = n
    log(f"  total={total:,} snp_rows={ns:,} rsid_matched={matched:,} written={written:,}  ({time.time()-t0:.0f}s)")
    log(f"  -> {outpath}")
    return written

if __name__ == "__main__":
    t_start = time.time()
    eur = load_eur()
    for name, (path, gtype, Nval, hg38) in TRAITS.items():
        process(name, path, gtype, Nval, hg38, eur, None)
    log(f"\nALL DONE in {time.time()-t_start:.0f}s")
