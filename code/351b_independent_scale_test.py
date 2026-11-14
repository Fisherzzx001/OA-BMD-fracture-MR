# -*- coding: utf-8 -*-
"""
code/351b_independent_scale_test.py
====================================
M1 -- Variant-level scale-reconstruction test against a GENUINELY INDEPENDENT
BMD GWAS: GEFOS total-body BMD (GCST005348; Medina-Gomez et al., AJHG 2018),
N = 56,284 European, DXA-derived total-body BMD.

WHY THIS FILE EXISTS
--------------------
code/351 --independent mode matches on the SNP *string* (rsID) and expects the
exposure harmonised file to carry A1/A2.  Two data realities block that path:
  (i)  The independent GEFOS harmonised release keys on `variant_id` = "chr:pos"
       (GRCh37) and carries NO rsID column.
  (ii) The exposure harmonised BMD_total__fracture.csv carries beta_e/se_e for
       only the 31 MR instruments and NO A1/A2 column.
Both files are GRCh37 (verified: GEFOS 1:22702231 == rs12568930 GRCh37 coord,
matching the exposure's BP).  So we match by (CHR, BP) and pull the exposure's
A1/A2 from the munged Qian file (BMD_total_N31873.sumstats.gz) by rsID.

TEST LOGIC (identical adjudication to code/351)
------------------------------------------------
The exposure beta_e is on the DEPOSITED scale.  The manuscript multiplies it by
sqrt(400000/31873) = 3.5426 to "reconstruct" the true scale.
  * median(beta_GEFOS / beta_exposure) ~ 1.0       -> 3.5426 is OVER-INFLATION (WRONG)
  * median(beta_GEFOS / beta_exposure) ~ 3.5426     -> reconstruction CONFIRMED (CORRECT)
  * neither                                                -> INCONCLUSIVE

USAGE
-----
  python3 code/351b_independent_scale_test.py
  (reads GEFOS file from code/output/genetics/independent_bmd/...build37.f.tsv.gz)
"""
import os, sys, gzip, csv, json, statistics

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HARMONISED = os.path.join(PROJ, "code", "output", "genetics", "harmonised")
MUNGE_DIR  = os.path.join(PROJ, "code", "output", "genetics", "ldsc_rg_verified", "munged")
INDEP_DIR  = os.path.join(PROJ, "code", "output", "genetics", "independent_bmd")
OUT_DIR    = os.path.join(PROJ, "code", "output", "genetics", "m1_audit")
OUT_JSON   = os.path.join(OUT_DIR, "m1_independent_gefs_tbmd_test.json")

EXP_HARM   = os.path.join(HARMONISED, "BMD_total__fracture.csv")
QIAN_MUNGE = os.path.join(MUNGE_DIR, "BMD_total_N31873.sumstats.gz")
GEFOS      = os.path.join(INDEP_DIR, "GCST005348_harmonised_build37.f.tsv.gz")

SCALE_FACTOR_CLAIMED = (400000.0 / 31873.0) ** 0.5   # 3.5426

def load_exposure():
    """Return rsID -> dict(CHR,BP,beta_e,se_e) from the 31-instrument harmonised file."""
    exp = {}
    with open(EXP_HARM, newline="") as fh:
        for row in csv.DictReader(fh):
            try:
                exp[row["SNP"]] = {
                    "chr": row["CHR"], "bp": row["BP"],
                    "beta_e": float(row["beta_e"]), "se_e": float(row["se_e"]),
                }
            except (KeyError, ValueError):
                continue
    return exp

def load_qian_alleles():
    """rsID -> (A1, A2) from the munged Qian exposure (GRCh37)."""
    alle = {}
    with gzip.open(QIAN_MUNGE, "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(hdr)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            snp = p[idx.get("SNP", 0)]
            if not snp.startswith("rs"):
                continue
            a1 = p[idx.get("A1", 1)].upper(); a2 = p[idx.get("A2", 2)].upper()
            if a1 and a2:
                alle[snp] = (a1, a2)
    return alle

def load_gefs_matches(exp_keys):
    """Stream GEFOS; keep only rows whose (chromosome, base_pair_location) is in exp_keys.
    exp_keys: set of 'chr:bp' strings. Returns dict key -> (effect_allele, other_allele, beta, se)."""
    matches = {}
    with gzip.open(GEFOS, "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(hdr)}
        i_ea = idx["effect_allele"]; i_oa = idx["other_allele"]
        i_b = idx["beta"]; i_se = idx["standard_error"]
        # GEFOS harmonised usually carries `chromosome`/`base_pair_location`; some
        # releases encode position only inside `variant_id` ("chr:pos"). Fall back.
        i_chr = idx.get("chromosome"); i_bp = idx.get("base_pair_location")
        i_vid = idx.get("variant_id")
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if i_chr is not None and i_bp is not None:
                key = f"{p[i_chr]}:{p[i_bp]}"
            elif i_vid is not None:
                key = p[i_vid].split(":")[:2]
                key = ":".join(key) if len(key) == 2 else None
            else:
                continue
            if not key or key not in exp_keys:
                continue
            try:
                b = float(p[i_b]); se = float(p[i_se])
            except ValueError:
                continue
            matches[key] = (p[i_ea].upper(), p[i_oa].upper(), b, se)
    return matches

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    if not (os.path.exists(EXP_HARM) and os.path.exists(QIAN_MUNGE) and os.path.exists(GEFOS)):
        print("[FATAL] missing input:",
              {k: os.path.exists(v) for k, v in
               [("EXP_HARM", EXP_HARM), ("QIAN_MUNGE", QIAN_MUNGE), ("GEFOS", GEFOS)]})
        sys.exit(2)

    exp = load_exposure()
    alle = load_qian_alleles()
    exp_keys = {f"{v['chr']}:{v['bp']}" for v in exp.values()}
    gef = load_gefs_matches(exp_keys)
    print(f"[load] exposure instruments = {len(exp)}; qian alleles = {len(alle)}; "
          f"GEFOS matches = {len(gef)}")

    ratios, zratios, detail = [], [], []
    n_ambig = 0
    for rsid, e in exp.items():
        key = f"{e['chr']}:{e['bp']}"
        g = gef.get(key)
        if not g:
            detail.append({"snp": rsid, "status": "no_gefs_match"}); continue
        ea, oa, gb, gse = g
        qa = alle.get(rsid)
        if not qa:
            detail.append({"snp": rsid, "status": "no_qian_allele"}); continue
        a1, a2 = qa
        # orient GEFOS beta to the exposure A1 (deposited effect allele) frame
        if ea == a1:
            gbo = gb
        elif ea == a2:
            gbo = -gb
        else:
            n_ambig += 1
            detail.append({"snp": rsid, "status": "ambiguous_allele",
                           "gefs_ea": ea, "qian_a1": a1, "qian_a2": a2})
            continue
        if abs(e["se_e"]) < 1e-12 or abs(gse) < 1e-12:
            continue
        r = gbo / e["beta_e"]
        zr = (gbo / gse) / (e["beta_e"] / e["se_e"])
        ratios.append(r); zratios.append(zr)
        detail.append({"snp": rsid, "beta_exp": round(e["beta_e"], 5),
                       "beta_gefs": round(gb, 5), "beta_gefs_oriented": round(gbo, 5),
                       "beta_ratio": round(r, 4)})

    n = len(ratios)
    med = statistics.median(ratios) if ratios else float("nan")
    medz = statistics.median(zratios) if zratios else float("nan")
    q1 = statistics.quantiles(ratios, n=4)[0] if n >= 4 else (min(ratios) if ratios else float("nan"))
    q3 = statistics.quantiles(ratios, n=4)[2] if n >= 4 else (max(ratios) if ratios else float("nan"))

    tol = 0.25 * SCALE_FACTOR_CLAIMED
    if abs(med - 1.0) <= 0.5:
        adj = "RECONSTRUCTION_OVER_INFLATED_WRONG"
        txt = (f"median(beta_GEFOS/beta_exp)={med:.3f} ~ 1.0 => the deposited exposure beta "
               "was already on the true scale; the x3.5426 reconstruction is a 3.5x "
               "over-inflation and must be removed.")
    elif abs(med - SCALE_FACTOR_CLAIMED) <= tol:
        adj = "RECONSTRUCTION_CONFIRMED_CORRECT"
        txt = (f"median(beta_GEFOS/beta_exp)={med:.3f} ~ {SCALE_FACTOR_CLAIMED:.3f} => the deposited "
               "exposure beta was compressed; x3.5426 is confirmed.")
    else:
        adj = "INCONCLUSIVE"
        txt = (f"median(beta_GEFOS/beta_exp)={med:.3f} is neither ~1.0 nor ~{SCALE_FACTOR_CLAIMED:.3f}; "
               "trait/site mismatch or residual scale issue -- investigate.")

    result = {
        "test": "exposure_beta_vs_independent_GEFOS_total_body_BMD",
        "independent_source": "GCST005348 (GEFOS Lifecourse total-body BMD, Medina-Gomez AJHG 2018), N=56284 EUR",
        "scale_factor_claimed": SCALE_FACTOR_CLAIMED,
        "n_exposure_instruments": len(exp),
        "n_gefs_matched": len(gef),
        "n_used": n,
        "n_ambiguous_allele": n_ambig,
        "median_beta_ratio": med,
        "q1_beta_ratio": q1, "q3_beta_ratio": q3,
        "median_Z_ratio": medz,
        "adjudication": adj,
        "adjudication_text": txt,
        "can_adjudicate_M1": adj != "INCONCLUSIVE",
        "caveat": ("Total-body BMD vs total-body BMD (trait-matched); GRCh37 position match; "
                   "31 MR instruments only -- a focused but directly relevant test of the "
                   "instruments that drive the MR."),
        "detail": detail,
    }
    with open(OUT_JSON, "w") as fh:
        json.dump(result, fh, indent=2)

    print("\n" + "#" * 78)
    print("HEADLINE ADJUDICATION:", adj)
    print(txt)
    print(f"  n_used={n}  median_beta_ratio={med:.4f}  (Q1={q1:.3f}, Q3={q3:.3f})  "
          f"median_Z_ratio={medz:.4f}")
    print("#" * 78)
    print("[output] JSON ->", OUT_JSON)

if __name__ == "__main__":
    main()
