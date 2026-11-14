# -*- coding: utf-8 -*-
"""
code/351_m1_variant_level_scale_test.py
=======================================
M1 — Variant-level scale-reconstruction test for the BMD exposure.

PURPOSE
-------
The manuscript reconstructs the BMD effect scale by multiplying the deposited
exposure beta/SE by sqrt(400000/31873) = 3.5426 (a "placeholder N = 400,000"
correction).  The intended independent anchor for this reconstruction was meant
to be a *separate* BMD GWAS: Qian 2025 (N = 31,873, file BMD_total_N31873) vs
UKB total-body DXA BMD / GCST90568455 (N = 31,986, file BMD_total_N31986).

CRITICAL PROVENANCE FINDING (verified 2026-11-13, see code/351 run output)
-------------------------------------------------------------------------
BMD_total_N31873 and BMD_total_N31986 are the SAME underlying summary
statistics.  Across 1,217,311 SNPs the A1/A2/Z columns are byte-identical;
only the N column differs (31873 vs 31986).  code/201_l8_bc_coloc.py documents
the source as "sumstats_total.sumstats.gz, true N = 31986 ... SAME exposure".
=> N31986 is NOT an independent study; it is the Qian exposure relabelled.
=> Comparing Qian(N31873) vs UKB(N31986) is CIRCULAR and CANNOT adjudicate M1.

This script:
  1. Loads both munged files, harmonises alleles, computes the per-SNP Z ratio.
  2. Runs a PROVENANCE / CIRCULARITY GUARD that detects identical Z and refuses
     to report a misleading "validation".
  3. Emits a NON_INDEPENDENT verdict with the supporting evidence.
  4. Documents the CORRECT test (requires a genuinely independent BMD GWAS that
     carries beta & SE) and the exact public accessions to obtain one.
  5. Optional --independent <file> mode: runs the real beta-ratio test against a
     user-supplied independent BMD GWAS (columns: SNP,A1,A2,beta,se).

WHY THE MUNGED FILES ALONE CANNOT TEST THE BETA-SCALE FACTOR
-----------------------------------------------------------
The munged files store only (SNP, A1, A2, Z, N).  Z = beta / SE is N-normalised,
so Z_Qian / Z_UKB ~ sqrt(31873/31986) ~ 1.0 for ANY two studies of the same
trait and ~same N.  The 3.5426 factor operates on beta/SE, NOT on Z.  Therefore
a Z-only comparison can never reveal (or refute) a 3.5426 beta-scale factor.
The decisive comparison requires raw beta & SE from BOTH an exposure file and a
genuinely independent BMD GWAS -- which is why --independent mode exists.

USAGE
-----
  python3 code/351_m1_variant_level_scale_test.py
        # default: Qian(N31873) vs UKB(N31986) with circularity guard

  python3 code/351_m1_variant_level_scale_test.py --independent <bmd.csv>
        # real variant-level test vs a genuinely independent BMD GWAS
        # (exposure beta taken from the harmonised BMD_total__fracture.csv)

  python3 code/351_m1_variant_level_scale_test.py --try-download
        # attempt to fetch ieu-a-980 (GEFOS femoral-neck BMD) via OpenGWAS
        # (requires network; degrades gracefully if unreachable)

Author: remediation pass (M1), 2026-11-13.
"""

import os
import sys
import gzip
import json
import csv
import argparse
import statistics

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MUNGE_DIR = os.path.join(PROJ, "code", "output", "genetics", "ldsc_rg_verified", "munged")
HARMONISED_DIR = os.path.join(PROJ, "code", "output", "genetics", "harmonised")
OUT_DIR = os.path.join(PROJ, "code", "output", "genetics", "m1_audit")
OUT_JSON = os.path.join(OUT_DIR, "m1_variant_level_scale_test.json")

QIAN = os.path.join(MUNGE_DIR, "BMD_total_N31873.sumstats.gz")   # exposure, true N=31873
UKB  = os.path.join(MUNGE_DIR, "BMD_total_N31986.sumstats.gz")   # labelled GCST90568455, N=31986

# The manuscript's claimed scale factor (placeholder N = 400,000 correction).
SCALE_FACTOR_CLAIMED = (400000.0 / 31873.0) ** 0.5   # = 3.5426

EXPOSURE_HARMONISED = os.path.join(HARMONISED_DIR, "BMD_total__fracture.csv")


# --------------------------------------------------------------------------- #
#  Loading
# --------------------------------------------------------------------------- #
def load_munge(path):
    """Return {snp: (A1, A2, Z)} from a ldsc-munged file (SNP A1 A2 Z N)."""
    data = {}
    skipped = 0
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        # tolerate space-separated as well
        if len(header) < 4:
            header = header[0].split()
        idx = {h: i for i, h in enumerate(header)}
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                parts = line.rstrip("\n").split()
            snp = parts[idx.get("SNP", 0)]
            a1 = parts[idx.get("A1", 1)].upper()
            a2 = parts[idx.get("A2", 2)].upper()
            try:
                z = float(parts[idx.get("Z", 3)])
            except (ValueError, IndexError):
                skipped += 1
                continue
            data[snp] = (a1, a2, z)
    return data, skipped


def load_independent_beta(path):
    """Load an independent BMD GWAS carrying SNP,A1,A2,beta,se (CSV/TSV)."""
    data = {}
    with open(path, newline="") as fh:
        sniff = fh.read(2048)
        fh.seek(0)
        delim = "\t" if "\t" in sniff else ","
        reader = csv.DictReader(fh, delimiter=delim)
        f_snp = reader.fieldnames[0]
        for row in reader:
            snp = row.get("SNP") or row.get(f_snp)
            a1 = (row.get("A1") or "").upper()
            a2 = (row.get("A2") or "").upper()
            try:
                b = float(row.get("beta"))
                se = float(row.get("se"))
            except (TypeError, ValueError):
                continue
            data[snp] = (a1, a2, b, se)
    return data


def load_exposure_beta_from_harmonised(path):
    """Exposure beta_e / se_e already harmonised (Qian deposited scale)."""
    data = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            snp = row["SNP"]
            try:
                b = float(row["beta_e"]); se = float(row["se_e"])
                a1 = (row.get("A1") or "").upper()
                a2 = (row.get("A2") or "").upper()
            except (KeyError, ValueError):
                continue
            data[snp] = (a1, a2, b, se)
    return data


# --------------------------------------------------------------------------- #
#  Harmonisation helper
# --------------------------------------------------------------------------- #
def harmonize_sign(z_or_b, a1_e, a2_e, a1_o, a2_o):
    """Flip sign of the 'other' effect if alleles are swapped (strand flip).

    Returns the other effect oriented to the exposure allele frame, or None if
    alleles are ambiguous (neither match nor swap)."""
    if a1_e == a1_o and a2_e == a2_o:
        return z_or_b, True          # concordant
    if a1_e == a2_o and a2_e == a1_o:
        return -z_or_b, True         # swapped -> flip
    return None, False               # ambiguous


# --------------------------------------------------------------------------- #
#  Statistics
# --------------------------------------------------------------------------- #
def pearson(xs, ys):
    n = len(xs)
    if n < 2:
        return float("nan")
    mx = sum(xs) / n
    my = sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx == 0 or syy == 0:
        return float("nan")
    return sxy / (sxx ** 0.5 * syy ** 0.5)


# --------------------------------------------------------------------------- #
#  Core 1: Qian vs UKB circularity guard
# --------------------------------------------------------------------------- #
def run_qian_vs_ukb():
    print("=" * 78)
    print("M1 variant-level test  —  Qian(N31873) vs UKB(N31986)")
    print("=" * 78)
    if not (os.path.exists(QIAN) and os.path.exists(UKB)):
        print("[FATAL] munged files not found:", QIAN, UKB)
        sys.exit(2)

    q, q_skip = load_munge(QIAN)
    u, u_skip = load_munge(UKB)
    print(f"[load] Qian  SNPs={len(q):,}  skipped(blank Z)={q_skip:,}")
    print(f"[load] UKB   SNPs={len(u):,}  skipped(blank Z)={u_skip:,}")

    common = set(q) & set(u)
    print(f"[overlap] common SNPs = {len(common):,}")

    zq, zu, ratios, exact = [], [], [], 0
    allele_concordant = 0
    for snp in common:
        a1q, a2q, zq0 = q[snp]
        a1u, a2u, zu0 = u[snp]
        zq_h, okq = harmonize_sign(zq0, a1q, a2q, a1u, a2u)
        zu_h, oku = harmonize_sign(zu0, a1u, a2u, a1q, a2q)
        if not (okq and oku):
            continue
        if a1q == a1u and a2q == a2u:
            allele_concordant += 1
        zq.append(zq_h); zu.append(zu_h)
        if abs(zu_h) > 1e-12:
            ratios.append(zq_h / zu_h)
        if abs(zq_h - zu_h) < 1e-9:
            exact += 1

    n = len(zq)
    r = pearson(zq, zu)
    med_ratio = statistics.median(ratios) if ratios else float("nan")
    exact_frac = exact / n if n else float("nan")
    same_set = (len(q) == len(u) == len(common))

    print(f"[metrics] n_harmonised = {n:,}")
    print(f"[metrics] Pearson r(Zq, Zu)      = {r:.6f}")
    print(f"[metrics] median Z-ratio (Q/U)   = {med_ratio:.6f}")
    print(f"[metrics] exact-Z-match fraction = {exact_frac:.6f}")
    print(f"[metrics] identical SNP set       = {same_set}")
    print(f"[metrics] allele-concordant frac  = {allele_concordant / n:.6f}" if n else "")

    # ---- circularity verdict ----
    circular = (exact_frac > 0.999) and (r > 0.9999) and same_set
    if circular:
        verdict = "CIRCULAR_NON_INDEPENDENT"
        explain = (
            "The two files are the SAME summary statistics (identical A1/A2/Z "
            "across all SNPs; only the N column differs). This is a relabelled "
            "copy of the exposure, NOT an independent BMD GWAS. A comparison "
            "between them is circular and CANNOT adjudicate the 3.5426 "
            "scale-reconstruction factor."
        )
    else:
        verdict = "INCONCLUSIVE_REVIEW_MANUALLY"
        explain = (
            "Z values differ beyond numerical noise; review the metrics above "
            "to decide independence. Note Z-only files still cannot test the "
            "beta-scale factor (see --independent mode)."
        )

    result = {
        "test": "Qian_N31873_vs_UKB_N31986_z_ratio",
        "verdict": verdict,
        "explanation": explain,
        "scale_factor_claimed": SCALE_FACTOR_CLAIMED,
        "n_qian": len(q),
        "n_ukb": len(u),
        "n_common": len(common),
        "n_harmonised": n,
        "pearson_r_Zq_Zu": r,
        "median_Z_ratio": med_ratio,
        "exact_Z_match_fraction": exact_frac,
        "identical_snp_set": same_set,
        "allele_concordant_fraction": (allele_concordant / n) if n else None,
        "skipped_blank_Z": {"qian": q_skip, "ukb": u_skip},
        "can_adjudicate_M1": False,
        "reason_cannot_adjudicate": (
            "Pair is circular (same source). Even if independent, Z-only files "
            "cannot reveal a beta-scale factor; the correct test needs beta&SE "
            "from a genuinely independent BMD GWAS (--independent mode)."
        ),
    }
    return result


# --------------------------------------------------------------------------- #
#  Core 2: correct beta-ratio test against a genuinely independent BMD GWAS
# --------------------------------------------------------------------------- #
def run_independent(indep_path):
    print("=" * 78)
    print("M1 variant-level test  —  exposure(beta_e) vs INDEPENDENT BMD GWAS")
    print("=" * 78)
    if not os.path.exists(EXPOSURE_HARMONISED):
        print("[FATAL] exposure harmonised file missing:", EXPOSURE_HARMONISED)
        sys.exit(2)
    if not os.path.exists(indep_path):
        print("[FATAL] independent BMD file missing:", indep_path)
        sys.exit(2)

    exp = load_exposure_beta_from_harmonised(EXPOSURE_HARMONISED)
    ind = load_independent_beta(indep_path)
    print(f"[load] exposure beta_e SNPs = {len(exp):,}")
    print(f"[load] independent BMD SNPs = {len(ind):,}")

    common = set(exp) & set(ind)
    print(f"[overlap] common SNPs = {len(common):,}")

    betaratios, zratios = [], []
    for snp in common:
        a1e, a2e, be, see = exp[snp]
        a1i, a2i, bi, sei = ind[snp]
        b_i, okb = harmonize_sign(bi, a1i, a2i, a1e, a2e)
        if not okb or abs(see) < 1e-12 or abs(sei) < 1e-12:
            continue
        z_e = be / see
        z_i = b_i / sei
        betaratios.append(b_i / be)
        zratios.append(z_i / z_e)

    n = len(betaratios)
    med_beta = statistics.median(betaratios) if betaratios else float("nan")
    med_z = statistics.median(zratios) if zratios else float("nan")

    # Adjudication:
    #   median(beta_indep / beta_exposure) ~ 1.0   -> 3.5426 is OVER-INFLATION (WRONG)
    #   median(beta_indep / beta_exposure) ~ 3.5426 -> reconstruction CORRECT
    tol = 0.25 * SCALE_FACTOR_CLAIMED
    if abs(med_beta - 1.0) <= 0.5:           # within 50% of 1 -> not 3.5x
        adj = "RECONSTRUCTION_OVER_INFLATED_WRONG"
        adj_txt = (
            f"median(beta_indep/beta_exp)={med_beta:.3f} ~ 1.0 => the deposited "
            "exposure beta was already on the true scale; the x3.5426 "
            "reconstruction is a 3.5x over-inflation and must be removed."
        )
    elif abs(med_beta - SCALE_FACTOR_CLAIMED) <= tol:
        adj = "RECONSTRUCTION_CONFIRMED_CORRECT"
        adj_txt = (
            f"median(beta_indep/beta_exp)={med_beta:.3f} ~ {SCALE_FACTOR_CLAIMED:.3f} "
            "=> the deposited exposure beta was compressed; x3.5426 is confirmed."
        )
    else:
        adj = "INCONCLUSIVE"
        adj_txt = (
            f"median(beta_indep/beta_exp)={med_beta:.3f} is neither ~1.0 nor "
            f"~{SCALE_FACTOR_CLAIMED:.3f}; sample/trait mismatch or residual "
            "scale issue -- investigate before concluding."
        )

    result = {
        "test": "exposure_beta_vs_independent_bmd_beta_ratio",
        "independent_file": indep_path,
        "scale_factor_claimed": SCALE_FACTOR_CLAIMED,
        "n_exposure": len(exp),
        "n_independent": len(ind),
        "n_common": len(common),
        "n_used": n,
        "median_beta_ratio": med_beta,
        "median_Z_ratio": med_z,
        "adjudication": adj,
        "adjudication_text": adj_txt,
        "can_adjudicate_M1": adj != "INCONCLUSIVE",
        "caveat": (
            "Trait/site mismatch (e.g. total-body vs femoral-neck BMD) adds "
            "genetic-correlation noise; interpret alongside the GEFOS aggregate "
            "OR benchmark (code/350)."
        ),
    }
    return result


# --------------------------------------------------------------------------- #
#  Optional download (graceful)
# --------------------------------------------------------------------------- #
def try_download():
    print("[try-download] attempting OpenGWAS ieu-a-980 (GEFOS femoral-neck BMD)...")
    try:
        import urllib.request
        url = "https://gwas.mrcieu.ac.uk/api/gwas/ieu-a-980"
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=25) as resp:
            body = resp.read(500)
        print("[try-download] response head:", body[:200])
        print("[try-download] NOTE: full fetch requires the 'ieugwasr' R package "
              "or 'gwas-summary' tool; automate separately where network allows.")
    except Exception as e:
        print(f"[try-download] UNREACHABLE in this environment: {e}")
        print("[try-download] Use these commands where network is available:")
        print("    # OpenGWAS / IEU (R): ieugwasr::associations('ieu-a-980') "
              "or gwas-summary ieu-a-980")
        print("    # GEFOS femoral neck (N=32735): ieu-a-980")
        print("    # GEFOS lumbar spine (N=28498): ieu-a-982")
        print("    # Total-body BMD (N=56284):     ebi-a-GCST005348")
        print("    # Real UKB total-body DXA (if distinct from Qian): GCST90568455")


# --------------------------------------------------------------------------- #
#  Main
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description="M1 variant-level scale test")
    ap.add_argument("--independent", help="path to a genuinely independent BMD GWAS (SNP,A1,A2,beta,se)")
    ap.add_argument("--try-download", action="store_true", help="attempt OpenGWAS fetch (needs network)")
    args = ap.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)
    results = {"scale_factor_claimed": SCALE_FACTOR_CLAIMED}

    if args.try_download:
        try_download()

    if args.independent:
        results["independent_test"] = run_independent(args.independent)
    else:
        results["qian_vs_ukb"] = run_qian_vs_ukb()

    with open(OUT_JSON, "w") as fh:
        json.dump(results, fh, indent=2)
    print("\n[output] JSON ->", OUT_JSON)

    # Human-readable headline
    print("\n" + "#" * 78)
    if "qian_vs_ukb" in results:
        r = results["qian_vs_ukb"]
        print("HEADLINE VERDICT:", r["verdict"])
        print(r["explanation"])
        print("=> M1 variant-level independent validation: NOT POSSIBLE with this pair.")
        print("   Correct path: obtain a genuinely independent BMD GWAS")
        print("   (ieu-a-980 / ieu-a-982 / ebi-a-GCST005348 / real GCST90568455)")
        print("   and re-run with --independent.")
    else:
        r = results["independent_test"]
        print("HEADLINE ADJUDICATION:", r["adjudication"])
        print(r["adjudication_text"])
    print("#" * 78)


if __name__ == "__main__":
    main()
