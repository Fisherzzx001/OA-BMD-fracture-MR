#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
205_run_ldsc_munge_rg.py  --  downstream munge + cross-trait LDSC rg launcher
=============================================================================
One-command driver for the four Hatzikotoulas-2025 site-stratified OA GWAS:
code/204 converts each GWAS to munge_sumstats.py input, munge standardises it,
then official ldsc estimates single-trait h2 and cross-trait rg against BMD.

WEIGHTS MODES
-------------
  --weights standard   (DEFAULT, publication-grade configuration)
      Broad HapMap3 LD Scores, locally available at OA_FLAGSHIP/07_reference/ldsc:
        --ref-ld-chr  .../eur_w_ld_chr_   (22 files, 1,290,028 SNPs)
        --w-ld-chr    .../eur_w_ld_chr_   (same HapMap3 SNP set as regression weights)
      Munge additionally uses --merge-alleles w_hm3.snplist (SNP/A1/A2,
      1,217,311 SNPs), which restricts every OA sumstats to the HapMap3 set
      and orients alleles consistently. This is the standard Broad/ldsc setup.

  --weights self       (earlier exploratory configuration, kept as a contrast)
      --ref-ld / --w-ld both = the self-computed single merged file
      G:/REF/1000G_EUR/EUR_ldscore/EUR.l2.ldscore.gz (8,550,110 SNPs,
      --ld-wind-kb 1000). Using the same genome-wide file as both predictor
      and regression weight is NOT publication-grade; results are sensitivity
      analyses only. --fallback-M applies to this mode only.

-----------------------------------------------------------------------------
METHOD / HONEST BOUNDARIES  (read before interpreting ANY number)
-----------------------------------------------------------------------------
1. BMD EXPOSURE SUMSTATS: THIS SCRIPT KEEPS THE EXPLORATORY BASELINE.
   sumstats_total.sumstats.gz comes from docs/07; its N column is a
   PLACEHOLDER (400,000) and it was munged WITHOUT --merge-alleles. rg is
   invariant to a global N rescaling (N cancels in
   rg = b_gencov/sqrt(b_h2_1*b_h2_2)), so the cross-trait rg below is valid;
   single-trait h2 is NOT invariant and MUST NOT be quoted from this baseline.
   Both caveats were RESOLVED on 2026-11-04 by code/207_bmd_verified_ldsc.py:
   BMD was re-munged with the verified N=31,986 (GCST90568455, total-body DXA)
   and --merge-alleles w_hm3.snplist. That run gives BMD h2 = 0.3479 (0.0355)
   and leaves rg unchanged to <=0.0002 (N=31,873 sensitivity identical).
   Use code/207 outputs (ldsc_rg_verified/COMPARISON.md) for any citable
   number; keep this directory as the auditable before/after control.

2. OA GWAS IS A MULTI-COHORT META with per-SNP-varying NCASES/NCONTROLS
   (36,480 distinct case counts in the first 200k knee rows). We emit
   N = NCASES + NCONTROLS PER ROW and pass --n-min 1, because munge's
   default n_min = quantile(0.9)/1.5 would otherwise systematically drop
   SNPs covered only by smaller cohorts (6,471/58,169 in the smoke test).

3. Z IS DERIVED FROM P (munge's p_to_z), NOT FROM SE. The GWAS
   standard_error column is "#NA" throughout -- harmless. OR is used only
   as the signed sumstat (null value 1) to set the sign of Z.

4. rsID HARMONISATION. OA GWAS use position-encoded variant IDs
   ("9_15867726_G_GC"); code/204 maps (CHR,BP) -> rsID via 1000G EUR.bim
   (GRCh37, verified 99.1% rsID overlap with the BMD sumstats). Palindromic
   and non-biallelic variants are dropped.

5. --merge-alleles writes ALL 1,217,311 HapMap3 rows with Z=NA for SNPs
   absent from the GWAS; ldsc drops those automatically (verified: a 9,623
   non-missing SNP smoke sumstats is read as exactly 9,623 SNPs).

6. --fallback-M (self mode only): if EUR.l2.M_5_50 is empty, run with an
   explicit --M 8550110 --not-M-5-50. h2 is unreliable there, rg ~ invariant.
   Standard mode has per-chromosome .l2.M_5_50 files and needs no fallback.

7. This is a POST-SUBMISSION / EXPLORATORY layer. It does NOT alter
   docs/08 (the submission-ready manuscript, OAC body 3999/4000) and does
   NOT change any reported beta/SE/p/OR/CI. A submission-ready methods
   paragraph and result table (candidate Table S29) are held in
   docs/Supplementary/Table_S30_ldsc_site_rg.md (S30 in docs/08; formerly S29)
   because the OAC body has only 1 word of headroom. Inserting it requires
   equal-word compression plus a re-run of gates 148/149b/182.

-----------------------------------------------------------------------------
USAGE
-----
  python3 code/205_run_ldsc_munge_rg.py --check
        Verify every prerequisite; print readiness. Does not run anything.

  python3 code/205_run_ldsc_munge_rg.py --smoke
        Small 200k-row end-to-end run (pipeline validation, ~3 min).

  python3 code/205_run_ldsc_munge_rg.py
        FULL run, standard HapMap3 weights: prep 4 OA GWAS -> munge (with
        --merge-alleles) -> h2 x5 -> BMD x OA rg (4 pairs) -> OA site-site
        rg (3 pairs). ~40-80 min (dominated by scanning 3.8 GB of GWAS gz).

  python3 code/205_run_ldsc_munge_rg.py --weights self
        Contrast run with the self-computed whole-genome EUR LD Scores.

  python3 code/205_run_ldsc_munge_rg.py --sites knee spine --no-h2
        Subset of OA sites, skip single-trait h2.

Outputs: code/output/genetics/ldsc_rg/  (logs, .sumstats.gz, .log, summary.json)
"""
import argparse
import json
import os
import subprocess
import sys
import time

HERE = "C:/Users/Administrator/WorkBuddy/2026-09-08-22-17-40/OA_B_PROJECT"
PYBIN = "C:/Users/Administrator/.workbuddy/binaries/python/versions/3.13.12/Scripts"
LDSC = os.path.join(PYBIN, "ldsc.py")
MUNGE = os.path.join(PYBIN, "munge_sumstats.py")
PREP = os.path.join(HERE, "code", "204_prep_oa_gwas_for_munge.py")

# --- standard Broad HapMap3 weights (publication-grade) ---
STD_DIR = ("C:/Users/Administrator/WorkBuddy/2026-09-08-20-17-24/"
           "OA_FLAGSHIP/07_reference/ldsc")
STD_PREFIX = STD_DIR + "/eur_w_ld_chr_"      # ldsc appends 1..22
STD_M = STD_DIR + "/eur_w_ld_chr.l2.M"        # 1,290,050
STD_M_5_50 = STD_DIR + "/eur_w_ld_chr.l2.M_5_50"  # 1,173,569
HM3_SNPLIST = STD_DIR + "/w_hm3.snplist"      # SNP A1 A2, 1,217,311 rows

# --- self-computed whole-genome EUR LD Scores (exploratory contrast) ---
REF_PREFIX = "G:/REF/1000G_EUR/EUR_ldscore/EUR"
REF_L2 = REF_PREFIX + ".l2.ldscore.gz"
REF_M = REF_PREFIX + ".l2.M"
REF_M_5_50 = REF_PREFIX + ".l2.M_5_50"
M_TOTAL_FALLBACK = 8550110   # exact SNP count in EUR.l2.ldscore.gz (verified)

BMD_SUMSTATS = STD_DIR + "/sumstats_total.sumstats.gz"

MUNGE_INPUT = os.path.join(HERE, "code", "output", "genetics", "ldsc_munge_input")
OUTDIR = os.path.join(HERE, "code", "output", "genetics", "ldsc_rg")
SITES = ["hip", "knee", "spine", "hand"]
GWAS_ACC = {"hip": "798", "knee": "800", "spine": "801", "hand": "803"}


def sh(cmd, logpath, cwd=HERE):
    """Run a command, stream output to a log file, return (rc, logpath)."""
    os.makedirs(os.path.dirname(logpath), exist_ok=True)
    with open(logpath, "w") as lf:
        lf.write("CMD: " + " ".join(cmd) + "\n\n")
        lf.flush()
        p = subprocess.run(cmd, cwd=cwd, stdout=lf, stderr=subprocess.STDOUT)
    return p.returncode, logpath


def size(p):
    return os.path.getsize(p) if os.path.exists(p) else 0


def check(weights):
    print("=== 205 · prerequisite check ===", flush=True)
    print(f"  weights mode              : {weights}"
          f" ({'HapMap3, publication-grade' if weights == 'standard' else 'self-computed, exploratory'})",
          flush=True)
    items = [
        ("ldsc.py", LDSC), ("munge_sumstats.py", MUNGE),
        ("204 prep script", PREP), ("BMD sumstats", BMD_SUMSTATS),
    ]
    for s in SITES:
        items.append((f"OA {s} GWAS",
                      f"G:/OA_B_PROJECT_cache/l7/GCST90566{GWAS_ACC[s]}.tsv.gz"))
    if weights == "standard":
        items += [("HapMap3 w_hm3.snplist", HM3_SNPLIST),
                  ("HapMap3 .l2.M", STD_M), ("HapMap3 .l2.M_5_50", STD_M_5_50)]
    else:
        items += [("EUR LD score", REF_L2),
                  ("EUR .l2.M", REF_M), ("EUR .l2.M_5_50", REF_M_5_50)]

    all_ok = True
    for name, path in items:
        ok = os.path.exists(path) and size(path) > 0
        if name in ("EUR .l2.M", "EUR .l2.M_5_50") and not ok:
            pass  # reported, but self mode has --fallback-M
        elif not ok:
            all_ok = False
        print(f"  [{'OK ' if ok else 'MISS'}] {name:22s} {size(path)/1e6:9.1f} MB  {path}",
              flush=True)

    if weights == "standard":
        l2 = [f"{STD_PREFIX}{i}.l2.ldscore.gz" for i in range(1, 23)]
        m5 = [f"{STD_PREFIX}{i}.l2.M_5_50" for i in range(1, 23)]
        n_l2 = sum(1 for p in l2 if size(p) > 0)
        n_m5 = sum(1 for p in m5 if size(p) > 0)
        print(f"  [{'OK ' if n_l2 == 22 else 'MISS'}] HapMap3 per-chr ldscore : {n_l2}/22",
              flush=True)
        print(f"  [{'OK ' if n_m5 == 22 else 'MISS'}] HapMap3 per-chr M_5_50 : {n_m5}/22",
              flush=True)
        weights_ready = (n_l2 == 22 and n_m5 == 22 and all_ok)
    else:
        weights_ready = size(REF_M) > 0 and size(REF_M_5_50) > 0
        print(f"\n  EUR M files ready       : {'YES' if weights_ready else 'NO'}", flush=True)
        if not weights_ready:
            print("     -> wait for code/203, or pass --fallback-M (self mode only).", flush=True)

    print(f"\n  weight files ready      : {'YES' if weights_ready else 'NO'}", flush=True)
    print(f"  all other inputs ready  : {'YES' if all_ok else 'NO'}", flush=True)
    print(f"  READY TO RUN            : {'YES' if (all_ok and weights_ready) else 'NO'}",
          flush=True)
    return all_ok, weights_ready


def run_prep(sites, limit):
    cmd = [sys.executable, "-u", PREP, "--sites"] + sites
    if limit:
        cmd += ["--limit", str(limit)]
    rc, log = sh(cmd, os.path.join(OUTDIR, "logs", "01_prep.log"))
    print(f"[prep] rc={rc} log={log}", flush=True)
    return rc


def run_munge(site, weights):
    src = os.path.join(MUNGE_INPUT, f"OA_{site}.munge_input.tsv.gz")
    out = os.path.join(OUTDIR, "munged", f"OA_{site}")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    cmd = [sys.executable, "-u", MUNGE,
           "--sumstats", src, "--out", out,
           "--signed-sumstats", "OR,1",
           "--a1", "A1", "--a2", "A2", "--snp", "SNP", "--p", "P", "--N-col", "N",
           "--n-min", "1", "--chunksize", "200000"]
    if weights == "standard":
        cmd += ["--merge-alleles", HM3_SNPLIST]
    rc, log = sh(cmd, os.path.join(OUTDIR, "logs", f"02_munge_{site}.log"))
    print(f"[munge {site}] rc={rc} log={log}", flush=True)
    return rc, out + ".sumstats.gz"


def ldsc_base(weights, ready, fallback):
    if weights == "standard":
        return ["--ref-ld-chr", STD_PREFIX, "--w-ld-chr", STD_PREFIX]
    base = ["--ref-ld", REF_PREFIX, "--w-ld", REF_PREFIX]
    if not ready and fallback:
        print("WARNING: M files empty -> using --M %d --not-M-5-50 "
              "(h2 unreliable; rg ~invariant)" % M_TOTAL_FALLBACK, flush=True)
        base += ["--M", str(M_TOTAL_FALLBACK), "--not-M-5-50"]
    elif not ready:
        raise RuntimeError("M files not ready; rerun code/203 or pass --fallback-M")
    return base


def run_h2(label, sumstats, base, tag):
    out = os.path.join(OUTDIR, "h2", tag)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    cmd = [sys.executable, "-u", LDSC, "--h2", sumstats] + base + ["--out", out]
    rc, log = sh(cmd, os.path.join(OUTDIR, "logs", f"03_h2_{tag}.log"))
    print(f"[h2 {label}] rc={rc} log={log}", flush=True)
    return rc


def run_rg(files, base, tag, note=""):
    """files: ordered list; ldsc pairs files[0] with each of files[1:]."""
    out = os.path.join(OUTDIR, "rg", tag)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    cmd = [sys.executable, "-u", LDSC, "--rg", ",".join(files)] + base + ["--out", out]
    rc, log = sh(cmd, os.path.join(OUTDIR, "logs", f"04_rg_{tag}.log"))
    print(f"[rg {tag}] rc={rc} {note} log={log}", flush=True)
    return rc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="verify prerequisites only")
    ap.add_argument("--smoke", action="store_true",
                    help="200k-row knee subset end-to-end (pipeline validation)")
    ap.add_argument("--weights", choices=["standard", "self"], default="standard",
                    help="LD Score reference / regression weights (default: standard HapMap3)")
    ap.add_argument("--sites", nargs="+", choices=SITES, default=None)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--no-h2", action="store_true")
    ap.add_argument("--fallback-M", action="store_true",
                    help="self mode only: run even if .l2.M_5_50 is empty")
    ap.add_argument("--skip-prep", action="store_true",
                    help="reuse existing munge_input files")
    args = ap.parse_args()

    os.makedirs(OUTDIR, exist_ok=True)
    inputs_ok, weights_ready = check(args.weights)
    if args.check:
        return 0 if (inputs_ok and weights_ready) else 1

    if not inputs_ok:
        print("ERROR: missing inputs (see above). Nothing run.", file=sys.stderr)
        return 2
    if args.weights == "standard" and not weights_ready:
        print("ERROR: HapMap3 weight files incomplete. Nothing run.", file=sys.stderr)
        return 2
    if args.weights == "self" and not weights_ready and not args.fallback_M:
        print("ERROR: EUR .l2.M / .l2.M_5_50 empty. Wait for code/203 or pass "
              "--fallback-M (self mode only). Nothing run.", file=sys.stderr)
        return 2
    if args.weights == "standard" and args.fallback_M:
        print("NOTE: --fallback-M is ignored in standard mode (per-chromosome "
              "M_5_50 files are present).", flush=True)

    sites = ["knee"] if args.smoke else (args.sites or SITES)
    limit = 200000 if args.smoke else args.limit
    base = ldsc_base(args.weights, weights_ready, args.fallback_M)
    summary = {
        "mode": "smoke" if args.smoke else "full",
        "weights_mode": args.weights,
        "ref_ld": (STD_PREFIX + "[1-22]") if args.weights == "standard" else REF_L2,
        "w_ld": (STD_PREFIX + "[1-22]") if args.weights == "standard" else REF_L2,
        "merge_alleles": HM3_SNPLIST if args.weights == "standard" else None,
        "sites": sites, "limit": limit,
        "weight_files_ready": weights_ready,
        "fallback_M_used": (args.weights == "self" and not weights_ready and args.fallback_M),
        "boundaries": [
            "BMD N=400000 is a placeholder: h2 NOT quotable; rg invariant to global N",
            "BMD sumstats was munged WITHOUT --merge-alleles (asymmetry, see docstring)",
            "RESOLVED by code/207: verified N=31986 + merge-alleles; see ldsc_rg_verified/",
            "results are exploratory / post-submission; docs/08 untouched",
        ],
        "started": time.strftime("%Y-%m-%d %H:%M:%S"), "steps": [],
    }
    t0 = time.time()

    if not args.skip_prep:
        rc = run_prep(sites, limit)
        summary["steps"].append({"step": "prep", "rc": rc})
        if rc != 0:
            print("ERROR: prep failed", file=sys.stderr)
            return 3
    else:
        print("[prep] skipped (--skip-prep)", flush=True)

    munged = {}
    for s in sites:
        rc, f = run_munge(s, args.weights)
        summary["steps"].append({"step": f"munge_{s}", "rc": rc, "file": f})
        if rc != 0 or not os.path.exists(f):
            print(f"ERROR: munge failed for {s}", file=sys.stderr)
            return 4
        munged[s] = f

    if not args.no_h2:
        rc = run_h2("BMD_total", BMD_SUMSTATS, base, "BMD_total")
        summary["steps"].append({"step": "h2_BMD_total", "rc": rc})
        for s in sites:
            rc = run_h2(s, munged[s], base, f"OA_{s}")
            summary["steps"].append({"step": f"h2_OA_{s}", "rc": rc})

    # Primary: BMD (p1) against each OA site -> one ldsc call, N-1 rg pairs
    files = [BMD_SUMSTATS] + [munged[s] for s in sites]
    rc = run_rg(files, base, "BMD_vs_OA_sites",
                note=f"({len(sites)} BMD x OA pairs)")
    summary["steps"].append({"step": "rg_BMD_vs_OA", "rc": rc})

    # Secondary: OA site vs site (knee = cartilage-dominant reference)
    if len(sites) >= 2:
        ref_site = "knee" if "knee" in sites else sites[0]
        others = [s for s in sites if s != ref_site]
        files = [munged[ref_site]] + [munged[s] for s in others]
        rc = run_rg(files, base, "OA_site_vs_site",
                    note=f"({len(others)} site-site pairs, p1={ref_site})")
        summary["steps"].append({"step": "rg_OA_site_vs_site", "rc": rc})

    summary["elapsed_sec"] = round(time.time() - t0, 1)
    summary["outdir"] = OUTDIR
    sp = os.path.join(OUTDIR, "summary.json")
    with open(sp, "w") as jf:
        json.dump(summary, jf, indent=2)
    print(f"\n[done] {summary['elapsed_sec']}s  summary={sp}", flush=True)
    if args.weights == "standard":
        print("REMINDER: standard HapMap3 weights used, but this BASELINE keeps "
              "BMD N=400000 (placeholder) -> h2 NOT quotable. For citable numbers "
              "run code/207_bmd_verified_ldsc.py (verified N=31986 + merge-alleles); "
              "see ldsc_rg_verified/COMPARISON.md.", flush=True)
    else:
        print("REMINDER: --weights self is EXPLORATORY ONLY (non-standard --w-ld); "
              "BMD N is a placeholder. See boundary notes in this script.", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
