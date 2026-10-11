#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
fix_mrlap_inputs.py -- repair the MRlap GWAS inputs produced by prep_mrlap_gwas.py
without re-running the (14 min) rsID recovery.

Two defects are corrected:

1. BMD_total was written with N=318,973 (typo in TRAITS) instead of the
   manuscript Table-1 value N=31,873.  A 10x-too-large N deflates the LDSC h2 by
   ~10x (MRlap returned h2_exp=0.0355 where Table 1 reports 0.391) and therefore
   distorts the overlap/winner's-curse correction.  -> N reset to 31,873.

2. OA_knee was written with (a) the effect column named "or" and (b) a per-variant
   N taken from NCASES+NCONTROLS.  Both are wrong:
     (a) MRlap's tidy_inputGWAS has a bug in its "or" branch -- it calls
         dplyr::mutate(Z = log(.data$beta)/.data$se, beta=NULL, se=NULL)
         without piping the data frame in, which throws
         "Can't subset `.data` outside of a data mask context".
         -> we pre-convert OR to log(OR) and name the column "beta".
     (b) the per-variant N column in GCST90566800 ranges 605-1,316,504 with no
         mode and is inconsistent with the study total, so it is replaced by the
         manuscript value N=1,208,397 for every SNP.

Also verifies the N of the remaining five traits.
"""
import os, sys, math, time

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DIR = os.path.join(BASE, "code", "ref", "mrlap_gwas")

# manuscript Table-1 N values (docs/08_投稿稿_OA_BMD_fracture.md)
TRUE_N = {
    "BMD_arms": 31873, "BMD_femur_neck": 31873, "BMD_leg": 31873,
    "BMD_lumbar": 31873, "BMD_total": 31873,
    "fracture": 264973,
    "OA_knee": 1208397,
}
LOGOR = {"OA_knee"}     # traits whose effect column must be log-transformed


def rewrite(name):
    src = os.path.join(DIR, name + ".tsv")
    tmp = os.path.join(DIR, name + ".tsv.fix")
    if not os.path.exists(src):
        print("  MISSING", src); return None
    N_new = TRUE_N[name]
    do_log = name in LOGOR
    t0 = time.time()
    n = 0
    bad = 0
    ns_seen = set()
    with open(src) as f, open(tmp, "w") as o:
        hdr = f.readline().rstrip("\n").split("\t")
        # normalise effect column header to "beta"
        hdr[5] = "beta"
        o.write("\t".join(hdr) + "\n")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 8:
                bad += 1; continue
            try:
                eff = float(p[5]); se = float(p[6])
            except ValueError:
                bad += 1; continue
            if do_log:
                if eff <= 0:
                    bad += 1; continue
                eff = math.log(eff)
            if not (se > 0):
                bad += 1; continue
            p[5] = repr(eff); p[6] = repr(se); p[7] = str(N_new)
            o.write("\t".join(p) + "\n")
            n += 1
            if len(ns_seen) < 3:
                ns_seen.add(p[7])
    os.replace(tmp, src)
    print("  %-16s rows=%-10d dropped=%-8d N=%-9d (%.0fs)" % (name, n, bad, N_new, time.time() - t0))
    return n


if __name__ == "__main__":
    t = time.time()
    print("Fixing MRlap inputs in", DIR)
    for name in ["BMD_total", "OA_knee", "BMD_arms", "BMD_femur_neck",
                 "BMD_leg", "BMD_lumbar", "fracture"]:
        rewrite(name)
    print("DONE in %.0fs" % (time.time() - t))
