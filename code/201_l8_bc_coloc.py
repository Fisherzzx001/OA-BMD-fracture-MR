#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
201_l8_bc_coloc.py  --  PATH C, F4-style Bayesian colocalisation extended to osteophyte-dominant sites
=================================================================================================
Replicates code/29's coloc_abf engine (pure numpy; priors identical: p1=p2=1e-4, p12=5e-6,
sd_prior quant=0.15 / cc=0.2, window=500 kb) but re-implements the file I/O in stdlib only
(no pandas) and -- crucially -- WINDOW-FILTERS both GWAS files so only SNPs within +/-500 kb of
the 31 instrument loci are retained (memory ~100x lighter, no GC thrash).

For each of the 31 BMD_total instrument loci we colocalise:
    EXPOSURE : BMD_total  (sumstats_total.sumstats.gz, true N=31986, quant)  -- SAME exposure as the
               already-published F4 coloc (knee OA 0/31, fracture 9/31) so numbers are comparable.
    OUTCOME  : knee  OA (GCST90566800)   -> validation vs published 0/31
             : spine OA (GCST90566801)   -> osteophyte-dominant
             : hand  OA (GCST90566803)   -> osteophyte-dominant

Reports PP.H4 per locus and counts PP.H4 >= 0.8.

Inputs:
  BMD   = OA_FLAGSHIP/07_reference/ldsc/sumstats_total.sumstats.gz
  KNEE  = G:/OA_B_PROJECT_cache/l7/GCST90566800.tsv.gz
  SPINE = G:/OA_B_PROJECT_cache/l7/GCST90566801.tsv.gz
  HAND  = G:/OA_B_PROJECT_cache/l7/GCST90566803.tsv.gz
  HARM  = code/output/genetics/harmonised/BMD_total__OA_knee.csv  (GRCh37 positions)
Output:
  code/output/genetics/L8_BC_coloc.json
"""
import os, gzip, math, json, csv, sys
from collections import defaultdict
import numpy as np

HERE = "C:/Users/Administrator/WorkBuddy/2026-09-08-22-17-40/OA_B_PROJECT"
L7   = "G:/OA_B_PROJECT_cache/l7"
BMD  = "C:/Users/Administrator/WorkBuddy/2026-09-08-20-17-24/OA_FLAGSHIP/07_reference/ldsc/sumstats_total.sumstats.gz"
KNEE = os.path.join(L7, "GCST90566800.tsv.gz")
SPINE= os.path.join(L7, "GCST90566801.tsv.gz")
HAND = os.path.join(L7, "GCST90566803.tsv.gz")
GEN  = os.path.join(HERE, "code/output/genetics")
HARM = os.path.join(GEN, "harmonised", "BMD_total__OA_knee.csv")
OUT  = os.path.join(GEN, "L8_BC_coloc.json")
WINDOW = 500_000
BMD_N = 31986.0                      # GEFOS total-femur true N (code/29 BMD_TRUE_N["total"])

# ---- coloc engine (copied verbatim from code/29; pure numpy) ----
SD_PRIOR = {"quant": 0.15, "cc": 0.2}
def logsumexp(a):
    a = np.asarray(a, dtype=float); m = float(np.max(a))
    if not np.isfinite(m): return m
    return m + math.log(float(np.sum(np.exp(a - m))))
def logdiff(a, b):
    if not np.isfinite(a): return -np.inf
    if b >= a: return -np.inf
    return a + math.log1p(-math.exp(b - a))
def var_beta(n, f, typ, s):
    f = np.clip(f, 1e-6, 1 - 1e-6)
    v = 1.0 / (2.0 * np.asarray(n, dtype=float) * f * (1.0 - f))
    if typ == "cc": v = v / (s * (1.0 - s))
    return v
def log_abf(z, v, sd_prior):
    w = float(sd_prior) ** 2; r = w / (w + np.asarray(v, dtype=float))
    return 0.5 * np.log1p(-r) + 0.5 * r * (np.asarray(z, dtype=float) ** 2)
def coloc_abf(z1, v1, z2, v2, typ1="quant", typ2="cc", p1=1e-4, p2=1e-4, p12=5e-6, sd1=None, sd2=None):
    sd1 = SD_PRIOR[typ1] if sd1 is None else sd1
    sd2 = SD_PRIOR[typ2] if sd2 is None else sd2
    l1 = log_abf(z1, v1, sd1); l2 = log_abf(z2, v2, sd2)
    L1 = logsumexp(l1); L2 = logsumexp(l2); L4 = logsumexp(l1 + l2)
    L3 = logdiff(L1 + L2, L4)
    lh = np.array([0.0, math.log(p1)+L1, math.log(p2)+L2, math.log(p1)+math.log(p2)+L3, math.log(p12)+L4])
    lh = lh - logsumexp(lh)
    return dict(PP0=np.exp(lh[0]),PP1=np.exp(lh[1]),PP2=np.exp(lh[2]),PP3=np.exp(lh[3]),PP4=np.exp(lh[4]))

COMP = {"A":"T","T":"A","C":"G","G":"C"}
def align(a1g, a2g, a1e, a2e, ze):
    a1g, a2g, a1e, a2e = str(a1g), str(a2g), str(a1e), str(a2e)
    if {a1g, a2g} == {a1e, a2e}: return ze if a1g == a1e else -ze
    ca1, ca2 = COMP.get(a1e,""), COMP.get(a2e,"")
    if {a1g, a2g} == {ca1, ca2}: return ze if a1g == ca2 else -ze
    return None

def load_loci():
    loci = []
    with open(HARM) as f:
        for r in csv.DictReader(f):
            loci.append((int(r["CHR"]), int(r["BP"])))
    return loci   # 31 entries, index = order

def make_windows(loci):
    wins = defaultdict(list)   # chr -> list of (lo, hi, idx)
    for i,(ch,bp) in enumerate(loci):
        wins[ch].append((bp-WINDOW, bp+WINDOW, i))
    return wins

def load_oa_windowed(path, wins):
    out = defaultdict(dict)    # idx -> {(ch,bp): fields}
    n=0
    with gzip.open(path,"rt") as f:
        h=f.readline().rstrip("\n").split("\t"); idx={x:i for i,x in enumerate(h)}
        for line in f:
            n+=1
            if n % 3000000 == 0:
                print(f"  [{os.path.basename(path)}] {n//1000000}M lines scanned", flush=True)
            p=line.rstrip("\n").split("\t")
            try: ch=int(p[idx["chromosome"]]); bp=int(p[idx["base_pair_location"]])
            except: continue
            for (lo,hi,i) in wins.get(ch,[]):
                if lo<=bp<=hi:
                    try:
                        orr=float(p[idx["odds_ratio"]]); cu=float(p[idx["ci_upper"]]); cl=float(p[idx["ci_lower"]])
                    except: break
                    if orr<=0 or cu<=0 or cl<=0: break
                    se=(math.log(cu)-math.log(cl))/3.92
                    if se<=0 or not math.isfinite(se): break
                    out[i][(ch,bp)] = dict(ea=p[idx["effect_allele"]].upper(), oa=p[idx["other_allele"]].upper(),
                        z=math.log(orr)/se, eaf=float(p[idx["effect_allele_frequency"]]),
                        ncase=float(p[idx["NCASES"]]), nctrl=float(p[idx["NCONTROLS"]]))
                    break
    return out

def load_bmd_windowed(path, wins):
    out = defaultdict(dict)
    n=0
    with gzip.open(path,"rt") as f:
        h=f.readline().rstrip("\n").split("\t"); idx={x:i for i,x in enumerate(h)}
        for line in f:
            n+=1
            if n % 1000000 == 0:
                print(f"  [BMD sumstats] {n//1000000}M lines scanned", flush=True)
            p=line.rstrip("\n").split("\t")
            try: ch=int(p[idx["CHR"]]); bp=int(p[idx["BP"]])
            except: continue
            for (lo,hi,i) in wins.get(ch,[]):
                if lo<=bp<=hi:
                    try: maf=float(p[idx["MAF"]]); z=float(p[idx["Z"]])
                    except: break
                    out[i][(ch,bp)] = dict(a1=p[idx["A1"]].upper(), a2=p[idx["A2"]].upper(), maf=maf, z=z)
                    break
    return out

def coloc_locus(bmd_d, oa_d):
    common = set(bmd_d) & set(oa_d)
    if len(common) < 20: return None
    z1=[];v1=[];z2=[];v2=[]
    for k in common:
        b=bmd_d[k]; o=oa_d[k]
        za = align(b["a1"], b["a2"], o["ea"], o["oa"], o["z"])
        if za is None: continue
        z1.append(b["z"]); v1.append(var_beta(BMD_N, b["maf"], "quant", None))
        nt=o["ncase"]+o["nctrl"]; s=o["ncase"]/nt if nt>0 else 0.1
        z2.append(za); v2.append(var_beta(nt, o["eaf"], "cc", s))
    if len(z1) < 20: return None
    z1=np.array(z1);v1=np.array(v1);z2=np.array(z2);v2=np.array(v2)
    res = coloc_abf(z1,v1,z2,v2,"quant","cc")
    k=int(np.argmax(z1**2))
    concordant = (np.sign(z1[k])==np.sign(z2[k]))
    return dict(n_common=len(common), PP4=res["PP4"],
                top_concordant=bool(concordant))

def run_outcome(name, path, wins, loci):
    if not os.path.exists(path):
        print(f"[skip] {name}: missing"); return None
    print(f"[load] {name} OA GWAS (window-filtered) ...", flush=True)
    oa_all = load_oa_windowed(path, wins)
    bmd_all = load_bmd_windowed(BMD, wins)
    print(f"  {name}: extracted OA window SNPs at {sum(1 for i in range(len(loci)) if i in oa_all)}/{len(loci)} loci", flush=True)
    out_loci=[]
    for i,(ch,bp) in enumerate(loci):
        oa_d = oa_all.get(i, {})
        bmd_d = bmd_all.get(i, {})
        r = coloc_locus(bmd_d, oa_d)
        if r is None:
            out_loci.append(dict(rs=None, chr=ch, bp=bp, status="insufficient_overlap", n_common=0, PP4=None))
        else:
            out_loci.append(dict(rs=None, chr=ch, bp=bp, status="ok",
                                n_common=r["n_common"], PP4=r["PP4"], top_concordant=r["top_concordant"]))
    n_pp4 = sum(1 for x in out_loci if x.get("PP4") is not None and x["PP4"]>=0.8)
    n_ok  = sum(1 for x in out_loci if x.get("status")=="ok")
    print(f"[coloc] {name}: {n_pp4}/{n_ok} loci PP.H4>=0.8", flush=True)
    return dict(name=name, n_loci_analyzed=n_ok, n_PP4_ge_0p8=n_pp4, loci=out_loci)

def main():
    loci = load_loci()
    wins = make_windows(loci)
    print(f"[load] {len(loci)} instrument loci; window +/-{WINDOW//1000} kb", flush=True)
    results={}
    results["knee"] = run_outcome("knee", KNEE, wins, loci)   # validation vs published 0/31
    results["spine"]= run_outcome("spine",SPINE,wins,loci)
    results["hand"] = run_outcome("hand", HAND, wins, loci)
    json.dump(results, open(OUT,"w"), indent=2)
    print(f"\n[done] wrote {OUT}", flush=True)
    print("\n=== COLOC SUMMARY (PP.H4 >= 0.8) ===", flush=True)
    for k in ("knee","spine","hand"):
        r=results[k]
        if r is None: print(f"  {k}: n/a", flush=True); continue
        print(f"  {k:6s}: {r['n_PP4_ge_0p8']}/{r['n_loci_analyzed']} loci PP.H4>=0.8", flush=True)

if __name__=="__main__":
    main()
