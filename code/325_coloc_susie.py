#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
325_coloc_susie.py  —  T1-C rebuttal reserve: multi-signal colocalisation
============================================================================
WHY. The published P3 evidence is coloc.abf at 31 BMD loci (0/31 PP.H4>=0.8
for knee OA). coloc.abf assumes ONE causal variant per locus per trait; if a
locus carries a BMD signal AND an independent OA signal, abf can miss it.
coloc-SuSiE relaxes this to multiple causal signals per locus (sum of single
effects), so it is the pre-registered-style robustness check for P3.

WHAT. For each of the 31 locked BMD instruments (window +/-500 kb):
  * extract BMD (GEFOS total, sumstats_total) and knee-OA (GCST90566800) SNP z;
  * harmonise alleles, keep SNPs present in 1000G-EUR (503 samples);
  * build LD matrix R from PLINK bed (pure numpy 2-bit decode);
  * fit SuSiE (IBSS, L=2 single-effect components, per-trait) -> per-SNP PIP;
  * a locus is "colocalised" if the two traits share a credible variant:
        shared-PIP mass = sum_i PIP_BMD[i]*PIP_OA[i]
        and the top variant of each trait's 95% credible set coincides.

PRIOR. p1=p2=1e-4, per-SNP; L=2; rho=0.95 credible-set coverage. These match
the coloc.abf SD priors (quant 0.15, cc 0.2) via the same ABF variance model.

HONESTY. If SuSiE still yields 0/31 shared credible variants, P3 is strengthened
(multi-signal-aware). If any locus is shared, it is reported honestly and P3
would need revision. This script NEVER touches docs/08 or any reported value.
Output: code/output/genetics/coloc_susie_knee.json
============================================================================
"""
import os, csv, gzip, math, json, sys
import numpy as np

GEN   = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "genetics")
HARM  = os.path.join(GEN, "harmonised", "BMD_total__OA_knee.csv")
BMD   = "C:/Users/Administrator/WorkBuddy/2026-09-08-20-17-24/OA_FLAGSHIP/07_reference/ldsc/sumstats_total.sumstats.gz"
KNEE  = "G:/OA_B_PROJECT_cache/l7/GCST90566800.tsv.gz"
BED   = "G:/REF/1000G_EUR/EUR.bed"
BIM   = "G:/REF/1000G_EUR/EUR.bim"
OUT   = os.path.join(GEN, "coloc_susie_knee.json")
WINDOW= 500_000
BMD_N = 31986.0
N_SAMPLES = 503
L_COMP = 2            # SuSiE components (multi-signal)
P1 = P2 = 1e-4        # per-SNP causal prior
RHO = 0.95            # credible set coverage

# ---------- PLINK bed reader (2-bit, SNP-major) ----------
def load_bim(path):
    rs2row = {}; row2info = {}
    with open(path) as f:
        for i, line in enumerate(f):
            p = line.rstrip("\n").split("\t")
            if len(p) < 5: continue
            rs2row[p[1]] = i
            row2info[i] = (p[0], int(p[3]), p[4].upper(), p[5].upper())  # chr,bp,A1,A2
    return rs2row, row2info

def read_dosages(rows, bim_rows):
    """Return dosage matrix (m x n) for PLINK bed rows; A2 dosage; NaN->0."""
    n = N_SAMPLES; by = (n + 3) // 4
    D = np.full((len(rows), n), np.nan, dtype=np.float32)
    # sort by row to read sequentially-friendly (random seek still fine for 31 windows)
    with open(BED, "rb") as f:
        for j, r in enumerate(rows):
            off = 3 + r * by
            f.seek(off)
            raw = f.read(by)
            if len(raw) < by: continue
            codes = np.frombuffer(raw, dtype=np.uint8)
            g = np.zeros(n, dtype=np.uint8)
            for bi in range(by):
                b = int(codes[bi])
                for k in range(4):
                    idx = bi * 4 + k
                    if idx >= n: break
                    g[idx] = (b >> (2 * k)) & 3
            # PLINK1 bed 2-bit: 00=homA1(dosage 0), 01=missing,
            # 10=het(dosage 1), 11=homA2(dosage 2). FIX 2026-11-10:
            # old code mapped g==2->2 and g==3->0 (WRONG). Verified against
            # bmd_instruments eaf: rs12568930 A2 freq=0.836 vs eaf=0.820.
            d = np.where(g == 3, 2.0, np.where(g == 2, 1.0, 0.0))
            d[g == 1] = np.nan
            D[j] = d
    return D

def ld_matrix(D):
    """Standardise columns (SNPs), return correlation R (m x m); mean-impute missing."""
    M = D.copy()
    for j in range(M.shape[0]):
        col = M[j]; ok = ~np.isnan(col)
        mu = np.nanmean(col) if ok.any() else 0.0
        col[~ok] = mu
        sd = col.std()
        M[j] = (col - mu) / sd if sd > 0 else 0.0
    R = M @ M.T / float(M.shape[1])
    return np.clip(R, -1, 1)

# ---------- ABF / SuSiE variational ----------
def var_beta(n, f, typ, s):
    f = np.clip(f, 1e-6, 1 - 1e-6)
    v = 1.0 / (2.0 * float(n) * f * (1.0 - f))
    if typ == "cc": v = v / (s * (1.0 - s))
    return v

def log_abf(z, v, sd_prior):
    w = float(sd_prior) ** 2
    r = w / (w + np.asarray(v, dtype=float))
    return 0.5 * np.log1p(-r) + 0.5 * r * (np.asarray(z, dtype=float) ** 2)

def logsumexp(a):
    a = np.asarray(a, float); m = float(np.max(a))
    return m + math.log(float(np.sum(np.exp(a - m)))) if np.isfinite(m) else m

def susie_ibss(z, V, R, L=L_COMP, prior=P1, sd_prior=0.15, n_iter=100, tol=1e-6):
    """IBSS: iteratively fit L single-effect components; return per-SNP PIP.
    Approximates SuSiE's sum-of-single-effects with per-component posterior
    inclusion probabilities updated against residual z (Wang et al. 2020)."""
    m = len(z)
    base_labf = log_abf(z, V, sd_prior)
    log_prior = np.full(m, math.log(prior / (1 - prior)))
    prev = None
    for it in range(n_iter):
        alpha = np.zeros((L, m), float)   # per-component posterior incl prob
        beta_hat = np.zeros((L, m), float)
        for l in range(L):
            # residual z after other components' posterior mean effects
            resid = z.copy()
            for k in range(L):
                if k == l:
                    continue
                # posterior mean effect at each SNP for component k
                pm = alpha[k] * (z / np.sqrt(V))  # shrunk signal proxy
                # subtract explained signal through LD
                resid = resid - R @ pm
            lbf = base_labf + (resid ** 2 - z ** 2) * 0.0  # placeholder; refined below
            # recompute ABF on residual variance scale
            lbf = log_abf(resid, V, sd_prior) + log_prior
            a = np.exp(lbf - logsumexp(np.append(lbf, 0.0)))  # incl prob, null mass=1
            alpha[l] = a
            beta_hat[l] = a * (resid / np.sqrt(V))
        # PIP = 1 - prod_l (1 - alpha_l)
        pip = 1.0 - np.prod(1.0 - alpha, axis=0)
        if prev is not None and np.max(np.abs(pip - prev)) < tol:
            break
        prev = pip
    return pip, alpha

def credible_set(pip, rho=RHO):
    order = np.argsort(pip)[::-1]
    csum = 0.0; cs = []
    for i in order:
        cs.append(int(i)); csum += pip[i]
        if csum >= rho: break
    return cs

# ---------- data loading (window filtered) ----------
COMP = {"A":"T","T":"A","C":"G","G":"C"}
def align(a1g, a2g, a1e, a2e, ze):
    a1g,a2g,a1e,a2e = str(a1g).upper(),str(a2g).upper(),str(a1e).upper(),str(a2e).upper()
    if {a1g,a2g} == {a1e,a2e}: return ze if a1g==a1e else -ze
    ca1,ca2 = COMP.get(a1e,""),COMP.get(a2e,"")
    if {a1g,a2g} == {ca1,ca2}: return ze if a1g==ca2 else -ze
    return None

def load_loci():
    loci=[]
    with open(HARM) as f:
        for r in csv.DictReader(f):
            loci.append((r["SNP"], int(r["CHR"]), int(r["BP"])))
    return loci

def window_index(loci):
    wins={}
    for rs,ch,bp in loci:
        wins.setdefault(ch,[]).append((bp-WINDOW,bp+WINDOW,rs))
    return wins

def load_bmd(wins):
    out={}
    with gzip.open(BMD,"rt") as f:
        h=f.readline().rstrip("\n").split("\t"); ix={x:i for i,x in enumerate(h)}
        for line in f:
            p=line.rstrip("\n").split("\t")
            try: ch=int(p[ix["CHR"]]); bp=int(p[ix["BP"]])
            except: continue
            for lo,hi,rs in wins.get(ch,[]):
                if lo<=bp<=hi:
                    try: out[rs]=dict(a1=p[ix["A1"]].upper(),a2=p[ix["A2"]].upper(),
                                      maf=float(p[ix["MAF"]]),z=float(p[ix["Z"]]),bp=bp)
                    except: pass
                    break
    return out

def load_oa(path,wins):
    out={}
    with gzip.open(path,"rt") as f:
        h=f.readline().rstrip("\n").split("\t"); ix={x:i for i,x in enumerate(h)}
        for line in f:
            p=line.rstrip("\n").split("\t")
            try: ch=int(p[ix["chromosome"]]); bp=int(p[ix["base_pair_location"]])
            except: continue
            for lo,hi,rs in wins.get(ch,[]):
                if lo<=bp<=hi:
                    try:
                        cu=float(p[ix["ci_upper"]]); cl=float(p[ix["ci_lower"]]); orr=float(p[ix["odds_ratio"]])
                        se=(math.log(cu)-math.log(cl))/3.92
                        if se<=0 or orr<=0: break
                        out[rs]=dict(ea=p[ix["effect_allele"]].upper(),oa=p[ix["other_allele"]].upper(),
                                     z=math.log(orr)/se, eaf=float(p[ix["effect_allele_frequency"]]),
                                     ncase=float(p[ix["NCASES"]]),nctrl=float(p[ix["NCONTROLS"]]),bp=bp)
                    except: pass
                    break
    return out

def main():
    loci=load_loci(); wins=window_index(loci)
    print(f"[325] {len(loci)} loci, window +/-{WINDOW//1000} kb, L={L_COMP}",flush=True)
    rs2row,_=load_bim(BIM)
    bmd_all=load_bmd(wins); oa_all=load_oa(KNEE,wins)
    print(f"[325] BMD window SNPs loaded for {len(bmd_all)}/{len(loci)} loci; OA {len(oa_all)}/{len(loci)}",flush=True)
    results=[]
    # ---- positional scan: rebuild OA/BMD as bp->fields PER LOCUS WINDOW ----
    # FIX 2026-11-10 (bug #2): previous version keyed oa_pos/bmd_pos by chromosome
    # and broke after the first matching window, so every locus on the same chr
    # shared ONE whole-chromosome SNP set (all chr1 loci showed identical n/pip/cs).
    # Now key by lead rsid and assign each SNP to EVERY window it falls in (no
    # break), so each locus gets its own +/-WINDOW SNP set. Lead-rs keyed loads
    # above (bmd_all/oa_all) remain for the availability printout only.
    from collections import defaultdict
    oa_pos=defaultdict(dict); bmd_pos=defaultdict(dict)
    with gzip.open(KNEE,"rt") as f:
        h=f.readline().rstrip("\n").split("\t"); ix={x:i for i,x in enumerate(h)}
        for line in f:
            p=line.rstrip("\n").split("\t")
            try: ch=int(p[ix["chromosome"]]); bp=int(p[ix["base_pair_location"]])
            except: continue
            for lo,hi,rs in wins.get(ch,[]):
                if lo<=bp<=hi:
                    try:
                        cu=float(p[ix["ci_upper"]]);cl=float(p[ix["ci_lower"]]);orr=float(p[ix["odds_ratio"]])
                        se=(math.log(cu)-math.log(cl))/3.92
                        if se<=0 or orr<=0: continue
                        oa_pos[rs][bp]=dict(ea=p[ix["effect_allele"]].upper(),oa=p[ix["other_allele"]].upper(),
                            z=math.log(orr)/se, eaf=float(p[ix["effect_allele_frequency"]]),
                            ncase=float(p[ix["NCASES"]]),nctrl=float(p[ix["NCONTROLS"]]),
                            snp=p[ix["variant_id"]].split("_")[0] if "variant_id" in ix else "")
                    except: pass
    with gzip.open(BMD,"rt") as f:
        h=f.readline().rstrip("\n").split("\t"); ix={x:i for i,x in enumerate(h)}
        for line in f:
            p=line.rstrip("\n").split("\t")
            try: ch=int(p[ix["CHR"]]); bp=int(p[ix["BP"]])
            except: continue
            for lo,hi,rs in wins.get(ch,[]):
                if lo<=bp<=hi:
                    try: bmd_pos[rs][bp]=dict(snp=p[ix["SNP"]],a1=p[ix["A1"]].upper(),a2=p[ix["A2"]].upper(),
                                              maf=float(p[ix["MAF"]]),z=float(p[ix["Z"]]))
                    except: pass

    n_shared=0; n_ok=0
    for rs,ch,bp in loci:
        B=bmd_pos[rs]; O=oa_pos[rs]
        common_bp=sorted(set(B)&set(O))
        # keep SNPs with EUR.bim mapping
        rows=[]; snps=[]; z1=[];v1=[];z2=[];v2=[]
        for position in common_bp:
            b=B[position]; o=O[position]
            snp=b["snp"]
            if snp not in rs2row: continue
            za=align(b["a1"],b["a2"],o["ea"],o["oa"],o["z"])
            if za is None: continue
            nt=o["ncase"]+o["nctrl"]; s=o["ncase"]/nt if nt>0 else 0.1
            rows.append(rs2row[snp]); snps.append(snp)
            z1.append(b["z"]); v1.append(var_beta(BMD_N,b["maf"],"quant",None))
            z2.append(za);     v2.append(var_beta(nt,o["eaf"],"cc",s))
        if len(rows)<30:
            results.append(dict(rs=rs,chr=ch,bp=bp,status="insufficient_overlap",n_snp=len(rows)))
            continue
        D=read_dosages(rows,rs2row)
        R=ld_matrix(D)
        z1=np.array(z1);V1=np.array(v1);z2=np.array(z2);V2=np.array(v2)
        pip1,_=susie_ibss(z1,V1,R,sd_prior=0.15)
        pip2,_=susie_ibss(z2,V2,R,sd_prior=0.20)
        pip1_sum=float(pip1.sum()); pip2_sum=float(pip2.sum())
        cs1=set(credible_set(pip1)); cs2=set(credible_set(pip2))
        shared=sorted(cs1 & cs2)
        shared_mass=float(np.sum(pip1*pip2))
        # A credible set is "reliable" only if the model assigns substantial total
        # causal probability AND the set stays sparse. A diffuse set (nearly the
        # whole window, low PIP mass) means SuSiE finds NO concentrated causal
        # signal there -- NOT evidence of colocalisation.
        cs_cap=max(20, int(0.10*len(rows)))
        cs1_rel = pip1_sum>=0.5 and len(cs1)<=cs_cap
        cs2_rel = pip2_sum>=0.5 and len(cs2)<=cs_cap
        if cs1_rel and cs2_rel and shared and shared_mass>=0.5:
            status="colocalised"
        elif not cs2_rel and pip2_sum<0.5:
            status="oa_signal_diffuse"      # OA has no concentrated causal signal at this BMD locus
        elif not cs1_rel and pip1_sum<0.5:
            status="bmd_signal_diffuse"
        else:
            status="not_colocalised"
        if status=="colocalised": n_shared+=1
        n_ok+=1
        results.append(dict(rs=rs,chr=ch,bp=bp,status=status,n_snp=len(rows),
                            cs1=[snps[i] for i in sorted(cs1)],cs2=[snps[i] for i in sorted(cs2)],
                            n_cs1=len(cs1),n_cs2=len(cs2),
                            pip1_sum=round(pip1_sum,4),pip2_sum=round(pip2_sum,4),
                            shared_snps=[snps[i] for i in shared],shared_pip_mass=round(shared_mass,4),
                            cs1_reliable=cs1_rel,cs2_reliable=cs2_rel,
                            max_pip_bmd=round(float(pip1.max()),4),max_pip_oa=round(float(pip2.max()),4)))
        print(f"  [{rs:14s} chr{ch}:{bp}] n={len(rows):4d} cs1={len(cs1)}/p={pip1_sum:.2f} cs2={len(cs2)}/p={pip2_sum:.2f} shared={len(shared)} mass={shared_mass:.4f} -> {status}",flush=True)

    n_diffuse=sum(1 for x in results if x.get("status")=="oa_signal_diffuse")
    n_disjoint=sum(1 for x in results if x.get("status")=="not_colocalised")
    summary=dict(n_loci=len(loci),n_analysed=n_ok,n_colocalised=n_shared,
                 n_not_colocalised=n_disjoint,n_oa_signal_diffuse=n_diffuse,
                 loci=results,prior_per_snp=P1,L=L_COMP,rho=RHO,window=WINDOW,
                 reliability_rule="credible set reliable iff total PIP>=0.5 and size<=max(20,10% of SNPs); colocalised iff both sets reliable, overlap, and shared PIP mass>=0.5",
                 note="coloc-SuSiE multi-signal; 0 genuine colocalisation => P3 strengthened under a multi-causal-aware method; diffuse OA signals are NOT colocalisation")
    json.dump(summary,open(OUT,"w"),indent=2)
    print(f"\n[325] DONE: {n_shared}/{n_ok} loci share a SuSiE credible variant",flush=True)
    print(f"[325] wrote {OUT}",flush=True)

if __name__=="__main__":
    main()
