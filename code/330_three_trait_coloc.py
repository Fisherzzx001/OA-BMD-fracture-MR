#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
330_three_trait_coloc.py  —  T2-C: three-trait (BMD/OA/fracture) colocalisation
==============================================================================
HyPrColoc (Foley et al. 2021) tests whether a shared causal variant underlies
>=2 traits. Full HyPrColoc uses a multi-signal branch-and-bound; with no R here
we implement the K-trait single-causal ABF generalisation of coloc (Wallace
2020, "coloc.abf for K traits"), the SAME unit-single-causal assumption as the
manuscript's published 0/31. This is a robustness *extension*, not a claim of
full HyPrColoc equivalence (honestly flagged).

For each of the 31 BMD loci (+/-500 kb), three traits:
  t1 BMD (quant, N=31,986), t2 knee OA (cc, GCST90566800),
  t3 any fracture (cc, CPBB Qian 2025, 35,192 cases/317,599 controls).
Per-SNP log-ABF for each trait; hypotheses over causal configuration:
  H0 none; H1 only BMD; H2 only OA; H3 only fracture;
  H4 ALL THREE share one causal SNP (the SOS "one set-point" prediction);
  H123 (any two share) reported as pairwise BMD-OA / BMD-FX / OA-FX.
Posterior of H4 (three-way) and of any-pair-shared per locus.
If H4 ~0 at all loci -> no single shared causal variant (consistent with SOS
being diffuse/polygenic and operating via sign opposition across MANY variants,
not one locus); reported as such, never forced.

Output: code/output/genetics/three_trait_coloc.json
"""
import os, csv, gzip, math, json
import numpy as np

ROOT=os.path.dirname(os.path.abspath(__file__))
GEN=os.path.join(ROOT,"output","genetics")
HARM=os.path.join(GEN,"harmonised","BMD_total__OA_knee.csv")
BMD="C:/Users/Administrator/WorkBuddy/2026-09-08-20-17-24/OA_FLAGSHIP/07_reference/ldsc/sumstats_total.sumstats.gz"
OA="G:/OA_B_PROJECT_cache/l7/GCST90566800.tsv.gz"
FX="G:/OA_B_PROJECT_cache/l7/CPBB_Fracture_any-type.tsv"
OUT=os.path.join(GEN,"three_trait_coloc.json")
WINDOW=500_000; BMD_N=31986.0
SD={"bmd":0.15,"oa":0.2,"fx":0.2}
P1=P2=P3=1e-4; P123=5e-6   # per-trait & shared priors (coloc convention)

def logsumexp(a):
    a=np.asarray(a,float); m=float(np.max(a))
    return m+math.log(float(np.sum(np.exp(a-m)))) if np.isfinite(m) else m

def var_beta(n,f,typ,s):
    f=np.clip(f,1e-6,1-1e-6); v=1.0/(2.0*float(n)*f*(1.0-f))
    if typ=="cc": v=v/(s*(1.0-s))
    return v

def log_abf(z,v,sd):
    w=sd**2; r=w/(w+np.asarray(v,float))
    return 0.5*np.log1p(-r)+0.5*r*(np.asarray(z,float)**2)

COMP={"A":"T","T":"A","C":"G","G":"C"}
def align(a1g,a2g,ea,oa,ze):
    a1g,a2g,ea,oa=map(lambda x:str(x).upper(),(a1g,a2g,ea,oa))
    if {a1g,a2g}=={ea,oa}: return ze if a1g==ea else -ze
    if {a1g,a2g}=={COMP.get(ea,""),COMP.get(oa,"")}: return ze if a1g==COMP.get(oa,"") else -ze
    return None

def load_loci():
    with open(HARM) as f:
        return [(int(r["CHR"]),int(r["BP"])) for r in csv.DictReader(f)]

def windows(loci):
    # per-chromosome list of (lo,hi,idx) so each SNP is assigned to its OWN locus
    w={}
    for i,(ch,bp) in enumerate(loci): w.setdefault(ch,[]).append((bp-WINDOW,bp+WINDOW,i))
    return w

def _in_window(ch,bp,w):
    """Return list of locus indices whose window contains (ch,bp)."""
    return [i for lo,hi,i in w.get(ch,[]) if lo<=bp<=hi]

def scan_bmd(w):
    from collections import defaultdict
    out=defaultdict(dict)
    with gzip.open(BMD,"rt") as f:
        h=f.readline().rstrip("\n").split("\t"); ix={x:i for i,x in enumerate(h)}
        for line in f:
            p=line.rstrip("\n").split("\t")
            try: ch=int(p[ix["CHR"]]);bp=int(p[ix["BP"]])
            except: continue
            for i in _in_window(ch,bp,w):
                try: out[i][bp]=dict(a1=p[ix["A1"]].upper(),a2=p[ix["A2"]].upper(),maf=float(p[ix["MAF"]]),z=float(p[ix["Z"]]))
                except: pass
    return out

def scan_oa(w):
    from collections import defaultdict
    out=defaultdict(dict)
    with gzip.open(OA,"rt") as f:
        h=f.readline().rstrip("\n").split("\t"); ix={x:i for i,x in enumerate(h)}
        for line in f:
            p=line.rstrip("\n").split("\t")
            try: ch=int(p[ix["chromosome"]]);bp=int(p[ix["base_pair_location"]])
            except: continue
            for i in _in_window(ch,bp,w):
                try:
                    cu=float(p[ix["ci_upper"]]);cl=float(p[ix["ci_lower"]]);orr=float(p[ix["odds_ratio"]])
                    se=(math.log(cu)-math.log(cl))/3.92
                    if se<=0 or orr<=0: continue
                    nt=float(p[ix["NCASES"]])+float(p[ix["NCONTROLS"]]); s=float(p[ix["NCASES"]])/nt
                    out[i][bp]=dict(ea=p[ix["effect_allele"]].upper(),oa=p[ix["other_allele"]].upper(),
                                    z=math.log(orr)/se,eaf=float(p[ix["effect_allele_frequency"]]),nt=nt,s=s)
                except: pass
    return out

def scan_fx(w):
    from collections import defaultdict
    out=defaultdict(dict)
    nc=35192.0; nt=nc+317599.0; s=nc/nt
    with open(FX) as f:
        h=f.readline().rstrip("\n").split("\t"); ix={x:i for i,x in enumerate(h)}
        for line in f:
            p=line.rstrip("\n").split("\t")
            try: ch=int(p[ix["chromosome"]]);bp=int(p[ix["base_pair_location"]])
            except: continue
            for i in _in_window(ch,bp,w):
                try:
                    beta=float(p[ix["beta"]]);se=float(p[ix["standard_error"]])
                    if se<=0: continue
                    out[i][bp]=dict(ea=p[ix["effect_allele"]].upper(),oa=p[ix["other_allele"]].upper(),
                                    z=beta/se,eaf=float(p[ix["effect_allele_frequency"]]),nt=nt,s=s)
                except: pass
    return out

def coloc3(z1,v1,z2,v2,z3,v3):
    """K=3 single-causal ABF; posterior over {none,only1,only2,only3,1&2,1&3,2&3,all3}."""
    l1=log_abf(z1,v1,SD["bmd"]); l2=log_abf(z2,v2,SD["oa"]); l3=log_abf(z3,v3,SD["fx"])
    L1=logsumexp(l1);L2=logsumexp(l2);L3=logsumexp(l3)
    L12=logsumexp(l1+l2);L13=logsumexp(l1+l3);L23=logsumexp(l2+l3)
    L123=logsumexp(l1+l2+l3)
    lh=np.array([0.0,
        math.log(P1)+L1, math.log(P2)+L2, math.log(P3)+L3,
        math.log(P123)+L12, math.log(P123)+L13, math.log(P123)+L23,
        math.log(P123)+L123])
    lh=lh-logsumexp(lh); pp=np.exp(lh)
    return dict(H0=pp[0],H1=pp[1],H2=pp[2],H3=pp[3],
                H12=pp[4],H13=pp[5],H23=pp[6],H123=pp[7],
                any_pair=float(pp[4]+pp[5]+pp[6]+pp[7]))

def main():
    loci=load_loci(); w=windows(loci)
    print(f"[330] {len(loci)} loci; scanning 3 traits...",flush=True)
    B=scan_bmd(w); print(f"  BMD window SNPs {len(B)}",flush=True)
    O=scan_oa(w);  print(f"  OA window SNPs {len(O)}",flush=True)
    F=scan_fx(w);  print(f"  Fracture window SNPs {len(F)}",flush=True)
    rows=[]
    for idx,(ch,bp) in enumerate(loci):
        Bd=B.get(idx,{});Od=O.get(idx,{});Fd=F.get(idx,{})
        common=sorted(set(Bd)&set(Od)&set(Fd))
        z1=[];v1=[];z2=[];v2=[];z3=[];v3=[]
        for k in common:
            b=Bd[k];o=Od[k];f=Fd[k]
            za_o=align(b["a1"],b["a2"],o["ea"],o["oa"],o["z"])
            za_f=align(b["a1"],b["a2"],f["ea"],f["oa"],f["z"])
            if za_o is None or za_f is None: continue
            z1.append(b["z"]);v1.append(var_beta(BMD_N,b["maf"],"quant",None))
            z2.append(za_o);v2.append(var_beta(o["nt"],o["eaf"],"cc",o["s"]))
            z3.append(za_f);v3.append(var_beta(f["nt"],f["eaf"],"cc",f["s"]))
        if len(z1)<20:
            rows.append(dict(chr=ch,bp=bp,n=len(z1),status="insufficient")); continue
        r=coloc3(np.array(z1),np.array(v1),np.array(z2),np.array(v2),np.array(z3),np.array(v3))
        r.update(chr=ch,bp=bp,n=len(z1),status="ok")
        # SNP-level sign pattern at the strongest BMD SNP.
        # SOS kernel: one BMD-raising allele is PROTECTIVE at the whole-bone
        # endpoint (lowers fracture => opposite sign to BMD) yet PATHOGENIC at
        # the joint (raises OA => SAME sign as BMD). So the correct pattern is
        # sign(z_BMD)==sign(z_OA) AND sign(z_BMD)!=sign(z_FX), comparing all at
        # the SAME allele (z2/z3 already oriented to BMD A1 by align()).
        # FIX 2026-11-10: old code required BMD!=OA AND BMD!=FX, which wrongly
        # flagged anti-SOS loci (BMD vs OA opposite) as SOS and mislabelled the
        # 3 genuine H123>=0.8 SOS loci as not-SOS. ABF (H123) is unaffected.
        k=int(np.argmax(np.array(z1)**2))
        r["top_bmd_z"]=round(float(z1[k]),3); r["top_oa_z_at_bmd"]=round(float(z2[k]),3); r["top_fx_z_at_bmd"]=round(float(z3[k]),3)
        r["sos_sign_pattern"]=bool(np.sign(z1[k])==np.sign(z2[k]) and np.sign(z1[k])!=np.sign(z3[k]))
        rows.append(r)
    ok=[r for r in rows if r["status"]=="ok"]
    n_h123=sum(1 for r in ok if r["H123"]>=0.8)
    n_anypair=sum(1 for r in ok if r["any_pair"]>=0.8)
    n_sos_sign=sum(1 for r in ok if r.get("sos_sign_pattern"))
    out=dict(design="T2-C: 3-trait coloc (BMD/knee-OA/any-fracture) at 31 BMD loci",
             method="K=3 single-causal ABF generalisation of coloc (Wallace 2020); "
                    "SAME unit-single-causal assumption as manuscript 0/31; "
                    "NOT full HyPrColoc multi-signal branch-and-bound (no R available)",
             priors=dict(p1=P1,p2=P2,p3=P3,p_shared=P123),
             n_loci=len(loci),n_analysed=len(ok),
             n_three_way_H4_ge_0p8=n_h123,n_any_pair_ge_0p8=n_anypair,
             n_loci_SOS_sign_pattern=n_sos_sign,
             interpretation=("No locus has a single shared causal variant for all three traits "
                             "=> consistent with SOS being polygenic/diffuse (one set-point across "
                             "MANY variants with sign opposition), not one pleiotropic locus."
                             if n_h123==0 else "Three-way shared locus detected; inspect sign pattern."),
             loci=rows)
    json.dump(out,open(OUT,"w",encoding="utf-8"),ensure_ascii=False,indent=2)
    print(f"\n[330] {n_h123}/{len(ok)} loci H4(3-way)>=0.8; {n_anypair}/{len(ok)} any-pair; "
          f"{n_sos_sign}/{len(ok)} loci show SOS sign pattern (BMD SAME sign as OA, OPPOSITE to FX)",flush=True)
    print(f"[write] {OUT}")

if __name__=="__main__":
    main()
