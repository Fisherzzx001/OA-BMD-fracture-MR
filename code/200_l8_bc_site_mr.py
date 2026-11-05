#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
200_l8_bc_site_mr.py  --  PATH B/C, MR contrast (the core of T4 "formation-side vs destruction-side")
=================================================================================================
Replicates the proven L7 two-sample MR pipeline (code/191 + G:/.../l7_f1_f2_mr.py) but extends it
from {knee,hip} to the SITE-STRATIFIED Hatzikotoulas-2025 OA GWAS:

    BMD_total (exposure, 31 instruments)  ->  knee OA  (cartilage/meniscus-DOMINANT)
                                         ->  spine OA  (OSTEOPHYTE-DOMINANT)
                                         ->  hand  OA  (OSTEOPHYTE-DOMINANT: Heberden nodes)

This is the public-data surrogate for F5 (osteophyte outcome) and the faithful execution of T4
(formation side = osteophyte-predominant sites  vs  destruction side = cartilage-predominant knee).

Method (identical to L7, for comparability):
  * 31 instruments = harmonised BMD_total__OA_knee.csv (beta_e = BMD exposure; beta_o = OA_knee on
    BMD exposure allele).
  * Orient each outcome to the BMD exposure allele using the KNEE reference (empirically recovered
    per SNP from sign agreement with harmonised beta_o) -> guarantees spine/hand share the SAME
    orientation as the already-published knee result.
  * IVW (weights 1/se_outcome^2), weighted median, MR-Egger. Ratio method: ratio = beta_oa_oriented
    / beta_bmd ; se via delta method.

Inputs (all local; knee/hip already cached, spine/hand downloaded by code/199):
  HARM = code/output/genetics/harmonised/BMD_total__OA_knee.csv
  KNEE = G:/OA_B_PROJECT_cache/l7/GCST90566800.tsv.gz   (cartilage-dominant reference)
  SPINE= G:/OA_B_PROJECT_cache/l7/GCST90566801.tsv.gz   (osteophyte-dominant)
  HAND = G:/OA_B_PROJECT_cache/l7/GCST90566803.tsv.gz   (osteophyte-dominant)
Output:
  code/output/genetics/L8_BC_site_mr.json
  code/output/genetics/L8_BC_matched_spine.tsv , _hand.tsv  (archival 31-SNP extracts)
"""
import os, csv, math, gzip, json, sys

HERE = "C:/Users/Administrator/WorkBuddy/2026-09-08-22-17-40/OA_B_PROJECT"
L7   = "G:/OA_B_PROJECT_cache/l7"
HARM = os.path.join(HERE, "code/output/genetics/harmonised/BMD_total__OA_knee.csv")
KNEE = os.path.join(L7, "GCST90566800.tsv.gz")
SPINE= os.path.join(L7, "GCST90566801.tsv.gz")
HAND = os.path.join(L7, "GCST90566803.tsv.gz")
GEN  = os.path.join(HERE, "code/output/genetics")
OUT  = os.path.join(GEN, "L8_BC_site_mr.json")

# hg38 positions for the 31 instruments (Ensembl GRCh38) -- identical to L7
# Position table is DERIVED from the harmonised file in main() (GRCh37 coordinates, which match
# the Hatzikotoulas GWAS Catalog base_pair_location and the GEFOS sumstats). The harmonised BP is
# authoritative for matching. A hardcoded GRCh38 dict would mismatch the GRCh37 GWAS.
HG38 = {}
COMP = {"A":"T","T":"A","C":"G","G":"C"}

def normal_cdf(x): return 0.5*(1+math.erf(x/math.sqrt(2)))

def parse_gwas(path, needed):
    out, missing = {}, set(needed)
    with gzip.open(path, "rt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        idx = {h:i for i,h in enumerate(hdr)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            try:
                ch = int(p[idx["chromosome"]]); ps = int(p[idx["base_pair_location"]])
            except:
                continue
            key = (ch, ps)
            if key in needed:
                out[key] = dict(
                    ea=p[idx["effect_allele"]].upper(), oa=p[idx["other_allele"]].upper(),
                    orr=p[idx["odds_ratio"]], cu=p[idx["ci_upper"]], cl=p[idx["ci_lower"]],
                    p=p[idx["p_value"]], eaf=p[idx["effect_allele_frequency"]],
                    ncase=p[idx["NCASES"]], nctrl=p[idx["NCONTROLS"]],
                    vid=p[idx["variant_id"]])
                missing.discard(key)
    return out, missing

def to_beta(row):
    try:
        orr=float(row["orr"]); cu=float(row["cu"]); cl=float(row["cl"])
    except: return None,None,False
    if orr<=0 or cu<=0 or cl<=0: return None,None,False
    b=math.log(orr); se=(math.log(cu)-math.log(cl))/3.92
    if se<=0 or not math.isfinite(se): return None,None,False
    return b,se,True

def orient(outcome_row, bmd_ea_for_snp, outcome_tag):
    """Return (beta_oriented, se_oriented) of the OA outcome on the BMD exposure allele."""
    b,se,ok = to_beta(outcome_row)
    if not ok: return None,None,False
    ea=outcome_row["ea"]; oa=outcome_row["oa"]
    if ea==bmd_ea_for_snp:       return b, se, True
    if oa==bmd_ea_for_snp:       return -b, se, True
    # strand complement
    if COMP.get(ea,"")==bmd_ea_for_snp: return b, se, True
    if COMP.get(oa,"")==bmd_ea_for_snp: return -b, se, True
    return None,None,False

def ivw(rows):
    num=sum(bh*be/se_o**2 for (_,be,_,bh,se_o,_,_) in rows if se_o>0)
    den=sum(be**2/se_o**2 for (_,be,_,bh,se_o,_,_) in rows if se_o>0)
    b=num/den; se=math.sqrt(1/den); z=b/se; p=2*(1-normal_cdf(abs(z)))
    return b,se,z,p
def wmedian(rows):
    rs=sorted(rows,key=lambda t:t[5]); w=[1/t[6]**2 for t in rs]; cw=sum(w); half=cw/2; acc=0
    for t,wk in zip(rs,w):
        acc+=wk
        if acc>=half: return t[5]
    return rs[len(rs)//2][5]
def mregger(rows):
    import numpy as np
    X=np.array([[1.0,be] for (_,be,_,bh,se_o,_,_) in rows if se_o>0])
    Y=np.array([bh for (_,be,_,bh,se_o,_,_) in rows if se_o>0])
    W=np.array([1/se_o**2 for (_,be,_,bh,se_o,_,_) in rows if se_o>0])
    WX=X*W[:,None]; betahat=np.linalg.solve(X.T@WX, X.T@(W*Y))
    cov=np.linalg.inv(X.T@WX)
    se=np.sqrt(np.diag(cov))
    return betahat[1],betahat[0],se[1],se[0]

def run_outcome(name, gwas, bmd_ea_map, inst):
    rows=[]
    for rs,r in inst.items():
        ch,ps=HG38[rs]
        be=float(r["beta_e"]); see=float(r["se_e"])
        if (ch,ps) not in gwas: continue
        bo,so,ok = orient(gwas[(ch,ps)], bmd_ea_map[rs], name)
        if not ok or abs(be)<1e-9: continue
        ratio=bo/be; se_r=abs(ratio)*math.sqrt((so/bo)**2+(see/be)**2)
        rows.append((rs,be,see,bo,so,ratio,se_r))
    if not rows:
        return None
    b,se,z,p=ivw(rows); wm=wmedian(rows); slope,inter,se_slope,se_inter=mregger(rows)
    from math import erf
    def ncdf(x): return 0.5*(1+erf(x/math.sqrt(2)))
    egger_inter_p = 2*(1-ncdf(abs(inter/se_inter))) if se_inter>0 else float('nan')
    return dict(name=name, n=len(rows),
                ivw_beta=b, ivw_se=se, ivw_z=z, ivw_p=p,
                or_per_sd=math.exp(b), ci=[math.exp(b-1.96*se),math.exp(b+1.96*se)],
                wm_ratio=wm, egger_slope=slope, egger_intercept=inter,
                egger_int_se=se_inter, egger_int_p=egger_inter_p,
                per_snp=[dict(rs=rs,be=be,bh=bo,ratio=ratio,se_r=se_r) for (rs,be,_,bo,_,ratio,se_r) in rows])

def main():
    if not (os.path.exists(KNEE) and os.path.exists(SPINE) and os.path.exists(HAND)):
        print("[FATAL] missing GWAS input(s). KNEE/SPINE/HAND must all exist:",
              {k:os.path.exists(v) for k,v in [("KNEE",KNEE),("SPINE",SPINE),("HAND",HAND)]})
        sys.exit(1)
    inst={}
    with open(HARM) as f:
        for r in csv.DictReader(f): inst[r["SNP"]]=r
    # positions from harmonised (GRCh37) -- authoritative for matching the GWAS
    for rs, r in inst.items():
        HG38[rs] = (int(r["CHR"]), int(r["BP"]))
    needed=set(HG38.values())
    print(f"[load] {len(inst)} instruments; {len(needed)} hg38 positions")

    knee,_ = parse_gwas(KNEE, needed)
    spine,_= parse_gwas(SPINE, needed)
    hand,_ = parse_gwas(HAND, needed)
    print(f"[parse] knee {len(knee)}/{len(needed)}  spine {len(spine)}/{len(needed)}  hand {len(hand)}/{len(needed)}")

    # ---- derive bmd_ea per SNP from KNEE (empirical, validated vs harmonised beta_o) ----
    bmd_ea_map={}; agree=0; tot=0; bad=[]
    for rs,r in inst.items():
        ch,ps=HG38[rs]
        bo=float(r["beta_o"]); 
        if (ch,ps) not in knee: 
            bmd_ea_map[rs]=r["exp_allele"] if "exp_allele" in r else None
            continue
        rk=knee[(ch,ps)]; bk,_,ok=to_beta(rk)
        if not ok or abs(bo)<1e-9:
            bmd_ea_map[rs]=rk["ea"]; continue
        same = (bk/bo)>0
        ea = rk["ea"] if same else rk["oa"]
        bmd_ea_map[rs]=ea
        pred = bk if same else -bk
        tot+=1
        if abs(pred-bo)<0.004: agree+=1
        else: bad.append((rs,round(pred,4),round(bo,4),round(bk,4)))
    print(f"[orient] knee-validation recover BMD EA: {agree}/{tot} within 0.004" + (f"; MISMATCH {bad[:6]}" if bad else ""))

    results={}
    results["knee"] = run_outcome("knee", knee, bmd_ea_map, inst)
    results["spine"]= run_outcome("spine",spine,bmd_ea_map, inst)
    results["hand"] = run_outcome("hand", hand, bmd_ea_map, inst)

    # ---- T4 contrast summary ----
    print("\n=== T4 SITE CONTRAST: BMD->OA (OR per SD, 95% CI) ===")
    for k in ("knee","spine","hand"):
        r=results[k]
        if r is None: print(f"  {k}: NO DATA"); continue
        tag={"knee":"cartilage-dominant","spine":"osteophyte-dominant","hand":"osteophyte-dominant"}[k]
        print(f"  {k:6s} [{tag:20s}] OR/SD={r['or_per_sd']:.4f} [{r['ci'][0]:.4f}-{r['ci'][1]:.4f}] "
              f"p={r['ivw_p']:.2e} nSNPs={r['n']}  WM={r['wm_ratio']:.4f} EggerInt={r['egger_intercept']:.4f}")

    # ---- write JSON ----
    json.dump(results, open(OUT,"w"), indent=2)
    # ---- archival matched.tsv (31-SNP extracts) for spine/hand ----
    for tag,gw in (("spine",spine),("hand",hand)):
        outp=os.path.join(GEN, f"L8_BC_matched_{tag}.tsv")
        with open(outp,"w") as f:
            f.write("rs\tchr\tbp\tE\tO\tOR\tSE\teaf\tp\tvariant_id\tci_u\tci_l\tNCASE\tNCTRL\n")
            for rs in inst:
                ch,ps=HG38[rs]
                if (ch,ps) not in gw: continue
                rw=gw[(ch,ps)]
                f.write("\t".join([rs,str(ch),str(ps),rw["ea"].lower(),rw["oa"].lower(),
                    rw["orr"],"#NA",rw["eaf"],rw["p"],rw["vid"],rw["cu"],rw["cl"],
                    rw["ncase"],rw["nctrl"]])+"\n")
        print(f"[write] {outp}")
    print(f"\n[done] wrote {OUT}")

if __name__=="__main__":
    main()
