#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
328_chai2026_independent_replication.py  —  T2-A: independent instrument-set replication
==========================================================================================
Chai RC et al. Nat Genet 2026;58:1530-1545 (PMID 42432248), 448,010 UKB eBMD
(= GWAS Catalog GCST90726625, 1,246 GWS associations). 本文主结果用的是
Qian 2025 DXA-aBMD 31 工具 (N=31,986)。Chai 是**大一个数量级的独立 BMD 工具来源**。

设计：
  1) 抓 GCST90726625 全部 1,246 lead (v2 associations API, size=200 × 7 页)；
  2) 解析坐标/β方向(unit increase/decrease)/CI->SE/EAF；
  3) 1000G EUR.bim 按 rsid 取 effect/other allele，按坐标取 LD；
  4) LD clumping: r2<0.001, 500 kb, 取 p 最小者 -> 独立工具集；
  5) 对 knee OA (GCST90566800, Hatzikotoulas 2025, 本地) 按坐标精确匹配 β/SE；
  6) 等位协调 -> IVW(固定/随机)/Egger/加权中位数/Q,I2；
  7) 与本文 31 工具的 OR 1.060-1.099 区间作**尺度外**方向比较（eBMD vs DXA 不同尺度，
     只比方向不比大小，与稿件 §4 对 Hartley 的处理一致）。

诚实：Chai eBMD 是 heel-ultrasound/量级不同表型、仅 EUR；这是独立工具的方向复制，
      不是精确量级复制。不改动 docs/08；不新增引用（refs 满额，留返修）。
输出：code/output/genetics/Chai2026_independent_replication.json
       code/output/genetics/Chai2026_leads_full.tsv
"""
import os, json, math, csv, gzip, time, struct, importlib.util
import numpy as np
import urllib.request

ROOT = os.path.dirname(os.path.abspath(__file__))
GEN  = os.path.join(ROOT, "output", "genetics")
CACHE= "G:/OA_B_PROJECT_cache/l7"
BIM  = "G:/REF/1000G_EUR/EUR.bim"
BED  = "G:/REF/1000G_EUR/EUR.bed"
OA   = os.path.join(CACHE, "GCST90566800.tsv.gz")
API  = "https://www.ebi.ac.uk/gwas/api/v2/studies/GCST90726625/associations"
LEADS_TSV = os.path.join(GEN, "Chai2026_leads_full.tsv")
OUT  = os.path.join(GEN, "Chai2026_independent_replication.json")
N_SAMPLES=503

_spec = importlib.util.spec_from_file_location("t323", os.path.join(ROOT, "323_f10_mvp_tabix_mr.py"))
t323 = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(t323)


def fetch_all():
    cache=os.path.join(GEN,"_chai_all_pages.json")
    if os.path.exists(cache):
        pages=json.load(open(cache,encoding="utf-8")); print(f"[328] cached {len(pages)} pages"); return pages
    pages=[]
    for p in range(7):
        url=f"{API}?size=200&page={p}"
        for attempt in range(3):
            try:
                req=urllib.request.Request(url,headers={"User-Agent":"OA-BMD/328"})
                with urllib.request.urlopen(req,timeout=60) as r:
                    d=json.loads(r.read().decode())
                pages.append(d); print(f"  page {p+1}/7: {len(d.get('_embedded',{}).get('associations',[]))} assocs",flush=True)
                break
            except Exception as ex:
                print(f"  [retry {attempt+1}] page {p}: {type(ex).__name__}"); time.sleep(3)
        else:
            print(f"[328] FATAL page {p} failed"); return None
    json.dump(pages,open(cache,"w",encoding="utf-8"))
    return pages


def parse_beta(beta_str, ci_str, pval, pexp, risk_freq):
    """'0.02 unit decrease' -> signed beta; CI '[0.014-0.026]' -> SE."""
    if not beta_str: return None
    parts=str(beta_str).split()
    try: mag=float(parts[0])
    except: return None
    sign=-1.0 if "decrease" in str(beta_str).lower() else 1.0
    b=sign*mag
    se=None
    if ci_str:
        import re
        m=re.findall(r"\d*\.?\d+",str(ci_str))   # CI bounds are positive; avoid '-' as sign
        if len(m)>=2:
            lo,hi=float(m[0]),float(m[1]); se=(hi-lo)/3.92 if hi>lo else None
    p=float(pval)*(10.0**int(pexp)) if pexp is not None else None
    eaf=None
    if risk_freq not in (None,""):
        try: eaf=float(risk_freq)
        except: pass
    return b,se,p,eaf


def parse_leads(pages):
    rows=[]
    for d in pages:
        for a in d.get("_embedded",{}).get("associations",[]):
            rs=a.get("riskAllele",[{}])[0].get("key")
            loc=a.get("locations",[None])[0]
            if not rs or not loc: continue
            try: ch,bp=loc.split(":")[0],int(loc.split(":")[1])
            except: continue
            pb=parse_beta(a.get("beta"),a.get("ci"),a.get("pValue"),a.get("pValueExponent"),a.get("riskFrequency"))
            if not pb: continue
            b,se,p,eaf=pb
            if se is None or se<=0: continue
            rows.append(dict(rs=rs,chr=ch,bp=bp,beta=b,se=se,p=p,eaf=eaf,
                             gene=",".join(a.get("mappedGenes",[])[:3])))
    # dedup by rs (keep most significant)
    best={}
    for r in rows:
        if r["rs"] not in best or (r["p"] or 1)<(best[r["rs"]]["p"] or 1): best[r["rs"]]=r
    rows=sorted(best.values(),key=lambda x:(x["p"] or 1))
    return rows


import bisect

def load_bim():
    rs2={}; by_chr={}
    with open(BIM) as f:
        for i,line in enumerate(f):
            p=line.rstrip("\n").split("\t")
            if len(p)<5: continue
            rs2[p[1]]=dict(row=i,chr=p[0],bp=int(p[3]),a1=p[4].upper(),a2=p[5].upper())
            by_chr.setdefault(p[0],[]).append((int(p[3]),p[1],i))
    for c in by_chr: by_chr[c].sort()          # sort by bp for bisect
    return rs2,by_chr

# module-level bed handle + dosage cache + per-chromosome bp array (memory-safe)
_BED_FH=None; _DOS_CACHE={}; _BP_CACHE={}
def dosage(row):
    # MEMORY FIX 2026-11-10: do NOT permanently cache standardised dosages.
    # 921 leads x ~1800 window SNPs = >1M float64[503] arrays = several GB -> OOM.
    # Keep the single reusable file handle; rely on OS page cache for the 126-byte
    # reads. Statistical cost is UNCHANGED (all leads, r2<0.001, +/-500 kb).
    global _BED_FH
    if _BED_FH is None: _BED_FH=open(BED,"rb")
    by=(N_SAMPLES+3)//4; off=3+row*by
    _BED_FH.seek(off); raw=_BED_FH.read(by)
    if len(raw)<by: return None
    # vectorised PLINK 2-bit decode (little-endian: low 2 bits = SNP0).
    # code 0=homA1(dosage0), 1=missing, 2=het(dosage1), 3=homA2(dosage2)
    bits=np.unpackbits(np.frombuffer(raw,dtype=np.uint8),bitorder="little")
    codes=bits[:N_SAMPLES*2].reshape(N_SAMPLES,2)
    # PLINK 2-bit: low 2 bits = sample0; code = bit_high*2 + bit_low
    # 00=homA1(0), 01=missing, 10=het(1), 11=homA2(2)
    g=codes[:,1].astype(np.uint8)*2+codes[:,0].astype(np.uint8)
    d=np.where(g==3,2.0,np.where(g==2,1.0,0.0)).astype(float); d[g==1]=np.nan
    mu=float(np.nanmean(d)) if np.isfinite(d).any() else 0.0
    if np.isnan(d).any(): d=np.where(np.isnan(d),mu,d)
    sd=float(d.std()); return (d-mu)/sd if sd>0 else None

def clump(leads, by_chr, r2_cut=0.001, kb=500):
    """Greedy LD clump on p-ordered leads; bisect window lookup, on-demand dosages."""
    keep=[]; excluded=set()
    for li,lead in enumerate(leads):
        if li and li%100==0:
            print(f"[328] clump {li}/{len(leads)} kept={len(keep)} excluded={len(excluded)}",flush=True)
        if lead["rs"] in excluded: continue
        keep.append(lead)
        # Chai lead coordinates come from GWAS Catalog (GRCh38); 1000G EUR is GRCh37.
        # Use the rsid-mapped 1000G row (rs2) for LD, not the catalog bp.
        r=rs2_global.get(lead["rs"])
        if not r: continue
        arr=by_chr.get(r["chr"],[])
        # cache the sorted bp numpy array ONCE per chromosome (avoids rebuilding
        # a ~1M-element python list for every lead -> MemoryError)
        c=r["chr"]
        if c not in _BP_CACHE:
            _BP_CACHE[c]=np.array([x[0] for x in arr],dtype=np.int64)
        bps=_BP_CACHE[c]
        lo=int(np.searchsorted(bps,r["bp"]-kb*1000,side="left"))
        hi=int(np.searchsorted(bps,r["bp"]+kb*1000,side="right"))
        ld=dosage(r["row"])
        if ld is None: continue
        for bp,rs,row in arr[lo:hi]:
            if rs==lead["rs"] or rs in excluded: continue
            other=dosage(row)
            if other is None: continue
            r2=float(np.dot(ld,other)/N_SAMPLES)**2
            if r2>=r2_cut: excluded.add(rs)
    return keep


def scan_oatarget(positions):
    """Scan knee OA sumstats once; return (chr,bp)->(beta_logor, se, ea, oa, eaf, ncase, nctrl)."""
    want=set(positions)
    out={}
    with gzip.open(OA,"rt") as f:
        h=f.readline().rstrip("\n").split("\t"); ix={x:i for i,x in enumerate(h)}
        n=0
        for line in f:
            n+=1
            if n%6_000_000==0: print(f"  [OA scan] {n//1_000_000}M",flush=True)
            p=line.rstrip("\n").split("\t")
            try: ch=p[ix["chromosome"]]; bp=int(p[ix["base_pair_location"]])
            except: continue
            if (ch,bp) not in want: continue
            try:
                orr=float(p[ix["odds_ratio"]]); cu=float(p[ix["ci_upper"]]); cl=float(p[ix["ci_lower"]])
                se=(math.log(cu)-math.log(cl))/3.92
                if se<=0 or orr<=0: continue
                out[(ch,bp)]=dict(b=math.log(orr),se=se,ea=p[ix["effect_allele"]].upper(),
                                  oa=p[ix["other_allele"]].upper(),eaf=float(p[ix["effect_allele_frequency"]]),
                                  ncase=float(p[ix["NCASES"]]),nctrl=float(p[ix["NCONTROLS"]]))
            except: pass
    return out


def comp_align(a1g,a2g,ea,oa,bo):
    if {ea,oa}!={a1g,a2g}: return None
    return bo if ea==a1g else -bo


def run_mr(recs, exposure_label):
    be=np.array([r["be"] for r in recs]); se_e=np.array([r["se_e"] for r in recs])
    bo=np.array([r["bo"] for r in recs]); se_o=np.array([r["se_o"] for r in recs])
    ratio=bo/be; var=se_o**2/be**2 + (bo**2/be**4)*se_e**2; w=1/var
    b=float(np.sum(w*ratio)/np.sum(w)); se=float(math.sqrt(1/np.sum(w)))
    Q=float(np.sum(w*(ratio-b)**2)); df=len(recs)-1; C=max(0,(Q-df)/np.sum(w)); wr=1/(var+C)
    br=float(np.sum(wr*ratio)/np.sum(wr)); ser=float(math.sqrt(1/np.sum(wr)))
    I2=max(0,(Q-df)/Q)*100 if Q>0 else 0
    # egger
    W=1/se_o**2; sw=W.sum(); sx=np.sum(W*be)/sw; sy=np.sum(W*bo)/sw
    sxx=np.sum(W*(be-sx)**2); sxy=np.sum(W*(be-sx)*(bo-sy)); slope=sxy/sxx; inter=sy-slope*sx
    se_i=math.sqrt(1/sw+sx**2/sxx)
    order=np.argsort(ratio); cum=0; wm=None
    for i in order:
        cum+=wr[i]
        if cum>=0.5*wr.sum(): wm=float(ratio[i]); break
    return dict(exposure=exposure_label,n=len(recs),
                ivw_fixed_OR=round(math.exp(b),4),ivw_fixed_p=float(math.erfc(abs(b/se)/math.sqrt(2))),
                ivw_random_OR=round(math.exp(br),4),ivw_random_p=float(math.erfc(abs(br/ser)/math.sqrt(2))),
                Q=round(Q,2),df=df,I2=round(I2,1),egger_intercept=round(inter,5),
                egger_intercept_p=float(math.erfc(abs(inter/se_i)/math.sqrt(2))),
                wmedian_OR=round(math.exp(wm),4) if wm is not None else None)


def verdict_for(mr):
    """Honest verdict: separate DIRECTION agreement from STATISTICAL robustness.
    Consistent with the manuscript convention that heterogeneity-dominated
    random-effects results are not read as significant causal estimates."""
    fixed_or=mr["ivw_fixed_OR"]; rand_or=mr["ivw_random_OR"]; rand_p=mr["ivw_random_p"]
    i2=mr["I2"]; egger_p=mr["egger_intercept_p"]
    direction_up = fixed_or>1.0 and rand_or>1.0
    direction_dn = fixed_or<1.0 and rand_or<1.0
    robust = rand_p<0.05 and i2<50
    het_dom = i2>=50
    if direction_up and robust:
        return ("DIRECTION-ROBUST REPLICATED: independent Chai 2026 eBMD instruments raise knee-OA risk "
                "under random effects (OR>1, p<0.05, I2<50), corroborating the DXA-aBMD 31-instrument "
                "result. Magnitudes are not comparable across eBMD/DXA scales.")
    if direction_up and het_dom:
        return ("DIRECTION-CONSISTENT, EFFECT-UNRELIABLE: independent Chai 2026 eBMD instruments give "
                "OR>1 under both fixed and random effects, matching the DIRECTION of the DXA-aBMD 31-"
                f"instrument result, BUT the random-effects estimate is not significant (OR={rand_or}, "
                f"p={rand_p:.3g}) and heterogeneity dominates (I2={i2}%). The fixed-effect significance "
                "is not read causally. This is directional corroboration, NOT an independent significant "
                "replication; eBMD/DXA scales are not magnitude-comparable.")
    if direction_up or direction_dn:
        return ("DIRECTION-ONLY: point estimates agree in sign but neither robust random-effects "
                "significance nor low heterogeneity was met; reported honestly, not a replication claim.")
    return ("DIRECTION-NOT-REPLICATED: Chai eBMD instruments do not give a consistent sign; reported honestly.")


def _reverdict():
    """Recompute verdict from the already-written JSON (avoids 10-min re-run)."""
    d=json.load(open(OUT,encoding="utf-8"))
    if "mr_knee" not in d:
        print("[328] no mr_knee in JSON; nothing to re-verdict"); return
    d["verdict"]=verdict_for(d["mr_knee"])
    d["verdict_basis"]="recomputed 2026-11-10 with heterogeneity-aware verdict_for(); ABF/MR numbers unchanged"
    json.dump(d,open(OUT,"w",encoding="utf-8"),ensure_ascii=False,indent=2)
    print("[328] revised verdict:",d["verdict"])


def main():
    global pos2_global, rs2_global
    pages=fetch_all()
    if not pages: return
    leads=parse_leads(pages)
    print(f"[328] parsed {len(leads)} unique Chai lead associations",flush=True)
    # write full leads
    with open(LEADS_TSV,"w",newline="") as f:
        wr=csv.writer(f,delimiter="\t"); wr.writerow(["rs","chr","bp","beta_eBMD","se","p","eaf","gene"])
        for r in leads: wr.writerow([r["rs"],r["chr"],r["bp"],f"{r['beta']:.5g}",f"{r['se']:.5g}",f"{r['p']:.3g}",r["eaf"],r["gene"]])
    rs2,by_chr=load_bim()
    rs2_global=rs2
    # map leads to 1000G GRCh37 coordinates + alleles (catalog locations are GRCh38!)
    mapped=[]
    for r in leads:
        m=rs2.get(r["rs"])
        if m: mapped.append(dict(r,a1=m["a1"],a2=m["a2"],row=m["row"],chr37=m["chr"],bp37=m["bp"]))
    print(f"[328] {len(mapped)}/{len(leads)} leads mapped to 1000G EUR (GRCh37)",flush=True)
    clumped=clump(mapped,by_chr,r2_cut=0.001,kb=500)
    print(f"[328] LD-clumped (r2<0.001,500kb) -> {len(clumped)} independent eBMD instruments",flush=True)
    # scan knee OA (GRCh37) using the rsid-mapped GRCh37 coordinates
    positions=[(c["chr37"],c["bp37"]) for c in clumped]
    oa=scan_oatarget(set(positions))
    print(f"[328] knee-OA matched {len(oa)}/{len(clumped)} independent instruments",flush=True)
    recs=[]
    for c in clumped:
        o=oa.get((c["chr37"],c["bp37"]))
        if not o: continue
        recs.append(dict(rs=c["rs"],se_e=c["se"],se_o=o["se"],bo=o["b"],
                         a1=c["a1"],a2=c["a2"],ea=o["ea"],oa=o["oa"],
                         chai_beta=c["beta"]))
    # Chai beta is the effect of the reported risk allele; orient it to 1000G A1.
    # Re-parse risk-allele label (rs-ALLELE) from the API to get its identity.
    ra={}
    for d in pages:
        for a in d.get("_embedded",{}).get("associations",[]):
            lab=a.get("riskAllele",[{}])[0].get("label","")
            rs=a.get("riskAllele",[{}])[0].get("key")
            if rs and "-" in lab: ra[rs]=lab.rsplit("-",1)[1].upper()
    final=[]
    for r in recs:
        risk=ra.get(r["rs"])
        a1,a2=r["a1"],r["a2"]
        if risk is None or {risk} & {a1,a2} == set(): continue
        comp={"A":"T","T":"A","C":"G","G":"C"}
        # Chai beta is effect of 'risk' allele. Orient to a1 (1000G A1).
        if risk==a1: be=r["chai_beta"]
        elif risk==a2: be=-r["chai_beta"]
        elif comp.get(risk)==a1: be=r["chai_beta"]
        elif comp.get(risk)==a2: be=-r["chai_beta"]
        else: continue
        # knee-OA beta is per its effect allele; orient to 1000G A1 to align with Chai be
        # (comp_align needs an outcome beta; OA dict stored log-OR under r["bo"] earlier)
        bo=comp_align(a1,a2,r["ea"],r["oa"],r.get("bo"))
        if bo is None or be==0: continue
        final.append(dict(rs=r["rs"],be=be,se_e=r["se_e"],bo=bo,se_o=r["se_o"]))
    print(f"[328] harmonised to knee OA: {len(final)} independent instruments",flush=True)
    out=dict(design="T2-A independent instrument-set replication: Chai 2026 eBMD leads -> knee OA",
             exposure="Chai 2026 Nat Genet eBMD (GCST90726625, n=448010 EUR)",
             outcome="knee OA GCST90566800 (Hatzikotoulas 2025)",
             n_catalog_assocs=sum(len(d.get("_embedded",{}).get("associations",[])) for d in pages),
             n_leads=len(leads),n_mapped_1000g=len(mapped),n_independent=len(clumped),
             n_harmonised_kneeOA=len(final),
             scale_note="eBMD (heel ultrasound rank-normal) vs DXA aBMD differ in scale; direction compared, magnitude not")
    if len(final)>=10:
        out["mr_knee"]=run_mr(final,"Chai eBMD -> knee OA")
        out["verdict"]=verdict_for(out["mr_knee"])
    else:
        out["verdict"]=f"UNDERPOWERED: only {len(final)} independent instruments harmonised; no MR reported."
    json.dump(out,open(OUT,"w",encoding="utf-8"),ensure_ascii=False,indent=2)
    print(f"\n[328] {out['verdict']}")
    if "mr_knee" in out: print(json.dumps(out["mr_knee"],indent=2))
    print(f"[write] {OUT}")

if __name__=="__main__":
    import sys
    if "--reverdict" in sys.argv:
        _reverdict()
    else:
        main()
