#!/usr/bin/env python3
"""L7 F2 (geometry arm): hip-shape-mode (HSM) genetic instruments -> knee OA and hip OA.

The SOS hypothesis reads areal BMD as a systemic osteogenic set-point. Areal BMD fuses
mass, size and geometry, so the BMD->OA effect could in principle be transmitted by bone
geometry rather than by mass/formation. Height indexes stature but not bone width or
geometry, so the manuscript leaves this split open ("resolved only by test F2").

HSM1-HSM10 (Faber et al. 2025, Hum Mol Genet 34:207-217; PMID 39574169) are statistical
hip shape modes derived from DXA imaging; they encode geometry only and contain no
density information. A null HSM->OA effect therefore identifies the BMD->OA effect as
mass/density-transmitted (formation side) rather than geometry-transmitted.

Exposure effects come from the GWAS Catalog association table (author-reported beta with
direction word); outcome effects are extracted from the Hatzikotoulas 2025 summary files.
"""
import csv, math, json, re, os

HARM = "C:/Users/Administrator/WorkBuddy/2026-09-08-22-17-40/OA_B_PROJECT/code/output/genetics/harmonised/BMD_total__OA_knee.csv"
D = "G:/OA_B_PROJECT_cache/l7"
OUT = os.path.join(D, "l7_f2_hsm.json")
Z975 = 1.959963985

def ncdf(x):
    return 0.5 * (1 + math.erf(x / math.sqrt(2)))

def p_two(b, se):
    return 2 * (1 - ncdf(abs(b / se)))

def inv_norm(p):
    """Acklam-style rational approximation of the normal quantile function."""
    a=[-3.969683028665376e+01,2.209460984245205e+02,-2.759285104469687e+02,1.383577518672690e+02,-3.066479806614716e+01,2.506628277459239e+00]
    b=[-5.447609879822406e+01,1.615858368580409e+02,-1.556989798598866e+02,6.680131188771972e+01,-1.328068155288572e+01]
    c=[-7.784894002430293e-03,-3.223964580411365e-01,-2.400758277161838e+00,-2.549732539343734e+00,4.374664141464968e+00,2.938163982698783e+00]
    d=[7.784695709041462e-03,3.224671290700398e-01,2.445134137142996e+00,3.754408661907416e+00]
    pl=0.02425
    p=min(max(p,1e-300),1-1e-16)
    if p<pl:
        q=math.sqrt(-2*math.log(p))
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    if p<=1-pl:
        q=p-0.5; r=q*q
        return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q/(((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)
    q=math.sqrt(-2*math.log(1-p))
    return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)

# ---------- BMD instruments (for overlap check) ----------
bmd_rs = set()
with open(HARM) as f:
    for row in csv.DictReader(f):
        bmd_rs.add(row["SNP"])

# ---------- HSM exposure effects ----------
HSM = {}
for fn in sorted(os.listdir(os.path.join(D, "hsm_assoc"))):
    if not fn.endswith(".tsv"):
        continue
    acc = fn.split(".")[0]
    with open(os.path.join(D, "hsm_assoc", fn), encoding="utf-8", errors="replace") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("STUDY ACCESSION") != acc:
                continue
            try:
                p = float(row["P-VALUE"])
                mag = float(row["OR or BETA"])
            except Exception:
                continue
            if p >= 5e-9 or mag == 0:
                continue
            sra = row["STRONGEST SNP-RISK ALLELE"]
            allele = sra.split("-", 1)[1] if "-" in sra else ""
            ci_txt = row.get("95% CI (TEXT)", "") or ""
            lo_s, hi_s = None, None
            m = re.search(r"\[([0-9.eE+-]+)\s*-\s*([0-9.eE+-]+)\]", ci_txt)
            if m:
                lo_s, hi_s = float(m.group(1)), float(m.group(2))
            dec = "decrease" in ci_txt.lower()
            beta = -mag if dec else mag
            # SE: from p-value (author p rounded to 1 s.f. -> z precision ~1-2%) ; CI as cross-check
            z = inv_norm(1 - p / 2)
            se_p = abs(beta) / z if z > 0 else float("nan")
            se_ci = (abs(hi_s - lo_s) / (2 * 1.959963985)) if (lo_s is not None and hi_s is not None) else float("nan")
            HSM.setdefault(acc, []).append(dict(
                rs=row["SNPS"], allele=allele, beta=beta, se=se_p, se_ci=se_ci, p=p,
                raf=float(row["RISK ALLELE FREQUENCY"]) if row.get("RISK ALLELE FREQUENCY") else None,
                gene=row.get("MAPPED_GENE", ""), trait=row.get("DISEASE/TRAIT", "")))
        #

# ---------- outcome effects ----------
def load_out(path):
    d = {}
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            rs = p[0]
            OR = float(p[5]); cu = float(p[10]); cl = float(p[11])
            d[rs] = dict(E=p[3], O=p[4], ln_or=math.log(OR),
                         se=(math.log(cu) - math.log(cl)) / 3.92,
                         eaf=float(p[7]), ncase=float(p[12]), nctrl=float(p[13]))
    return d

hip = load_out(os.path.join(D, "hsm_hip_out.tsv"))
knee = load_out(os.path.join(D, "hsm_knee_out.tsv"))

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
PALIN = {frozenset(("A", "T")), frozenset(("C", "G"))}

def orient(out, allele):
    """Align the outcome effect to the HSM effect allele.
    Returns (beta_o, se_o, flipped) or (None, None, None) if the alleles cannot be aligned.
    Palindromic (AT/CG) SNPs are rejected because the strand is not resolvable from
    forward-strand alleles alone."""
    a = allele.upper()
    e = out["E"].upper(); o = out["O"].upper()
    if PALIN and frozenset((e, o)) in PALIN:
        return None, None, None
    if e == a:
        return out["ln_or"], out["se"], False
    if o == a:
        return -out["ln_or"], out["se"], True
    if COMP.get(e) == a:
        # the outcome file is on the opposite strand: the exposure allele is the
        # complement of the outcome's effect allele -> effect reverses sign
        return -out["ln_or"], out["se"], True
    if COMP.get(o) == a:
        return out["ln_or"], out["se"], False
    return None, None, None

def ivw(rows):
    num = sum(r["beta_o"] * r["beta_e"] / r["se_o"]**2 for r in rows)
    den = sum(r["beta_e"]**2 / r["se_o"]**2 for r in rows)
    b = num / den; se = 1 / math.sqrt(den)
    return b, se

def egger(rows):
    xs=[r["beta_e"] for r in rows]; ys=[r["beta_o"] for r in rows]; ws=[1/r["se_o"]**2 for r in rows]
    sw=sum(ws); swx=sum(w*x for w,x in zip(ws,xs)); swxx=sum(w*x*x for w,x in zip(ws,xs))
    swy=sum(w*y for w,y in zip(ws,ys)); swxy=sum(w*x*y for w,x,y in zip(ws,xs,ys))
    dd=sw*swxx-swx*swx
    if dd==0: return (float("nan"),)*4
    b=(sw*swxy-swx*swy)/dd; ic=(swxx*swy-swx*swxy)/dd
    return b, math.sqrt(sw/dd), ic, math.sqrt(swxx/dd)

def wmedian(rows):
    items=sorted(rows, key=lambda r: r["beta_o"]/r["beta_e"])
    wts=[(r["beta_e"]/r["se_o"])**2 for r in items]; tot=sum(wts); cum=0.0
    for r,w in zip(items,wts):
        cum+=w
        if cum>=tot/2: return r["beta_o"]/r["beta_e"]
    return items[-1]["beta_o"]/items[-1]["beta_e"]

def run_mode(acc, out, label):
    rows=[]; skipped=[]
    for x in HSM.get(acc, []):
        if x["rs"] not in out:
            skipped.append(x["rs"]); continue
        b_o, s_o, ok = orient(out[x["rs"]], x["allele"])
        if b_o is None:
            skipped.append(x["rs"]); continue
        rows.append(dict(rs=x["rs"], beta_e=x["beta"], se_e=x["se"], beta_o=b_o, se_o=s_o))
    if len(rows) < 3:
        return dict(acc=acc, trait=HSM[acc][0]["trait"] if HSM.get(acc) else "", n=len(rows),
                    note="insufficient instruments", skipped=skipped)
    b, se = ivw(rows)
    Q = sum((r["beta_o"] - b*r["beta_e"])**2/r["se_o"]**2 for r in rows)
    df = len(rows)-1
    I2 = max(0.0,(Q-df)/Q)*100 if Q>0 else 0.0
    se_r = se*math.sqrt(max(1.0, Q/df)) if df>0 else se
    be, seb, ic, seic = egger(rows)
    Fs = [(r["beta_e"]/r["se_e"])**2 for r in rows]
    return dict(acc=acc, trait=HSM[acc][0]["trait"], n=len(rows), skipped=skipped,
                ivw=dict(beta=b, se=se, p=p_two(b,se), OR=math.exp(b),
                         lo=math.exp(b-Z975*se), hi=math.exp(b+Z975*se),
                         p_random=p_two(b,se_r), se_random=se_r),
                Q=Q, Q_df=df, I2=I2,
                egger=dict(beta=be, se=seb, p=p_two(be,seb), intercept=ic, p_int=p_two(ic,seic)),
                wmedian=dict(beta=wmedian(rows), OR=math.exp(wmedian(rows))),
                F_mean=sum(Fs)/len(Fs), F_min=min(Fs),
                overlap_bmd=[r["rs"] for r in rows if r["rs"] in bmd_rs])

print("="*100)
print("F2 GEOMETRY ARM :  hip shape modes (pure geometry, no density) -> OA")
print("="*100)
summary=[]
for acc in sorted(HSM):
    h = run_mode(acc, hip, "hip")
    k = run_mode(acc, knee, "knee")
    summary.append((acc, h, k))
    print(f"\n{acc}  {h.get('trait','')}   [n_hip={h['n']}  n_knee={k['n']}]")
    for nm, r in (("HIP ", h), ("KNEE", k)):
        if "ivw" not in r:
            print(f"   {nm}: {r.get('note')}"); continue
        i=r["ivw"]
        print(f"   {nm}: IVW beta={i['beta']:+.4f} SE={i['se']:.4f} p={i['p']:.3g} "
              f"OR={i['OR']:.3f} [{i['lo']:.3f},{i['hi']:.3f}]  p_rand={i['p_random']:.3g} "
              f"Q={r['Q']:.1f}({r['Q_df']})p={p_two(r['Q']-r['Q_df'], math.sqrt(2*r['Q_df'])):.2f} I2={r['I2']:.0f}% "
              f"Egger b={r['egger']['beta']:+.4f}(p={r['egger']['p']:.2g}) int p={r['egger']['p_int']:.2g} "
              f"WM={r['wmedian']['beta']:+.4f} F_mean={r['F_mean']:.0f}")
        if r["overlap_bmd"]:
            print(f"        overlap with BMD instruments: {r['overlap_bmd']}")

# ---------- aggregate across modes ----------
print("\n"+"="*100)
print("F2 AGGREGATE (meta-analysis across the 10 HSM modes, log-odds scale)")
print("="*100)
agg={}
for nm, key in (("hip",1),("knee",2)):
    bs=[]; ses=[]
    for acc,h,k in summary:
        r = h if key==1 else k
        if "ivw" not in r: continue
        bs.append(r["ivw"]["beta"]); ses.append(r["ivw"]["se"])
    w=[1/s**2 for s in ses]
    b=sum(wi*bi for wi,bi in zip(w,bs))/sum(w)
    se=1/math.sqrt(sum(w))
    # heterogeneity across modes
    Q=sum(wi*(bi-b)**2 for wi,bi in zip(w,bs)); df=len(bs)-1
    I2=max(0.0,(Q-df)/Q)*100 if Q>0 else 0
    # Bonferroni across 10 modes
    bonf=0.05/10
    agg[nm]=dict(n_modes=len(bs), beta=b, se=se, p=p_two(b,se), OR=math.exp(b),
                 lo=math.exp(b-Z975*se), hi=math.exp(b+Z975*se),
                 Q=Q, Q_df=df, I2=I2, bonferroni_alpha=bonf,
                 max_abs_beta=max(abs(x) for x in bs), min_p=min(p_two(x,s) for x,s in zip(bs,ses)),
                 any_significant_bonf=any(p_two(x,s)<bonf for x,s in zip(bs,ses)))
    a=agg[nm]
    print(f"{nm:5s}: pooled beta={b:+.4f} SE={se:.4f} p={a['p']:.3g} OR={a['OR']:.3f} [{a['lo']:.3f},{a['hi']:.3f}]")
    print(f"        across-mode heterogeneity Q={Q:.2f} ({df} df) I2={I2:.0f}%; "
          f"max |beta_mode| = {a['max_abs_beta']:.4f}; smallest per-mode p = {a['min_p']:.3g}; "
          f"any mode significant at Bonferroni alpha=0.005: {a['any_significant_bonf']}")

json.dump(dict(modes=[dict(acc=a, hip=h, knee=k) for a,h,k in summary], aggregate=agg),
          open(OUT,"w"), indent=2)
print("\nSaved", OUT)
