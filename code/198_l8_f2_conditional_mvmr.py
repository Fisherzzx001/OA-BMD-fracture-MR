#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
198_l8_f2_conditional_mvmr.py  —  F2 geometry arm UPGRADED to a formal conditional test
===================================================================================
Upgrades the F2 "hip-shape-mode (HSM) -> OA" univariable screen into a BMD + HSM
multivariable (conditional) MR. For each HSM mode k we model

        OA = theta_BMD|HSM_k * BMD  +  theta_HSM_k|BMD * HSM_k   (+ Egger intercept)

using the HSM lead-variant instruments, with BMD_total looked up at those same SNPs
(GEFOS tidy_femur_total.tsv). All three columns are oriented to the OA effect allele.

Why this is the formal F2 test:
  * theta_BMD|HSM_k  answers "does BMD's causal effect on OA survive after accounting
    for hip geometry?"  -> if it stays significant & near the marginal BMD->OA estimate,
    geometry is NOT the mediator (BMD acts via mass/formation = SOS 'formation side').
  * theta_HSM_k|BMD  answers "does hip shape still matter beyond BMD?"

Estimators: IVW-MVMR (Sanderson & Windmeijer 2015) + Egger-MVMR (intercept = horizontal
pleiotropy). Only the OUTCOME SE is used for weighting, so HSM's p-value-derived SE is
irrelevant to estimation. Robustness: per-mode Cochran-Q heterogeneity and leave-one-out.

Inputs (all pre-existing; no new downloads):
  OA outcome (hip/knee):  G:/OA_B_PROJECT_cache/l7/hsm_hip_out.tsv , hsm_knee_out.tsv
  BMD exposure:           .../OA_FLAGSHIP/07_reference/ldsc/tidy_femur_total.tsv
  HSM exposure (per mode): .../OA_B_PROJECT/code/output/genetics/L7_hsm_associations/*.tsv
Output:
  code/output/genetics/L8_F2_conditional_mvmr.json
"""
import os, csv, math, json, glob
import numpy as np

ROOT = "C:/Users/Administrator/WorkBuddy"
PROJ = ROOT + "/2026-09-08-22-17-40/OA_B_PROJECT"
GEN  = PROJ + "/code/output/genetics"
HSM_DIR = GEN + "/L7_hsm_associations"
BMD_TIDY = ROOT + "/2026-09-08-20-17-24/OA_FLAGSHIP/07_reference/ldsc/tidy_femur_total.tsv"
OA_HIP = "G:/OA_B_PROJECT_cache/l7/hsm_hip_out.tsv"
OA_KNEE = "G:/OA_B_PROJECT_cache/l7/hsm_knee_out.tsv"
OUT = GEN + "/L8_F2_conditional_mvmr.json"
Z975 = 1.959963985

COMP = {"A":"T","T":"A","C":"G","G":"C"}
PALIN = {frozenset(("A","T")), frozenset(("C","G"))}

def ncdf(x): return 0.5*(1+math.erf(x/math.sqrt(2)))
def p_two(b,s): return 2*(1-ncdf(abs(b/s)))
def inv_norm(p):
    a=[-3.969683028665376e+01,2.209460984245205e+02,-2.759285104469687e+02,1.383577518672690e+02,-3.066479806614716e+01,2.506628277459239e+00]
    b=[-5.447609879822406e+01,1.615858368580409e+02,-1.556989798598866e+02,6.680131188771972e+01,-1.328068155288572e+01]
    c=[-7.784894002430293e-03,-3.223964580411365e-01,-2.400758277161838e+00,-2.549732539343734e+00,4.374664141464968e+00,2.938163982698783e+00]
    d=[7.784695709041462e-03,3.224671290700398e-01,2.445134137142996e+00,3.754408661907416e+00]
    p=min(max(p,1e-300),1-1e-16)
    if p<0.02425:
        q=math.sqrt(-2*math.log(p)); return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    if p<=1-0.02425:
        q=p-0.5; r=q*q; return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q/(((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)
    q=math.sqrt(-2*math.log(1-p)); return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)

# ---------------------------------------------------------------- load OA outcome
def load_oa(path):
    d={}
    with open(path) as f:
        for line in f:
            p=line.rstrip("\n").split("\t")
            if len(p)<14: continue
            rs=p[0]; E=p[3].upper(); O=p[4].upper()
            try:
                OR=float(p[5]); cu=float(p[10]); cl=float(p[11]); P=float(p[8])
            except ValueError:
                continue
            if OR<=0 or cu<=0 or cl<=0: continue
            ln_or=math.log(OR); lo=math.log(cl); hi=math.log(cu)
            se=(hi-lo)/(2*Z975) if hi>lo else abs(ln_or)/inv_norm(1-P/2)
            d[rs]=dict(E=E,O=O,beta_o=ln_or,se_o=se,P=P,ncase=float(p[12]),nctrl=float(p[13]))
    return d

# ---------------------------------------------------------------- load BMD (GEFOS total)
def load_bmd(path):
    d={}
    with open(path) as f:
        r=csv.DictReader(f,delimiter="\t")
        for row in r:
            d[row["SNP"].strip()]={"A1":row["A1"].upper(),"BETA":float(row["BETA"]),"SE":float(row["SE"])}
    return d

# ---------------------------------------------------------------- orient exp beta to OA E allele
def orient_exp(beta, exp_allele, E, O):
    """Return beta oriented to OA effect allele E (None if unresolvable)."""
    a=exp_allele.upper()
    if PALIN and frozenset((E,O)) in PALIN:   # OA itself palindromic -> cannot resolve
        return None
    if a==E: return beta
    if a==O: return -beta
    if COMP.get(a)==E: return -beta            # exp allele is complement of E -> strand flip
    if COMP.get(a)==O: return beta
    return None

# ---------------------------------------------------------------- MVMR estimators
def ivw_mvmr(X, y, w):
    """X: (n,2) [bmd, hsm]; y: (n,); w: (n,) weights=1/se_o^2. returns theta, se, Q, df."""
    W=np.diag(w); XtW=X.T@W; M=XtW@X
    try:
        Minv=np.linalg.inv(M)
    except np.linalg.LinAlgError:
        return (np.nan,)*2, np.nan, len(y)-2
    theta=Minv@(XtW@y)
    se=np.sqrt(np.diag(Minv))
    resid=y-X@theta
    Q=float(resid.T@(W@resid))
    return theta, se, Q, len(y)-2

def egger_mvmr(X, y, w):
    Xe=np.column_stack([np.ones(len(y)), X])   # [1, bmd, hsm]
    W=np.diag(w); XtW=Xe.T@W; M=XtW@Xe
    try: Minv=np.linalg.inv(M)
    except np.linalg.LinAlgError: return (np.nan,)*3,(np.nan,)*3
    theta=Minv@(XtW@y); se=np.sqrt(np.diag(Minv))
    return theta, se   # [intercept, bmd|hsm, hsm|bmd]

def cond_F(X):
    """Conditional F for column 0 (BMD) given column 1 (HSM)."""
    bmd=X[:,0]; hsm=X[:,1]; n=len(bmd)
    if n<3: return np.nan
    r=(np.corrcoef(bmd,hsm)[0,1])**2
    k=2
    F=((r)/(k-1))/((1-r)/(n-k)) if (1-r)>0 else np.inf
    return F

# ---------------------------------------------------------------- load everything
oa_hip=load_oa(OA_HIP); oa_knee=load_oa(OA_KNEE)
bmd=load_bmd(BMD_TIDY)
print(f"OA hip SNPs={len(oa_hip)}  OA knee SNPs={len(oa_knee)}  BMD SNPs={len(bmd)}")

def build_mode(acc, oa):
    """Build MVMR dataset for one HSM mode against one OA outcome."""
    fn=os.path.join(HSM_DIR, acc+".tsv")
    rows=[]
    with open(fn, encoding="utf-8", errors="replace") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            try: p=float(row["P-VALUE"]); mag=float(row["OR or BETA"])
            except: continue
            if p>=5e-9 or mag==0: continue
            sra=row["STRONGEST SNP-RISK ALLELE"]; allele=sra.split("-",1)[1] if "-" in sra else ""
            ci=row.get("95% CI (TEXT)","") or ""
            dec="decrease" in ci.lower()
            beta_hsm=-mag if dec else mag
            rs=row["SNPS"].strip()
            if rs not in oa or rs not in bmd: continue
            E=oa[rs]["E"]; O=oa[rs]["O"]
            bo=orient_exp(beta_hsm, allele, E, O)
            bb=orient_exp(bmd[rs]["BETA"], bmd[rs]["A1"], E, O)
            if bo is None or bb is None: continue
            rows.append(dict(rs=rs, bmd=bb, hsm=bo, o=oa[rs]["beta_o"], se=oa[rs]["se_o"],
                             E=E, O=O, ra_hsm=allele, A1_bmd=bmd[rs]["A1"]))
    return rows

def analyse(rows):
    if len(rows)<4:
        return dict(n=len(rows), note="insufficient instruments for stable MVMR (need >=4)")
    X=np.array([[r["bmd"], r["hsm"]] for r in rows])
    y=np.array([r["o"] for r in rows]); w=1.0/np.array([r["se"] for r in rows])**2
    theta,se,Q,df=ivw_mvmr(X,y,w)
    theta_e,se_e=egger_mvmr(X,y,w)
    # heterogeneity
    I2=max(0.0,(Q-df)/Q)*100 if Q>0 else 0.0
    p_Q=2*(1-ncdf(math.sqrt(Q)))  # rough upper-tail chi2(df) via norm approx for large df
    # leave-one-out on theta_BMD|HSM
    loo=[]
    for i in range(len(rows)):
        m=np.ones(len(rows),bool); m[i]=False
        t,s,Q2,df2=ivw_mvmr(X[m],y[m],w[m])
        loo.append(t[0])
    # marginal contrasts (same SNPs)
    bmd_only=np.array([r["bmd"] for r in rows]); hsm_only=np.array([r["hsm"] for r in rows])
    def univ(exp):
        num=sum(exp*r["o"]/r["se"]**2 for exp,r in zip(exp,rows))
        den=sum(e*e/r["se"]**2 for e,r in zip(exp,rows))
        b=num/den; s=1/math.sqrt(den); return b,s
    b_u_bmd,s_u_bmd=univ(bmd_only)
    b_u_hsm,s_u_hsm=univ(hsm_only)
    cf=cond_F(X)
    return dict(
        n=len(rows),
        theta_BMD_given_HSM=float(theta[0]), se_BMD=float(se[0]),
        p_BMD=p_two(theta[0],se[0]),
        theta_HSM_given_BMD=float(theta[1]), se_HSM=float(se[1]),
        p_HSM=p_two(theta[1],se[1]),
        egger_intercept=float(theta_e[0]), se_intercept=float(se_e[0]), p_intercept=p_two(theta_e[0],se_e[0]),
        egger_BMD=float(theta_e[1]), egger_HSM=float(theta_e[2]),
        Q=float(Q), Q_df=df, I2=float(I2),
        cond_F_BMD_given_HSM=float(cf),
        univ_BMD=float(b_u_bmd), univ_BMD_se=float(s_u_bmd), univ_BMD_p=p_two(b_u_bmd,s_u_bmd),
        univ_HSM=float(b_u_hsm), univ_HSM_se=float(s_u_hsm), univ_HSM_p=p_two(b_u_hsm,s_u_hsm),
        loo_BMD_given_HSM=[float(x) for x in loo],
        loo_min=min(loo), loo_max=max(loo),
        snps=[r["rs"] for r in rows])

# ---------------------------------------------------------------- run per mode
modes=["GCST90482704","GCST90482705","GCST90482706","GCST90482707","GCST90482708",
       "GCST90482709","GCST90482710","GCST90482711","GCST90482712","GCST90482713"]
mode_names={m:"HSM%d"%(i+1) for i,m in enumerate(modes)}

result={"hip":{}, "knee":{}, "meta":{}}
for outcome,oa in (("hip",oa_hip),("knee",oa_knee)):
    print("\n"+"="*100); print("OUTCOME =",outcome.upper(),"OA"); print("="*100)
    agg_BMD=[]; agg_HSM=[]
    for m in modes:
        rows=build_mode(m,oa)
        a=analyse(rows)
        result[outcome][m]=a
        if "theta_BMD_given_HSM" in a:
            agg_BMD.append((a["theta_BMD_given_HSM"],a["se_BMD"],a["n"]))
            agg_HSM.append((a["theta_HSM_given_BMD"],a["se_HSM"],a["n"]))
            b=a
            print(f"{mode_names[m]:5s}: n={b['n']:2d}  "
                  f"theta_BMD|HSM={b['theta_BMD_given_HSM']:+.4f}(SE {b['se_BMD']:.4f},p={b['p_BMD']:.2g})  "
                  f"theta_HSM|BMD={b['theta_HSM_given_BMD']:+.4f}(p={b['p_HSM']:.2g})  "
                  f"I2={b['I2']:.0f}%  condF={b['cond_F_BMD_given_HSM']:.1f}  "
                  f"univBMD={b['univ_BMD']:+.4f}(p={b['univ_BMD_p']:.2g})")
            print(f"        Egger int p={b['p_intercept']:.2g}  LLOO_BMD|HSM=[{b['loo_min']:+.3f},{b['loo_max']:+.3f}]")
        else:
            print(f"{mode_names[m]:5s}: {a.get('note')}")
    # fixed-effect meta across modes
    def meta(agg):
        wts=[1/s**2 for _,s,_ in agg]; b=sum(wi*bi for wi,(bi,si,_) in zip(wts,agg))/sum(wts)
        se=1/math.sqrt(sum(wts)); return b,se
    if agg_BMD:
        mb,sb=meta(agg_BMD); result["meta"][outcome+"_BMD_given_HSM"]=dict(
            theta=float(mb),se=float(sb),p=p_two(mb,sb),OR=math.exp(mb),lo=math.exp(mb-Z975*sb),hi=math.exp(mb+Z975*sb),n_modes=len(agg_BMD))
        mh,sh=meta(agg_HSM); result["meta"][outcome+"_HSM_given_BMD"]=dict(
            theta=float(mh),se=float(sh),p=p_two(mh,sh),OR=math.exp(mh),lo=math.exp(mh-Z975*sh),hi=math.exp(mh+Z975*sh),n_modes=len(agg_HSM))
        print(f"\n  META theta_BMD|HSM = {mb:+.4f} (SE {sb:.4f}, p={p_two(mb,sb):.2g}, OR={math.exp(mb):.3f} [{math.exp(mb-Z975*sb):.3f},{math.exp(mb+Z975*sb):.3f}])")
        print(f"  META theta_HSM|BMD = {mh:+.4f} (p={p_two(mh,sh):.2g})")

json.dump(result, open(OUT,"w"), indent=2)
print("\nSaved", OUT)
