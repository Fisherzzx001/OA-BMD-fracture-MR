#!/usr/bin/env python3
"""Fresh, independent re-derivation audit for OA/BMD/fracture manuscript.
Recomputes: N-back from tophits SEs; Lee-2011 liability h2 (z = pdf height);
MR OR + CI from table3_corrected.csv; correction-factor arithmetic; GEFOS gap.
Pure stdlib only (no scipy) using an erf-based normal CDF/PDF.
"""
import csv, math, os

def erf(x):
    # Abramowitz-Stegun 7.1.26
    s = 1 if x >= 0 else -1
    x = abs(x)
    t = 1/(1+0.3275911*x)
    y = 1 - (((((1.061405429*t -1.453152027)*t)+1.421413741)*t -0.284496736)*t +0.254829592)*t*math.exp(-x*x)
    return s*y

def Phi(x): return 0.5*(1+erf(x/math.sqrt(2)))
def phi(x): return math.exp(-0.5*x*x)/math.sqrt(2*math.pi)
def inv_Phi(p):
    # Acklam's rational approximation
    a=[-3.969683028665376e+01,2.209460984245205e+02,-2.759285104469687e+02,1.383577518672690e+02,-3.066479806614716e+01,2.506628277459239e+00]
    b=[-5.447609879822406e+01,1.615858368580409e+02,-1.556989798598866e+02,6.680131188771972e+01,-1.328068155288572e+01]
    c=[-7.784894002430293e-03,-3.223964580411365e-01,-2.400758277161838e+00,-2.549732539343734e+00,4.374664141464968e+00,2.938163982698783e+00]
    d=[7.784695709041462e-03,3.224671290700398e-01,2.445134137142996e+00,3.754408661907416e+00]
    plow, phigh = 0.02425, 1-0.02425
    if p < plow:
        q=math.sqrt(-2*math.log(p)); return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    if p > phigh:
        q=math.sqrt(-2*math.log(1-p)); return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    q=p-0.5; r=q*q
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q/(((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)

def liab_h2(h2_obs, K, P=None):
    """Lee 2011 Eq 23. z = pdf height at threshold T=inv_Phi(1-K)."""
    T = inv_Phi(1-K)
    z = phi(T)
    if P is None: P = K
    return h2_obs * (K*(1-K))/(z*z) * (K*(1-K))/(P*(1-P))

def n_back_se(se, eaf):
    """N implied by SE for quantitative trait: SE = 1/sqrt(N*2*EAF*(1-EAF))."""
    return 1.0/(se*se*2*eaf*(1-eaf))

def n_back_binary(se, eaf, K):
    """Effective N for binary: SE = 1/sqrt(Neff*2*EAF*(1-EAF)); Neff = N*K*(1-K)."""
    Neff = 1.0/(se*se*2*eaf*(1-eaf))
    return Neff/(K*(1-K)), Neff

print("="*70)
print("1. N-BACK FROM TOPHITS (raw data, independent)")
print("="*70)
base = "code/output/genetics"
for f, trait, kind in [("tophits_Total_BMD.csv","BMD_total","quant"),
                       ("tophits_fracture_any-type.csv","fracture","binary")]:
    rows=[]
    with open(os.path.join(base,f)) as fh:
        for r in csv.DictReader(fh):
            rows.append(r)
    ns=[]; neffs=[]
    for r in rows:
        se=float(r["standard_error"]); eaf=float(r["effect_allele_frequency"])
        z=float(r["beta"])/se
        if kind=="quant":
            n=n_back_se(se,eaf); ns.append(n)
        else:
            n,neff=n_back_binary(se,eaf,0.143); ns.append(n); neffs.append(neff)
    ns_s=sorted(ns)
    med=ns_s[len(ns_s)//2]
    print(f"  {trait}: n={len(ns)} SNPs; median N-back={med:,.0f}; range [{min(ns):,.0f}, {max(ns):,.0f}]")
    if neffs:
        print(f"      median Neff-back={ns_s and sorted(neffs)[len(neffs)//2]:,.0f}; implied N*K*(1-K)")
    # verify one p-value from Z
    r0=rows[0]; z0=float(r0["beta"])/float(r0["standard_error"])
    p_from_z=2*(1-Phi(abs(z0)))
    print(f"      top SNP {r0['rs_id']}: beta={r0['beta']} SE={r0['standard_error']} Z={z0:.2f} p(file)={r0['p_value']} p(from Z)={p_from_z:.2e}")

print()
print("="*70)
print("2. LEE-2011 LIABILITY h2 (z = pdf HEIGHT at threshold), independent")
print("="*70)
for name,h2,K in [("OA_knee",0.0354,0.12),("OA_selfrep",0.0223,0.083),("fracture",0.0258,0.143)]:
    h2l=liab_h2(h2,K)
    print(f"  {name}: h2_obs={h2} K={K} -> h2_liab={h2l:.4f}  (manuscript claims per trait)")

print()
print("="*70)
print("3. MR OR + 95% CI from table3_corrected.csv (strict clump)")
print("="*70)
with open(os.path.join(base,"table3_corrected.csv")) as fh:
    for r in csv.DictReader(fh):
        b=float(r["beta"]); se=float(r["se"])
        # OR per 1-SD LOWER exposure = exp(+b) for protective negative b
        lor_low=b  # log OR per 1 SD lower = b (since b is per 1 SD higher, negative)
        lo=(b-1.96*se); hi=(b+1.96*se)
        # per 1 SD lower BMD -> fracture: OR = exp(|b|) with CI exp(|hi|..|lo|) appropriately
        # b negative => per SD lower OR = exp(-b); CI = exp(-hi)..exp(-lo)
        if b<0:
            or_low=math.exp(-b); ci_lo=math.exp(-hi); ci_hi=math.exp(-lo)
            tag="per 1-SD LOWER"
        else:
            or_low=math.exp(b); ci_lo=math.exp(lo); ci_hi=math.exp(hi); tag="per 1-SD HIGHER"
        print(f"  {r['exp']:>14} -> {r['out']:<12} k={r['k']:>2} b={b:+.4f} OR{tag}={or_low:.3f} (95% CI {ci_lo:.3f}-{ci_hi:.3f}) p={float(r['p']):.1e}")

print()
print("="*70)
print("4. CORRECTION FACTOR ARITHMETIC")
print("="*70)
N_used=400000; N_true_cited=31873; N_true_tophits=29342
f_cited=math.sqrt(N_used/N_true_cited); f_top=math.sqrt(N_used/N_true_tophits)
print(f"  cited:  sqrt(400000/31873) = {f_cited:.4f}  -> divide factor 1/{f_cited:.4f} = {1/f_cited:.4f}")
print(f"  tophits: sqrt(400000/29342) = {f_top:.4f} -> divide factor 1/{f_top:.4f} = {1/f_top:.4f}")
print(f"  -> using tophits N would make corrected BMD->fracture b = -1.194 * {1/f_top:.4f} = {-1.194/f_top:.4f} (vs cited -0.337)")
print(f"  -> relative difference in corrected b: {(1/f_top)/(1/f_cited)-1:+.1%}")

print()
print("="*70)
print("5. GEFOS GAP (Trajanoska 2018 OR = 1.55, 95% CI 1.48-1.63)")
print("="*70)
for label,b in [("strict femoral-neck b=-0.35066 -> OR",math.exp(0.35066)),
                ("1-Mb femoral-neck b=-0.3654 -> OR",math.exp(0.3654)),
                ("total-body b=-0.33713 -> OR",math.exp(0.33713))]:
    gap=(1.55-b)/1.55*100
    print(f"  {label} {b:.3f}: gap vs GEFOS 1.55 = {gap:.1f}%  (overlap CI? {1.48 < b < 1.63})")
print()
print("DONE")
