"""Independent re-derivation / audit of OA_BMD_fracture key statistics (no scipy)."""
import numpy as np, math

def phi(x):  # standard normal pdf
    return math.exp(-0.5*x*x)/math.sqrt(2*math.pi)

print("="*70)
print("1. TRUE N FROM SE  (quant trait: SE = 1/sqrt(N * 2pq))")
print("="*70)
b_se, se, eaf = -0.1722, 0.0094, 0.7391
pq = 2*eaf*(1-eaf)
N_bmd = 1.0/(se**2 * pq)
print(f"BMD_total lead rs3779381: beta={b_se} se={se} eaf={eaf} 2pq={pq:.4f}")
print(f"  => N_implied = {N_bmd:,.0f}   (manuscript claims 31,873)")
N2 = 1.0/(0.0092**2 * 2*0.722*(1-0.722))
print(f"BMD_total rs3801387: => N_implied = {N2:,.0f}")

print("\n1b. FRACTURE N from SE (logistic: SE = 1/sqrt(N*K*(1-K)*2pq))")
fb, fse, feaf = -0.0792, 0.0084, 0.3693
neff = 1.0/(fse**2 * 2*feaf*(1-feaf))
print(f"fracture lead rs62476338: beta={fb} se={fse} eaf={feaf}")
print(f"  => N*K*(1-K) = {neff:,.0f}  (manuscript measured 31,605)")
for K in [0.10, 0.1429, 0.15]:
    print(f"     if K={K}: N_implied = {neff/(K*(1-K)):,.0f}")
print(f"  GEFOS K=0.1429 -> N = {neff/(0.1429*0.8571):,.0f}  (manuscript claims 264,973)")

print("\n"+"="*70)
print("2. LIABILITY h2  (Lee 2011 Eq23): h2_liab = h2_obs * K(1-K)/z^2 ; z=phi(Phi^-1(1-K))")
print("="*70)
# precomputed Phi^-1(1-K):
inv = {0.12:1.174986, 0.0831:1.38303, 0.1429:1.0660}
for name, h2, K, claim in [("OA_knee",0.0354452350956293,0.12,0.094),
                           ("OA_selfrep",0.02230416663074558,0.0831,0.073),
                           ("fracture",0.025767681247168977,0.1429,0.062)]:
    t = inv[K]; z = phi(t); fac = K*(1-K)/(z**2); hl = h2*fac
    print(f"{name:11s} h2_obs={h2:.5f} K={K} z={z:.4f} fac={fac:.4f} h2_liab={hl:.4f} (ms {claim}) {'OK' if abs(hl-claim)<0.005 else 'MISMATCH'}")

print("\n"+"="*70)
print("3. BMD MR CORRECTION FACTOR  beta_corrected = beta_placeholder * sqrt(N_true/N_used)")
print("="*70)
c = math.sqrt(31873/400000)
print(f"correction factor = {c:.4f}  (placeholder beta is {1/c:.3f}x too LARGE; corrected smaller)")
for exp, bu in [("BMD_total->fracture",-1.194319491050181),("BMD_femur_neck->fracture",-1.2422370063075012)]:
    bc = bu*c; OR = math.exp(-bc)
    print(f"{exp}: unc={bu:.4f} corrected={bc:.4f} OR(per SD lower)={OR:.3f}")
for bu, lab in [(-1.194319491050181,"total"),(-1.2422370063075012,"femurneck")]:
    OR = math.exp(-bu*c)
    print(f"  {lab}: OR={OR:.3f} vs GEFOS 1.55 -> {(1.55-OR)/1.55*100:.1f}% lower")

print("\n"+"="*70)
print("4. REVERSE MR (BMD outcome) should scale x 1/c = x{:.2f}".format(1/c))
print("="*70)
print(f"OA_knee->BMD_total uncorrected(full)=0.0300 -> corrected={0.0300*(1/c):.4f}")
print("  Table3/§3.5 report +0.048 (uncorrected OUTLIER-REMOVED) -> NOT on corrected scale")
print("  => forward MR corrected, reverse MR uncorrected: INCONSISTENT within Table 3.")
