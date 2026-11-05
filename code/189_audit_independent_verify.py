#!/usr/bin/env python3
# code/189 — INDEPENDENT re-verification of headline MR statistics.
# Route: manuscript Table 2/3  <->  code/output/genetics/AUDIT_mr_recomputed.csv
#        (produced by code/13_mr_analysis.py)  <->  closed-form re-derivation here.
# This is the "trust-but-verify" artefact: it recomputes OR/CI and power FROM THE
# deposited CSV and checks them against the numbers printed in docs/08, instead of
# only checking that docs/08 equals docs/08.
import csv, math
from statistics import NormalDist

def power_from_z(z, alpha=0.05):
    # two-sided achieved power for observed |z|
    za = NormalDist().inv_cdf(1 - alpha/2)
    return NormalDist().cdf(abs(z) - za) + NormalDist().cdf(-abs(z) - za)

rows = {}
with open("code/output/genetics/AUDIT_mr_recomputed.csv") as f:
    for r in csv.DictReader(f):
        rows[f"{r['exposure']}->{r['outcome']}"] = r

# manuscript Table 2 (IVW fixed) expected: pair -> (beta, se, p_exp, OR_exp, lo, hi)
EXP = {
 "BMD_total->fracture":      (-0.337, 0.022, 3.7e-51, 1.401, 1.341, 1.464),
 "BMD_femur_neck->fracture": (-0.351, 0.030, 7.7e-31, 1.420, 1.338, 1.507),
 "BMD_total->OA_knee":       ( 0.078, 0.013, 3.9e-10, 1.081, 1.055, 1.108),
 "BMD_lumbar->OA_knee":      ( 0.095, 0.016, 7.7e-9,  1.099, 1.064, 1.135),
 "BMD_arms->OA_knee":        ( 0.075, 0.016, 1.7e-6,  1.078, 1.045, 1.112),
 "BMD_leg->OA_knee":         ( 0.060, 0.014, 2.3e-5,  1.062, 1.033, 1.092),
 "BMD_femur_neck->OA_knee":  ( 0.058, 0.017, 6.4e-4,  1.060, 1.025, 1.095),
 "OA_knee->fracture":        ( 0.027, 0.032, 0.395,   0.973, 0.914, 1.036),
 "OA_knee->BMD_total":       ( 0.031, 0.032, 0.326,   None,  None,  None),
 "OA_knee->BMD_femur_neck":  ( 0.042, 0.032, 0.190,   None,  None,  None),
}

print(f"{'pair':28} {'beta':>9} {'ORcalc':>8} {'ORms':>7} {'CIcalc':>18} {'CImatch':>7} {'pow':>5} {'powMS':>6}")
bad = 0
for pair,(eb,ese,ep,eo,elo,ehi) in EXP.items():
    r = rows[pair]
    b=float(r["beta"]); se=float(r["se"]); p=float(r["p"])
    # closed-form OR for binary outcome = exp(|beta|) on log-odds scale
    orc = math.exp(abs(b))
    lo = math.exp(abs(b)-1.96*se); hi = math.exp(abs(b)+1.96*se)
    z = abs(float(r["z"]))
    pw = power_from_z(z)
    # manuscript power (Table 3) keyed separately
    orm = eo if eo else float('nan')
    ok_or = (eo is None) or (abs(orc-eo)/eo < 0.01)
    ok_ci = (elo is None) or (abs(lo-elo)/elo < 0.01 and abs(hi-ehi)/ehi < 0.01)
    ok_beta = abs(b-eb) < 1e-3
    ok_p = (p/ep < 3) and (ep/p < 3)  # within 3x (p tiny)
    flag = "OK" if (ok_or and ok_ci and ok_beta and ok_p) else "MISMATCH"
    if flag=="MISMATCH": bad+=1
    pwms = {"BMD_total->fracture":1.00,"BMD_femur_neck->fracture":1.00,"BMD_total->OA_knee":1.00,
            "BMD_lumbar->OA_knee":1.00,"BMD_arms->OA_knee":1.00,"BMD_leg->OA_knee":0.99,
            "BMD_femur_neck->OA_knee":0.93,"OA_knee->fracture":0.14,"OA_knee->BMD_total":0.17,
            "OA_knee->BMD_femur_neck":0.26}[pair]
    print(f"{pair:28} {b:9.4f} {orc:8.3f} {orm:7.3f} [{lo:.3f},{hi:.3f}] {str(ok_ci):>7} {pw:5.3f} {pwms:6.2f}  {flag}")

# I2 cross-check vs Table 2
print("\nI2 (manuscript -> csv):")
I2EXP={"BMD_total->fracture":78,"BMD_femur_neck->fracture":83,"BMD_total->OA_knee":81,
       "BMD_lumbar->OA_knee":83,"BMD_arms->OA_knee":86,"BMD_leg->OA_knee":85,
       "BMD_femur_neck->OA_knee":68,"OA_knee->fracture":49,"OA_knee->BMD_total":49,
       "OA_knee->BMD_femur_neck":30}
for pair,ei in I2EXP.items():
    ci2=int(round(float(rows[pair]["I2"])*100))
    print(f"  {pair:28} ms={ei:3d}%  csv={ci2:3d}%  {'OK' if abs(ci2-ei)<=1 else 'DIFF'}")

print(f"\nMISMATCH count (OR/CI/beta/p): {bad}")
print("Source CSV: code/output/genetics/AUDIT_mr_recomputed.csv  (upstream: code/13_mr_analysis.py)")
