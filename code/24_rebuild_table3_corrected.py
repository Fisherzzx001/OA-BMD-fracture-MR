"""
code/24_rebuild_table3_corrected.py
Rebuild Table 3 of the manuscript on the CORRECTED, STRICT-clump effect scale.

Principles (re-derived, scale-invariance verified):
- Forward MR (BMD exposure, placeholder N=400,000): reconstructed beta is
  underestimated by sqrt(31873/400000)=0.2823; IVW estimate is inflated 3.5426x.
  Correction = multiply IVW beta, SE, Egger beta, weighted median by 0.2823.
  p, Egger-intercept p, I2, LOO are SCALE-INVARIANT -> keep as-is.
- Reverse MR (BMD outcome, placeholder N): outcome beta underestimated 0.2823x
  -> IVW estimate underestimated 0.2823x -> correction = multiply by 3.5426
  (beta, SE, Egger beta, weighted median all x3.5426; p invariant).
- OA_knee->fracture: both OA_knee (true per-SNP N) and fracture (GEFOS true N)
  are on true scale -> NO correction; use strict values directly.

Sources:
- mr_corrected_N.csv  : STRICT clump, CORRECTED IVW beta/se/p/OR (authoritative core)
- mr_results.csv      : 1-Mb clump, UNCORRECTED full stats (reliable for
                         Egger-intercept p, I2, LOO; scale-correct Egger beta/wmed)
Random-effect p is recomputed via a DerSimonian-Laird style approximation from
fixed SE and I2 (verified to reproduce the known forward random p ~1e-12).

Outputs: code/output/genetics/table3_corrected.csv  + printed markdown table.
"""
import csv, math

OUT = "code/output/genetics"
FWD = 0.282307   # forward correction factor (sqrt(31873/400000))
REV = 3.54261    # reverse correction factor (1/FWD)

def z_to_p(z):
    # two-sided normal p
    return 2.0 * (1.0 - 0.5 * (1.0 + math.erf(abs(z) / math.sqrt(2.0))))

def dl_random_p(beta, se_fixed, i2):
    # DerSimonian-Laird approximation: tau2 = I2/(1-I2) * typical within-study var
    if i2 <= 0:
        i2 = 1e-6
    s2 = se_fixed ** 2
    tau2 = (i2 / (1.0 - i2)) * s2
    se_re = math.sqrt(s2 + tau2)
    z = beta / se_re
    return z_to_p(z), se_re

# load corrected strict core
core = {}
with open(f"{OUT}/mr_corrected_N.csv") as f:
    for r in csv.DictReader(f):
        core[(r["exposure"], r["outcome"])] = r

# load uncorrected full (reliable invariant + scale-correctable)
full = {}
with open(f"{OUT}/mr_results.csv") as f:
    for r in csv.DictReader(f):
        full[(r["exposure"], r["outcome"])] = r

# pair order for Table 3
order = [
    ("BMD_total", "fracture"),
    ("BMD_femur_neck", "fracture"),
    ("BMD_total", "OA_knee"),
    ("BMD_femur_neck", "OA_knee"),
    ("BMD_lumbar", "OA_knee"),
    ("BMD_arms", "OA_knee"),
    ("BMD_leg", "OA_knee"),
    ("OA_knee", "fracture"),
    ("OA_knee", "BMD_total"),
    ("OA_knee", "BMD_femur_neck"),
]

rows = []
for (exp, out) in order:
    c = core[(exp, out)]
    fl = full[(exp, out)]
    k = int(c["K"])
    beta = float(c["beta"]); se = float(c["se"]); p = float(c["p"])
    # Egger intercept p, I2, LOO: SCALE-INVARIANT, take from uncorrected full
    egger_int_p = float(fl["egger_intercept_p"])
    i2 = float(fl["I2"])
    loo = int(fl["loo_sign_flips"])
    # Egger beta + weighted median: scale-correct
    if exp.startswith("BMD"):
        f = FWD
    elif out.startswith("BMD"):
        f = REV
    else:
        f = 1.0  # OA_knee -> fracture, both true scale
    egger_beta = float(fl["egger"]) * f
    wmed = float(fl["wmed"]) * f
    # random-effect p
    p_re, _ = dl_random_p(beta, se, i2)
    rows.append({
        "exp": exp, "out": out, "k": k, "beta": beta, "se": se, "p": p,
        "p_re": p_re, "egger_beta": egger_beta, "egger_int_p": egger_int_p,
        "wmed": wmed, "i2": i2, "loo": loo,
    })

# write csv
fields = ["exp","out","k","beta","se","p","p_re","egger_beta","egger_int_p","wmed","i2","loo"]
with open(f"{OUT}/table3_corrected.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields); w.writeheader()
    for r in rows: w.writerow(r)

# pretty print
def fmtp(p):
    if p < 1e-3:
        e = int(math.floor(math.log10(p)))
        return f"{p:.1e}".replace("e","×10^").replace("^-0","^-").replace("^-","^-")
    return f"{p:.3f}"

print("| Exposure -> Outcome | n_snp | IVW beta | SE | p | IVW p(RE) | Egger int p | Wmedian | I2 | LOO |")
print("|---|---|---|---|---|---|---|---|---|---|")
for r in rows:
    print(f"| {r['exp']} -> {r['out']} | {r['k']} | {r['beta']:+.4f} | {r['se']:.4f} | "
          f"{fmtp(r['p'])} | {fmtp(r['p_re'])} | {r['egger_int_p']:.3f} | {r['wmed']:+.4f} | "
          f"{r['i2']*100:.0f}% | {r['loo']} |")
print("\nWrote code/output/genetics/table3_corrected.csv")
