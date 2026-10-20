"""Independent IVW recomputation (does NOT import the project MR module).
Goal: verify the BMD_total -> OA_knee headline effect from the validated
harmonised CSV, and compare against (a) the manuscript Table 2 and
(b) the stale mr_results.csv artifact.
"""
import csv, math
from statistics import NormalDist

PATH = "code/output/genetics/harmonised/BMD_total__OA_knee.csv"

rows = []
with open(PATH) as f:
    r = csv.DictReader(f)
    for row in r:
        rows.append(row)

be = [float(x["beta_e"]) for x in rows]
se_e = [float(x["se_e"]) for x in rows]
bo = [float(x["beta_o"]) for x in rows]
so = [float(x["se_o"]) for x in rows]
n = len(rows)

# IVW (inverse-variance weighted on outcome): beta = sum(be*bo/so^2) / sum(be^2/so^2)
w = [be[i]**2 / so[i]**2 for i in range(n)]
beta = sum(be[i]*bo[i]/so[i]**2 for i in range(n)) / sum(w)
se = 1.0 / math.sqrt(sum(w))
z = beta / se
p = 2 * (1 - NormalDist().cdf(abs(z)))
or_per_sd = math.exp(beta)
lo = math.exp(beta - 1.96*se)
hi = math.exp(beta + 1.96*se)

print(f"n_snp                = {n}")
print(f"IVW beta (per SD)    = {beta:.5f}")
print(f"IVW SE              = {se:.5f}")
print(f"z                    = {z:.3f}")
print(f"two-sided p          = {p:.2e}")
print(f"OR per SD higher BMD = {or_per_sd:.3f} (95% CI {lo:.3f}-{hi:.3f})")
print()
print("COMPARISON")
print(f"  manuscript Table2  : beta=+0.023, SE=0.004, p=3.9e-10, OR 1.023 (1.016-1.031)")
print(f"  mr_results.csv     : beta=+0.0819, SE=0.0131, p=3.9e-10, OR 1.085")
print(f"  this independent    : beta={beta:.5f}, SE={se:.5f}, p={p:.2e}, OR {or_per_sd:.3f}")
match_ms = abs(beta - 0.023) < 0.005 and p < 1e-8
match_csv = abs(beta - 0.0819) < 0.01
print()
print(f"  -> matches MANUSCRIPT scale? {match_ms}")
print(f"  -> matches mr_results.csv scale? {match_csv}")
print(f"  external GEFOS BMD->fracture OR = 1.55 (manuscript control OR 1.40-1.42)")
print(f"     mr_results.csv control OR = exp(1.195)=3.30  -> INCONSISTENT with GEFOS")
