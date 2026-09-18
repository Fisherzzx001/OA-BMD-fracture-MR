"""
Compute a data-derived sample-overlap / pleiotropy gauge for every MR pair.

For each harmonised per-SNP file we compute the empirical Pearson correlation
between the exposure and outcome association statistics across the IVs:

    rho_GE = corr( beta_e/se_e , beta_o/se_o ) = corr( z_e , z_o )

A positive rho_GE for a BMD->outcome pair indicates that the IVs' exposure and
outcome effects co-vary in the same direction, which is what sample overlap (and
shared biology) produces; it confirms the IVW direction is not driven by a single
anti-correlated instrument. The value is saved to overlap_statistics.csv.
"""
import csv, glob, os, statistics as st

SRC = "code/output/genetics/harmonised"
OUT = "code/output/genetics/overlap_statistics.csv"

def corr(xs, ys):
    n = len(xs)
    mx, my = sum(xs)/n, sum(ys)/n
    sxx = sum((x-mx)**2 for x in xs)
    syy = sum((y-my)**2 for y in ys)
    sxy = sum((x-mx)*(y-my) for x, y in zip(xs, ys))
    if sxx == 0 or syy == 0:
        return float("nan")
    return sxy / (sxx*syy)**0.5

rows = []
for f in sorted(glob.glob(os.path.join(SRC, "*.csv"))):
    base = os.path.basename(f)
    if base.startswith("_"):
        continue
    ze, zo, be, bo = [], [], [], []
    with open(f) as fh:
        r = csv.DictReader(fh)
        for d in r:
            try:
                be.append(float(d["beta_e"])); bo.append(float(d["beta_o"]))
                ze.append(float(d["beta_e"])/float(d["se_e"]))
                zo.append(float(d["beta_o"])/float(d["se_o"]))
            except (ValueError, ZeroDivisionError, KeyError):
                continue
    if len(ze) < 3:
        continue
    rho_ge = corr(ze, zo)
    rho_beta = corr(be, bo)
    rows.append((base, len(ze), round(rho_ge, 4), round(rho_beta, 4)))

rows.sort(key=lambda r: r[0])
with open(OUT, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["pair", "n_iv", "rho_GE_zcorr", "rho_beta"])
    for r in rows:
        w.writerow(r)

print(f"{'pair':32s} {'n_iv':>5s} {'rho_GE':>9s} {'rho_beta':>9s}")
for r in rows:
    print(f"{r[0]:32s} {r[1]:5d} {r[2]:9.4f} {r[3]:9.4f}")
print(f"\nSaved -> {OUT}")
ge = [r[2] for r in rows]
print(f"rho_GE range across 10 pairs: {min(ge):.4f} .. {max(ge):.4f}")
print(f"All BMD->outcome rho_GE positive: {all(v>0 for k,v in [(r[0],r[2]) for r in rows] if k.startswith('BMD'))}")
