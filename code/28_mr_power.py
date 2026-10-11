"""
Statistical power / minimum detectable effect (MDE) for the ten two-sample MR pairs.

Rationale: the reverse (OA -> BMD) and OA -> fracture estimates are reported as
"null", but a non-significant IVW estimate is only informative if the analysis
could have detected an effect of a plausible size. We therefore report, for every
pair, the minimum |beta| detectable at 80% power with a two-sided alpha = 0.05,
and the achieved power for the observed point estimate.

MDE          = ( z_{1-alpha/2} + z_{power} ) * SE = (1.959964 + 0.841621) * SE
power(beta)  = Phi(|beta|/SE - 1.959964) + Phi(-|beta|/SE - 1.959964)

The reading of |beta|/SE is unaffected by the N-scale correction (beta and SE
scale identically), so power conclusions are N-scale-invariant.
"""
import csv, math

SRC = "code/output/genetics/mr_corrected_N.csv"
OUT = "code/output/genetics/mr_power.csv"

Z_A = 1.959963984540054   # two-sided alpha = 0.05
Z_P = 0.8416212335729143  # 80% power


def phi(x):
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


rows = []
with open(SRC) as fh:
    for d in csv.DictReader(fh):
        beta = float(d["beta"]); se = float(d["se"]); K = int(d["K"])
        mde = (Z_A + Z_P) * se
        ratio = abs(beta) / se
        power = phi(ratio - Z_A) + phi(-ratio - Z_A)
        or_bound = math.exp(mde) if d["outcome"] == "fracture" else ""
        rows.append({
            "exposure": d["exposure"], "outcome": d["outcome"],
            "K": K, "beta": round(beta, 4), "se": round(se, 4),
            "abs_z": round(ratio, 2),
            "MDE_beta_80pct": round(mde, 4),
            "power_observed": round(power, 3),
            "OR_bound_fracture": round(or_bound, 3) if or_bound != "" else "",
            "verdict": "powered" if power >= 0.80 else "underpowered (uninformative null)"
            if abs(beta) < mde else "borderline",
        })

with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    for r in rows:
        w.writerow(r)

hdr = f"{'pair':34s} {'K':>3s} {'beta':>8s} {'SE':>8s} {'|z|':>6s} {'MDE':>8s} {'power':>7s} {'OR_bnd':>8s}"
print(hdr); print("-" * len(hdr))
for r in rows:
    print(f"{r['exposure']+'->'+r['outcome']:34s} {r['K']:3d} {r['beta']:8.4f} "
          f"{r['se']:8.4f} {r['abs_z']:6.2f} {r['MDE_beta_80pct']:8.4f} "
          f"{r['power_observed']:7.3f} {str(r['OR_bound_fracture']):>8s}")
print(f"\nSaved -> {OUT}")
nulls = [r for r in rows if r["power_observed"] < 0.80]
print(f"Underpowered (null) pairs: {len(nulls)} -> "
      f"power range {min(r['power_observed'] for r in nulls):.3f}-{max(r['power_observed'] for r in nulls):.3f}")
sig = [r for r in rows if r["power_observed"] >= 0.80]
print(f"Powered pairs: {len(sig)} -> min power {min(r['power_observed'] for r in sig):.3f}")
