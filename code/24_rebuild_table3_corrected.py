"""
code/24_rebuild_table3_corrected.py
Rebuild Table 3 of the manuscript on the STRICT-clump, N-corrected effect scale.

Principle (verified against the rendered manuscript Table 3 and AUDIT_mr_recomputed.csv):
- AUDIT_mr_recomputed.csv is the authoritative strict-clump, N-corrected recomputation
  (IVW fixed/random, MR-Egger slope + intercept, weighted median, Cochran I2, LOO sign
  flips) that MATCHES the manuscript Table 3 on every one of the ten pairs. It is used here
  as the single source so that the derived table3_corrected.csv can never drift from the
  strict set the manuscript reports. A previous version mixed in the 1-Mb/loose ivw and
  Egger columns from mr_results.csv, which produced stale loose-set Egger-intercept and
  random-effect p-values (the audit's "loose vs strict" P0/P1 flags).
- Forward MR beta/se/p are on the N-corrected scale (BMD N = 31,873; §3.7); p, Egger-intercept
  p, I2 and LOO sign flips are N-scale-invariant, so no further rescaling is required.

Output: code/output/genetics/table3_corrected.csv + printed markdown table.
"""
import csv
import math
import os

import pandas as pd

OUT = "code/output/genetics"
AUD = os.path.join(OUT, "AUDIT_mr_recomputed.csv")

# pair order for Table 3
PAIRS = [
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


def z_to_p(z):
    return 2.0 * (1.0 - 0.5 * (1.0 + math.erf(abs(z) / math.sqrt(2.0))))


def fmtp(p):
    if p < 1e-3:
        e = int(math.floor(math.log10(p)))
        return f"{p:.1e}".replace("e", "×10^").replace("^-0", "^-").replace("^-", "^-")
    return f"{p:.3f}"


def main():
    aud = pd.read_csv(AUD).set_index("pair")

    rows = []
    for (exp, out) in PAIRS:
        r = aud.loc[f"{exp}__{out}"]
        rows.append(dict(
            exp=exp, out=out, k=int(r.K),
            beta=float(r.beta), se=float(r.se), p=float(r.p),
            p_re=float(r.p_re),
            egger_beta=float(r.egger_beta),
            egger_int_p=float(r.egger_int_p),
            wmed=float(r.wmed),
            i2=float(r.I2),
            loo=int(r.loo_flips),
        ))

    fields = ["exp", "out", "k", "beta", "se", "p", "p_re",
              "egger_beta", "egger_int_p", "wmed", "i2", "loo"]
    with open(f"{OUT}/table3_corrected.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print("| Exposure -> Outcome | n_snp | IVW beta | SE | p | IVW p(RE) | Egger int p | Wmedian | I2 | LOO |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        print(f"| {r['exp']} -> {r['out']} | {r['k']} | {r['beta']:+.4f} | {r['se']:.4f} | "
              f"{fmtp(r['p'])} | {fmtp(r['p_re'])} | {r['egger_int_p']:.3f} | {r['wmed']:+.4f} | "
              f"{r['i2'] * 100:.0f}% | {r['loo']} |")
    print("\nWrote code/output/genetics/table3_corrected.csv "
          "(source: AUDIT_mr_recomputed.csv — strict clump + N-corrected)")


if __name__ == "__main__":
    main()
