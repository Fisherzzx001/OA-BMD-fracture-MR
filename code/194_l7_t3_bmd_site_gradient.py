#!/usr/bin/env python3
"""L7 T3 (BMD-site axis) / F9 formalisation: is the BMD->OA effect ordered by
weight-bearing status of the BMD measurement site?

The manuscript reports the site-level ordering clause of P2 as under-powered
(arms - femoral neck = +0.017, one-sided p = 0.228) from a single pairwise contrast.
Here the full ten-site frame of the five prespecified DXA sites is used:

  non-weight-bearing : arms
  weight-bearing     : total body, lumbar spine, legs, femoral neck

Two contrasts are pre-specified:
  (a) arms - femoral neck  (the manuscript's original pairwise contrast, reproduced)
  (b) arms - mean(weight-bearing sites)  (the systematic contrast the ordering clause needs)

Reported with the manuscript's multiplicative random-effects convention
(SE_random = SE_fixed * sqrt(max(1, Q/df))) and with a Cochran's Q test of
between-site homogeneity (F6 already reported Q = 6.92, 4 df, p = 0.14, tau2 = 0).
"""
import csv, math, os

BASE = "C:/Users/Administrator/WorkBuddy/2026-09-08-22-17-40/OA_B_PROJECT/code/output/genetics/harmonised"
Z975 = 1.959963985

SITES = [("total body", "BMD_total__OA_knee.csv", "weight-bearing"),
         ("lumbar spine", "BMD_lumbar__OA_knee.csv", "weight-bearing"),
         ("arms", "BMD_arms__OA_knee.csv", "non-weight-bearing"),
         ("legs", "BMD_leg__OA_knee.csv", "weight-bearing"),
         ("femoral neck", "BMD_femur_neck__OA_knee.csv", "weight-bearing")]

def ncdf(x):
    return 0.5 * (1 + math.erf(x / math.sqrt(2)))

def p_two(b, se):
    return 2 * (1 - ncdf(abs(b / se)))

def analyse(path):
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            rows.append(dict(beta_e=float(r["beta_e"]), se_e=float(r["se_e"]),
                             beta_o=float(r["beta_o"]), se_o=float(r["se_o"])))
    n = len(rows)
    num = sum(r["beta_o"] * r["beta_e"] / r["se_o"]**2 for r in rows)
    den = sum(r["beta_e"]**2 / r["se_o"]**2 for r in rows)
    b = num / den; se = 1 / math.sqrt(den)
    Q = sum((r["beta_o"] - b * r["beta_e"])**2 / r["se_o"]**2 for r in rows)
    df = n - 1
    I2 = max(0.0, (Q - df) / Q) * 100 if Q > 0 else 0
    se_r = se * math.sqrt(max(1.0, Q / df))
    flips = 0
    for i in range(n):
        sub = [r for j, r in enumerate(rows) if j != i]
        nl = sum(r["beta_o"] * r["beta_e"] / r["se_o"]**2 for r in sub)
        dl = sum(r["beta_e"]**2 / r["se_o"]**2 for r in sub)
        if (nl / dl) * b <= 0:
            flips += 1
    return dict(n=n, beta=b, se=se, p=p_two(b, se), OR=math.exp(b),
                lo=math.exp(b - Z975 * se), hi=math.exp(b + Z975 * se),
                p_random=p_two(b, se_r), se_random=se_r, Q=Q, Q_df=df, I2=I2, loo_flips=flips)

print("=" * 104)
print("T3 / F9  BMD-SITE GRADIENT (five prespecified DXA sites -> knee OA)")
print("=" * 104)
res = {}
for name, fn, cls in SITES:
    p = os.path.join(BASE, fn)
    if not os.path.exists(p):
        print(f"{name:14s} MISSING {fn}"); continue
    a = analyse(p)
    res[name] = dict(cls=cls, **a)
    print(f"{name:14s} [{cls:18s}] n={a['n']:2d} IVW beta={a['beta']:+.4f} SE={a['se']:.4f} p={a['p']:.3g} "
          f"OR={a['OR']:.3f} [{a['lo']:.3f},{a['hi']:.3f}] p_rand={a['p_random']:.3g} "
          f"Q={a['Q']:.1f}({a['Q_df']}) I2={a['I2']:.0f}% LOO flips={a['loo_flips']}")

print()
arms = res.get("arms"); fnk = res.get("femoral neck")
wb = [v for k, v in res.items() if v["cls"] == "weight-bearing"]
nwb = [v for k, v in res.items() if v["cls"] == "non-weight-bearing"]

# (a) manuscript's original contrast
if arms and fnk:
    d = arms["beta"] - fnk["beta"]
    se = math.sqrt(arms["se"]**2 + fnk["se"]**2)
    print(f"(a) arms - femoral neck            = {d:+.4f}  SE={se:.4f}  z={d/se:+.2f}  "
          f"two-sided p={p_two(d,se):.3g}  one-sided (arms>FN) p={1-ncdf(d/se):.3g}")
    print(f"    manuscript reports +0.017, one-sided p = 0.228")
    print()

# (b) systematic contrast: arms vs mean of weight-bearing sites
def pool(vs):
    bs = [v["beta"] for v in vs]; ses = [v["se"] for v in vs]
    w = [1 / s**2 for s in ses]
    b = sum(wi * bi for wi, bi in zip(w, bs)) / sum(w)
    return b, 1 / math.sqrt(sum(w))
if wb and nwb:
    bw, sew = pool(wb); ba, sea = pool(nwb)
    d = ba - bw
    se = math.sqrt(sea**2 + sew**2)
    print(f"(b) arms - mean(weight-bearing)     = {d:+.4f}  SE={se:.4f}  z={d/se:+.2f}  "
          f"two-sided p={p_two(d,se):.3g}  one-sided (arms>WB) p={1-ncdf(d/se):.3g}")
    print(f"    arms beta={ba:+.4f} (SE {sea:.4f}); weight-bearing pooled beta={bw:+.4f} (SE {sew:.4f})")

# between-site heterogeneity (Cochran Q across the five site estimates)
bs = [v["beta"] for v in res.values()]; ses = [v["se"] for v in res.values()]
w = [1 / s**2 for s in ses]
b = sum(wi * bi for wi, bi in zip(w, bs)) / sum(w)
Q = sum(wi * (bi - b) ** 2 for wi, bi in zip(w, bs)); df = len(bs) - 1
I2 = max(0.0, (Q - df) / Q) * 100 if Q > 0 else 0
print(f"\nbetween-site heterogeneity (Cochran Q on the five site estimates): "
      f"Q={Q:.2f} ({df} df) p={1-ncdf((Q-df)/math.sqrt(2*df)):.3g} I2={I2:.0f}%")
print(f"    manuscript F6 reports Q = 6.92, 4 df, p = 0.14, tau2 = 0")
print(f"\nordering of site estimates (descending): " +
      ", ".join(f"{k} {v['beta']:+.3f}" for k, v in sorted(res.items(), key=lambda kv: -kv[1]['beta'])))
print(f"    weight-bearing sites occupy ranks: " +
      ", ".join(str(i + 1) for i, (k, v) in enumerate(sorted(res.items(), key=lambda kv: -kv[1]['beta'])) if v["cls"] == "weight-bearing"))
print(f"    non-weight-bearing (arms) rank: " +
      str([i + 1 for i, (k, v) in enumerate(sorted(res.items(), key=lambda kv: -kv[1]['beta'])) if v["cls"] == "non-weight-bearing"][0]))
