#!/usr/bin/env python3
"""L7 T3: joint-site gradient. BMD (total body, one fixed instrument set) -> OA at
weight-bearing vs non-weight-bearing joints.

P2's discriminating clause is that the BMD->OA effect is systemic (a set-point that is
read at non-weight-bearing sites) rather than local mechanical loading (which would
inflate weight-bearing joints only). With the SAME 31 BMD instruments applied to every
joint, the joint-site gradient is a clean test: a positive effect at thumb/hand/finger
excludes local mechanical loading; the ordering of effect sizes across the six sites
tests the P2 ordering clause.

Joint-site outcomes come from Hatzikotoulas 2025 (Nature 641:1217-1224; PMID 40205036),
mixed-ancestry meta-analyses, all GRCh37.
"""
import gzip, json, math, os, subprocess, sys

D = "G:/OA_B_PROJECT_cache/l7"
INST = json.load(open(os.path.join(D, "bmd_instruments.json")))
INST = INST["inst"]
Z975 = 1.959963985

# label, file (in joints/ or l7 root), weight-bearing class
SITES = [
    ("hip",        os.path.join(D, "GCST90566798.tsv.gz"),              "weight-bearing"),
    ("knee",       os.path.join(D, "knee_matched.tsv"),                 "weight-bearing"),
    ("spine",      os.path.join(D, "joints/GCST90566801.tsv.gz"),       "weight-bearing"),
    ("thumb",      os.path.join(D, "joints/GCST90566803.tsv.gz"),       "non-weight-bearing"),
    ("hand",       os.path.join(D, "joints/GCST90566797.tsv.gz"),       "non-weight-bearing"),
    ("finger",     os.path.join(D, "joints/GCST90566796.tsv.gz"),       "non-weight-bearing"),
]

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}

def ncdf(x):
    return 0.5 * (1 + math.erf(x / math.sqrt(2)))

def p_two(b, se):
    return 2 * (1 - ncdf(abs(b / se)))

def extract(path):
    """Return {rs: (E,O,ln_or,se,eaf)} for the 31 BMD instruments at a joint site."""
    out = {}
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
        # harmonised matched files carry the rsid in column 0 and the same column layout
        idx = {c: i for i, c in enumerate(header)}
        rs_col = 0 if "SNP" in header else 0
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 14:
                continue
            rs = p[rs_col]
            if rs not in INST:
                continue
            if p[6] == "#NA":                      # standard_error missing -> derive from CI
                cu, cl = float(p[10]), float(p[11])
                se = (math.log(cu) - math.log(cl)) / 3.92
            else:
                se = float(p[6])
            out[rs] = dict(E=p[3], O=p[4], ln_or=math.log(float(p[5])), se=se, eaf=float(p[7]))
    return out

def orient(o, allele):
    a = allele.upper(); e = o["E"].upper(); oth = o["O"].upper()
    if e == a:
        return o["ln_or"], o["se"], False
    if oth == a:
        return -o["ln_or"], o["se"], True
    if COMP.get(e) == a:
        return -o["ln_or"], o["se"], True
    if COMP.get(oth) == a:
        return o["ln_or"], o["se"], False
    return None, None, None

def analyse(rows):
    n = len(rows)
    num = sum(r["beta_o"] * r["beta_e"] / r["se_o"]**2 for r in rows)
    den = sum(r["beta_e"]**2 / r["se_o"]**2 for r in rows)
    b = num / den; se = 1 / math.sqrt(den)
    Q = sum((r["beta_o"] - b * r["beta_e"])**2 / r["se_o"]**2 for r in rows)
    df = n - 1
    I2 = max(0.0, (Q - df) / Q) * 100 if Q > 0 else 0.0
    se_r = se * math.sqrt(max(1.0, Q / df)) if df > 0 else se
    xs = [r["beta_e"] for r in rows]; ys = [r["beta_o"] for r in rows]; ws = [1 / r["se_o"]**2 for r in rows]
    sw = sum(ws); swx = sum(w * x for w, x in zip(ws, xs)); swxx = sum(w * x * x for w, x in zip(ws, xs))
    swy = sum(w * y for w, y in zip(ws, ys)); swxy = sum(w * x * y for w, x, y in zip(ws, xs, ys))
    dd = sw * swxx - swx * swx
    if dd != 0:
        be = (sw * swxy - swx * swy) / dd; ic = (swxx * swy - swx * swxy) / dd
        sebe = math.sqrt(sw / dd); seic = math.sqrt(swxx / dd)
    else:
        be = ic = sebe = seic = float("nan")
    # leave-one-out sign flips
    flips = 0
    for i in range(n):
        sub = [r for j, r in enumerate(rows) if j != i]
        nl = sum(r["beta_o"] * r["beta_e"] / r["se_o"]**2 for r in sub)
        dl = sum(r["beta_e"]**2 / r["se_o"]**2 for r in sub)
        if (nl / dl) * b <= 0:
            flips += 1
    return dict(n=n, beta=b, se=se, p=p_two(b, se), OR=math.exp(b),
                lo=math.exp(b - Z975 * se), hi=math.exp(b + Z975 * se),
                p_random=p_two(b, se_r), se_random=se_r, Q=Q, Q_df=df, I2=I2,
                egger=dict(beta=be, se=sebe, p=p_two(be, sebe), intercept=ic, p_int=p_two(ic, seic)),
                loo_flips=flips)

results = {}
print("=" * 108)
print("T3  JOINT-SITE GRADIENT :  BMD(total body, same 31 instruments) -> OA at 6 joint sites")
print("=" * 108)
for label, path, cls in SITES:
    if not os.path.exists(path):
        print(f"{label:8s} [{cls:18s}] MISSING FILE {path}")
        continue
    site = extract(path)
    rows = []
    for rs, h in INST.items():
        if rs not in site:
            continue
        bo, so, _ = orient(site[rs], h["exp_allele"])
        if bo is None:
            continue
        rows.append(dict(rs=rs, beta_e=h["beta_e"], se_e=h["se_e"], beta_o=bo, se_o=so))
    if len(rows) < 5:
        print(f"{label:8s} [{cls:18s}] only {len(rows)} instruments -> skipped")
        continue
    a = analyse(rows)
    results[label] = dict(cls=cls, **a)
    print(f"{label:8s} [{cls:18s}] n={a['n']:2d}  IVW beta={a['beta']:+.4f} SE={a['se']:.4f} "
          f"p={a['p']:.3g}  OR={a['OR']:.3f} [{a['lo']:.3f},{a['hi']:.3f}]  p_rand={a['p_random']:.3g}  "
          f"Q={a['Q']:.1f}({a['Q_df']}) I2={a['I2']:.0f}%  Egger={a['egger']['beta']:+.4f}(p={a['egger']['p']:.2g}) "
          f"int p={a['egger']['p_int']:.2g}  LOO flips={a['loo_flips']}")

# ---------- gradient contrast ----------
print("\n" + "=" * 108)
print("T3 GRADIENT CONTRAST  (weight-bearing vs non-weight-bearing joints)")
print("=" * 108)
def pool(labels):
    bs = []; ses = []
    for l in labels:
        if l in results:
            bs.append(results[l]["beta"]); ses.append(results[l]["se"])
    if not bs:
        return None
    w = [1 / s**2 for s in ses]
    b = sum(wi * bi for wi, bi in zip(w, bs)) / sum(w)
    se = 1 / math.sqrt(sum(w))
    return dict(labels=labels, n=len(bs), beta=b, se=se, p=p_two(b, se),
                OR=math.exp(b), lo=math.exp(b - Z975 * se), hi=math.exp(b + Z975 * se))

wb = [l for l, _, c in SITES if c == "weight-bearing" and l in results]
nwb = [l for l, _, c in SITES if c == "non-weight-bearing" and l in results]
pw = pool(wb); pn = pool(nwb)
out = {}
if pw and pn:
    d = pn["beta"] - pw["beta"]
    se_ind = math.sqrt(pn["se"]**2 + pw["se"]**2)   # conservative: assumes independence
    print(f"weight-bearing     ({', '.join(wb)}): pooled beta={pw['beta']:+.4f} SE={pw['se']:.4f} "
          f"p={pw['p']:.3g} OR={pw['OR']:.3f} [{pw['lo']:.3f},{pw['hi']:.3f}]")
    print(f"non-weight-bearing ({', '.join(nwb)}): pooled beta={pn['beta']:+.4f} SE={pn['se']:.4f} "
          f"p={pn['p']:.3g} OR={pn['OR']:.3f} [{pn['lo']:.3f},{pn['hi']:.3f}]")
    print(f"contrast (non-WB - WB) = {d:+.4f}  SE(conservative indep)={se_ind:.4f}  z={d/se_ind:+.2f}  "
          f"two-sided p={p_two(d, se_ind):.3g}")
    print(f"one-sided H1: non-WB > WB (P2 ordering clause) -> p = {1-ncdf(d/se_ind):.3g}")
    out["contrast"] = dict(diff=d, se=se_ind, z=d / se_ind, p_two=p_two(d, se_ind),
                            p_one_greater=1 - ncdf(d / se_ind))
    out["weight_bearing"] = pw; out["non_weight_bearing"] = pn
    # heterogeneity across all sites
    bs = [results[l]["beta"] for l in results]; ses = [results[l]["se"] for l in results]
    w = [1 / s**2 for s in ses]
    b = sum(wi * bi for wi, bi in zip(w, bs)) / sum(w)
    Q = sum(wi * (bi - b) ** 2 for wi, bi in zip(w, bs)); df = len(bs) - 1
    I2 = max(0.0, (Q - df) / Q) * 100 if Q > 0 else 0.0
    print(f"across-site heterogeneity: Q={Q:.2f} ({df} df) p={1-ncdf((Q-df)/math.sqrt(2*df)):.3g} I2={I2:.0f}%")
    out["across_site"] = dict(Q=Q, df=df, I2=I2, beta=b)

json.dump(dict(sites=results, **out), open(os.path.join(D, "l7_t3_gradient.json"), "w"), indent=2)
print("\nSaved", os.path.join(D, "l7_t3_gradient.json"))
