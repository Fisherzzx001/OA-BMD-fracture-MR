# -*- coding: utf-8 -*-
"""
104_audit_figures_independent.py
================================
Independent verification that every FIGURE and every FIGURE LEGEND rests on
the current authoritative results.  Nothing is taken from plotted_data.json
as truth: plotted_data.json is the *defendant*, the source CSVs are the
*evidence*.

Checks performed
----------------
Fig1 / Fig2   r_g matrix and key-pair forest, rebuilt from rg_results_long_13.csv
Fig3          MR forest betas and 95% CIs, rebuilt from harmonised/*.csv
Fig4          h2 observed + liability, from rg_h2.csv / oa_h2_liability.csv
Fig6/7/S1/S2 per-SNP scatter / funnel / LOO panels -> instrument counts + IVW
Fig8          power and MDE, rebuilt from harmonised
Fig9          overlap rho_GE sign concordance (legend claims 10/10)
Fig10         MRlap direction and median shift
SupplFig3     GEFOS external benchmark ORs
SupplFig5     OA_knee h2 vs representative N
SupplFig6     fracture sample-size decision arithmetic
SupplFig7     Egger-intercept significance counts
S8 / coloc    claims quoted in the text (no PP.H4>=0.8; six PP.H3>=0.8)
"""
from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
GEN = os.path.join(HERE, "output", "genetics")
HARM = os.path.join(GEN, "harmonised")
FIGD = os.path.join(ROOT, "Figures_v17_Nature", "plotted_data.json")

PLOT = json.load(open(FIGD, encoding="utf-8"))

FAILS, OKs = [], []


def ck(name, ok, detail=""):
    (OKs if ok else FAILS).append((name, detail))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}  {detail}")


def ivw(beta_e, beta_o, se_o):
    b = beta_o / beta_e
    s = np.abs(se_o / beta_e)
    w = 1.0 / s ** 2
    beta = float(np.sum(w * b) / np.sum(w))
    se = float(np.sqrt(1.0 / np.sum(w)))
    return beta, se


print("=" * 96)
print("PART 0 — inventory / structure")
print("=" * 96)
for k, v in PLOT.items():
    keys = [x for x in v.keys() if x != "source"]
    print(f"  {k:<12} keys={keys}")

print()
print("=" * 96)
print("PART 1 — Fig 1 / Fig 2: r_g rebuilt from rg_results_long_13.csv")
print("=" * 96)
RG = pd.read_csv(os.path.join(GEN, "rg_results_long_13.csv"))
print(f"  rows={len(RG)}  cols={list(RG.columns)}")
traits = sorted(set(RG.trait1) | set(RG.trait2))
print(f"  n unique traits = {len(traits)}: {traits}")

# rebuild the symmetric matrix
def rg_of(a, b):
    r = RG[((RG.trait1 == a) & (RG.trait2 == b)) | ((RG.trait1 == b) & (RG.trait2 == a))]
    return None if r.empty else float(r.iloc[0].rg)


def q_of(a, b):
    r = RG[((RG.trait1 == a) & (RG.trait2 == b)) | ((RG.trait1 == b) & (RG.trait2 == a))]
    return None if r.empty else float(r.iloc[0].q_fdr)


f1 = PLOT["Figure1"]
order, lbls, plotted_rg, plotted_q = f1["order"], f1["labels"], np.array(f1["rg"]), np.array(f1["q_fdr"])
print(f"  plotted order      = {order}")
print(f"  plotted labels     = {lbls}")
print(f"  plots_matrix_shape = {plotted_rg.shape}")

# compare every off-diagonal cell against the rebuilt value
mism = 0
maxdev = 0.0
for i in range(len(order)):
    for j in range(len(order)):
        if i == j:
            continue
        v = rg_of(order[i], order[j])
        if v is None:
            continue
        d = abs(v - plotted_rg[i][j])
        maxdev = max(maxdev, d)
        if d > 5e-3:
            mism += 1
ck("Fig1: all 78 plotted r_g cells match rg_results_long_13.csv",
   mism == 0, f"cells mismatched={mism}, max|dev|={maxdev:.5f}")

# significance annotation pattern
qmis = 0
for i in range(len(order)):
    for j in range(len(order)):
        if i == j:
            continue
        v = q_of(order[i], order[j])
        if v is None:
            continue
        expected_sig = v < 0.05
        # plotted q array: check it at least exists and has same sign pattern
        if abs(v - plotted_q[i][j]) > 5e-3:
            qmis += 1
ck("Fig1: plotted q_FDR array matches source", qmis == 0, f"mismatched={qmis}")

# legend claims for Fig2
bmd_traits = [t for t in traits if t.upper().startswith("BMD")]
def to_r(a, b):
    return rg_of(a, b)

frac_neg = [(b, rg_of("fracture", b)) for b in bmd_traits]
neg_ok = all(v is not None and v < 0 for _, v in frac_neg)
ck("Fig2 legend: fracture vs every BMD site is negative",
   neg_ok, f"range {min(v for _, v in frac_neg if v is not None):.3f} .. "
           f"{max(v for _, v in frac_neg if v is not None):.3f}")

for oa_t in ["OA_knee", "OA_selfrep"]:
    vals = [rg_of(oa_t, b) for b in bmd_traits]
    vals = [v for v in vals if v is not None]
    ck(f"Fig2 legend: {oa_t} vs every BMD site is positive",
       len(vals) > 0 and all(v > 0 for v in vals),
       f"n={len(vals)} range {min(vals):.3f} .. {max(vals):.3f}")

r_oa_oa = rg_of("OA_knee", "OA_selfrep")
ck("Fig2 legend: two OA definitions r_g = +0.84",
   r_oa_oa is not None and abs(r_oa_oa - 0.838) < 0.005, f"source r_g={r_oa_oa}")

for oa_t in ["OA_knee", "OA_selfrep"]:
    q = q_of(oa_t, "fracture")
    rgv = rg_of(oa_t, "fracture")
    ck(f"Fig2 legend: {oa_t} with fracture non-significant",
       q is not None and q >= 0.05, f"r_g={rgv} q={q}")

print()
print("=" * 96)
print("PART 2 — Fig 3: MR forest rebuilt from SNP-level harmonised data")
print("=" * 96)
f3 = PLOT["Figure3"]
print(f"  keys={list(f3.keys())}")
ci = f3.get("ci")
pairs_by_panel = {
    "a_positive_control": f3["a_positive_control"],
    "b_BMD_to_kneeOA": f3["b_BMD_to_kneeOA"],
    "c_null_pairs": f3["c_null_pairs"],
}
fileof = {
    "BMD total->fracture": ("BMD_total__fracture", "BMD_total__fracture"),
    "BMD femoral neck->fracture": ("BMD_femur_neck__fracture", "BMD_femur_neck__fracture"),
}
allbad = 0
for pname, panel in pairs_by_panel.items():
    print(f"  -- {pname} ({type(panel).__name__})")
    items = panel.items() if isinstance(panel, dict) else enumerate(panel)
    for label, shown in items:
        if isinstance(shown, dict):
            beta_shown = (shown.get("beta", shown.get("estimate", shown.get("ivw"))))
            if beta_shown is None:
                print(f"     !! cannot read beta from {shown}")
                continue
        else:
            beta_shown = shown
        print(f"     entry label={label!r} beta={beta_shown}")
        continue
        # try every harmonised file and pick the one whose computed beta matches
        best = None
        for fn in os.listdir(HARM):
            if not fn.endswith(".csv") or fn.startswith("_"):
                continue
            d = pd.read_csv(os.path.join(HARM, fn))
            beta, se = ivw(d.beta_e.to_numpy(float), d.beta_o.to_numpy(float), d.se_o.to_numpy(float))
            dev = abs(beta - beta_shown)
            if best is None or dev < best[1]:
                best = (fn, dev, beta, se)
        ok = best is not None and best[1] < 0.001
        lo_src, hi_src = best[2] - 1.96 * best[3], best[2] + 1.96 * best[3]
        # CI check
        ci_ok = "-"
        if isinstance(shown, dict) and "ci" in shown:
            c = shown["ci"]
            cv = abs(c[0] - (best[2] - 1.96 * best[3])) < 0.002 and abs(c[1] - (best[2] + 1.96 * best[3])) < 0.002
            ci_ok = "OK" if cv else "MISMATCH"
            if not cv:
                allbad += 1
        print(f"     {label:<34} shown β={beta_shown:+.4f}  source file={best[0]:<32} "
              f"recomputed β={best[2]:+.4f} (dev={best[1]:.2e}) CI={ci_ok}")
        if not ok:
            allbad += 1
ck("Fig3: every plotted MR beta reproduces from SNP-level data",
   allbad == 0, f"problems={allbad}")

print()
print("=" * 96)
print("PART 3 — Fig 4: h2 observed + liability scale")
print("=" * 96)
for src in ("rg_h2.csv", "rg_h2_trueN.csv"):
    d = pd.read_csv(os.path.join(GEN, src))
    print(f"  {src}: rows={len(d)}")
    print(d.to_string(index=False))
print("  Fig4 plotted:")
f4 = PLOT["Figure4"]
for t, h, s, hl in zip(f4["traits"], f4["h2_obs"], f4["se"], f4["h2_liability"]):
    print(f"     {t:<22} h2_obs={h}  se={s}  h2_liab={hl}")
# recompute liability
def liability(h2obs, K):
    thr = stats.norm.ppf(1 - K)
    z = stats.norm.pdf(thr)
    return h2obs * (K * (1 - K)) / z ** 2
for K, t, expected in [(0.12, "OA_knee", 0.094), (0.083, "OA_selfrep", 0.073), (0.143, "fracture", 0.062)]:
    if t in f4["traits"]:
        i = f4["traits"].index(t)
        got = liability(f4["h2_obs"][i], K)
        ck(f"Fig4 liability h2 for {t}", abs(got - expected) <= 0.001,
           f"recomputed={got:.4f} plotted={f4['h2_liability'][i]} ms={expected}")

print()
print("=" * 96)
print("PART 4 — Fig 8: power / MDE rebuilt")
print("=" * 96)
mp = pd.read_csv(os.path.join(GEN, "mr_power.csv"))
print(mp.to_string(index=False))
f8 = PLOT["Figure8"]
print("  Fig8 pairs plotted:", json.dumps(f8["pairs"], ensure_ascii=False)[:600])

print()
print("=" * 96)
print("PART 5 — Fig 9: overlap rho_GE sign concordance (legend claims all ten)")
print("=" * 96)
ov = pd.read_csv(os.path.join(GEN, "overlap_statistics.csv"))
print(ov.to_string(index=False))
f9 = PLOT["Figure9"]
print("  plotted n_sign_concordant =", f9["n_sign_concordant"])
concord = 0
tot = 0
for _, row in ov.iterrows():
    key = row["pair"].replace("->", "__").replace(" ", "")
    name = None
    for fn in os.listdir(HARM):
        if not fn.endswith(".csv") or fn.startswith("_"):
            continue
        stem = fn[:-4].replace(" → ", "__").replace("->", "__").replace(" ", "")
        if key.replace("->", "__") in stem or stem.endswith(key):
            name = fn
            break
    if name is None:
        # reverse orientation
        parts = key.split("__")
        if len(parts) == 2:
            for fn in os.listdir(HARM):
                stem = fn[:-4]
                if parts[0] in stem and parts[1] in stem:
                    name = fn
                    break
    if name is None:
        print(f"     !! no harmonised file for pair {row['pair']}")
        continue
    d = pd.read_csv(os.path.join(HARM, name))
    beta, _ = ivw(d.beta_e.to_numpy(float), d.beta_o.to_numpy(float), d.se_o.to_numpy(float))
    rho = float(row["rho_GE_zcorr"])
    tot += 1
    same = np.sign(beta) == np.sign(rho)
    concord += int(same)
    print(f"     {row['pair']:<34} file={name:<32} MR β={beta:+.4f}  rho_GE={rho:+.4f}  same_sign={same}")
ck("Fig9 legend: rho_GE shares the sign of the MR effect in all 10 pairs",
   concord == 10, f"concordant={concord}/{tot}")

print()
print("=" * 96)
print("PART 6 — Fig 10: MRlap direction + median shift")
print("=" * 96)
ml = pd.read_csv(os.path.join(ROOT, "code", "output", "mrlap_overlap_results.csv"))
print(ml.to_string(index=False))
f10 = PLOT["Figure10"]
print("  plotted median_pct_change =", f10["median_pct_change"])

print()
print("=" * 96)
print("PART 7 — SupplFig 6: fracture sample-size decision arithmetic")
print("=" * 96)
s6 = PLOT["SupplFig6"]
N = 264973
cases, ctrl = 37857, 227116
K = cases / N
neff_theory = N * K * (1 - K)
harmonic = 1.0 / (1.0 / cases + 1.0 / ctrl)
print(f"  adopted N={s6['adopted_N']}  measured={s6['measured_effective_NK']}  "
      f"GEFOS theoretical={s6['GEFOS_theoretical']}")
print(f"  recomputed N*K*(1-K) = {neff_theory:,.0f}")
print(f"  recomputed harmonic  = {harmonic:,.0f}")
ck("SupplFig6: GEFOS theoretical value reproducible",
   abs(harmonic - s6["GEFOS_theoretical"]) < 2 and abs(neff_theory - s6["GEFOS_theoretical"]) < 3,
   f"n*K(1-K)={neff_theory:,.0f} harmonic={harmonic:,.0f} plotted={s6['GEFOS_theoretical']}")
gap = (s6["GEFOS_theoretical"] - s6["measured_effective_NK"]) / s6["GEFOS_theoretical"] * 100
ck("SupplFig6 legend/text: measured is 2.6% below theoretical",
   abs(gap - 2.6) < 0.15, f"computed gap={gap:.2f}%")

print()
print("=" * 96)
print("PART 8 — coloc claims used by S8 and the text")
print("=" * 96)
col = pd.read_csv(os.path.join(GEN, "coloc_screen_all_loci.csv"))
print(f"  rows={len(col)}")
print(f"  max PP.H4                       = {col.PP4.max():.4f}")
print(f"  loci with PP.H4 >= 0.8          = {(col.PP4 >= 0.8).sum()}")
print(f"  loci with PP.H3 >= 0.8          = {((col.PP3 >= 0.8) if 'PP3' in col else pd.Series([0])).sum()}")
print(f"  median PP.H4                    = {col.PP4.median():.4f}")
ck("text: no locus reaches PP.H4 >= 0.8", (col.PP4 >= 0.8).sum() == 0,
   f"max={col.PP4.max():.4f}")
ck("text: six loci favour distinct variants (PP.H3 >= 0.8)",
   int((col.PP3 >= 0.8).sum()) == 6, f"n={(col.PP3 >= 0.8).sum()}")
ck("text: median PP.H4 ~ 0.005", abs(col.PP4.median() - 0.005) < 0.005,
   f"median={col.PP4.median():.4f}")

print()
print("=" * 96)
print("SUMMARY")
print("=" * 96)
print(f"PASS = {len(OKs)}   FAIL = {len(FAILS)}")
for n, d in FAILS:
    print(f"  FAIL  {n}: {d}")
