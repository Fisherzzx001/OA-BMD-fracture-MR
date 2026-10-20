# -*- coding: utf-8 -*-
"""104c — corrected re-checks for Fig4 (liability dict), Fig9 (pair key),
and Fig10 (direction + median shift)."""
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
PLOT = json.load(open(os.path.join(ROOT, "Figures_v17_Nature", "plotted_data.json"), encoding="utf-8"))
FAILS = []


def ck(n, ok, d=""):
    if not ok:
        FAILS.append((n, d))
    print(f"  [{'PASS' if ok else 'FAIL'}] {n}  {d}")


print("=" * 96)
print("PART 3 (fixed) — Fig 4: h2_liability is a dict; verify conversion")
print("=" * 96)
f4 = PLOT["Figure4"]
print("  traits      :", f4["traits"])
print("  h2_obs      :", f4["h2_obs"])
print("  h2_liability:", f4["h2_liability"])
hl = f4["h2_liability"]
if isinstance(hl, dict):
    for t, v in hl.items():
        print(f"     liability {t:<24} = {v}")
    # recompute from the plotted observed h2 for the matching trait
    def liab(h2obs, K):
        z = stats.norm.pdf(stats.norm.ppf(1 - K))
        return h2obs * K * (1 - K) / z ** 2
    for K, t, exp in ((0.12, "OA_knee", 0.094), (0.083, "OA_selfrep", 0.073), (0.143, "fracture", 0.062)):
        got = liab(hl.get(t, float("nan")), K) if False else None
        # map plotted display name -> h2_obs
        idx = None
        disp = {"OA_knee": "OA (knee)", "OA_selfrep": "OA (self-reported)", "fracture": "Any fracture"}.get(t)
        if disp and disp in f4["traits"]:
            idx = f4["traits"].index(disp)
        if idx is None:
            print(f"     !! cannot locate {t} in traits")
            continue
        recomputed = liab(f4["h2_obs"][idx], K)
        plotted = hl.get(t)
        ck(f"Fig4 liability {t}",
           plotted is not None and abs(recomputed - plotted) <= 0.001,
           f"h2_obs={f4['h2_obs'][idx]:.5f} recomputed={recomputed:.4f} plotted={plotted} ms={exp}")
        ck(f"Fig4 liability {t} matches manuscript Table 1",
           plotted is not None and abs(plotted - exp) <= 0.001, f"plotted={plotted} ms={exp}")
else:
    print("  !! h2_liability is not a dict:", type(hl))

print()
print("=" * 96)
print("PART 5 (fixed) — Fig 9: rho_GE sign concordance")
print("=" * 96)
ov = pd.read_csv(os.path.join(GEN, "overlap_statistics.csv"))
TRUTH = {}
for f in os.listdir(HARM):
    if not f.endswith(".csv") or f.startswith("_"):
        continue
    d = pd.read_csv(os.path.join(HARM, f))
    b = d.beta_o.to_numpy(float) / d.beta_e.to_numpy(float)
    s = np.abs(d.se_o.to_numpy(float) / d.beta_e.to_numpy(float))
    w = 1.0 / s ** 2
    TRUTH[f[:-4]] = (float(np.sum(w * b) / np.sum(w)), float(np.sqrt(1.0 / np.sum(w))))

conc = 0
for _, r in ov.iterrows():
    key = str(r["pair"]).replace(" ", "")
    if key.endswith(".csv"):
        key = key[:-4]
    key = key.replace("->", "__")
    name = key if key in TRUTH else None
    if name is None:
        a, b = (key.split("__") + [""])[:2]
        for k in TRUTH:
            if a in k and b in k:
                name = k
                break
    beta = TRUTH[name][0] if name else float("nan")
    rho = float(r["rho_GE_zcorr"])
    same = np.sign(beta) == np.sign(rho)
    conc += int(same)
    print(f"  {str(r['pair']):<38} -> {str(name):<28} MR β={beta:+.4f}  rho_GE={rho:+.4f}  same={same}")
ck("Fig9 legend: rho_GE shares the sign of the MR effect in all 10 pairs",
   conc == 10, f"concordant={conc}/10")

print()
print("=" * 96)
print("PART 6 (fixed) — Fig 10: MRlap direction + median shift")
print("=" * 96)
ml = pd.read_csv(os.path.join(ROOT, "code", "output", "mrlap_overlap_results.csv"))
pct, away = [], 0
for _, r in ml.iterrows():
    o, c = float(r.observed_effect), float(r.corrected_effect)
    p = 100.0 * (abs(c) - abs(o)) / abs(o)
    pct.append(p)
    away += int(abs(c) > abs(o))
    print(f"  {r.exposure:<15}->{r.outcome:<16} observed={o:+.6f} corrected={c:+.6f} "
          f"|moved away from null| by {p:+.1f}%  away={abs(c) > abs(o)}")
med = float(np.median(pct))
f10 = PLOT["Figure10"]
print(f"  computed median shift = {med:.2f}%   plotted median_pct_change = {f10['median_pct_change']:.2f}%")
ck("Fig10 legend: every pair moves away from the null", away == len(ml), f"{away}/{len(ml)}")
ck("Fig10 legend: median shift ~ +20.5%", abs(med - 20.5) < 0.6, f"computed={med:.2f}% plotted={f10['median_pct_change']:.2f}%")
ck("Fig10: n_difference_significant matches source",
   int((ml["p_difference"] < 0.05).sum()) == int(f10["n_difference_significant"]),
   f"source={(ml['p_difference'] < 0.05).sum()} plotted={f10['n_difference_significant']}")

print()
print("=" * 96)
print("SUMMARY  FAIL =", len(FAILS))
for n, d in FAILS:
    print("   FAIL:", n, "|", d)
