# -*- coding: utf-8 -*-
"""
104b — remaining figure/legend checks (Fig3, panels 6/7/S1/S2, Fig8-10,
SupplFig3/5/6/7, coloc).  Each part is wrapped so one failure cannot abort
the run.  Structures are printed first so the comparison can be trusted.
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
PLOT = json.load(open(os.path.join(ROOT, "Figures_v17_Nature", "plotted_data.json"), encoding="utf-8"))

FAILS = []


def ck(name, ok, detail=""):
    if not ok:
        FAILS.append((name, detail))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}  {detail}")


def ivw(d):
    b = d.beta_o.to_numpy(float) / d.beta_e.to_numpy(float)
    s = (d.se_o.to_numpy(float) / d.beta_e.to_numpy(float)).__abs__()
    w = 1.0 / s ** 2
    return float(np.sum(w * b) / np.sum(w)), float(np.sqrt(1.0 / np.sum(w)))


HARM_FILES = [f for f in os.listdir(HARM) if f.endswith(".csv") and not f.startswith("_")]
TRUTH = {}
for f in HARM_FILES:
    d = pd.read_csv(os.path.join(HARM, f))
    TRUTH[f[:-4]] = ivw(d)     # (beta, se)

print("=" * 96)
print("AUTHORITATIVE MR ESTIMATES (recomputed from SNP-level files)")
print("=" * 96)
for k, (b, s) in sorted(TRUTH.items()):
    print(f"  {k:<26} beta={b:+.6f}  se={s:.6f}  95%CI=[{b-1.96*s:+.4f},{b+1.96*s:+.4f}]")


def match_any(beta):
    """find the harmonised file whose recomputed beta equals a plotted value"""
    best = min(TRUTH.items(), key=lambda kv: abs(kv[1][0] - beta))
    return best[0], abs(best[1][0] - beta), best[1]


def walk(node, path=""):
    """yield (path, number) for every numeric leaf inside a nested structure"""
    if isinstance(node, dict):
        for k, v in node.items():
            yield from walk(v, f"{path}.{k}" if path else str(k))
    elif isinstance(node, list):
        for i, v in enumerate(node):
            yield from walk(v, f"{path}[{i}]")
    elif isinstance(node, (int, float)) and not isinstance(node, bool):
        yield path, float(node)


print()
print("=" * 96)
print("PART 2 — Fig 3: every plotted MR number must exist in the authoritative set")
print("=" * 96)
for fig in ("Figure3", "Figure6", "Figure7", "SupplFig1", "SupplFig2"):
    node = PLOT.get(fig)
    if node is None:
        print(f"  {fig}: ABSENT from plotted_data.json")
        continue
    leaves = list(walk({k: v for k, v in node.items() if k != "source"}))
    print(f"\n  -- {fig}: {len(leaves)} numeric leaves")
    # for each leaf that looks like an MR beta, try to match a known pair
    unmatched, checked = [], 0
    for p, val in leaves:
        if any(tag in p.lower() for tag in ("beta", "ivw", "est", "slope")):
            name, dev, _ = match_any(val)
            checked += 1
            if dev > 1.5e-3:
                unmatched.append((p, val, name, dev))
    if checked:
        ck(f"{fig}: plotted MR betas all traceable to SNP-level data",
           not unmatched, f"checked={checked} unmatched={len(unmatched)}")
        for p, val, name, dev in unmatched[:12]:
            print(f"        !! {p} = {val:+.6f}  closest={name} dev={dev:.2e}")
    else:
        print(f"     (no beta-tagged leaves; sample leaves: {leaves[:6]})")

print()
print("=" * 96)
print("PART 3 — Fig 4 / Table 1: h2 and liability")
print("=" * 96)
f4 = PLOT["Figure4"]
try:
    for t, h, s, hl in zip(f4["traits"], f4["h2_obs"], f4["se"], f4["h2_liability"]):
        print(f"  {t:<24} h2_obs={h}  se={s}  h2_liab={hl}")
    def liab(h2obs, K):
        z = stats.norm.pdf(stats.norm.ppf(1 - K))
        return h2obs * K * (1 - K) / z ** 2
    for K, t, exp in ((0.12, "OA_knee", 0.094), (0.083, "OA_selfrep", 0.073), (0.143, "fracture", 0.062)):
        try:
            i = f4["traits"].index(t)
            got = liab(f4["h2_obs"][i], K)
            ck(f"Fig4 liability {t}", abs(got - exp) <= 0.001,
               f"recomp={got:.4f} plotted={f4['h2_liability'][i]} ms={exp}")
        except ValueError:
            print(f"  !! trait {t} not in Fig4 traits list")
except Exception as e:
    print("  !! Fig4 check error:", e)

print()
print("=" * 96)
print("PART 4 — Fig 8: power / MDE vs Table 4")
print("=" * 96)
try:
    mp = pd.read_csv(os.path.join(GEN, "mr_power.csv"))
    bad = 0
    for _, r in mp.iterrows():
        key = f"{r.exposure}__{r.outcome}"
        if key in TRUTH:
            b, s = TRUTH[key]
            z = abs(b / s)
            mde = 2.8 * s
            powr = stats.norm.cdf(abs(z) - 1.96)
            dz = abs(z - float(r.abs_z)); dm = abs(mde - float(r.MDE_beta_80pct))
            dp = abs(powr - float(r.power_observed))
            flag = "OK" if (dz < 0.02 and dm < 0.002 and dp < 0.02) else "MISMATCH"
            if flag != "OK":
                bad += 1
            print(f"  {key:<28} z: {float(r.abs_z):>6.2f} vs {z:>6.2f} | "
                  f"MDE: {float(r.MDE_beta_80pct):.4f} vs {mde:.4f} | "
                  f"power: {float(r.power_observed):.2f} vs {powr:.2f}  [{flag}]")
    ck("Fig8: power table reproduces from SNP-level data", bad == 0, f"mismatched={bad}")
except Exception as e:
    print("  !! ", e)

print()
print("=" * 96)
print("PART 5 — Fig 9: rho_GE sign concordance (legend: all ten)")
print("=" * 96)
try:
    ov = pd.read_csv(os.path.join(GEN, "overlap_statistics.csv"))
    conc = 0
    for _, r in ov.iterrows():
        key = str(r["pair"]).replace("->", "__").replace(" ", "")
        name = None
        for k in TRUTH:
            if k.replace(" ", "") == key or k.replace(" ", "") == "__".join(reversed(key.split("__"))):
                name = k
                break
        if name is None:
            for k in TRUTH:
                a, b = key.split("__") if "__" in key else (key, "")
                if a in k and b in k:
                    name = k
                    break
        beta = TRUTH[name][0] if name else float("nan")
        rho = float(r["rho_GE_zcorr"])
        same = np.sign(beta) == np.sign(rho)
        conc += int(same)
        print(f"  {str(r['pair']):<36} file={str(name):<28} MR β={beta:+.4f} rho_GE={rho:+.4f} same={same}")
    ck("Fig9 legend: rho_GE shares MR-effect sign in all 10 pairs", conc == 10, f"concordant={conc}/10")
except Exception as e:
    print("  !! ", e)

print()
print("=" * 96)
print("PART 6 — Fig 10: MRlap direction and median shift")
print("=" * 96)
try:
    ml = pd.read_csv(os.path.join(ROOT, "code", "output", "mrlap_overlap_results.csv"))
    print(ml.to_string(index=False))
    f10 = PLOT["Figure10"]
    print("  plotted median_pct_change:", f10["median_pct_change"],
          "| n_difference_significant:", f10["n_difference_significant"])
except Exception as e:
    print("  !! ", e)

print()
print("=" * 96)
print("PART 7 — SupplFig 6: fracture N decision")
print("=" * 96)
s6 = PLOT["SupplFig6"]
N, cases, ctrl = 264973, 37857, 227116
K = cases / N
neff = N * K * (1 - K)
harm = 1.0 / (1.0 / cases + 1.0 / ctrl)
print(f"  plotted: adopted_N={s6['adopted_N']} measured={s6['measured_effective_NK']} "
      f"theoretical={s6['GEFOS_theoretical']} IQR={s6['IQR']}")
print(f"  recomputed N*K*(1-K)={neff:,.0f}  harmonic={harm:,.0f}")
ck("SupplFig6: GEFOS theoretical reproducible",
   abs(neff - s6["GEFOS_theoretical"]) < 3, f"{neff:,.0f} vs {s6['GEFOS_theoretical']}")
gap = (s6["GEFOS_theoretical"] - s6["measured_effective_NK"]) / s6["GEFOS_theoretical"] * 100
ck("SupplFig6: '2.6% below' arithmetic", abs(gap - 2.6) < 0.15, f"gap={gap:.2f}%")

print()
print("=" * 96)
print("PART 8 — coloc claims (S8 + text)")
print("=" * 96)
try:
    col = pd.read_csv(os.path.join(GEN, "coloc_screen_all_loci.csv"))
    print(f"  rows={len(col)} maxPP.H4={col.PP4.max():.4f} median={col.PP4.median():.4f} "
          f"PP.H4>=0.8: {(col.PP4>=0.8).sum()}  PP.H3>=0.8: {(col.PP3>=0.8).sum()}")
    ck("text: no locus PP.H4>=0.8", (col.PP4 >= 0.8).sum() == 0, f"max={col.PP4.max():.4f}")
    ck("text: six loci PP.H3>=0.8", int((col.PP3 >= 0.8).sum()) == 6, f"n={(col.PP3>=0.8).sum()}")
    ck("text: median PP.H4 ~0.005", abs(col.PP4.median() - 0.005) < 0.005, f"{col.PP4.median():.4f}")
except Exception as e:
    print("  !! ", e)

print()
print("=" * 96)
print("SUMMARY  FAIL =", len(FAILS))
for n, d in FAILS:
    print("   FAIL:", n, "|", d)
