#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""EUR-only knee-OA MR (reviewer R2-2): BMD (EUR instruments, Qian 2025) -> FinnGen EUR-only knee OA.

Replicates code/13 harmonise + MR (IVW / MR-Egger / weighted median).
Exposure betas come from the already-harmonised main-MR files (beta_e, se_e, on
reference A1); alleles (A1,A2) come from the munged BMD sumstats so the FinnGen
outcome can be aligned to the exposure by name (EA==OA flip). Outcome = FinnGen
knee OA (EUR-only), pre-harmonised to BMD instruments (42 SNPs, chr:pos key).
Read-only on inputs; writes a JSON + prints a table.
"""
import os, gzip, json, math
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
HARM = os.path.join(HERE, "output", "genetics", "harmonised")
BMD_MUNGED = os.path.join(HERE, "output", "genetics", "ldsc_rg_verified", "munged", "BMD_total_N31873.sumstats.gz")
FINNGEN = r"C:\Users\Administrator\WorkBuddy\2026-09-08-20-17-24\OA_FLAGSHIP\02_data\11_bmd_snp_cache\finngen_knee.csv.gz"
EUR_BIM = r"G:\REF\1000G_EUR\EUR.bim"

SITES = {
    "BMD_total": "BMD_total__OA_knee.csv",
    "BMD_femur_neck": "BMD_femur_neck__OA_knee.csv",
    "BMD_lumbar": "BMD_lumbar__OA_knee.csv",
    "BMD_arms": "BMD_arms__OA_knee.csv",
    "BMD_leg": "BMD_leg__OA_knee.csv",
}

def load_bmd_alleles():
    """SNP -> (A1, A2) from munged BMD sumstats (reference A1/A2)."""
    d = {}
    with gzip.open(BMD_MUNGED, "rt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        si, a1i, a2i = hdr.index("SNP"), hdr.index("A1"), hdr.index("A2")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) <= a2i:
                continue
            d[p[si]] = (p[a1i].upper(), p[a2i].upper())
    return d

def load_finngen():
    """key=chr:pos -> dict with EA,OA,beta,se,eaf,p ; map to rsID via EUR.bim."""
    bim = {}
    with open(EUR_BIM) as f:
        for line in f:
            c = line.split("\t")
            if len(c) < 6:
                continue
            bim[(int(c[0]), int(c[3]))] = (c[1], c[4].upper(), c[5].upper())  # (rsID, A1, A2)
    fg = {}
    with gzip.open(FINNGEN, "rt") as f:
        hdr = f.readline().rstrip("\n").split(",")
        ki = hdr.index("key"); eai = hdr.index("EA_o"); oai = hdr.index("OA_o")
        bi = hdr.index("beta_o"); si = hdr.index("se_o"); ei = hdr.index("eaf_o"); pi = hdr.index("p_o")
        for line in f:
            p = line.rstrip("\n").split(",")
            if len(p) <= pi:
                continue
            try:
                ch, pos = p[ki].split(":")
                rs = bim.get((int(ch), int(pos)))
            except Exception:
                continue
            if rs is None:
                continue
            fg[rs[0]] = dict(EA=p[eai].upper(), OA=p[oai].upper(),
                             beta=float(p[bi]), se=float(p[si]), eaf=float(p[ei]), p=float(p[pi]))
    return fg

def harmonise(inst, fg):
    """inst: DataFrame with SNP,beta_e,se_e ; fg: rsID->outcome dict. Align by name."""
    rows = []
    for _, r in inst.iterrows():
        q = fg.get(r.SNP)
        if q is None:
            continue
        a1e, a2e = BMD_ALLELES.get(r.SNP, (None, None))
        a1o, a2o = q["EA"], q["OA"]
        if a1e is None:
            continue
        if a1e == a1o and a2e == a2o:
            f = 1.0
        elif a1e == a2o and a2e == a1o:
            f = -1.0
        else:
            continue
        rows.append(dict(SNP=r.SNP, beta_e=r.beta_e, se_e=r.se_e,
                         beta_o=f * q["beta"], se_o=q["se"], flip=f))
    return pd.DataFrame(rows)

def _ivw(be, se_e, bo, so):
    w = 1.0 / (so ** 2)              # outcome-variance weighted (standard IVW)
    b = float(np.sum(w * be * bo) / np.sum(w * be ** 2))
    se = 1.0 / np.sqrt(np.sum(w * be ** 2))
    return b, se

def _egger(be, bo, so):
    w = 1.0 / (so ** 2)
    x, y = be, bo
    n = len(x)
    mx = np.sum(w * x) / np.sum(w); my = np.sum(w * y) / np.sum(w)
    sxx = np.sum(w * (x - mx) ** 2)
    b = np.sum(w * (x - mx) * (y - my)) / sxx
    a = my - b * mx
    resid = y - (a + b * x)
    s2 = np.sum(w * resid ** 2) / max(n - 2, 1)
    se_b = np.sqrt(s2 / sxx)
    se_a = np.sqrt(s2 * (1.0 / np.sum(w) + mx ** 2 / sxx))
    return b, se_b, a, se_a

def _wmed(be, bo, so):
    """Weighted median of ratio estimates gamma_i = bo_i/be_i, weight w_i = be_i^2/so_i^2."""
    g = bo / be
    w = (be ** 2) / (so ** 2)
    idx = np.argsort(g)
    g = g[idx]; w = w[idx]
    cw = np.cumsum(w); tot = cw[-1]
    m = cw[-1] / 2.0
    # find median point
    k = np.searchsorted(cw, m)
    k = min(k, len(g) - 1)
    return float(g[k]), float(w[k]), float(tot)

def main():
    global BMD_ALLELES
    BMD_ALLELES = load_bmd_alleles()
    fg = load_finngen()
    print(f"FinnGen EUR-knee-OA SNPs mapped to rsID: {len(fg)}")
    out = {}
    print(f"\n{'site':16s} {'n':>3s} {'IVW_beta':>10s} {'se':>9s} {'p':>10s} {'OR/SD':>8s} {'Egger_int_p':>12s} {'WMed_beta':>11s}")
    for site, fn in SITES.items():
        inst = pd.read_csv(os.path.join(HARM, fn))
        h = harmonise(inst, fg)
        if len(h) < 3:
            print(f"{site:16s} {len(h):3d}  <3 SNPs, skipped"); continue
        be, se_e = h.beta_e.values, h.se_e.values
        bo, so = h.beta_o.values, h.se_o.values
        b, se = _ivw(be, se_e, bo, so)
        p = 2 * (1 - stats.norm.cdf(abs(b / se)))
        orsd = math.exp(b); lo = math.exp(b - 1.96 * se); hi = math.exp(b + 1.96 * se)
        eb, seb, ea, sea = _egger(be, bo, so)
        ep = 2 * (1 - stats.norm.cdf(abs(ea / sea)))
        wmb, wmw, wmt = _wmed(be, bo, so)
        out[site] = dict(n=int(len(h)), ivw_beta=round(b, 5), ivw_se=round(se, 5),
                         ivw_p=float(f"{p:.2e}"), or_sd=round(orsd, 4),
                         or_lo=round(lo, 4), or_hi=round(hi, 4),
                         egger_slope=round(eb, 5), egger_int=round(ea, 5), egger_int_p=float(f"{ep:.2e}"),
                         wmed_beta=round(wmb, 5))
        print(f"{site:16s} {len(h):3d} {b:10.5f} {se:9.5f} {p:10.2e} {orsd:8.4f} {ep:12.2e} {wmb:11.5f}")
    with open(os.path.join(HERE, "output", "genetics", "eur_only_knee_oa_mr.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("\nSaved: code/output/genetics/eur_only_knee_oa_mr.json")

if __name__ == "__main__":
    main()
