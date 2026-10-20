#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""临时探查：确认 OA_knee munged 文件的 N / MAF，以及真实 GWAS 文件可用性。"""
import os
import pandas as pd

LDSC = r"C:\Users\Administrator\WorkBuddy\2026-09-08-20-17-24\OA_FLAGSHIP\07_reference\ldsc"

for nm in ["sumstats_oa_knee.sumstats.gz", "sumstats_total.sumstats.gz",
           "sumstats_femur_neck.sumstats.gz"]:
    p = os.path.join(LDSC, nm)
    if not os.path.exists(p):
        print(nm, "MISSING")
        continue
    d = pd.read_csv(p, sep=r"\s+")
    print("=" * 60)
    print(nm, "cols=", list(d.columns), "shape=", d.shape)
    print("N  :", d.N.min(), d.N.max(), "median", d.N.median(),
          "nunique", d.N.nunique())
    if "MAF" in d.columns:
        print("MAF:", round(d.MAF.min(), 4), round(d.MAF.max(), 4))

# harmonised 里 OA 结局的 se 反推 implied N
H = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "genetics", "harmonised")
for f in ["BMD_total__OA_knee.csv", "BMD_arms__OA_knee.csv",
          "BMD_total__fracture.csv", "OA_knee__BMD_total.csv"]:
    p = os.path.join(H, f)
    if not os.path.exists(p):
        continue
    d = pd.read_csv(p)
    print("-" * 60)
    print(f, "n=", len(d), "cols=", list(d.columns))
    if "eaf" in d.columns and "se_o" in d.columns:
        imp = 1.0 / (d.se_o.values ** 2 * 2 * d.eaf.values * (1 - d.eaf.values))
        print("  implied N (from se_o): median %.0f  min %.0f  max %.0f"
              % (pd.Series(imp).median(), imp.min(), imp.max()))
    if "se_e" in d.columns:
        imp = 1.0 / (d.se_e.values ** 2 * 2 * d.eaf.values * (1 - d.eaf.values))
        print("  implied N (from se_e): median %.0f  min %.0f  max %.0f"
              % (pd.Series(imp).median(), imp.min(), imp.max()))
