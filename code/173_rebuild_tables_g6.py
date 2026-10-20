#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
code/173_rebuild_tables_g6.py — 从 G6 校正后的 CSV 重新生成稿件 4 张表的 markdown 行
输出到 code/output/genetics/_g6_tables.md 供人工核验后贴回 docs/08。
溯源: mr_results_corrected.csv (G6) + wm_bootstrap_ci.csv (审计 wmed) + mvmr_height_main.csv + mvmr_bmi_main.csv
F = 3.3887 (OA 结局 K(1-K)=0.087082)
"""
import os, numpy as np, pandas as pd
from scipy.stats import norm

G = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "genetics")
F = 3.3887

mr = pd.read_csv(os.path.join(G, "mr_results_corrected.csv"))
wm = pd.read_csv(os.path.join(G, "wm_bootstrap_ci.csv")).set_index("pair")
h  = pd.read_csv(os.path.join(G, "mvmr_height_main.csv"))
b  = pd.read_csv(os.path.join(G, "mvmr_bmi_main.csv"))

SUP = {"0":"⁰","1":"¹","2":"²","3":"³","4":"⁴","5":"⁵","6":"⁶","7":"⁷","8":"⁸","9":"⁹","-":"⁻"}
def sup(s): return "".join(SUP[c] for c in s)
def pstr(p):
    if p < 1e-3:
        e = int(np.ceil(-np.log10(p)))
        m = p*(10**e)
        return f"{m:.1f}×10{sup('-'+str(e))}"
    return f"{p:.3f}"
def num(x, d=3):
    return f"{x:+.3f}".replace("-","−")
def or_ci(beta, se, protective):
    b = -beta if protective else beta
    lo = np.exp(b - 1.96*se); hi = np.exp(b + 1.96*se)
    return f"{np.exp(b):.3f} ({lo:.3f}–{hi:.3f})"

order = [("BMD_total","fracture","**BMD total → fracture** (positive control)"),
         ("BMD_femur_neck","fracture","BMD femoral neck → fracture"),
         ("BMD_total","OA_knee","**BMD total → OA knee**"),
         ("BMD_lumbar","OA_knee","Lumbar spine → OA knee"),
         ("BMD_arms","OA_knee","Arms → OA knee"),
         ("BMD_leg","OA_knee","Legs → OA knee"),
         ("BMD_femur_neck","OA_knee","Femoral neck → OA knee"),
         ("OA_knee","fracture","OA knee → fracture"),
         ("OA_knee","BMD_total","OA knee → BMD total"),
         ("OA_knee","BMD_femur_neck","OA knee → BMD femoral neck")]

# ---------------- Table 2 ----------------
print("## Table 2 (G6 校正)\n")
print("| Exposure → Outcome | n_snp | IVW β (fixed) | SE | p | IVW p (random) | OR (95% CI) | MR-Egger β (SE) | Egger p | Egger intercept p | Weighted median (95% CI) | I² | LOO sign flips |")
print("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---")
for exp,out,lab in order:
    r = mr[(mr.exposure==exp)&(mr.outcome==out)].iloc[0]
    wb = wm.loc[f"{exp}__{out}"]
    if out == "OA_knee":
        s = F
    elif exp == "OA_knee":
        s = 1.0/F
    else:
        s = 1.0
    wm_beta, wlo, whi = wb.audited_wmed*s, wb.wmed_boot_lo*s, wb.wmed_boot_hi*s
    is_binary = (out in ("OA_knee","fracture"))
    orc = or_ci(r.ivw, r.ivw_se, protective=(out=="fracture")) if is_binary else "—"
    print(f"| {lab} | {int(r.n_snp)} | {num(r.ivw)} | {r.ivw_se:.3f} | {pstr(r.ivw_p)} | "
          f"{pstr(r.ivw_re_p)} | {orc} | {num(r.egger)} ({r.egger_se:.3f}) | {pstr(r.egger_p)} | "
          f"{pstr(r.egger_intercept_p)} | {num(wm_beta)} ({num(wlo)} to {num(whi)}) | {r.I2*100:.0f}% | {int(r.loo_sign_flips)} |")

# ---------------- Table 3 ----------------
print("\n## Table 3 (G6 校正)\n")
print("| Exposure → Outcome | K | IVW β | SE | abs z | MDE (80% power) | Achieved power | Interpretation |")
print("| --- | --- | --- | --- | --- | --- | --- | --- |")
for exp,out,lab in order:
    r = mr[(mr.exposure==exp)&(mr.outcome==out)].iloc[0]
    z = abs(r.ivw/r.ivw_se); mde = 2.8*r.ivw_se
    power = norm.cdf(z-1.96) + norm.cdf(-z-1.96)
    interp = "powered" if power>=0.8 else "underpowered"
    print(f"| {lab} | {int(r.n_snp)} | {num(r.ivw)} | {r.ivw_se:.3f} | {z:.2f} | {mde:.3f} | {power:.2f} | {interp} |")

# ---------------- Table S9 (height) ----------------
print("\n## Table S9 (G6 校正, n_h=100/EUR)\n")
print("| Exposure → Outcome | n_BMD (of primary) | n_H used | β_uni (Table 2) | MVMR β_BMD | Atten. % vs Table 2 | θ_H | θ_H p | Fc(BMD) | Fc(H) | Egger p | I² |")
print("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
surv = {"BMD_total":50,"BMD_femur_neck":49,"BMD_lumbar":50,"BMD_arms":50,"BMD_leg":49}
for _,r in h[h.outcome=="OA_knee"].iterrows():
    bu = float(mr[(mr.exposure==r.exposure)&(mr.outcome=="OA_knee")].ivw.iloc[0])  # Table 2 β (已 G6)
    att = 100*(bu - r.mvmr_beta)/bu
    prim = int(mr[(mr.exposure==r.exposure)&(mr.outcome=="OA_knee")].n_snp.iloc[0])
    print(f"| {r.exposure} → OA knee | {int(r.n_bmd)} (of {prim}) | {surv[r.exposure]} | {num(bu)} | {num(r.mvmr_beta)} | {att:.1f} | {num(r.uni_h_beta)} | {pstr(r.uni_h_p)} | {r.Fc_bmd:.1f} | {r.Fc_h:.1f} | {pstr(r.egger_int_p)} | {r.I2*100:.2f} |")
fc = h[(h.outcome=="fracture")].iloc[0]
prf = int(mr[(mr.exposure=="BMD_total")&(mr.outcome=="fracture")].n_snp.iloc[0])
bfx = float(mr[(mr.exposure=="BMD_total")&(mr.outcome=="fracture")].ivw.iloc[0])
print(f"| **BMD total → fracture (positive control)** | {int(fc.n_bmd)} (of {prf}) | 50 | {num(bfx)} | {num(fc.mvmr_beta)} | {100*(bfx-fc.mvmr_beta)/bfx:.1f} | {num(fc.uni_h_beta)} | {pstr(fc.uni_h_p)} | {fc.Fc_bmd:.1f} | {fc.Fc_h:.1f} | {pstr(fc.egger_int_p)} | {fc.I2*100:.2f} |")

# ---------------- Table S16 (BMI) ----------------
print("\n## Table S16 (G6 校正, n_bmi=100)\n")
print("| Exposure → Outcome | n_BMD | n_BMI | β_uni (BMD) | MVMR β_BMD | MVMR p_BMD | θ_BMI | θ_BMI p | Atten. % | Fc(BMD) | Fc(BMI) | Egger p | I² |")
print("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---")
for _,r in b[(b.outcome=="OA_knee")&(b.n_bmi_config==100)].iterrows():
    print(f"| {r.exposure} → OA knee | {int(r.n_bmd_iv)} | {int(r.n_bmi_iv)} | {num(r.uni_beta_bmd)} | {num(r.mvmr_beta_bmd)} | {pstr(r.mvmr_p_bmd)} | {num(r.mvmr_beta_bmi)} | {pstr(r.mvmr_p_bmi)} | {r.atten_vs_table2:.1f} | {r.F_bmd:.1f} | {r.F_bmi:.1f} | {pstr(r.egger_int_p)} | {r.I2*100:.2f} |")
