#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""172_g6_manuscript_values.py — 生成稿件回填所需的全部 G6 后数值（只读，不写结果文件）"""
import os
import sys
import importlib
import numpy as np
import pandas as pd
from scipy.stats import norm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
m21 = importlib.import_module("21_mr_corrected_N")

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "output", "genetics")
C = float(m21.F_OA_KNEE)

mr = pd.read_csv(os.path.join(OUT, "mr_results_corrected.csv"))
pw = pd.read_csv(os.path.join(OUT, "mr_power.csv"))
wm = pd.read_csv(os.path.join(OUT, "wm_recomputed_vs_tabulated.csv"))
pr = pd.read_csv(os.path.join(OUT, "mr_presso_results_corrected.csv"))
mh = pd.read_csv(os.path.join(OUT, "mvmr_height_main.csv"))
mb = pd.read_csv(os.path.join(OUT, "mvmr_bmi_main.csv"))
h2 = pd.read_csv(os.path.join(OUT, "oa_h2_liability.csv"))

WM = {r.pair: r for r in wm.itertuples()}
PR = {(r.exposure, r.outcome): r for r in pr.itertuples()}
PW = {(r.exposure, r.outcome): r for r in pw.itertuples()}


def sci(p):
    if p >= 0.001:
        return "%.3f" % p
    return "%.1f×10⁻%d" % (p * 10 ** int(-np.floor(np.log10(p))), int(-np.floor(np.log10(p))))


print("=" * 92)
print("TABLE 2  (n_snp | β | SE | p | p_RE | OR (95% CI) | Egger β (SE) | Egger p | "
      "int p | wmed (95% CI) | I² | LOO)")
print("=" * 92)
for r in mr.itertuples():
    k = f"{r.exposure}__{r.outcome}"
    w = WM.get(k)
    wm_s = ("%+.3f (%+.3f to %+.3f)" % (w.wm_recomp, w.wm_lo, w.wm_hi)) if w else "—"
    if r.outcome == "OA_knee":
        ors = "%.3f (%.3f–%.3f)" % (np.exp(r.ivw), np.exp(r.ivw - 1.96 * r.ivw_se),
                                    np.exp(r.ivw + 1.96 * r.ivw_se))
    elif r.outcome == "fracture" and r.exposure.startswith("BMD"):
        ors = "%.2f (%.2f–%.2f)" % (np.exp(-r.ivw), np.exp(-r.ivw - 1.96 * r.ivw_se),
                                    np.exp(-r.ivw + 1.96 * r.ivw_se))
    else:
        ors = "—"
    print("| %-16s → %-15s | %d | %+.3f | %.3f | %s | %s | %s | %+.3f (%.3f) | %s | %.3f | %s | %.0f%% | %d |"
          % (r.exposure, r.outcome, r.n_snp, r.ivw, r.ivw_se, sci(r.ivw_p), sci(r.ivw_re_p),
             ors, r.egger, r.egger_se, sci(r.egger_p), r.egger_intercept_p, wm_s,
             r.I2 * 100, r.loo_sign_flips))

print()
print("=" * 92)
print("TABLE 3  (K | β | SE | abs z | MDE | power)")
print("=" * 92)
for r in mr.itertuples():
    p = PW.get((r.exposure, r.outcome))
    print("| %-16s → %-15s | %d | %+.3f | %.3f | %.2f | %.3f | %.2f |"
          % (r.exposure, r.outcome, r.n_snp, r.ivw, r.ivw_se, abs(r.ivw / r.ivw_se),
             2.8 * r.ivw_se, p.power_observed if p else np.nan))

print()
print("=" * 92)
print("关键派生量")
print("=" * 92)
oas = mr[mr.outcome == "OA_knee"]
print("摘要 OR 区间：%.3f–%.3f ；CI 下界最小 %.3f，上界最大 %.3f"
      % (np.exp(oas.ivw.min()), np.exp(oas.ivw.max()),
         np.exp(oas.ivw - 1.96 * oas.ivw_se).min(), np.exp(oas.ivw + 1.96 * oas.ivw_se).max()))
print("与 Hartley 膝 OA log-odds 0.3365 的倍差：%.1f–%.1f 倍"
      % (0.3365 / oas.ivw.max(), 0.3365 / oas.ivw.min()))
rev = mr[mr.exposure == "OA_knee"]
print("反向 MDE 区间：%.3f–%.3f  (排除 |β|≥ %.3f)"
      % ((2.8 * rev.ivw_se).min(), (2.8 * rev.ivw_se).max(), (2.8 * rev.ivw_se).max()))
a = mr[(mr.exposure == "BMD_arms") & (mr.outcome == "OA_knee")].iloc[0]
f = mr[(mr.exposure == "BMD_femur_neck") & (mr.outcome == "OA_knee")].iloc[0]
d = a.ivw - f.ivw
sd = np.sqrt(a.ivw_se ** 2 + f.ivw_se ** 2)
print("T3a 部位梯度 arms−femoral neck：%+.5f (SE %.5f)，单侧 p=%.4f"
      % (d, sd, norm.sf(d / sd)))
for k, lab in ((("OA_knee", "BMD_total"), "total-body BMD"),
               (("OA_knee", "fracture"), "fracture")):
    p = PR[k]
    print("PRESSO 离群剔除后 %s：%+.3f (完整 %+.3f, 离群 %d)"
          % (lab, p.b_ivw_out, p.b_ivw_full, p.n_outlier))

print()
print("=" * 92)
print("MVMR — 身高（Table S9 主配置 n_h=100 / EUR）")
print("=" * 92)
for r in mh[mh.spec.str.contains("n_h=100")].itertuples():
    print("| %-16s → %-9s | n_BMD %d | β_uni %+.4f | MVMR β_BMD %+.4f | Atten %.1f%% | "
          "θ_H %+.4f | θ_H p %.1e | Fc(BMD) %.1f | Fc(H) %.1f | Egger p %.2f | I² %.2f |"
          % (r.exposure, r.outcome, r.n_bmd, r.uni_beta, r.mvmr_beta, r.atten_vs_full,
             r.mvmr_h, r.mvmr_h_p, r.Fc_bmd if hasattr(r, "Fc_bmd") else r.Fc_bmd,
             r.Fc_h, r.egger_int_p, r.I2))

print()
print("=" * 92)
print("MVMR — BMI（Table S16 主配置 n_bmi=100）")
print("=" * 92)
mbp = mb[mb.n_bmi_config == 100] if "n_bmi_config" in mb.columns else mb
for r in mbp.itertuples():
    print("| %-16s → %-9s | β_uni %+.4f | MVMR β_BMD %+.4f | Atten %.1f%% | θ_BMI %+.3f | "
          "p_BMD %.2e | Fc(BMD) %.1f | Fc(BMI) %.1f | Egger int p %.3f |"
          % (r.exposure, r.outcome, r.uni_beta_bmd, r.mvmr_beta_bmd, r.attenuation_pct,
             r.mvmr_beta_bmi, r.mvmr_p_bmd, r.Fcond_bmd, r.Fcond_bmi, r.egger_int_p))

print()
print("=" * 92)
print("Table 1 脚注：liability-scale h²")
print("=" * 92)
print(h2.to_string(index=False))
