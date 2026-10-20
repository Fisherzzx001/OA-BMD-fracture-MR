# -*- coding: utf-8 -*-
"""Table 2 全单元格核验：稿件 docs/08 Table 2 vs mr_results_corrected.csv
（IVW/SE/p/p_RE/OR/Egger/int_p/I²/LOO）+ wm_recomputed_vs_tabulated.csv（wmed+CI）。
输出每个单元格的稿件值、真源值、判定。"""
import csv
import math
import re

ROOT = "code/output/genetics"


def load(p):
    with open(p, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


mr = {(r["exposure"], r["outcome"]): r
      for r in load(f"{ROOT}/mr_results_corrected.csv")}
wm = {r["pair"]: r for r in load(f"{ROOT}/wm_recomputed_vs_tabulated.csv")}

# 稿件 Table 2 数据行（手工转录，避免解析 markdown 的转义字符）
MAN = [
    # (exposure, outcome, n_snp, ivw, se, p, p_re, OR, egger_b, egger_se, egger_p, int_p, wmed, wm_lo, wm_hi, I2, loo)
    ("BMD_total", "fracture", 31, -0.337, 0.022, 3.7e-51, 2.4e-12, (1.401, 1.341, 1.464),
     -0.337, 0.049, 4.0e-12, 0.558, -0.393, -0.511, -0.205, 78, 0),
    ("BMD_femur_neck", "fracture", 25, -0.351, 0.030, 7.7e-31, 2.3e-06, (1.420, 1.338, 1.507),
     -0.343, 0.068, 4.9e-07, 0.018, -0.365, -0.414, -0.129, 83, 0),
    ("BMD_total", "OA_knee", 31, 0.078, 0.013, 3.9e-10, 0.006, (1.081, 1.055, 1.108),
     0.078, 0.029, 0.007, 0.810, 0.065, 0.027, 0.177, 81, 0),
    ("BMD_lumbar", "OA_knee", 21, 0.095, 0.016, 7.7e-09, 0.017, (1.099, 1.064, 1.135),
     0.089, 0.040, 0.024, 0.295, 0.115, 0.014, 0.189, 83, 0),
    ("BMD_arms", "OA_knee", 16, 0.075, 0.016, 1.7e-06, 0.070, (1.078, 1.045, 1.112),
     0.075, 0.043, 0.080, 0.982, 0.026, 0.019, 0.217, 86, 0),
    ("BMD_leg", "OA_knee", 24, 0.060, 0.014, 2.3e-05, 0.098, (1.062, 1.033, 1.092),
     0.058, 0.038, 0.128, 0.815, 0.026, -0.003, 0.132, 85, 0),
    ("BMD_femur_neck", "OA_knee", 25, 0.058, 0.017, 6.4e-04, 0.055, (1.060, 1.025, 1.095),
     0.058, 0.031, 0.059, 0.875, 0.052, -0.000, 0.109, 68, 0),
    ("OA_knee", "fracture", 69, 0.027, 0.032, 0.395, 0.542, (0.973, 0.914, 1.036),
     0.031, 0.045, 0.490, 0.476, 0.019, -0.112, 0.126, 49, 0),
    ("OA_knee", "BMD_total", 69, 0.031, 0.032, 0.326, 0.481, None,
     0.036, 0.045, 0.420, 0.376, 0.026, -0.032, 0.073, 49, 0),
    ("OA_knee", "BMD_femur_neck", 69, 0.042, 0.032, 0.190, 0.272, None,
     0.048, 0.038, 0.206, 0.169, 0.023, -0.063, 0.122, 30, 0),
]

Z = 1.959963985
n_bad = 0
for (e, o, k, ivw, se, p, p_re, orr, eb, ese, ep, ip, wmd, wlo, whi, i2, loo) in MAN:
    r = mr[(e, o)]
    w = wm[f"{e}__{o}"]
    checks = []
    checks.append(("n_snp", k, int(r["n_snp"])))
    checks.append(("ivw", ivw, float(r["ivw"])))
    checks.append(("se", se, float(r["ivw_se"])))
    checks.append(("p", p, float(r["ivw_p"])))
    checks.append(("p_re", p_re, float(r["ivw_re_p"])))
    checks.append(("egger_b", eb, float(r["egger"])))
    checks.append(("egger_se", ese, float(r["egger_se"])))
    checks.append(("egger_p", ep, float(r["egger_p"])))
    checks.append(("int_p", ip, float(r["egger_intercept_p"])))
    checks.append(("I2", i2, round(float(r["I2"]) * 100)))
    checks.append(("loo", loo, int(r["loo_sign_flips"])))
    checks.append(("wmed", wmd, float(w["wm_recomp"])))
    checks.append(("wm_lo", wlo, float(w["wm_lo"])))
    checks.append(("wm_hi", whi, float(w["wm_hi"])))
    bad = []
    for name, man_v, true_v in checks:
        tol = 0.0006 if name in ("ivw", "se", "egger_b", "egger_se", "wmed", "wm_lo", "wm_hi") else (
            0.06 if name in ("I2",) else 0)
        if name in ("p", "p_re", "egger_p", "int_p"):
            ok = abs(man_v - true_v) <= 0.06 * max(abs(true_v), 1e-300) or man_v == true_v
        elif name in ("n_snp", "loo", "I2"):
            ok = man_v == true_v
        else:
            ok = abs(man_v - true_v) <= tol
        if not ok:
            bad.append(f"{name}: 稿件={man_v} 真源={true_v}")
    # OR 列
    if orr is not None:
        b, s = float(r["ivw"]), float(r["ivw_se"])
        if o == "fracture":
            true_or = (math.exp(-b), math.exp(-b - Z * s), math.exp(-b + Z * s))
        else:
            true_or = (math.exp(b), math.exp(b - Z * s), math.exp(b + Z * s))
        for i, (mv, tv) in enumerate(zip(orr, true_or)):
            if abs(mv - tv) > 0.0015:
                bad.append(f"OR[{i}]: 稿件={mv} 真源={tv:.4f}")
    status = "✓" if not bad else "✗ " + "; ".join(bad)
    if bad:
        n_bad += 1
    print(f"{e} → {o}: {status}")
print(f"\n{n_bad}/10 行存在不一致")
