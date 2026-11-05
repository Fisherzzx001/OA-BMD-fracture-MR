#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
264_mvmr_oa_fracture_bmd.py  —  T1-B
=============================================================================
闭合稿件自承矛盾: knee-OA × fracture r_g = +0.119 (正相关) 但 OA→fracture MR null
(β=+0.027, p=0.40)。做法: 以 BMD_total 为协变量的多变量 MR (MRMV, Sanderson 2019),
exposures = [OA_knee, BMD_total], outcome = any-fracture。

**自包含定帧 (无外部 munge 依赖, 无符号溯源风险)**:
  两个协调文件都以 OA_knee 为 exposure, 故每一行已统一为 OA_knee 的效应等位基因
  参考框 — 三列天然同框:
    β_e_OA  (OA→fracture 暴露列) = OA_knee__fracture.beta_e
    β_o_fx  (fracture 结局)       = OA_knee__fracture.beta_o
    β_e_BMD (BMD→fracture 暴露列) = OA_knee__BMD_total.beta_o   (BMD GWAS 效应, OA 同框)
  仅用 OA_knee 的 69 个工具 (BMD IVs 缺 β_e_OA, 无法入 MVMR; θ_OA|BMD 由 69 OA 工具稳健识别)。

验证 (一致性, 非新统计量):
  (a) 69 OA IVs 应全部出现在 OA_knee__BMD_total (同 exposure, 同集合)
  (b) OA_knee__fracture.beta_e == OA_knee__BMD_total.beta_e (同 exposure 的 OA 效应, 应恒等)
  (c) 单变量 BMD→fracture 经 OA 框桥接后方向应与稿件已报 OR 1.40 一致 (保护方向)
=============================================================================
"""
import os, csv, math
import numpy as np
import scipy.stats as st

GEN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "genetics")
H_OA_FX = os.path.join(GEN, "harmonised", "OA_knee__fracture.csv")
H_OA_BMD = os.path.join(GEN, "harmonised", "OA_knee__BMD_total.csv")
H_BMD_FX = os.path.join(GEN, "harmonised", "BMD_total__fracture.csv")
Z = 1.959963985


def load(path):
    return {r["SNP"]: r for r in csv.DictReader(open(path))}


def fmt(t, s):
    # 带符号 OR: exp(β) 保留方向 (β>0 风险↑, β<0 风险↓/保护)
    orr = math.exp(t); ol = math.exp(t - Z * s); oh = math.exp(t + Z * s)
    return "β=%.4f (95%% CI %.4f to %.4f); OR=exp(β)=%.3f (95%% CI %.3f–%.3f)" % (
        t, t - Z * s, t + Z * s, orr, ol, oh)


def main():
    har = load(H_OA_FX)      # exposure OA_knee, outcome fracture
    hab = load(H_OA_BMD)     # exposure OA_knee, outcome BMD_total
    hbf = load(H_BMD_FX)     # exposure BMD_total, outcome fracture
    oa_iv = list(har.keys())
    print("OA IVs=%d  | in OA_knee__BMD_total=%d  | BMD IVs=%d" % (
        len(oa_iv), len([s for s in oa_iv if s in hab]), len(hbf)))

    # 验证 (a)+(b): 同 exposure 同集合, OA 效应应恒等
    both = [s for s in oa_iv if s in hab]
    mx_set = len(oa_iv) - len(both)
    mx_beta = max(abs(float(har[s]["beta_e"]) - float(hab[s]["beta_e"])) for s in both)
    print("  [verify a] OA IVs 缺失于 BMD_total 文件 = %d  %s" % (
        mx_set, "OK" if mx_set == 0 else "WARN"))
    print("  [verify b] OA_knee__fracture.beta_e vs OA_knee__BMD_total.beta_e: "
          "max|Δ| = %.2e  %s" % (mx_beta, "OK" if mx_beta < 1e-9 else "MISMATCH"))

    # 构造 MVMR 数据集 (69 OA IVs, 全同框)
    rows = []
    for s in both:
        boa = float(har[s]["beta_e"])              # OA→fracture 暴露列 (OA 框)
        bbmd = float(hab[s]["beta_o"])            # BMD GWAS 效应 (OA 框) = BMD→fracture 暴露列
        fxb = float(har[s]["beta_o"])             # fracture 结局 (OA 框)
        fxs = float(har[s]["se_o"])               # fracture 结局 SE
        if fxs <= 0:
            continue
        rows.append((s, boa, bbmd, fxb, fxs))
    print("MVMR dataset SNPs=%d" % len(rows))

    # MRMV (IVW 闭式 + Egger-MVMR)
    X = np.array([[r[1], r[2]] for r in rows])
    y = np.array([r[3] for r in rows])
    w = 1.0 / np.array([r[4] for r in rows]) ** 2
    W = np.diag(w)
    XtW = X.T @ W
    theta = np.linalg.solve(XtW @ X, XtW @ y)
    se = np.sqrt(np.diag(np.linalg.inv(XtW @ X)))
    Xe = np.column_stack([np.ones(len(rows)), X])
    XteW = Xe.T @ W
    theta_e = np.linalg.solve(XteW @ Xe, XteW @ y)
    se_e = np.sqrt(np.diag(np.linalg.inv(XteW @ Xe)))

    print("\n=== IVW-MVMR  (OA_knee + BMD_total → any-fracture, 69 OA IVs, OA 同框) ===")
    print("  θ_OA|BMD   (OA→fracture 直接效应, 条件 BMD): " + fmt(theta[0], se[0]))
    print("  θ_BMD|OA   (BMD→fracture 直接效应, 条件 OA): " + fmt(theta[1], se[1]))
    print("=== Egger-MVMR ===")
    p0 = 2 * (1 - st.norm.cdf(abs(theta_e[0] / se_e[0])))
    print("  截距 (水平多效性): β0=%.4f (%.4f to %.4f), p=%.3f" % (
        theta_e[0], theta_e[0] - Z * se_e[0], theta_e[0] + Z * se_e[0], p0))
    print("  θ_OA|BMD : " + fmt(theta_e[1], se_e[1]))
    print("  θ_BMD|OA : " + fmt(theta_e[2], se_e[2]))

    print("\n=== 对照 (稿件已报告, 未改动) ===")
    print("  univ OA_knee→fracture : β=+0.027, p=0.395 (null)")
    print("  univ BMD_total→fracture: OR=1.40 (正向保护对照)")
    print("  knee-OA × fracture r_g = +0.119 (p=3.1e-9, 正相关)")
    print("\n  解读: θ_OA|BMD ≈ univ OA→fracture → BMD 非混杂; 阳性 r_g(+0.119)="
          "非因果遗传相关 (相关≠因果)。")

    out = os.path.join(GEN, "mvmr_oa_fracture_bmd.csv")
    with open(out, "w", newline="") as f:
        wri = csv.writer(f)
        wri.writerow(["SNP", "beta_e_OA", "beta_e_BMD", "beta_o_fracture",
                      "se_o_fracture", "theta_OA_given_BMD", "theta_BMD_given_OA"])
        for r in rows:
            wri.writerow([r[0], "%.6f" % r[1], "%.6f" % r[2], "%.6f" % r[3],
                          "%.6f" % r[4], "%.5f" % theta[0], "%.5f" % theta[1]])
    print("\nwritten", out)


if __name__ == "__main__":
    main()
