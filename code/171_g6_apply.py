#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
171_g6_apply.py — G6 落地：OA 结局的 case-control（有效样本量）尺度修正
=============================================================================
背景（审计 docs/123 §2）
----------------------
code/13 load_munged 对二分类结局 OA_knee 用**总**样本量 N 做 Z 线性重构：

    se = 1/sqrt(N · 2 f (1-f)),  beta = Z · se

但二分类结局的 log-OR 满足 se = 1/sqrt(N_eff · 2 f (1-f))，N_eff = N·K(1-K)。
故重构出的效应落在线性概率尺度，比真实 log-OR 小 1/sqrt(K(1-K)) 倍。
骨折结局因为原始文件直接给 log-OR 的 beta/se 而幸免 —— 同一条流水线对两个
结局用了不一致的尺度处理。

修正因子（实测，非假设）
------------------------
code/_probe_g6b.py：39 个 SNP 的真实 log-OR 与 munged Z 逐点比对
  |Z_true / Z_munged| = 1.000000 (std 5.3e-14)  —— 确认同源
  |beta_true| = |Z_munged| × se_true (误差 < 1e-16)
  → K(1-K) = N_eff_true / N_munged，逐 SNP 中位 0.087082 (IQR 0.0835–0.0881)
  → K = 0.0964 ； 因子 C = 1/sqrt(K(1-K)) = 3.3882

已在 code/21_mr_corrected_N.py 中定义为 F_OA_KNEE 并由 code/26 施加到
harmonised/*.csv（2026-10-20 重跑，10/10 一致性校验通过）。

本脚本做什么
------------
1) 以 harmonised（已修正、严格 clump）为唯一真源，重算
   mr_results_corrected.csv（此前是 1-Mb clumping 的旧源，K=79 而非 69 —— G15）
   mr_presso_results_corrected.csv
2) 对其余涉 OA 的产物按「OA 为结局 ×C / OA 为暴露 ÷C」缩放 β 与 SE 类列，
   p 值、z、Q、I²、F、Fc、衰减百分比一概不动（尺度不变量）
3) 重算 OA_knee 的 liability-scale h²（K 由 0.12 改为实测 0.0964）
4) 输出 G6_change_table.csv（每个受影响数字的新旧对照，供稿件逐项回填）

用法：python 171_g6_apply.py
"""
import os
import sys
import importlib
import numpy as np
import pandas as pd
from scipy.stats import norm, chi2

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
m13 = importlib.import_module("13_mr_analysis")
m21 = importlib.import_module("21_mr_corrected_N")
m26 = importlib.import_module("26_export_harmonised_snp")

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "output", "genetics")
HARM = os.path.join(OUT, "harmonised")
DOCS_SUP = os.path.join(os.path.dirname(ROOT), "docs", "Supplementary")

C = float(m21.F_OA_KNEE)                       # 3.3882 ...
K_NEW = float(m21.K_OA_KNEE)                   # 0.0964
PAIRS = m13.PAIRS

CHANGES = []


def _log(fname, col, old, new, note=""):
    CHANGES.append(dict(file=fname, field=col, old=old, new=new, note=note))


def oa_factor(exp, out):
    """OA_knee 为结局 → ×C；OA_knee 为暴露 → ÷C；不涉及 OA → 1。"""
    if out == "OA_knee":
        return C
    if exp == "OA_knee":
        return 1.0 / C
    return 1.0


# ------------------------------------------------------------------ 1. MR 主表
def rebuild_mr_results():
    """从 harmonised（严格 clump + G6 修正）重算完整统计量。"""
    rows = []
    for exp, out, note in PAIRS:
        fp = os.path.join(HARM, f"{exp}__{out}.csv")
        if not os.path.exists(fp):
            continue
        h = pd.read_csv(fp)
        res = m13.run_mr(h.beta_e.values, h.se_e.values,
                         h.beta_o.values, h.se_o.values, snps=h.SNP.tolist())
        rows.append(dict(exposure=exp, outcome=out, note=note, **{
            k: v for k, v in res.items() if k != "snps"}))
    df = pd.DataFrame(rows)
    # 列顺序与既有 mr_results_corrected.csv 保持一致（下游兼容）
    cols = ["exposure", "outcome", "note", "n_snp", "ivw", "ivw_se", "ivw_z", "ivw_p",
            "egger", "egger_se", "egger_p", "egger_intercept", "egger_intercept_se",
            "egger_intercept_p", "wmed", "Q", "Q_df", "Q_p", "I2", "ivw_re_se",
            "ivw_re_p", "n_outlier", "ivw_radial", "ivw_radial_se", "loo_min",
            "loo_max", "loo_sign_flips"]
    df = df[[c for c in cols if c in df.columns]]
    df.to_csv(os.path.join(OUT, "mr_results_corrected.csv"), index=False)
    print("✅ mr_results_corrected.csv  (K=%s)" % list(df.n_snp))
    return df


def rebuild_presso():
    rows = []
    for exp, out, note in PAIRS:
        fp = os.path.join(HARM, f"{exp}__{out}.csv")
        if not os.path.exists(fp):
            continue
        h = pd.read_csv(fp)
        be = h.beta_e.values.astype(float)
        se_e = h.se_e.values.astype(float)
        bo = h.beta_o.values.astype(float)
        so = h.se_o.values.astype(float)

        def _q(b, m=None):
            m = np.ones(len(be), bool) if m is None else m
            bb = m26.ivw(be[m], bo[m], so[m])
            sig = np.sqrt(so[m] ** 2 + (bb * se_e[m]) ** 2)
            return float(np.sum(((bo[m] - bb * be[m]) / sig) ** 2))

        p_i, outl = m26.presso_per_snp(be, se_e, bo, so)
        q_full = _q(None)
        keep = ~outl
        rows.append(dict(
            exposure=exp, outcome=out, note=note, K=len(be), Q_full=q_full,
            n_outlier=int(outl.sum()),
            b_ivw_full=m26.ivw(be, bo, so),
            b_ivw_out=m26.ivw(be[keep], bo[keep], so[keep]) if keep.sum() >= 3 else np.nan,
            Q_out=_q(keep) if keep.sum() >= 3 else np.nan,
            chi2_distortion=float(q_full - (_q(keep) if keep.sum() >= 3 else q_full)),
            p_distortion=float(chi2.sf(max(q_full - (_q(keep) if keep.sum() >= 3 else q_full), 0.0), 1)),
        ))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "mr_presso_results_corrected.csv"), index=False)
    print("✅ mr_presso_results_corrected.csv")
    return df


# ------------------------------------------------------------------ 2. 通用缩放
SCALE_SPECS = [
    # (文件名, 主键对 (exp列, out列), 需要缩放的列)
    ("table3_corrected.csv", ("exp", "out"),
     ["beta", "se", "egger_beta", "wmed"]),
    ("mr_power.csv", ("exposure", "outcome"),
     ["beta", "se", "MDE_beta_80pct"]),
    ("wm_recomputed_vs_tabulated.csv", (None, None), None),   # 特殊：按 pair 名
    ("wm_bootstrap_ci.csv", (None, None), None),
    ("mvmr_height_main.csv", ("exposure", "outcome"),
     ["uni_beta", "uni_se", "uni_bmd_only", "uni_bmd_only_se",
      "mvmr_beta", "mvmr_se", "mvmr_se_re", "mvmr_h", "mvmr_h_se",
      "uni_h_beta", "uni_h_se", "egger_int", "egger_int_se",
      "egger_beta", "egger_se", "trim_beta", "trim_se",
      "ref_uni_beta", "ref_uni_se"]),
    ("mvmr_height_n_h1000.csv", ("exposure", "outcome"),
     ["uni_beta", "uni_se", "uni_bmd_only", "uni_bmd_only_se",
      "mvmr_beta", "mvmr_se", "mvmr_se_re", "mvmr_h", "mvmr_h_se",
      "uni_h_beta", "uni_h_se", "egger_int", "egger_int_se",
      "egger_beta", "egger_se", "trim_beta", "trim_se",
      "ref_uni_beta", "ref_uni_se"]),
    ("mvmr_height_grid.csv", ("exposure", "outcome"),
     ["uni_beta", "uni_se", "uni_bmd_only", "uni_bmd_only_se",
      "mvmr_beta", "mvmr_se", "mvmr_se_re", "mvmr_h", "mvmr_h_se",
      "uni_h_beta", "uni_h_se", "egger_int", "egger_int_se",
      "egger_beta", "egger_se", "trim_beta", "trim_se",
      "ref_uni_beta", "ref_uni_se"]),   # G6 补丁：ref 列曾漏缩放（2026-10-20）
    ("mvmr_bmi_main.csv", ("exposure", "outcome"),
     ["uni_beta", "uni_se", "mvmr_beta", "mvmr_se", "mvmr_se_re",
      "mvmr_bmi", "mvmr_bmi_se", "egger_int", "egger_int_se",
      "egger_beta", "egger_se"]),
    ("mvmr_unify_baselines.csv", ("exposure", "outcome"),
     ["uni_beta", "uni_se", "mvmr_beta", "mvmr_se", "mvmr_h", "mvmr_h_se"]),
]


def scale_generic(fname, expcol, outcol, cols):
    fp = os.path.join(OUT, fname)
    if not os.path.exists(fp):
        print("  ⚠ 缺失", fname)
        return
    df = pd.read_csv(fp)
    cols = [c for c in cols if c in df.columns]
    if not cols:
        print("  ⚠ 无可缩放列", fname)
        return
    for i, r in df.iterrows():
        f = oa_factor(r[expcol], r[outcol])
        if f == 1.0:
            continue
        for c in cols:
            old = df.at[i, c]
            if pd.isna(old):
                continue
            df.at[i, c] = old * f
            _log(fname, f"{r[expcol]}→{r[outcol]}:{c}", old, old * f, f"×{f:.4f}")
    df.to_csv(fp, index=False)
    print("✅ %s (%d 列缩放)" % (fname, len(cols)))


def scale_wm(fname):
    """wm_* 文件：pair 名形如 '<exp>__<out>'。"""
    fp = os.path.join(OUT, fname)
    if not os.path.exists(fp):
        print("  ⚠ 缺失", fname)
        return
    df = pd.read_csv(fp)
    keycol = "pair" if "pair" in df.columns else df.columns[0]
    cols = [c for c in df.columns
            if c in ("ivw_recomp", "se_recomp", "egger_recomp", "egger_se",
                     "wm_recomp", "wm_lo", "wm_hi", "wm_tab")]
    hit = 0
    for i, r in df.iterrows():
        key = str(r[keycol])
        if "__" not in key:
            continue
        exp, out = key.split("__", 1)
        f = oa_factor(exp, out)
        if f == 1.0:
            continue
        hit += 1
        for c in cols:
            old = df.at[i, c]
            if pd.isna(old):
                continue
            df.at[i, c] = old * f
            _log(fname, f"{key}:{c}", old, old * f, f"×{f:.4f}")
    df.to_csv(fp, index=False)
    print("✅ %s (%d 行缩放)" % (fname, hit))


# ------------------------------------------------------------------ 3. h²_liab
def recompute_liability():
    fp = os.path.join(OUT, "oa_h2_liability.csv")
    if not os.path.exists(fp):
        return
    df = pd.read_csv(fp)
    phi = norm.pdf(norm.ppf(1 - np.array(df.K_primary)))
    old = df.h2_liab_primary.copy()
    new = df.h2_obs * (df.K_primary * (1 - df.K_primary)) / phi ** 2
    df["h2_liab_primary"] = new
    # OA_knee 的 K 换成实测值
    m = df.trait == "OA_knee"
    if m.any():
        k = K_NEW
        ph = norm.pdf(norm.ppf(1 - k))
        df.loc[m, "K_primary"] = k
        df.loc[m, "h2_liab_primary"] = float(df.loc[m, "h2_obs"].iloc[0]) * k * (1 - k) / ph ** 2
        _log("oa_h2_liability.csv", "OA_knee K", 0.12, k, "实测校准")
        _log("oa_h2_liability.csv", "OA_knee h2_liab",
             float(old[m].iloc[0]), float(df.loc[m, "h2_liab_primary"].iloc[0]), "")
    df.to_csv(fp, index=False)
    print("✅ oa_h2_liability.csv  OA_knee h²_liab=%.4f (K=%.4f)"
          % (float(df.loc[m, "h2_liab_primary"].iloc[0]), K_NEW))


# ------------------------------------------------------------------ 4. S16
def scale_s16():
    for fname in ("Table_S16_mvmr_bmi.csv",):
        fp = os.path.join(DOCS_SUP, fname)
        if not os.path.exists(fp):
            print("  ⚠ 缺失", fname)
            continue
        df = pd.read_csv(fp)
        cols = [c for c in df.columns if c not in
                ("exposure", "outcome", "K", "n_bmd", "n_bmi", "atten_pct",
                 "p", "p_re", "egger_p", "Fc_bmd", "Fc_bmi", "I2")]
        # 只缩放数值列中的 beta/se 型
        cols = [c for c in cols if pd.api.types.is_numeric_dtype(df[c])
                and any(t in c.lower() for t in ("beta", "se", "theta", "ci"))]
        print("   S16 缩放列:", cols)
        for i, r in df.iterrows():
            f = oa_factor(r.get("exposure"), r.get("outcome"))
            if f == 1.0:
                continue
            for c in cols:
                old = df.at[i, c]
                if pd.isna(old):
                    continue
                df.at[i, c] = old * f
                _log(fname, f"{r.get('exposure')}→{r.get('outcome')}:{c}", old, old * f, "")
        df.to_csv(fp, index=False)
        print("✅", fname)


# ------------------------------------------------------------------ main
def main():
    print("=" * 78)
    print("G6 落地：OA_knee case-control 尺度修正   C = 1/sqrt(K(1-K)) = %.4f  (K=%.4f)"
          % (C, K_NEW))
    print("=" * 78)

    print("\n▶ 1) 从 harmonised（严格 clump + 已修正）重算 MR 主表")
    mr = rebuild_mr_results()
    print(mr[["exposure", "outcome", "n_snp", "ivw", "ivw_se", "ivw_p",
              "ivw_re_p", "I2"]].to_string(index=False))

    print("\n▶ 2) 重算 MR-PRESSO（离群剔除后 IVW）")
    pr = rebuild_presso()
    print(pr[["exposure", "outcome", "K", "n_outlier", "b_ivw_full",
              "b_ivw_out", "p_distortion"]].to_string(index=False))

    print("\n▶ 3) 缩放其余涉 OA 的产物")
    for fname, (ec, oc), cols in SCALE_SPECS:
        if cols is None:
            scale_wm(fname)
        else:
            scale_generic(fname, ec, oc, cols)

    print("\n▶ 4) liability-scale h² 重算（K 0.12 → 实测 %.4f）" % K_NEW)
    recompute_liability()

    print("\n▶ 5) Supplementary 表")
    scale_s16()

    ch = pd.DataFrame(CHANGES)
    ch.to_csv(os.path.join(OUT, "G6_change_table.csv"), index=False)
    print("\n📁 变更对照表 -> %s (%d 项)" % (os.path.join(OUT, "G6_change_table.csv"), len(ch)))

    # ---- 关键数字打印（供稿件回填）----
    print("\n" + "=" * 78)
    print("稿件 Table 2 / Table 3 回填值")
    print("=" * 78)
    for r in mr.itertuples():
        if r.outcome == "OA_knee":
            lo, hi = np.exp(r.ivw - 1.96 * r.ivw_se), np.exp(r.ivw + 1.96 * r.ivw_se)
            print("  %-16s→OA_knee  K=%d  β=%+.3f (SE %.3f)  p=%.2e  OR=%.3f (%.3f–%.3f)"
                  % (r.exposure, r.n_snp, r.ivw, r.ivw_se, r.ivw_p, np.exp(r.ivw), lo, hi))
        elif r.exposure == "OA_knee":
            print("  OA_knee→%-16s K=%d  β=%+.3f (SE %.3f)  p=%.3f  MDE=%.3f"
                  % (r.outcome, r.n_snp, r.ivw, r.ivw_se, r.ivw_p, 2.8 * r.ivw_se))
        else:
            print("  %-16s→fracture  K=%d  β=%+.3f (SE %.3f)  p=%.2e  OR/SD_lower=%.3f"
                  % (r.exposure, r.n_snp, r.ivw, r.ivw_se, r.ivw_p, np.exp(-r.ivw)))


if __name__ == "__main__":
    main()
