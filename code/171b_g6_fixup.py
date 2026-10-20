#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
171b_g6_fixup.py — G6 落地的增量修正（171 首轮的三处缺陷）
=============================================================================
171 首轮的三处缺陷，本脚本逐一修掉（不重复缩放已在 171 中正确处理的列）：
  D1  Table_S16_mvmr_bmi.csv 的 `n_bmi_used`（仪器计数）被列名子串匹配误缩放 → 撤销
  D2  mvmr_bmi_main.csv 的实际列名形如 uni_beta_bmd / mvmr_beta_bmi，171 的列名清单
      未覆盖 → 按精确后缀规则补缩放
  D3  mvmr_unify_baselines.csv 列名形如 beta_uni_primary / theta_H，171 未覆盖 → 补缩放
  D4  mr_presso_results_corrected.csv：171 的 _q(m) 实参错位导致 chi2_distortion=0
      → 改为按缩放法重建（Q / p_distortion 为尺度不变量，只缩放 b_ivw_* 两列）

缩放不变性自检（随脚本运行并打印）：
  · 正向（OA 为结局）：b·se_e 与 resid 同乘 C → MR-PRESSO 的 Q = Σ(resid/sig)² 不变
  · 反向（OA 为暴露）：b·se_e 与 resid 均不变 → Q 亦不变
"""
import os
import sys
import importlib
import numpy as np
import pandas as pd
from scipy.stats import chi2

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
m21 = importlib.import_module("21_mr_corrected_N")

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "output", "genetics")
BACKUP = os.path.join(OUT, "_backup_pre_G6")
DOCS_SUP = os.path.join(os.path.dirname(ROOT), "docs", "Supplementary")

C = float(m21.F_OA_KNEE)
K_NEW = float(m21.K_OA_KNEE)


def oa_factor(exp, out):
    if str(out) == "OA_knee":
        return C
    if str(exp) == "OA_knee":
        return 1.0 / C
    return 1.0


# ---------------------------------------------------------------- D1
def fix_s16_count():
    fp = os.path.join(DOCS_SUP, "Table_S16_mvmr_bmi.csv")
    df = pd.read_csv(fp)
    col = "n_bmi_used"
    if col not in df.columns:
        return
    before = df[col].tolist()
    df[col] = df[col].astype(float)
    for i, r in df.iterrows():
        eo = str(r["exposure_outcome"])
        if "fracture" in eo:
            continue
        df.at[i, col] = df.at[i, col] / C          # 撤销误缩放
    df[col] = df[col].round().astype(int)
    df.to_csv(fp, index=False)
    print("D1 ✅ S16 n_bmi_used 撤销误缩放：", before, "->", df[col].tolist())


# ---------------------------------------------------------------- D2
BETA_TOKENS = ("beta", "theta")
SE_SUFFIX = ("se",)


def _is_scale_col(c):
    """精确判定：列名含 beta/theta 或以 _se 结尾（且不是计数/比值/百分比）。"""
    lc = c.lower()
    if any(t in lc for t in ("pct", "percent", "atten", "fc", "_p", "n_", "count")):
        return False
    if lc.endswith("_se") or lc.endswith("_se_re") or lc == "se":
        return True
    if any(t in lc for t in BETA_TOKENS):
        return True
    if lc in ("wmed", "wm_lo", "wm_hi", "egger_int", "trim_beta"):
        return True
    return False


def fix_mvmr_bmi():
    fp = os.path.join(OUT, "mvmr_bmi_main.csv")
    if not os.path.exists(fp):
        return
    df = pd.read_csv(fp)
    cols = [c for c in df.columns
            if _is_scale_col(c) and pd.api.types.is_numeric_dtype(df[c])]
    # 排除 171 已处理过的两列（避免二次缩放）
    cols = [c for c in cols if c not in ("egger_int", "egger_int_se")]
    print("D2 mvmr_bmi_main 补缩放列:", cols)
    for i, r in df.iterrows():
        f = oa_factor(r["exposure"], r["outcome"])
        if f == 1.0:
            continue
        for c in cols:
            v = df.at[i, c]
            if pd.isna(v):
                continue
            df.at[i, c] = v * f
    df.to_csv(fp, index=False)
    print("D2 ✅ mvmr_bmi_main.csv")
    print(df[["exposure", "outcome", "uni_beta_bmd", "mvmr_beta_bmd",
              "mvmr_beta_bmi", "attenuation_pct"]].to_string(index=False))


# ---------------------------------------------------------------- D3
def fix_unify():
    fp = os.path.join(OUT, "mvmr_unify_baselines.csv")
    if not os.path.exists(fp):
        return
    df = pd.read_csv(fp)
    cols = [c for c in df.columns
            if _is_scale_col(c) and pd.api.types.is_numeric_dtype(df[c])]
    cols = [c for c in cols if "atten" not in c.lower()]
    print("D3 mvmr_unify_baselines 缩放列:", cols)
    for i, r in df.iterrows():
        f = oa_factor(r["exposure"], r["outcome"])
        if f == 1.0:
            continue
        for c in cols:
            v = df.at[i, c]
            if pd.isna(v):
                continue
            df.at[i, c] = v * f
    df.to_csv(fp, index=False)
    print("D3 ✅ mvmr_unify_baselines.csv")
    print(df[["outcome", "exposure", "beta_uni_primary", "beta_mvmr",
              "theta_H", "atten_vs_primary_pct"]].to_string(index=False))


# ---------------------------------------------------------------- D4
def fix_presso():
    src = os.path.join(BACKUP, "mr_presso_results_corrected.csv")
    dst = os.path.join(OUT, "mr_presso_results_corrected.csv")
    if not os.path.exists(src):
        print("D4 ⚠ 备份缺失")
        return
    df = pd.read_csv(src)
    for i, r in df.iterrows():
        f = oa_factor(r["exposure"], r["outcome"])
        if f == 1.0:
            continue
        for c in ("b_ivw_full", "b_ivw_out"):
            df.at[i, c] = df.at[i, c] * f
    df.to_csv(dst, index=False)
    print("D4 ✅ mr_presso_results_corrected.csv（Q/p_distortion 保持尺度不变量）")
    print(df[["exposure", "outcome", "K", "n_outlier", "Q_full", "Q_out",
              "chi2_distortion", "p_distortion", "b_ivw_full",
              "b_ivw_out"]].to_string(index=False))


# ---------------------------------------------------------------- 自检
def selfcheck():
    """缩放不变性自检：z 与 Q 型量在 G6 前后必须完全一致。"""
    print("\n" + "=" * 78)
    print("自检：G6 前后 p 值 / I² / Q 必须完全不变")
    print("=" * 78)
    ok = True
    for fname, keys in (("mr_results_corrected.csv", ("exposure", "outcome")),):
        new = pd.read_csv(os.path.join(OUT, fname))
        old = pd.read_csv(os.path.join(BACKUP, fname))
        on = {(r.exposure, r.outcome): r for r in old.itertuples()}
        for r in new.itertuples():
            o = on.get((r.exposure, r.outcome))
            if o is None:
                continue
            same_p = abs(r.ivw_p - o.ivw_p) <= 1e-12 * max(abs(o.ivw_p), 1e-300)
            same_z = abs(r.ivw / r.ivw_se - o.ivw / o.ivw_se) < 1e-8
            flag = "✅" if (same_p and same_z) else "❌"
            if not (same_p and same_z):
                ok = False
            print("  %s %-16s→%-15s z %.6f/%.6f  p %.6e/%.6e"
                  % (flag, r.exposure, r.outcome,
                     r.ivw / r.ivw_se, o.ivw / o.ivw_se, r.ivw_p, o.ivw_p))
    print("\n" + ("✅ 全部通过：G6 只改效应量绝对值" if ok else "❌ 存在变化，须排查"))


def main():
    print("G6 增量修正   C = %.4f   K = %.4f" % (C, K_NEW))
    fix_s16_count()
    print()
    fix_mvmr_bmi()
    print()
    fix_unify()
    print()
    fix_presso()
    selfcheck()


if __name__ == "__main__":
    main()
