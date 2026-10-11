#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUDIT_MR.py —— 全面复盘审计：从逐 SNP harmonised 数据独立重算全部 MR 统计量。

原则：不 import 任何原分析脚本（13/17/21/24），全部公式在此重新实现，
      与交付文件（table3_corrected / mr_corrected_N）逐值对照。

重算内容
--------
  1. beta 尺度反推：beta = z/sqrt(2pqN)  -> 反解 N，判断该性状 β 的标度是否正确
  2. IVW（固定效应，加权最小二乘过原点）
  3. IVW 随机效应 SE（DerSimonian-Laird 式：se_re = se_fe * sqrt(max(1, Q/df))）
  4. Cochran Q / I²
  5. MR-Egger（过原点逆方差加权回归，误差在两轴，用标准 Egger 回归）
  6. Egger 截距及其 p
  7. 加权中位数（ penalised / simple weighted median ）
  8. 留一法（LOO）范围与符号翻转计数
  9. MR-PRESSO 全局检验的 Q（残差平方和）与 outlier 计数（按 residual SS > chi2 阈值）
"""
import os
import sys
import math
import numpy as np
import pandas as pd
from scipy.stats import norm, chi2

ROOT = os.path.dirname(os.path.abspath(__file__))
HARM = os.path.join(ROOT, "output", "genetics", "harmonised")
OUTG = os.path.join(ROOT, "output", "genetics")

# 真实 N（外部元数据，用于尺度反推）
N_TRUE = {
    "OA_knee": 1_208_397,
    "OA_selfrep": 462_933,
    "fracture": 264_973,
}
N_BMD_TRUE = 31_873.0
N_BMD_PLACEHOLDER = 400_000.0


# ---------------------------------------------------------------- 基础统计
def ivw_fixed(be, se_e, bo, se_o):
    """IVW 固定效应。ratio_i = bo/be, w_i = 1/se_ratio^2, se_ratio 用 delta 法。"""
    ratio = bo / be
    se_ratio = se_o / np.abs(be)
    w = 1.0 / se_ratio ** 2
    b = np.sum(w * ratio) / np.sum(w)
    se = np.sqrt(1.0 / np.sum(w))
    z = b / se
    p = 2 * norm.sf(abs(z))
    return b, se, z, p


def cochran_q(be, se_e, bo, se_o, b):
    ratio = bo / be
    se_ratio = se_o / np.abs(be)
    w = 1.0 / se_ratio ** 2
    q = float(np.sum(w * (ratio - b) ** 2))
    df = len(ratio) - 1
    p = chi2.sf(q, df) if df > 0 else np.nan
    i2 = max(0.0, (q - df) / q) if q > 0 else 0.0
    return q, df, p, i2


def ivw_random(be, se_e, bo, se_o, q, df):
    """随机效应：se_re = se_fe * sqrt(max(1, Q/df))。"""
    b, se_fe, _, _ = ivw_fixed(be, se_e, bo, se_o)
    if not np.isfinite(q) or df <= 0:
        return se_fe, np.nan
    mult = math.sqrt(max(1.0, q / df))
    se_re = se_fe * mult
    z = b / se_re
    p = 2 * norm.sf(abs(z))
    return se_re, p


def mr_egger(be, se_e, bo, se_o):
    """
    标准 MR-Egger (Bowden 2015)：
        以 beta_o 为因变量、beta_e 为自变量做加权线性回归，
        权重 w_i = 1/se_o^2（Egger 用第二轴的误差定权）。
    返回 斜率、斜率SE、斜率p、截距、截距SE、截距p。
    """
    w = 1.0 / se_o ** 2
    x, y = be, bo
    Sw = np.sum(w)
    xbar = np.sum(w * x) / Sw
    ybar = np.sum(w * y) / Sw
    Sxx = np.sum(w * (x - xbar) ** 2)
    b1 = np.sum(w * (x - xbar) * (y - ybar)) / Sxx
    b0 = ybar - b1 * xbar
    # 残差方差（不除以 Sw-2 的自由度校正时用 psi 估计）
    resid = y - (b0 + b1 * x)
    psi = np.sum(w * resid ** 2) / (len(x) - 2)
    se_b1 = np.sqrt(psi / Sxx)
    se_b0 = np.sqrt(psi / Sw + psi * xbar ** 2 / Sxx)
    p1 = 2 * norm.sf(abs(b1 / se_b1))
    p0 = 2 * norm.sf(abs(b0 / se_b0))
    return b1, se_b1, p1, b0, se_b0, p0


def weighted_median(be, se_e, bo, se_o):
    """简单加权中位数（Bowden 2016 的 weighted median, 非 penalised 版）。"""
    ratio = bo / be
    se_ratio = se_o / np.abs(be)
    w = 1.0 / se_ratio ** 2
    order = np.argsort(ratio)
    r, ww = ratio[order], w[order]
    cum = np.cumsum(ww) - 0.5 * ww
    tot = np.sum(ww)
    idx = np.searchsorted(cum / tot, 0.5)
    idx = min(idx, len(r) - 1)
    return float(r[idx])


def loo(be, se_e, bo, se_o):
    n = len(be)
    bs = []
    for i in range(n):
        m = np.ones(n, bool)
        m[i] = False
        b, _, _, _ = ivw_fixed(be[m], se_e[m], bo[m], se_o[m])
        bs.append(b)
    bs = np.array(bs)
    full, _, _, _ = ivw_fixed(be, se_e, bo, se_o)
    flips = int(np.sum(np.sign(bs) != np.sign(full)))
    return float(bs.min()), float(bs.max()), flips, bs


def presso_global(be, se_e, bo, se_o):
    """
    MR-PRESSO 全局检验：Q = sum( (bo - b*be)^2 / se_o^2 )，df = K-1，
    用卡方（而非置换）近似报 p；outlier 按 |resid|/se_o > 阈值 判定。
    """
    b, _, _, _ = ivw_fixed(be, se_e, bo, se_o)
    resid = (bo - b * be) / se_o
    q = float(np.sum(resid ** 2))
    df = len(be) - 1
    p = chi2.sf(q, df) if df > 0 else np.nan
    return q, df, p, np.abs(resid)


# ---------------------------------------------------------------- 尺度反推
def infer_scale(beta, pval, eaf, label):
    """由 beta = z/sqrt(2pq*N) 反解 N（若 beta 为标准 per-allele 效应）。"""
    if not np.isfinite(pval) or pval <= 0 or pval >= 1:
        return np.nan
    z = abs(norm.isf(pval / 2.0))
    if not np.isfinite(z) or beta == 0:
        return np.nan
    pq2 = 2.0 * eaf * (1.0 - eaf)
    if pq2 <= 0:
        return np.nan
    return (z / abs(beta)) ** 2 / pq2


def main():
    print("=" * 92)
    print("全面审计 A：从逐 SNP harmonised 数据独立重算全部 MR 统计量")
    print("=" * 92)

    files = sorted(f for f in os.listdir(HARM)
                   if f.endswith(".csv") and not f.startswith("_"))
    rows = []
    scale_rows = []

    for fn in files:
        base = fn[:-4]
        exp, out = base.split("__")
        d = pd.read_csv(os.path.join(HARM, fn))
        n = len(d)
        be = d.beta_e.values.astype(float)
        se_e = d.se_e.values.astype(float)
        bo = d.beta_o.values.astype(float)
        se_o = d.se_o.values.astype(float)
        eaf = d.eaf.values.astype(float)

        # ---- 尺度反推（对暴露与结局各取 |beta| 最大的几个 SNP，稳健）----
        for which, bb, pp in (("exposure", be, d.p_e.values.astype(float)),
                              ("outcome", bo, d.p_o.values.astype(float))):
            ns = []
            for i in np.argsort(-np.abs(bb))[:5]:
                v = infer_scale(bb[i], pp[i], eaf[i], which)
                if np.isfinite(v):
                    ns.append(v)
            if ns:
                scale_rows.append(dict(
                    pair=base, axis=which,
                    trait=exp if which == "exposure" else out,
                    N_implied_median=float(np.median(ns)),
                    N_implied_min=float(np.min(ns)),
                    N_implied_max=float(np.max(ns))))

        # ---- 统计量重算 ----
        b, se, z, p = ivw_fixed(be, se_e, bo, se_o)
        q, qdf, qp, i2 = cochran_q(be, se_e, bo, se_o, b)
        se_re, p_re = ivw_random(be, se_e, bo, se_o, q, qdf)
        e_b, e_se, e_p, e_i, e_ise, e_ip = mr_egger(be, se_e, bo, se_o)
        wm = weighted_median(be, se_e, bo, se_o)
        lo_min, lo_max, flips, _ = loo(be, se_e, bo, se_o)
        pq, pdf_, pp_, absres = presso_global(be, se_e, bo, se_o)

        rows.append(dict(
            pair=base, exposure=exp, outcome=out, K=n,
            beta=b, se=se, z=z, p=p,
            Q=q, Q_df=qdf, Q_p=qp, I2=i2,
            se_re=se_re, p_re=p_re,
            egger_beta=e_b, egger_se=e_se, egger_p=e_p,
            egger_int=e_i, egger_int_se=e_ise, egger_int_p=e_ip,
            wmed=wm,
            loo_min=lo_min, loo_max=lo_max, loo_flips=flips,
            presso_Q=pq, presso_Q_p=pp_,
        ))

        print("\n【%s】 K=%d" % (base, n))
        print("  IVW        β=%+.6f  SE=%.6f  z=%+.3f  p=%.4g" % (b, se, z, p))
        print("  IVW-RE     SE=%.6f  p=%.4g" % (se_re, p_re))
        print("  Cochran    Q=%.3f df=%d p=%.4g  I²=%.4f" % (q, qdf, qp, i2))
        print("  MR-Egger   β=%+.6f SE=%.6f p=%.4g | 截距=%+.6g SE=%.4g p=%.4g"
              % (e_b, e_se, e_p, e_i, e_ise, e_ip))
        print("  W-median   β=%+.6f" % wm)
        print("  LOO        [%.6f, %.6f] 符号翻转=%d" % (lo_min, lo_max, flips))
        print("  PRESSO-Q   Q=%.3f p=%.4g (卡方近似)" % (pq, pp_))

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUTG, "AUDIT_mr_recomputed.csv"), index=False)

    sd = pd.DataFrame(scale_rows)
    sd.to_csv(os.path.join(OUTG, "AUDIT_scale_inference.csv"), index=False)

    print("\n" + "=" * 92)
    print("B. β 尺度反推（beta = z/sqrt(2pq·N) 反解 N；判断每个性状 β 的真实标度）")
    print("=" * 92)
    for r in sd.itertuples():
        print("  %-28s %-9s %-16s N≈%10.0f  [%10.0f – %10.0f]"
              % (r.pair, r.axis, r.trait, r.N_implied_median,
                 r.N_implied_min, r.N_implied_max))

    # ---- 与交付表对照 ----
    print("\n" + "=" * 92)
    print("C. 与交付文件 table3_corrected.csv / mr_corrected_N.csv 逐值对照")
    print("=" * 92)
    t3 = pd.read_csv(os.path.join(OUTG, "table3_corrected.csv"))
    mcn = pd.read_csv(os.path.join(OUTG, "mr_corrected_N.csv"))
    name_map = {f[:-4].split("__")[0] + "|" + f[:-4].split("__")[1]: f[:-4]
                for f in files}
    print("\n%-28s %-4s | %12s %12s %10s | %12s %12s"
          % ("pair", "K", "audit_β", "table3_β", "Δ", "audit_p", "table3_p"))
    worst_b, worst_p = 0.0, 0.0
    for r in t3.itertuples():
        key = "%s|%s" % (r.exp, r.out)
        fn = name_map.get(key)
        if fn is None:
            print("  !! 找不到 harmonised 文件: %s" % key)
            continue
        a = df[df.pair == fn].iloc[0]
        db = abs(a.beta - r.beta)
        dp = abs(a.p - r.p)
        worst_b = max(worst_b, db)
        worst_p = max(worst_p, dp)
        flag = "" if (db < 1e-9 and a.K == r.k) else "   <<< 不一致"
        print("%-28s K=%-3d/%-3d | %12.6f %12.6f %10.2e | %12.4g %12.4g%s"
              % (fn, a.K, r.k, a.beta, r.beta, db, a.p, r.p, flag))
    print("\n  最大 Δβ = %.3e   最大 Δp = %.3e" % (worst_b, worst_p))

    # ---- p_re 交叉核对 ----
    print("\n" + "=" * 92)
    print("D. 随机效应 p 值交叉核对（table3 的 p_re vs 独立重算；另对比 mr_results_corrected）")
    print("=" * 92)
    mrc_path = os.path.join(OUTG, "mr_results_corrected.csv")
    mrc = pd.read_csv(mrc_path) if os.path.exists(mrc_path) else None
    print("%-28s %12s %12s %12s" % ("pair", "audit_p_re", "table3_p_re", "mrc_p_re"))
    for r in t3.itertuples():
        key = "%s|%s" % (r.exp, r.out)
        fn = name_map[key]
        a = df[df.pair == fn].iloc[0]
        mrcv = np.nan
        if mrc is not None:
            s = mrc[(mrc.exposure == r.exp) & (mrc.outcome == r.out)]
            if len(s):
                mrcv = float(s.iloc[0]["ivw_re_p"])
        print("%-28s %12.4g %12.4g %12.4g" % (fn, a.p_re, r.p_re, mrcv))

    # ---- Egger 截距显著性（阳性对照关键） ----
    print("\n" + "=" * 92)
    print("E. MR-Egger 截距（方向性多效性）逐对核对")
    print("=" * 92)
    print("%-28s %12s %12s %12s" % ("pair", "audit_int_p", "table3_int_p", "一致?"))
    for r in t3.itertuples():
        fn = name_map["%s|%s" % (r.exp, r.out)]
        a = df[df.pair == fn].iloc[0]
        ok = abs(a.egger_int_p - r.egger_int_p) < 1e-6
        print("%-28s %12.4g %12.4g %12s"
              % (fn, a.egger_int_p, r.egger_int_p, "OK" if ok else "<<< 不一致"))

    print("\n✅ 审计输出：AUDIT_mr_recomputed.csv / AUDIT_scale_inference.csv")


if __name__ == "__main__":
    main()
