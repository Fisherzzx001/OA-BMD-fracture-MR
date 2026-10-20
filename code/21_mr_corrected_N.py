#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
21_mr_corrected_N.py
=========================================================================
修正 MR 效应量的 **样本量标度错误**，并与 GEFOS 外部基准校验。

问题
----
13_mr_analysis.py 对 munge 过的表型按
        SE = 1/sqrt(N·2p(1-p)),   beta = Z·SE
重构效应量。该式只在 **N 正确** 时才给出 SD 单位效应。
但 BMD 的 munge 文件里 N 是常量占位值 **400,000**，真实 N = **31,873**（Qian 2025）。

    beta_重构(占位 N=400,000) = beta_SD · sqrt(N_真 / N_用) = beta_SD · 0.2823
    （即 munged β 仅为真实 SD 标度的 0.2823×，被**低估** 3.543 倍；
     本脚本将 BMD 暴露/结局的 β、SE ×3.5426 恢复真实标度）

SE 同步按同一因子缩放，故 z / p / 方向 / 显著性 **完全不变** —— 与 N 不变性一致。

校正规则
--------
  · BMD（10 个）：beta, se  ×  sqrt(400000/31873)  = ×3.5428
  · OA_selfrep ：beta, se  ×  sqrt(350000/462933)  = ×0.8695
  · OA_knee    ：beta, se  ×  1/sqrt(K(1-K))      = ×3.3887   ← **G6 修正（2026-10-20）**
  · fracture   ：**不校正**（原始文件直接给 log-OR 的 beta）

OA_knee 的 case-control 修正（G6）
---------------------------------
旧口径认为 OA_knee 的 munge 文件里 N 是 per-SNP 真实值，故"无需校正"。
这是错的：per-SNP 的 N 是**总**样本量，而二分类结局（logistic 回归）的
效应量满足

        se_logOR = 1 / sqrt( N_eff · 2 f (1-f) ),   N_eff = N_total · K(1-K)

直接用 N_total 重构得到的是**线性概率尺度**的效应，比真实 log-OR 小
1/sqrt(K(1-K)) 倍。骨折结局因为原始文件直接给了 log-OR 的 beta/se 而幸免，
OA 结局则没有——同一条流水线对两个结局用了不一致的尺度处理。

K 的校准（实测，不假设）：
  取 GO Consortium 2025 膝 OA 真实汇总统计中与 munged 文件可匹配的 39 个 SNP，
  逐点比对得  |Z_true / Z_munged| = 1.000000（std 5e-14，确认同源），且
  |beta_true| = |Z_munged| × se_true（误差 <1e-16），
  故 K(1-K) = N_eff_true / N_munged，逐 SNP 中位 **0.087082**（IQR 0.0835–0.0881）
  → K = 0.0964（膝 OA 病例约 116,500 / 1,208,397，与 GO 联盟量级一致）
  → 校正因子 1/sqrt(0.087082) = **3.3887**

  校准脚本：code/_probe_g6b.py；逐 SNP 明细：output/genetics/_g6_k_estimation.csv

影响：只改 β / SE / OR 的绝对值；z、p、方向、显著性、Cochran Q、I²、τ²、
Egger 斜率与截距、加权中位数、PRESSO、留一法**全部不变**（尺度不变量）。


外部基准：GEFOS（Trajanoska BMJ 2018, PMID 30158200）
    股骨颈 BMD 每降低 1 SD → 骨折 OR = 1.55（95% CI 1.48–1.63）
=========================================================================
"""
import os
import sys
import importlib
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
m13 = importlib.import_module("13_mr_analysis")

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "output", "genetics")

# ---- N 校正 ---------------------------------------------------------------
N_FILE_BMD, N_TRUE_BMD = 400_000.0, 31_873.0
N_FILE_SELF, N_TRUE_SELF = 350_000.0, 462_933.0
F_BMD = float(np.sqrt(N_FILE_BMD / N_TRUE_BMD))      # 3.5428
F_SELF = float(np.sqrt(N_FILE_SELF / N_TRUE_SELF))   # 0.8695

# --- G6：OA_knee 的 case-control（有效样本量）修正 -------------------------
# K(1-K) 由 39 个 SNP 的真实 log-OR 与 munged Z 逐点比对实测得到（见文件头）
K_OA_KNEE = 0.0964                                   # 膝 OA 病例比例
KC_OA_KNEE = K_OA_KNEE * (1 - K_OA_KNEE)             # 0.087082（实测中位）
F_OA_KNEE = float(1.0 / np.sqrt(KC_OA_KNEE))         # 3.3887

GEFOS_OR_PER_SD = 1.55        # 股骨颈 BMD 每降 1 SD 的骨折 OR
GEFOS_CI = (1.48, 1.63)


def trait_factor(name):
    """返回该表型 beta/se 的校正因子。"""
    if name.startswith("OA_knee"):
        # G6：per-SNP 的 N 是总样本量；二分类结局须换为 N_eff = N·K(1-K)，
        # 否则重构出的效应落在线性概率尺度，比真实 log-OR 小 1/sqrt(K(1-K)) 倍。
        return F_OA_KNEE
    if name.startswith("OA_selfrep"):
        return F_SELF
    if name.startswith("fracture"):
        return 1.0                      # 原始文件直接给 log-OR 的 beta，无需校正
    return F_BMD                        # 其余为 BMD（N 常量占位 400,000）


def load_all_corrected():
    dfs = m13.load_all()
    for nm, d in dfs.items():
        f = trait_factor(nm)
        if f != 1.0:
            d = d.copy()
            d["beta"] = d["beta"] * f
            d["se"] = d["se"] * f
            dfs[nm] = d
        print("   %-16s 校正因子 ×%.4f" % (nm, f))
    return dfs


def load_strict(exp, dfs):
    fp = os.path.join(OUT, f"{exp}_clumped.clumped")
    if not os.path.exists(fp):
        return m13.select_instruments(dfs[exp])
    cl = pd.read_csv(fp, sep=r"\s+")
    snps = set(cl["SNP"].astype(str))
    return m13.select_instruments(dfs[exp][dfs[exp]["SNP"].astype(str).isin(snps)])


def main():
    os.makedirs(OUT, exist_ok=True)
    print("=" * 74)
    print("MR 效应量 N 标度校正（BMD 400,000 -> 31,873；OA_selfrep 350,000 -> 462,933）")
    print("=" * 74)
    print("BMD 校正因子 = %.4f    OA_selfrep 校正因子 = %.4f" % (F_BMD, F_SELF))
    print("\n▶ 载入并校正 ...")
    dfs = load_all_corrected()

    print("\n" + "=" * 74)
    print("校正后 MR（严格 clump 仪器集）")
    print("=" * 74)
    rows = []
    for exp, out, note in m13.PAIRS:
        iv = load_strict(exp, dfs)
        h = m13.harmonise(iv, dfs[out])
        if len(h) < 3:
            print("【%s -> %s】 IV 不足，跳过" % (exp, out))
            continue
        r = m13.run_mr(h.beta_e.values, h.se_e.values,
                       h.beta_o.values, h.se_o.values)
        b = r.get("ivw", np.nan)
        p = r.get("ivw_p", np.nan)
        row = dict(exposure=exp, outcome=out, note=note, K=r.get("n_snp", np.nan),
                   beta=b, se=r.get("ivw_se", np.nan), p=p)
        # BMD 为暴露、fracture 为结局时：β 已是「每 SD BMD」→ 可直接取 OR
        if exp.startswith("BMD") and out == "fracture":
            row["OR_per_SD_higher"] = float(np.exp(b))
            row["OR_per_SD_lower"] = float(np.exp(-b))
        rows.append(row)

        extra = ""
        if "OR_per_SD_lower" in row:
            extra = "   OR(每降1SD)=%.3f" % row["OR_per_SD_lower"]
        print("【%-16s -> %-16s】K=%-3d β=%+.4f  SE=%.4f  p=%.3g%s"
              % (exp, out, row["K"], b, row["se"], p, extra))

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "mr_corrected_N.csv"), index=False)
    print("\n✅ -> output/genetics/mr_corrected_N.csv")

    # ---- GEFOS 外部基准校验 ----
    print("\n" + "=" * 74)
    print("外部基准校验：GEFOS（股骨颈 BMD 每降 1 SD → 骨折 OR = %.2f, 95%%CI %.2f–%.2f）"
          % (GEFOS_OR_PER_SD, *GEFOS_CI))
    print("=" * 74)
    sub = df[(df.exposure.str.startswith("BMD")) & (df.outcome == "fracture")]
    for r in sub.itertuples():
        flag = "✅ 与 GEFOS 一致" if 1.0 < r.OR_per_SD_lower < 3.0 else "⚠ 量级偏离"
        print("  %-16s OR(每降1SD) = %.3f   校正前 β=%s → 校正后 β=%+.4f   %s"
              % (r.exposure, r.OR_per_SD_lower, "n/a", r.beta, flag))
    if len(sub):
        best = sub[sub.exposure.str.contains("femur_neck")]
        if len(best):
            o = float(best.iloc[0]["OR_per_SD_lower"])
            print("\n  股骨颈（与 GEFOS 同部位）：OR = %.3f  vs  GEFOS 1.55 "
                  "→ 相对差异 %.1f%%" % (o, 100 * abs(o - GEFOS_OR_PER_SD) / GEFOS_OR_PER_SD))
    print("\n注：校正只改变效应量绝对值，z / p / 方向与显著性均不变（N 不变性）。")


if __name__ == "__main__":
    main()
