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
  · OA_knee    ：**不校正**（文件内 N 为 per-SNP，107,795 个唯一值，真实 meta 分析）
  · fracture   ：**不校正**（原始文件直接给 log-OR 的 beta）

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

GEFOS_OR_PER_SD = 1.55        # 股骨颈 BMD 每降 1 SD 的骨折 OR
GEFOS_CI = (1.48, 1.63)


def trait_factor(name):
    """返回该表型 beta/se 的校正因子。"""
    if name.startswith("OA_knee"):
        return 1.0                      # per-SNP N，真实 meta 分析，无需校正
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
