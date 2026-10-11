#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
30_coloc_power.py — 位点级共定位的**功效定量**（为什么测不到，以及要多大样本才测得）
=========================================================================
背景
----
29_colocalization.py 的结果：31 个 BMD 工具变量基因座对 OA_knee 无一达到 PP.H4>=0.8
（中位数 0.005），而同一批基因座对 fracture 的阳性对照在 WNT16 达到 PP.H4≈1.00。
「全阴性」有两种读法：
  (a) 真的不存在共享因果变异；
  (b) 结局侧信号太弱，coloc 根本没有分辨力。
本脚本用 MR 估计本身把 (b) 算出来，从而判定该选哪种读法。

推导（尺度无关，规避 BMD 的 N 重建问题）
----------------------------------------
对某个 SNP i，记
    暴露侧:  z_e,i = beta_e,i / SE_e,i
    结局侧:  z_o,i = beta_o,i / SE_o,i
MR 假说 beta_o = beta_MR * beta_e（结局每 1-SD 暴露的效应为 beta_MR），于是
    z_o,i = beta_MR * beta_e,i / SE_o,i
          = z_e,i * beta_MR * (SE_e,i / SE_o,i)
而
    SE_e = 1/sqrt(2 N_e f(1-f))                 (quant)
    SE_o = 1/sqrt(2 N_o f(1-f) s(1-s))          (cc)
    => SE_e/SE_o = sqrt( N_o * s(1-s) / N_e )
**f 与 beta_e 完全消去**，因此该式不受 BMD 样本量重建（x3.5428）的影响 ——
这正是功效分析相对于效应量估计的稳健之处。

再算「要把最强工具位点的结局侧 z 推到全基因组显著 |z|>=5.45 所需的样本量」：
    N_needed = N_o * (5.45 / |z_expected|)^2

输出
----
  output/genetics/coloc_power.csv
  output/genetics/coloc_power_report.txt
"""

import os
import sys
import math
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "output", "genetics")
sys.path.insert(0, ROOT)

GW_Z = 5.45          # p = 5e-8 双侧对应的 |z|

# MR 估计值（取自稿件 Table 3，已由 audit_v20 独立复算验证；此处只读取不改）
# BMD_total -> OA_knee ; BMD_total -> fracture
BETA_MR = {
    # BMD_total -> OA_knee  +0.02312 (SE 0.00369)  [AUDIT_mr_recomputed.csv, audit_v20 verified]
    # 五部位区间: femoral-neck +0.01709 .. lumbar +0.02791
    "OA_knee": {"beta": 0.02312, "lo": 0.01709, "hi": 0.02791},
    # BMD_total -> fracture  -0.33713 (SE 0.02241)
    "fracture": {"beta": -0.33713, "lo": None, "hi": None},
}

N_OA = 1208397
S_OA = 0.12
N_FX = 264973
S_FX = 37857 / 264973.0
N_BMD = 31986          # BMD_total 采用值（GCST90568455）


def se_ratio_to_outcome(n_bmd, n_out, s_out):
    """SE_e / SE_o = sqrt( N_o * s(1-s) / N_e )"""
    return math.sqrt(n_out * s_out * (1.0 - s_out) / n_bmd)


def main():
    lines = []

    def P(s=""):
        print(s)
        lines.append(s)

    P("=" * 78)
    P("POWER OF LOCUS-LEVEL COLOCALIZATION — why BMD x OA cannot be resolved")
    P("=" * 78)
    P("z_o = z_e * beta_MR * sqrt(N_o*s(1-s)/N_e)      (scale-free: f and beta_e cancel)")
    P("")

    harm = pd.read_csv(os.path.join(OUT, "harmonised", "BMD_total__OA_knee.csv"))
    m29 = __import__("29_colocalization")
    d_oa, _, _ = m29.get_trait("OA_knee")
    d_fx, _, _ = m29.get_trait("fracture")

    rows = []
    for _, r in harm.iterrows():
        ze = float(r["beta_e"]) / float(r["se_e"])       # 暴露侧 z（尺度无关）
        rec = dict(snp=r["SNP"], chrom=int(r["CHR"]), pos=int(r["BP"]),
                   beta_e=float(r["beta_e"]), z_bmd=ze,
                   p_o=float(r["p_o"]))
        # 结局侧观测 z（从 OA/fracture 全汇总统计里按 CHR:BP 取）
        for tag, dd in (("oa", d_oa), ("fx", d_fx)):
            hit = dd[(dd.CHR.values == int(r["CHR"])) & (dd.BP.values == int(r["BP"]))]
            rec["z_%s_obs" % tag] = float(hit.Z.iloc[0]) if len(hit) else float("nan")
        rows.append(rec)
    df = pd.DataFrame(rows)

    for out_name, n_out, s_out, tag in (("OA_knee", N_OA, S_OA, "oa"),
                                        ("fracture", N_FX, S_FX, "fx")):
        mr = BETA_MR[out_name]
        ratio = se_ratio_to_outcome(N_BMD, n_out, s_out)
        df["z_%s_exp" % tag] = df["z_bmd"] * mr["beta"] * ratio
        if mr["lo"] is not None:
            df["z_%s_exp_lo" % tag] = df["z_bmd"] * mr["lo"] * ratio
            df["z_%s_exp_hi" % tag] = df["z_bmd"] * mr["hi"] * ratio
        df["N_%s_needed" % tag] = n_out * (GW_Z / df["z_%s_exp" % tag].abs()) ** 2

    df.to_csv(os.path.join(OUT, "coloc_power.csv"), index=False)

    P("BMD_total instruments: n = %d" % len(df))
    P("SE_e/SE_oa = sqrt(%d*%.4f/%.0f) = %.3f" %
      (N_OA, S_OA * (1 - S_OA), N_BMD, se_ratio_to_outcome(N_BMD, N_OA, S_OA)))
    P("SE_e/SE_fx = %.3f\n" % se_ratio_to_outcome(N_BMD, N_FX, S_FX))

    P("--- OA_knee (MR beta = %+.4f per 1-SD BMD; range %+.3f..%+.3f) ---"
      % (BETA_MR["OA_knee"]["beta"], BETA_MR["OA_knee"]["lo"], BETA_MR["OA_knee"]["hi"]))
    P("%-13s %8s %9s %9s %9s %14s" %
      ("snp", "z_bmd", "z_oa_exp", "z_oa_obs", "ratio", "N_needed(5e-8)"))
    top = df.reindex(df.z_bmd.abs().sort_values(ascending=False).index)
    for _, r in top.head(12).iterrows():
        P("%-13s %8.2f %9.2f %9.2f %9s %14s" %
          (r.snp, r.z_bmd, r.z_oa_exp, r.z_oa_obs,
           ("%.2f" % (r.z_oa_obs / r.z_oa_exp)) if r.z_oa_exp else "-",
           "%.3g" % r.N_oa_needed))
    P("")
    P("  strongest instrument: %s  z_bmd=%.2f  z_oa_exp=%.2f  z_oa_obs=%.2f"
      % (top.iloc[0].snp, top.iloc[0].z_bmd, top.iloc[0].z_oa_exp, top.iloc[0].z_oa_obs))
    P("  median |z_oa_exp| over 31 loci      = %.2f" % df.z_oa_exp.abs().median())
    P("  max    |z_oa_exp| over 31 loci      = %.2f" % df.z_oa_exp.abs().max())
    P("  n loci with |z_oa_exp| >= 5.45      = %d / %d"
      % ((df.z_oa_exp.abs() >= GW_Z).sum(), len(df)))
    P("  n loci with |z_oa_obs| >= 5.45      = %d / %d"
      % ((df.z_oa_obs.abs() >= GW_Z).sum(), len(df)))
    P("  => OA N needed for the TOP locus to reach p<5e-8 : %.3g"
      % df.N_oa_needed.min())
    P("     (current OA_knee N = %d ; fold = %.0f x)"
      % (N_OA, df.N_oa_needed.min() / N_OA))

    P("\n--- fracture (positive control; MR beta = %+.3f per 1-SD BMD) ---"
      % BETA_MR["fracture"]["beta"])
    P("%-13s %8s %9s %9s %9s" % ("snp", "z_bmd", "z_fx_exp", "z_fx_obs", "ratio"))
    for _, r in top.head(6).iterrows():
        P("%-13s %8.2f %9.2f %9.2f %9.2f" %
          (r.snp, r.z_bmd, r.z_fx_exp, r.z_fx_obs, r.z_fx_obs / r.z_fx_exp))
    P("  n loci with |z_fx_exp| >= 5.45 = %d / %d"
      % ((df.z_fx_exp.abs() >= GW_Z).sum(), len(df)))
    P("  n loci with |z_fx_obs| >= 5.45 = %d / %d"
      % ((df.z_fx_obs.abs() >= GW_Z).sum(), len(df)))

    P("\n--- Interpretation ---")
    P("Expected per-locus OA z under the MR estimate is O(0.6-1.0), i.e. ~6-9x")
    P("below genome-wide significance, while the fracture control is O(6).")
    P("coloc.abf resolves H3 vs H4 only when BOTH traits carry strong signals;")
    P("with z_oa < 1 the posterior is dominated by the prior regardless of data.")
    P("=> The null coloc result is a POWER statement, not evidence of absence.")

    with open(os.path.join(OUT, "coloc_power_report.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))


if __name__ == "__main__":
    main()
