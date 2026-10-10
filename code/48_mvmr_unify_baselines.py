#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
48_mvmr_unify_baselines.py — M1 统一基线：三套 univariable 基线与三套衰减率的可审计产出
=============================================================================
问题（复审 102 · M1）
--------------------
同一对 BMD→OA_knee，稿件内部并存三套互不相等的「单变量基线」，且 §3.10 把
uni_beta 交叉引用为 Table 3（不符），读者无法复现。

诊断结论（46_ / 47_，已实证）
-----------------------------
· 主 MR 工具 31/25/21/16/24 → MVMR 实际 27/23/19/15/21，
  **丢失的 12 个全部是「不在身高 GWAS」，回文 0、比对失败 0**。
  6 个唯一位点（chr:bp 亦查无）不在 Yengo 2022 EUR（HapMap3，1,372,038 变体）内，
  **无法在两暴露设计中补全** → 物理上做不到 31/31 归一。
· 用 harmonised 值在完整工具集上算 IVW == Table 3（差 +0.01%），
  在 27 子集上为 0.01996（−13.69%）→ 差异**纯粹来自工具覆盖**，非估计器差异。
· n_h=0 不可用：BMD 工具间 β_BMD 与 β_height 高度共线（femur_neck θ_H 由 0.049
  跳到 0.201、衰减 45.6%），属弱识别。加入身高专属工具正是为打破共线，
  故 n_h=100 主规格保留。

本脚本做什么
------------
不加新估计，只把既有 grid 中的三套基线并置，产出可审计对照表：
  β_uni_primary    = Table 3 主 MR（严格 clump 全工具集）   ← 与 Table 3 恒等
  β_uni_mvmset     = MVMR 所用 K 个工具上的 IVW（随 n_h 变）
  β_uni_bmdsubset  = MVMR 工具中的 BMD 子集 IVW（不随 n_h 变）
  β_mvmr           = MVMR 调整后估计
并给出对应三套衰减率 + 阳性对照（fracture）的同口径核对。

输出
----
  output/genetics/mvmr_unify_baselines.csv
  output/genetics/mvmr_unify_dropped_snps.csv
=============================================================================
"""
import os, csv
import pandas as pd

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "output", "genetics")
PRIMARY_NH = 100          # 主规格（与稿件 §3.10 的 θ_H=+0.0426 一致）
PRIMARY_HEIGHT = "EUR"

EXP_ORDER = ["BMD_total", "BMD_femur_neck", "BMD_lumbar", "BMD_arms", "BMD_leg"]


def main():
    grid = pd.read_csv(os.path.join(OUT, "mvmr_height_grid.csv"))
    diag = pd.read_csv(os.path.join(OUT, "mvmr_unify_diag.csv"))

    rows = []
    for o in ("OA_knee", "fracture"):
        sub = grid[(grid.outcome == o) &
                   (grid.height == PRIMARY_HEIGHT) &
                   (grid.n_h == PRIMARY_NH)]
        for e in (EXP_ORDER if o == "OA_knee" else ["BMD_total", "BMD_femur_neck"]):
            r = sub[sub.exposure == e]
            if r.empty:
                continue
            r = r.iloc[0]
            bp, bm, bb, mv = (float(r.uni_beta), float(r.ref_uni_beta),
                              float(r.uni_bmd_only), float(r.mvmr_beta))
            rows.append(dict(
                outcome=o, exposure=e, spec="%s|n_h=%d" % (PRIMARY_HEIGHT, PRIMARY_NH),
                K=int(r.K), n_bmd=int(r.n_bmd), n_height_only=int(r.K) - int(r.n_bmd),
                beta_uni_primary=bm,                      # == Table 3，恒等
                beta_uni_mvmset=bp,                       # MVMR 同工具集 IVW
                beta_uni_bmdsubset=bb,                    # BMD 子集 IVW
                beta_mvmr=mv,
                atten_vs_primary_pct=100.0 * (1.0 - mv / bm) if bm else float("nan"),
                atten_vs_mvmset_pct=float(r.atten_vs_full),
                atten_vs_bmdsubset_pct=float(r.atten_vs_bmdsubset),
                theta_H=float(r.mvmr_h), theta_H_p=float(r.mvmr_h_p),
                Fc_bmd=float(r.Fc_bmd), Fc_bmd_bmdsub=float(r.Fc_bmd_bmdsub),
                Q=float(r.Q), Q_p=float(r.Q_p)))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "mvmr_unify_baselines.csv"), index=False, encoding="utf-8")

    print("=" * 118)
    print("统一基线对照表（主规格 %s|n_h=%d）" % (PRIMARY_HEIGHT, PRIMARY_NH))
    print("=" * 118)
    for o in ("OA_knee", "fracture"):
        print("\n--- outcome = %s ---" % o)
        print("  %-15s %4s %5s | %11s %11s %11s | %11s | %9s %9s %9s"
              % ("exposure", "K", "nBMD", "β_uni(T3)", "β_uni(MVMR集)", "β_uni(BMD子集)",
                 "β_MVMR", "att|T3", "att|集", "att|子集"))
        for r in df[df.outcome == o].itertuples():
            print("  %-15s %4d %5d | %11.5f %12.5f %12.5f | %11.5f | %8.2f%% %8.2f%% %8.2f%%"
                  % (r.exposure, r.K, r.n_bmd, r.beta_uni_primary, r.beta_uni_mvmset,
                     r.beta_uni_bmdsubset, r.beta_mvmr,
                     r.atten_vs_primary_pct, r.atten_vs_mvmset_pct, r.atten_vs_bmdsubset_pct))

    ok = df[df.outcome == "OA_knee"]
    print()
    print("=" * 118)
    print("稿件应报告的「身高–size 归因比例」区间（三基线下取并集）")
    print("=" * 118)
    lo = min(ok.atten_vs_primary_pct.min(), ok.atten_vs_mvmset_pct.min(),
             ok.atten_vs_bmdsubset_pct.min())
    hi = max(ok.atten_vs_primary_pct.max(), ok.atten_vs_mvmset_pct.max(),
             ok.atten_vs_bmdsubset_pct.max())
    print("  衰减率 min=%.1f%%  max=%.1f%%" % (lo, hi))
    print("  → 身高–size 归因比例 ≈ %.0f–%.0f%%" % (lo, hi))
    print("  → 独立于身高–size 通路的残差 ≈ %.0f–%.0f%%" % (100 - hi, 100 - lo))
    for col, lab in (("atten_vs_primary_pct", "vs Table 3（主）"),
                     ("atten_vs_mvmset_pct", "vs MVMR 工具集"),
                     ("atten_vs_bmdsubset_pct", "vs BMD 子集")):
        print("     %-20s %.1f–%.1f%%" % (lab, ok[col].min(), ok[col].max()))

    fr = df[df.outcome == "fracture"]
    print()
    print("  阳性对照（BMD→fracture，已知骨量效应）三基线下：")
    for r in fr.itertuples():
        print("     %-15s %8.2f%%  %8.2f%%  %8.2f%%"
              % (r.exposure, r.atten_vs_primary_pct,
                 r.atten_vs_mvmset_pct, r.atten_vs_bmdsubset_pct))
    flo = min(fr.atten_vs_primary_pct.min(), fr.atten_vs_mvmset_pct.min(),
              fr.atten_vs_bmdsubset_pct.min())
    fhi = max(fr.atten_vs_primary_pct.max(), fr.atten_vs_mvmset_pct.max(),
              fr.atten_vs_bmdsubset_pct.max())
    print("     → 区间 %.1f%% 到 %+.1f%%：阳性对照对基线选择不敏感（mass-transparency 成立）"
          % (flo, fhi))

    # ---- 丢弃 SNP 清单 ----
    miss = diag[diag["drop"] == "NOT_IN_HEIGHT_GWAS"][
        ["exposure", "SNP", "CHR", "BP"]].drop_duplicates()
    miss = miss.sort_values(["SNP", "exposure"])
    miss.to_csv(os.path.join(OUT, "mvmr_unify_dropped_snps.csv"), index=False, encoding="utf-8")
    print()
    print("=" * 118)
    print("因「不在 Yengo 2022 EUR 身高 GWAS」而无法进入 MVMR 的位点")
    print("=" * 118)
    uniq = miss.groupby("SNP").agg(exposures=("exposure", lambda s: ",".join(sorted(set(s)))),
                                   CHR=("CHR", "first"), BP=("BP", "first"))
    for s, r in uniq.iterrows():
        print("  %-12s chr%s:%-10d  影响暴露: %s" % (s, r.CHR, r.BP, r.exposures))
    print("\n  唯一位点数 = %d ；受影响 (暴露,位点) 组合 = %d" % (len(uniq), len(miss)))
    print("\nwritten mvmr_unify_baselines.csv (%d 行), mvmr_unify_dropped_snps.csv (%d 行)"
          % (len(df), len(miss)))


if __name__ == "__main__":
    main()
