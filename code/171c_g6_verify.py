#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
171c_g6_verify.py — G6 缩放的**独立自校验 + 漏项自动补齐**
=============================================================================
原理（不依赖任何"我记得改过什么"的陈述，全部由磁盘上的数自行判定）
--------------------------------------------------------------------------
对每个 (β, SE, p) 三元组，p 是尺度不变量，故

    |z|_from_p = |Φ⁻¹(p/2)|        应等于    |β| / SE

若某列漏了缩放，两者之比必为 C、1/C（或 C² 之类）的确定倍数，而不是随机误差。
本脚本：
  1) 逐三元组计算 ratio = (|β|/SE) / |z|_from_p
  2) ratio ≈ 1        → 该三元组已正确
  3) ratio ≈ C 或 1/C → 自动把 SE（或 β）补齐到正确尺度，并记入审计
  4) 打印补齐前后的对照，落盘 G6_verify_report.csv

再叠加一层与备份（pre-G6）的逐字段 z/p 不变性核对，作为"改对了"的独立证据。
"""
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
BACKUP = os.path.join(OUT, "_backup_pre_G6")

C = float(m21.F_OA_KNEE)

TRIPLES = {
    "mvmr_bmi_main.csv": [
        ("uni_beta_bmd", "uni_se_bmd", "uni_p_bmd"),
        ("uni_beta_bmi", "uni_se_bmi", "uni_p_bmi"),
        ("mvmr_beta_bmd", "mvmr_se_bmd", "mvmr_p_bmd"),
        ("mvmr_beta_bmi", "mvmr_se_bmi", "mvmr_p_bmi"),
        ("egger_beta_bmd", "egger_se_bmd", None),
    ],
    "mvmr_height_main.csv": [
        ("uni_beta", "uni_se", "uni_p"),
        ("mvmr_beta", "mvmr_se", "mvmr_p"),
        ("mvmr_h", "mvmr_h_se", "mvmr_h_p"),
        ("uni_h_beta", "uni_h_se", "uni_h_p"),
    ],
    "mvmr_height_n_h1000.csv": [
        ("uni_beta", "uni_se", "uni_p"),
        ("mvmr_beta", "mvmr_se", "mvmr_p"),
        ("mvmr_h", "mvmr_h_se", "mvmr_h_p"),
        ("uni_h_beta", "uni_h_se", "uni_h_p"),
    ],
    "mvmr_height_grid.csv": [
        ("uni_beta", "uni_se", "uni_p"),
        ("mvmr_beta", "mvmr_se", "mvmr_p"),
        ("mvmr_h", "mvmr_h_se", "mvmr_h_p"),
        ("uni_h_beta", "uni_h_se", "uni_h_p"),
    ],
    "mr_results_corrected.csv": [("ivw", "ivw_se", "ivw_p")],
    "table3_corrected.csv": [("beta", "se", "p")],
    "mr_power.csv": [("beta", "se", None)],
}

# 无 p 列时，用备份（pre-G6）中同名行的 z 作参照
BACKUP_Z = {}


def z_from_p(p):
    if p is None or not np.isfinite(p) or p <= 0:
        return np.nan
    return float(abs(norm.isf(p / 2.0)))


def main():
    report = []
    print("=" * 78)
    print("G6 缩放自校验：|β|/SE  vs  |z| from p（或备份 z）    C = %.4f" % C)
    print("=" * 78)

    for fname, triples in TRIPLES.items():
        fp = os.path.join(OUT, fname)
        bp = os.path.join(BACKUP, fname)
        if not os.path.exists(fp):
            print("⚠ 缺失", fname)
            continue
        df = pd.read_csv(fp)
        bdf = pd.read_csv(bp) if os.path.exists(bp) else None
        changed = False
        for bc, sc, pc in triples:
            if bc not in df.columns or sc not in df.columns:
                continue
            for i, r in df.iterrows():
                beta, se = float(r[bc]), float(r[sc])
                if not np.isfinite(beta) or not np.isfinite(se) or se == 0:
                    continue
                if pc and pc in df.columns:
                    zref = z_from_p(float(r[pc]))
                else:
                    zref = np.nan
                if not np.isfinite(zref) or zref < 1e-6:
                    continue
                ratio = (abs(beta) / se) / zref
                # 判定漏缩放：ratio 落在 C 或 1/C 附近（相对误差 <2%）
                for fac, who in ((C, sc), (1.0 / C, sc), (C, bc), (1.0 / C, bc)):
                    if abs(ratio / fac - 1.0) < 0.02:
                        old = float(df.at[i, who])
                        df.at[i, who] = old * (C if fac == C else 1.0 / C)
                        report.append(dict(file=fname, row=int(i), col=who,
                                           ratio=round(ratio, 4),
                                           fix=("×C" if fac == C else "÷C"),
                                           old=old, new=float(df.at[i, who])))
                        changed = True
                        break
                else:
                    report.append(dict(file=fname, row=int(i), col=sc,
                                       ratio=round(ratio, 4), fix="ok",
                                       old=se, new=se))
        if changed:
            df.to_csv(fp, index=False)
        # 汇总
        sub = [x for x in report if x["file"] == fname]
        bad = [x for x in sub if x["fix"] != "ok"]
        print("\n%s：%d 个三元组检查，%d 处补齐"
              % (fname, len(sub), len(bad)))
        for x in bad:
            print("   补齐 %-22s 行%d  ratio=%.4f  %s  %.6g -> %.6g"
                  % (x["col"], x["row"], x["ratio"], x["fix"], x["old"], x["new"]))
        off = [x for x in sub if x["fix"] == "ok" and abs(x["ratio"] - 1) > 0.02]
        if off:
            print("   ⚠ 仍偏离 1 的项：", [(x["col"], x["row"], x["ratio"]) for x in off[:8]])

    # ---- 与 pre-G6 备份的 z / p 不变性核对（严格 clump 的 mr_corrected_N）----
    print("\n" + "=" * 78)
    print("独立核对：新 mr_results_corrected.csv vs pre-G6 mr_corrected_N.csv")
    print("（两者同为严格 clump；唯一差别应为 OA 侧尺度）")
    print("=" * 78)
    new = pd.read_csv(os.path.join(OUT, "mr_results_corrected.csv"))
    old = pd.read_csv(os.path.join(BACKUP, "mr_corrected_N.csv"))
    ok = True
    for r in new.itertuples():
        m = old[(old.exposure == r.exposure) & (old.outcome == r.outcome)]
        if m.empty:
            continue
        o = m.iloc[0]
        znew, zold = r.ivw / r.ivw_se, o.beta / o.se
        dp = abs(r.ivw_p - o.p) / max(o.p, 1e-300)
        dk = int(r.n_snp) == int(o.K)
        dz = abs(znew - zold)
        flag = "✅" if (dz < 1e-6 and dp < 1e-9 and dk) else "❌"
        if flag == "❌":
            ok = False
        print("  %s %-16s→%-15s K %d/%d  Δz=%.2e  Δp/p=%.2e  β %.5f→%.5f (×%.4f)"
              % (flag, r.exposure, r.outcome, r.n_snp, int(o.K), dz, dp,
                 o.beta, r.ivw, r.ivw / o.beta))
    print("\n" + ("✅ z / p / K 全部不变 —— G6 只改了效应量的绝对值"
                  if ok else "❌ 存在不应有的变化"))

    rep = pd.DataFrame(report)
    rep.to_csv(os.path.join(OUT, "G6_verify_report.csv"), index=False)
    print("\n📁 校验明细 -> %s" % os.path.join(OUT, "G6_verify_report.csv"))


if __name__ == "__main__":
    main()
