#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
34_mvmr_subgroup.py — 不依赖 MVMR 假设的独立佐证：按「工具是否也是身高位点」分组
============================================================================
MVMR 的结论依赖「仪器在控制另一暴露后仍有效」等假设。为给出一个不依赖这些假设的
独立检验，本脚本把每个 BMD 部位的工具变量按**该 SNP 是否同时是身高 GWS 位点**分组，
分别在组内做单变量 IVW，并比较两组估计。

逻辑（可证伪的预测）
--------------------
若 BMD→OA 的效应主要由骨尺寸（骨架大）驱动，则「同时是身高位点」的那部分工具
应该给出**更大**的 BMD→OA 估计；若效应来自骨量本身，则两组估计应无系统差异。

方法
----
· 分组：height-GWS 组（p_h < 5×10⁻⁸）vs height-清洁组（p_h > 1×10⁻³）；
  （5e-8 ≤ p_h ≤ 1e-3 的中间组单独列出但不进入比较，避免误分）
· 组内 IVW（与手稿 13_ 一致的 ratio 形式）；
· 组间差异检验：
     z = (β_A − β_B)/√(SE_A² + SE_B²)          （两独立估计的 Wald 检验）
     Q_diff = Σ_g (β_g − β_pooled)²/SE_g², df=1 （Cochran 异质性）
· 同时报告每组内的平均 F 与工具数。
============================================================================
"""
import os, math, csv
import numpy as np
import pandas as pd
from scipy.stats import chi2

import importlib
m33 = importlib.import_module("33_mvmr_height_v2")

OUT = m33.OUT
load_munge, load_height = m33.load_munge, m33.load_height
ivw_uni = m33.ivw_uni
p_from_z = m33.p_from_z


def run(hfn, hlab):
    dh = load_height(hfn)
    rows = []
    for (e, o, me, mo) in m33.PAIRS:
        if o != "OA_knee":
            continue
        pc = os.path.join(OUT, "harmonised", m33.HN[(e, o)])
        d, meta = m33.build(pc, me, mo, dh, 0)      # n_h=0 → 仅 BMD 工具
        d = d[d.source == "BMD"].copy()
        g = d.p_h < 5e-8
        clean = d.p_h > 1e-3
        mid = (~g) & (~clean)
        res = {"exposure": e, "outcome": o, "height": hlab, "n_total": len(d)}
        for name, msk in (("height_GWS", g.values), ("height_clean", clean.values),
                          ("intermediate", mid.values)):
            sub = d[msk]
            if len(sub) >= 3:
                b, s = ivw_uni(sub.beta_e.values, sub.beta_o.values, sub.se_o.values)
                F = float(np.mean((sub.beta_e.values / sub.se_e.values) ** 2))
                res["n_" + name] = len(sub)
                res["beta_" + name] = b
                res["se_" + name] = s
                res["p_" + name] = p_from_z(b / s)
                res["F_" + name] = F
                res["mean_abszh_" + name] = float(np.mean(np.abs(sub.beta_h.values / sub.se_h.values)))
            else:
                for k in ("n", "beta", "se", "p", "F", "mean_abszh"):
                    res[k + "_" + name if k != "n" else "n_" + name] = np.nan
                res["n_" + name] = len(sub)
        if res.get("n_height_GWS", 0) >= 3 and res.get("n_height_clean", 0) >= 3:
            bA, sA = res["beta_height_GWS"], res["se_height_GWS"]
            bB, sB = res["beta_height_clean"], res["se_height_clean"]
            res["diff_z"] = (bA - bB) / math.sqrt(sA ** 2 + sB ** 2)
            res["diff_p"] = p_from_z(res["diff_z"])
            pooled = (bA / sA ** 2 + bB / sB ** 2) / (1 / sA ** 2 + 1 / sB ** 2)
            res["Q_diff"] = (bA - pooled) ** 2 / sA ** 2 + (bB - pooled) ** 2 / sB ** 2
            res["Q_diff_p"] = float(chi2.sf(res["Q_diff"], 1))
            res["ratio_GWS_over_clean"] = bA / bB if abs(bB) > 1e-12 else np.nan
        rows.append(res)
        print("  %-15s n=%2d | 身高GWS组 n=%2d β=%+.4f(p=%.2e, 平均|z_h|=%.1f) | "
              "身高清洁组 n=%2d β=%+.4f(p=%.2e) | 差 z=%+.2f (p=%.3f)"
              % (e, len(d), res.get("n_height_GWS", 0), res.get("beta_height_GWS", float('nan')),
                 res.get("p_height_GWS", float('nan')), res.get("mean_abszh_height_GWS", float('nan')),
                 res.get("n_height_clean", 0), res.get("beta_height_clean", float('nan')),
                 res.get("p_height_clean", float('nan')),
                 res.get("diff_z", float('nan')), res.get("diff_p", float('nan'))))
    return pd.DataFrame(rows)


def main():
    dfs = [run(m33.HEIGHT_EUR, "EUR"), run(m33.HEIGHT_NOUKB, "EUR_noUKB")]
    pd.concat(dfs).to_csv(os.path.join(OUT, "mvmr_height_subgroup.csv"),
                          index=False, encoding="utf-8")
    print("\nwritten mvmr_height_subgroup.csv")


if __name__ == "__main__":
    main()
