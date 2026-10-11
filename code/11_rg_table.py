#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
11_rg_table.py — rg 结果的多重检验校正 + 发表级结果表

校正方法：
  · Bonferroni：α = 0.05 / n_tests（保守）
  · BH-FDR（Benjamini-Hochberg）：控制假发现率（推荐主用）
产出：
  rg_results_long.csv  —— 长表：trait1, trait2, rg, se, z, p, p_bonf, q_fdr, sig
  rg_blocks_summary.csv —— 三块结构汇总（OA×OA / OA×BMD / BMD×BMD）
"""
import os
import sys
import importlib
import numpy as np
import pandas as pd
from scipy.stats import norm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
m09 = importlib.import_module("09_rg_matrix")

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "genetics")


def bh_fdr(p):
    """Benjamini-Hochberg FDR q 值。"""
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty(n)
    out[order] = q
    return out


def stars(p, q):
    if p < 0.001 and q < 0.001: return "***"
    if q < 0.01: return "**"
    if q < 0.05: return "*"
    return ""


def main():
    rg = pd.read_csv(os.path.join(OUT, "rg_matrix.csv"), index_col=0)
    se = pd.read_csv(os.path.join(OUT, "rg_se_matrix.csv"), index_col=0)
    pv = pd.read_csv(os.path.join(OUT, "rg_p_matrix.csv"), index_col=0)

    names = list(rg.index)
    rows = []
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = names[i], names[j]
            r, s, p = rg.loc[a, b], se.loc[a, b], pv.loc[a, b]
            if not (np.isfinite(r) and np.isfinite(p)):
                continue
            rows.append(dict(trait1=a, trait2=b, rg=r, se=s,
                             z=r / s if s > 0 else np.nan, p=p))
    df = pd.DataFrame(rows)
    n = len(df)
    df["p_bonf"] = np.minimum(df.p * n, 1.0)
    df["q_fdr"] = bh_fdr(df.p.values)
    df["sig"] = [stars(p, q) for p, q in zip(df.p, df.q_fdr)]

    # 区块标签
    def block(t):
        return "OA" if t.startswith("OA_") else "BMD"
    df["block"] = [f"{block(a)}×{block(b)}" for a, b in zip(df.trait1, df.trait2)]

    df = df.sort_values(["block", "p"]).reset_index(drop=True)
    df.to_csv(os.path.join(OUT, "rg_results_long.csv"), index=False)

    # ---- 区块汇总
    summ = (df.groupby("block")
              .agg(n_pairs=("rg", "size"),
                   rg_min=("rg", "min"), rg_max=("rg", "max"),
                   rg_mean=("rg", "mean"),
                   p_min=("p", "min"), n_sig_q05=("q_fdr", lambda x: (x < 0.05).sum()))
              .reset_index())
    summ.to_csv(os.path.join(OUT, "rg_blocks_summary.csv"), index=False)

    print(f"n_tests = {n}  →  Bonferroni α = {0.05/n:.2e}")
    print("\n=== 区块汇总 ===")
    print(summ.to_string(index=False))
    print("\n=== OA×BMD 全部 22 对（核心结果） ===")
    oa = df[df.block == "OA×BMD"].copy()
    print(oa[["trait1", "trait2", "rg", "se", "p", "q_fdr", "sig"]].to_string(index=False))
    print(f"\n=== 多重校正后仍显著（q_fdr<0.05）: "
          f"{(df.q_fdr < 0.05).sum()} / {n} 对 ===")


if __name__ == "__main__":
    main()
