#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
169_rg_long13_rebuild.py — 重建 13 表型 rg 长表（rg_results_long_13.csv）

背景（2026-10-20，docs/123 §8 G8 审计）：
  rg_results_long_13.csv 原由 9/16 的临时脚本生成，其 12 条 fracture 行取自
  code/10 旧产物；而 code/10 存在截距解包 bug——compute_h2 返回
  (hsq, h2_SE, intercept, intercept_SE)，旧代码 `h2, i2, _, _ = compute_h2(...)`
  把第 2 个返回值（h² 的 SE，如 0.00095）误当作截距传入 gencov 权重
  （i_h1/i_h2），截距项缩小约 1000 倍，fracture 行 rg/SE/p 全部被污染。
  code/09 于 10/20 重跑（PHENOS 补 fracture）得到正确的 13×13 矩阵；
  code/10 修复后重跑与 code/09 逐位一致（max|Δ|=1.04e-16）。本脚本按旧
  schema 从当前（正确）矩阵重建长表，使 Figure 1(FDR)/2 与全部审计引用一致。

  与 9/16 旧表比对结论：66 条 12×12 行逐位一致（旧表本就取自正确矩阵），
  12 条 fracture 行更新为修复值；p_bonf/q_fdr 均在 78 检验上重算（同旧口径）。

schema（与旧表一致，无 z/sig 列）：trait1,trait2,rg,se,p,p_bonf,q_fdr,block
块顺序：BMD×BMD, FRACTURE×BMD, OA×BMD, OA×FRACTURE, OA×OA；块内按 p 升序。
"""
import os
import numpy as np
import pandas as pd

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "genetics")
BAK = os.path.join(OUT, "_pre_g8_backup", "rg_results_long_13_9月16_bug.csv")


def bh_fdr(p):
    """Benjamini-Hochberg FDR q 值（与 code/11_rg_table.py 同实现）。"""
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


def block_of(a, b):
    """块标签：12×12 对用 OA/BMD 交叉（同 code/11），fracture 对用大写 FRACTURE（同 9/16 旧表）。"""
    if a == "fracture" or b == "fracture":
        other = b if a == "fracture" else a
        return "FRACTURE×BMD" if other.startswith("BMD") else "OA×FRACTURE"
    fa = "OA" if a.startswith("OA_") else "BMD"
    fb = "OA" if b.startswith("OA_") else "BMD"
    return f"{fa}×{fb}"


def main():
    rg = pd.read_csv(os.path.join(OUT, "rg_matrix.csv"), index_col=0)
    se = pd.read_csv(os.path.join(OUT, "rg_se_matrix.csv"), index_col=0)
    pv = pd.read_csv(os.path.join(OUT, "rg_p_matrix.csv"), index_col=0)
    names = list(rg.index)
    assert "fracture" in names, "rg_matrix.csv 应已含 fracture（code/09 G8 修复后为 13×13）"
    assert len(names) == 13, f"expected 13 phenotypes, got {len(names)}"

    rows = []
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = names[i], names[j]
            r, s, p = float(rg.loc[a, b]), float(se.loc[a, b]), float(pv.loc[a, b])
            if not (np.isfinite(r) and np.isfinite(p)):
                continue
            rows.append(dict(trait1=a, trait2=b, rg=r, se=s, p=p,
                             block=block_of(a, b)))
    df = pd.DataFrame(rows)
    n = len(df)
    assert n == 78, f"expected 78 pairs, got {n}"

    df["p_bonf"] = np.minimum(df.p * n, 1.0)
    df["q_fdr"] = bh_fdr(df.p.values)

    rank = {"BMD×BMD": 0, "FRACTURE×BMD": 1, "OA×BMD": 2,
            "OA×FRACTURE": 3, "OA×OA": 4}
    df["_r"] = df.block.map(rank)
    df = (df.sort_values(["_r", "p"])
            .drop(columns="_r")
            .reset_index(drop=True))

    out_path = os.path.join(OUT, "rg_results_long_13.csv")
    df.to_csv(out_path, index=False)
    print(f"✅ 已重建 {out_path}（{len(df)} 行）")
    print(f"   Bonferroni α = {0.05 / n:.3e}; BH-FDR q<0.05: {(df.q_fdr < 0.05).sum()}/{n}")

    # ---------------- 与 9/16 旧表逐行比对 ----------------
    if os.path.exists(BAK):
        old = pd.read_csv(BAK)
        key = ["trait1", "trait2"]
        m = df.merge(old, on=key, suffixes=("_new", "_old"), how="outer",
                     indicator=True)
        assert (m._merge == "both").all(), "新旧表键集合不一致"
        is_fx = (m.trait1 == "fracture") | (m.trait2 == "fracture")
        for col in ("rg", "se", "p"):
            d = (m[f"{col}_new"] - m[f"{col}_old"]).abs()
            print(f"   {col}: 12×12 max|Δ|={d[~is_fx].max():.2e}  "
                  f"fracture max|Δ|={d[is_fx].max():.2e}")
        print(f"   fracture 行逐条（新值 vs 9/16 污染值）:")
        for _, r in m[is_fx].iterrows():
            print(f"     {r.trait1:10s} {r.trait2:16s} "
                  f"rg {r.rg_old:+.4f} → {r.rg_new:+.4f}   "
                  f"p {r.p_old:.3e} → {r.p_new:.3e}")


if __name__ == "__main__":
    main()
