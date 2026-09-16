#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
10_fracture_row.py — 把真实骨折 GWAS 接入 12 表型 rg 矩阵（第 13 个表型）

难点：fracture_anytype.tsv 的 rs_id 100% 为 `chr:pos_ea_oa` 格式（无 dbSNP rs），
      无法与 LD 参考按 rsid 匹配。
解法：LD 参考同时含 CHR/BP/SNP -> 建 (chr,pos) -> rsid 桥接表，
      把骨折 GWAS 先映射到 rsid，再与其余 12 个表型走同一套等位基因对齐。

用法：
  python 10_fracture_row.py                 # N=400,000（默认）
  python 10_fracture_row.py --sens-n        # N 敏感性：100k/200k/400k/800k
"""
import os
import sys
import argparse
import importlib
import numpy as np
import pandas as pd
from scipy.stats import norm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
m09 = importlib.import_module("09_rg_matrix")
m08 = importlib.import_module("08_genetic_correlation")

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(ROOT, "output", "genetics")
FRACTURE = (r"C:\Users\Administrator\WorkBuddy\2026-09-09-14-17-27\Direction1_pipeline"
            r"\Direction1_pipeline\data\raw\fracture_anytype.tsv")
DEFAULT_N = 400_000


def load_ld_bridge():
    """LD 参考 -> DataFrame[chr,pos,SNP,L2]，并给出 (chr,pos)->SNP 桥接。"""
    frames = []
    for c in range(1, 23):
        p = os.path.join(m09.FLAGSHIP_LDSC, f"eur_w_ld_chr_{c}.l2.ldscore.gz")
        d = pd.read_csv(p, sep=r"\s+", usecols=["CHR", "BP", "SNP", "L2"])
        frames.append(d)
    df = pd.concat(frames, ignore_index=True)
    df = df.dropna(subset=["L2"])
    df = df[df.L2 > 0]
    return df


def load_fracture(path=FRACTURE):
    """读 harmonised 骨折 GWAS，Z = beta/se。"""
    df = pd.read_csv(
        path, sep="\t",
        usecols=["chromosome", "base_pair_location", "effect_allele",
                 "other_allele", "beta", "standard_error",
                 "effect_allele_frequency", "p_value"],
        dtype={"chromosome": "int32", "base_pair_location": "int64",
               "effect_allele": "category", "other_allele": "category",
               "beta": "float64", "standard_error": "float64"},
    )
    df = df.dropna(subset=["beta", "standard_error"])
    df = df[df.standard_error > 0]
    df["Z"] = df.beta / df.standard_error
    df["A1"] = df.effect_allele.astype(str).str.upper()
    df["A2"] = df.other_allele.astype(str).str.upper()
    return df[["chromosome", "base_pair_location", "A1", "A2", "Z",
               "effect_allele_frequency"]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sens-n", action="store_true", help="做 N 敏感性分析")
    ap.add_argument("--blocks", type=int, default=100)
    args = ap.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)
    print("▶ 建 LD 桥接 (chr,pos) -> rsid")
    ld_df = load_ld_bridge()
    print(f"  LD 参考 {len(ld_df):,} SNP")

    print("▶ 读骨折 GWAS")
    fx = load_fracture()
    print(f"  骨折 GWAS {len(fx):,} 变异")

    print("▶ 按 (chr,pos) 桥接到 rsid")
    merged = fx.merge(
        ld_df[["CHR", "BP", "SNP", "L2"]],
        left_on=["chromosome", "base_pair_location"],
        right_on=["CHR", "BP"], how="inner",
    )
    print(f"  桥接成功 {len(merged):,} / {len(fx):,} "
          f"({100*len(merged)/len(fx):.1f}%)  ← 低则说明基因组版本不匹配")
    merged = merged.drop_duplicates(subset=["SNP"], keep="first")
    fx_ss = merged.set_index("SNP")[["A1", "A2", "Z", "effect_allele_frequency"]]
    fx_ss = fx_ss.rename(columns={"effect_allele_frequency": "MAF"})
    fx_ss["N"] = float(DEFAULT_N)          # harmonise 需要 N 列

    # 逐染色体诊断：版本不匹配会表现为某些染色体匹配率骤降
    per = (merged.groupby("CHR").size()
                 .reindex(range(1, 23), fill_value=0)
                 .rename("matched"))
    ref_per = (ld_df.groupby("CHR").size()
                    .reindex(range(1, 23), fill_value=0).rename("in_ld_ref"))
    chk = pd.concat([ref_per, per], axis=1)
    chk["pct"] = (100 * chk.matched / chk.in_ld_ref).round(1)
    print("  逐染色体匹配率（LD 参考内占比 %）：")
    print("    " + "  ".join(f"c{c}={v}" for c, v in chk.pct.items()))

    # ---- 与 12 表型求共同交集 + 等位基因对齐（以 OA_knee 为基准）
    print("▶ 与 12 表型对齐")
    dfs = {nm: m09.load_sumstats(fn) for nm, fn in m09.PHENOS.items()}
    dfs["fracture"] = fx_ss
    ld = ld_df.set_index("SNP")["L2"].to_dict()
    Zdf, ell, Ns, diag = m09.harmonise(dfs, ld)
    Ns["fracture"] = DEFAULT_N
    print(f"  最终共同 SNP = {diag['n_final']:,}")
    print(f"  fracture z 翻转数 = {diag['n_flipped'].get('fracture', 'NA')}")

    M = m09.load_M()
    names = list(Zdf.columns)

    def run_one(N_fx):
        h2, i2, _, _ = m08.compute_h2(Zdf["fracture"].values, ell, N_fx, M,
                                      3, args.blocks)
        rows = []
        for nm in names:
            if nm == "fracture":
                continue
            h1, ih1, _, _ = m08.compute_h2(Zdf[nm].values, ell, Ns[nm], M,
                                           3, args.blocks)
            gc, gse, ig, _ = m08.compute_gencov(
                Zdf["fracture"].values, Zdf[nm].values, ell, N_fx, Ns[nm], M,
                h1=h2, h2=h1, i_h1=i2, i_h2=ih1, n_iter=3, n_blocks=args.blocks)
            den = np.sqrt(max(h2, 1e-12) * max(h1, 1e-12))
            rg = gc / den
            se = abs(gse / den) if np.isfinite(gse) else np.nan
            p = 2 * norm.sf(abs(rg / se)) if (np.isfinite(se) and se > 0) else np.nan
            rows.append(dict(trait=nm, h2_trait=h1, rg=rg, se=se, p=p,
                             n_snp=len(ell)))
        return h2, pd.DataFrame(rows)

    print(f"▶ LDSC（fracture N={DEFAULT_N:,}）")
    h2f, tab = run_one(DEFAULT_N)
    print(f"  fracture h2 = {h2f:.4f}")
    for _, r in tab.sort_values("rg").iterrows():
        print(f"    fracture ~ {r.trait:18s} rg={r.rg:+.3f} (SE {r.se:.3f}) P={r.p:.2e}")
    tab["h2_fracture"] = h2f
    tab.to_csv(os.path.join(OUT_DIR, "rg_fracture.csv"), index=False)

    if args.sens_n:
        print("▶ N 敏感性（验证 rg 对 N 不敏感）")
        out = []
        for N in (100_000, 200_000, 400_000, 800_000):
            _, t = run_one(N)
            out.append(dict(N=N, **{r.trait: r.rg for _, r in t.iterrows()}))
            print(f"    N={N:,}: " + " ".join(
                f"{r.trait.replace('BMD_','')}={r.rg:+.3f}" for _, r in t.iterrows()))
        pd.DataFrame(out).to_csv(os.path.join(OUT_DIR, "rg_fracture_N_sens.csv"),
                                 index=False)
    print(f"✅ 输出至 {OUT_DIR}")


if __name__ == "__main__":
    main()
