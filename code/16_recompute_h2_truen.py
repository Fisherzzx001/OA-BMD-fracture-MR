#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
16_recompute_h2_truen.py — 用真实样本量 N 重算 SNP-h²（Item 1 关键修正）

背景
----
原 rg_h2.csv 中 10 个 BMD 位点与骨折的 N 列被填为占位值 400,000，
但该值与观测信号强度矛盾（见 docs/08 稿 §4.5 / §5）。
经核实：
  · 10 个 BMD 部位来自 Qian et al. 2025 (Genomics Proteomics & Bioinformatics,
    PMID 41206123) —— UKB 中 >30,000 名欧洲人的 11 个 DXA-BMD 位点 GWAS。
    真实 N ≈ 32,000（DXA 亚组）。
  · OA_knee 的 N=1,208,397 与 OA_selfrep 的 N=350,000 取自汇总统计文件，为真实值。
  · 骨折 GWAS（Neale-lab UKB 格式）真实 N 待定（λ_gc≈1.056 提示为数万级，
    非 400,000）；本脚本以 400,000 占位并明确标注。

关键性质
--------
rg = gencov / sqrt(h1·h2) 与 N、M 无关（已实证：N 变 6.7× 三位小数不变）。
因此只重算 h²（N 依赖项），rg 矩阵保持不变，无需重跑。

输出
----
output/genetics/rg_h2_trueN.csv  (phenotype, N_used, h2, h2_se, intercept, intercept_se)
并打印旧值(占位 N) -> 新值(真实 N) 的校正因子。
"""
import os, sys, importlib, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
m09 = importlib.import_module("09_rg_matrix")
m08 = importlib.import_module("08_genetic_correlation")

OUT = os.path.join(m09.ROOT, "output", "genetics")

# 真实 N（BMD×10 已确认 Qian 2025 ≈ 32,000；其余沿用文件真实值；骨折占位待核）
TRUE_N = {
    "OA_knee": 1_208_397,
    "OA_selfrep": 350_000,
    "BMD_total": 32_000, "BMD_femur_neck": 32_000, "BMD_femur_total": 32_000,
    "BMD_lumbar": 32_000, "BMD_spine": 32_000, "BMD_arms": 32_000,
    "BMD_leg": 32_000, "BMD_pelvis": 32_000, "BMD_rib": 32_000,
    "BMD_trunk": 32_000,
}
FRACTURE_N = 400_000  # 占位；真实 N 待 Neale-lab manifest 确认

def main():
    ld = m09.load_ld_ref()
    M = m09.load_M()
    print(f"LD 参考 SNP={len(ld):,}  M={M:,.0f}")

    # ---- 12 表型在共同 SNP 集上重算 h²（与 rg 矩阵同一 SNP 集，仅改 N）----
    dfs = {nm: m09.load_sumstats(m09.PHENOS[nm]) for nm in m09.PHENOS}
    Zdf, ell, Ns, diag = m09.harmonise(dfs, ld)
    print(f"共同 SNP 集={len(Zdf):,}  翻转={diag['n_flipped']}")

    rows = []
    for nm in Zdf.columns:
        N = TRUE_N.get(nm, Ns[nm])
        r = m08.compute_h2(Zdf[nm].values, ell, N, M, n_iter=3, n_blocks=100)
        rows.append((nm, N, r[0], r[1], r[2], r[3]))
        print(f"  {nm:16s} N={N:>9,.0f}  h2={r[0]:.4f} (SE {r[1]:.4f}) intercept={r[2]:.3f}")

    df = pd.DataFrame(rows, columns=["phenotype", "N_used", "h2", "h2_se",
                                     "intercept", "intercept_se"])

    # ---- 骨折 h²（单独；rg 已在 rg_matrix_13.csv，N 无关）----
    try:
        import importlib as _il
        m10 = _il.import_module("10_fracture_row")
        bridge = m10.load_ld_bridge()
        pos2snp = {(int(r.CHR), int(r.BP)): r.SNP for r in bridge.itertuples()}
        fx = m10.load_fracture()
        fx["poskey"] = list(zip(fx.chromosome.astype(int),
                                fx.base_pair_location.astype(int)))
        fx = fx[fx.poskey.isin(pos2snp)].copy()
        fx["SNP"] = fx.poskey.map(pos2snp)
        snps = [s for s in fx.SNP if s in ld]
        fsub = fx.set_index("SNP").loc[snps]
        keep = [not m09._is_palindromic(a, b) for a, b in
                zip(fsub.effect_allele, fsub.other_allele)]
        fsub = fsub[keep]
        z = fsub.Z.values.astype(float)
        ellf = np.array([ld[s] for s in fsub.index], dtype=float)
        rf = m08.compute_h2(z, ellf, FRACTURE_N, M, n_iter=3, n_blocks=100)
        df.loc[len(df)] = ["fracture_any", FRACTURE_N, rf[0], rf[1], rf[2], rf[3]]
        print(f"  {'fracture_any':16s} N={FRACTURE_N:>9,.0f}  h2={rf[0]:.4f} "
              f"(SE {rf[1]:.4f}) intercept={rf[2]:.3f}  [占位 N 待核]")
    except Exception as e:
        print(f"  ⚠ 骨折 h² 重算跳过: {e}")

    df.to_csv(os.path.join(OUT, "rg_h2_trueN.csv"), index=False)

    # ---- 与旧值(占位 N)对比，打印校正因子 ----
    old = pd.read_csv(os.path.join(OUT, "rg_h2.csv"))
    old_idx = old.set_index("phenotype")["h2"].to_dict()
    print("\n=== 校正因子 (旧 h² × N_ratio) ===")
    for _, r in df.iterrows():
        nm = r.phenotype
        if nm in old_idx and nm.startswith("BMD"):
            ratio = 400_000 / r.N_used
            print(f"  {nm:16s} 旧 h2={old_idx[nm]:.4f} × {ratio:5.2f} = 新 {r.h2:.4f}")
    print(f"\n✅ 输出 -> {os.path.join(OUT, 'rg_h2_trueN.csv')}")


if __name__ == "__main__":
    main()
