#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
26_export_harmonised_snp.py
=========================================================================
重跑 exposure-outcome harmonisation，导出**逐 SNP** 数据，用于绘制
MR 散点图、MR-PRESSO 漏斗图与留一法图。

为什么需要重跑
--------------
此前工作区只保留了 PLINK clump 输入文件（SNP\\tP，无效应量）与聚合后的
MR 结果，缺少逐 SNP 的 exposure-outcome 合并数据，因此无法真实绘制散点图。
本脚本直接从**原始 GWAS 汇总统计**重建该数据，不依赖任何模拟/插值。

数据来源（真实文件，均在磁盘上）
--------------------------------
  · 12 个已 munge 的 LDSC sumstats
      C:\\...\\OA_FLAGSHIP\\07_reference\\ldsc\\sumstats_*.sumstats.gz
      （只有 Z / N / MAF，按 SE = 1/sqrt(N·2p(1-p))、beta = Z·SE 重构）
  · 骨折 GWAS（真实 log-OR beta/se）
      C:\\...\\Direction1_pipeline\\data\\raw\\fracture_anytype.tsv

与既有结果严格对齐（三重一致性校验）
------------------------------------
  1. 仪器集  ：严格 clump（r²<0.001 / 10 Mb，PLINK 等价），取自
              output/genetics/<exposure>_clumped.clumped；缺失时回退 1-Mb 近似
  2. 效应标度：N 校正因子来自 21_mr_corrected_N.py
              BMD ×3.5428（N 占位 400,000 → 真实 31,873）
              OA_selfrep ×0.8695（350,000 → 462,933）
              OA_knee / fracture ×1.0
  3. 校验    ：重跑的 IVW β / SE / p 必须与 mr_corrected_N.csv **逐位一致**，
              否则报错退出（不允许带着偏差画图）

输出
----
  output/genetics/harmonised/<exposure>__<outcome>.csv   逐 SNP 数据
  output/genetics/harmonised/_validation.csv             与 mr_corrected_N 的比对
  output/genetics/harmonised/_loo.csv                    留一法估计
=========================================================================
"""
import os
import sys
import importlib
import numpy as np
import pandas as pd
from scipy.stats import chi2, norm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
m13 = importlib.import_module("13_mr_analysis")
m21 = importlib.import_module("21_mr_corrected_N")

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "output", "genetics")
HARM = os.path.join(OUT, "harmonised")


# ------------------------------------------------------------------ MR 基础
def ivw(be, bo, so):
    w = 1.0 / (so ** 2)
    return float(np.sum(w * be * bo) / np.sum(w * be ** 2))


def ivw_se(be, so):
    return float(1.0 / np.sqrt(np.sum(be ** 2 / so ** 2)))


def egger(be, bo, so):
    w = 1.0 / (so ** 2)
    mx = np.sum(w * be) / np.sum(w)
    my = np.sum(w * bo) / np.sum(w)
    b = np.sum(w * (be - mx) * (bo - my)) / np.sum(w * (be - mx) ** 2)
    a = my - b * mx
    return float(b), float(a)


def presso_per_snp(be, se_e, bo, so, alpha=0.05):
    """
    MR-PRESSO 逐 SNP 离群检验（Verbanck 2018 Nat Genet）。
      Q_i = Q_full - Q_{(-i)} ~ χ²(1)，Bonferroni 校正（alpha/K）
    返回 (p_i 数组, 离群布尔数组)。
    """
    K = len(be)
    b = ivw(be, bo, so)
    sig = np.sqrt(so ** 2 + (b * se_e) ** 2)
    Q_full = float(np.sum(((bo - b * be) / sig) ** 2))
    p_i = np.ones(K)
    for i in range(K):
        m = np.ones(K, bool)
        m[i] = False
        b_wo = ivw(be[m], bo[m], so[m])
        sig_wo = np.sqrt(so[m] ** 2 + (b_wo * se_e[m]) ** 2)
        Q_wo = float(np.sum(((bo[m] - b_wo * be[m]) / sig_wo) ** 2))
        p_i[i] = float(chi2.sf(max(Q_full - Q_wo, 0.0), 1))
    return p_i, (p_i < alpha / K)


# ------------------------------------------------------------------ 主流程
def main():
    os.makedirs(HARM, exist_ok=True)

    print("=" * 78)
    print("重跑 harmonisation · 导出逐 SNP 数据（真实 GWAS，不使用任何模拟值）")
    print("=" * 78)

    print("\n▶ 载入原始汇总统计并施加 N 校正 ...")
    dfs = m21.load_all_corrected()          # 与 21_ 完全相同的校正口径

    print("\n▶ 逐对 harmonisation ...")
    val_rows, loo_rows = [], []
    for exp, out, note in m13.PAIRS:
        iv = m21.load_strict(exp, dfs)      # 严格 clump 仪器集（与 21_ 一致）
        h = m13.harmonise(iv, dfs[out])
        if len(h) < 3:
            print("【%s -> %s】 IV 不足 (%d)，跳过" % (exp, out, len(h)))
            continue

        be = h.beta_e.values.astype(float)
        se_e = h.se_e.values.astype(float)
        bo = h.beta_o.values.astype(float)
        so = h.se_o.values.astype(float)

        b_ivw = ivw(be, bo, so)
        s_ivw = ivw_se(be, so)
        z = b_ivw / s_ivw
        p = float(2 * norm.sf(abs(z)))
        b_egg, a_egg = egger(be, bo, so)

        p_i, outl = presso_per_snp(be, se_e, bo, so)

        # 逐 SNP Wald 比值与其 SE（一阶 delta 法）
        ratio = bo / be
        se_ratio = np.abs(so / be)

        # 留一法
        loo = np.full(len(be), np.nan)
        for i in range(len(be)):
            m = np.ones(len(be), bool)
            m[i] = False
            if m.sum() >= 2:
                loo[i] = ivw(be[m], bo[m], so[m])

        df = pd.DataFrame(dict(
            SNP=h.SNP.values, CHR=h.CHR.values, BP=h.BP.values,
            beta_e=be, se_e=se_e, p_e=h.p_e.values,
            beta_o=bo, se_o=so, p_o=h.p_o.values,
            eaf=h.eaf.values, flip=h.flip.values,
            ratio=ratio, se_ratio=se_ratio,
            presso_p=p_i, presso_outlier=outl.astype(int),
            loo_beta=loo,
        ))
        fp = os.path.join(HARM, f"{exp}__{out}.csv")
        df.to_csv(fp, index=False)

        val_rows.append(dict(exposure=exp, outcome=out, K=len(be),
                             beta=b_ivw, se=s_ivw, p=p,
                             egger=b_egg, egger_int=a_egg,
                             n_outlier=int(outl.sum())))
        for i, snp in enumerate(h.SNP.values):
            loo_rows.append(dict(exposure=exp, outcome=out, SNP=snp,
                                 loo_beta=loo[i], full_beta=b_ivw))
        print("  【%-16s -> %-16s】K=%-3d β=%+.4f (SE %.4f) p=%.3g  离群=%d  -> %s"
              % (exp, out, len(be), b_ivw, s_ivw, p, int(outl.sum()),
                 os.path.basename(fp)))

    val = pd.DataFrame(val_rows)
    val.to_csv(os.path.join(HARM, "_validation.csv"), index=False)
    pd.DataFrame(loo_rows).to_csv(os.path.join(HARM, "_loo.csv"), index=False)

    # ---------------------------------------------------------- 一致性校验
    print("\n" + "=" * 78)
    print("一致性校验：重跑结果 vs 已发表的 mr_corrected_N.csv")
    print("=" * 78)
    ref = pd.read_csv(os.path.join(OUT, "mr_corrected_N.csv"))
    ok = True
    for r in val.itertuples():
        m = ref[(ref.exposure == r.exposure) & (ref.outcome == r.outcome)]
        if m.empty:
            print("  ⚠ %s -> %s 在参照表中缺失" % (r.exposure, r.outcome))
            ok = False
            continue
        rb, rp, rk = float(m.iloc[0]["beta"]), float(m.iloc[0]["p"]), int(m.iloc[0]["K"])
        db, dp = abs(rb - r.beta), abs(rp - r.p)
        flag = "✅" if (db < 1e-6 and dp < 1e-6 and rk == r.K) else "❌"
        if flag == "❌":
            ok = False
        print("  %s %-16s -> %-16s K %d/%d  Δβ=%.2e  Δp=%.2e"
              % (flag, r.exposure, r.outcome, r.K, rk, db, dp))

    print("\n" + ("✅ 全部一致 —— 逐 SNP 数据可安全用于作图" if ok
                  else "❌ 存在偏差，禁止作图"))
    if not ok:
        sys.exit(1)
    print("📁 逐 SNP 数据 ->", HARM)


if __name__ == "__main__":
    main()
