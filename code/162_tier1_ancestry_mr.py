#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
code/162_tier1_ancestry_mr.py — Tier-1 F5 电池之 T2：祖先敏感性（东亚 EAS OA）
(零新样本；暴露=10 区域 BMD [EUR 工具变量]，结局=非欧 OA GWAS [EAS])

数据：Hatzikotoulas 2025 Nature 多祖先 OA 的东亚(EAS)分层汇总统计（GCST，EBI GWAS）。
  解析：EffectSize=OR(每等位基因) -> beta_o=ln(OR)；se_o=|ln(OR)|/Z, Z=norm.ppf(1-P/2)。
伪证逻辑：BMD→OA 方向应在非欧人群中仍为正（方向一致）；若反转则祖先特异性生物学推翻 SOS。
重要局限：EUR 工具变量在 EAS 的 LD 结构未重估（正式分析应换 EAS LD 重 clump）；
  本脚本为敏感性检验，方向一致性结论稳健，绝对效应量受跨祖先 LD 差异影响。
用法：将 EAS OA GWAMA/tsv 放到 data/gwas_stratified/EAS_OA.tsv.gz 后运行。
  文件需含列：Chromosome,Position,EffectAllele,AlternateAllele,EffectSize,Pvalue
  （若为 GO Consortium gwama 格式含 Male./Female. 列，则取合并 EffectSize/Pvalue 列）。
"""
import os, sys, importlib, json
import numpy as np
import pandas as pd
from scipy.stats import norm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
m09 = importlib.import_module("09_rg_matrix")
m13 = importlib.import_module("13_mr_analysis")

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "genetics")
os.makedirs(OUT, exist_ok=True)
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "gwas_stratified")
BMD_SITES = [k for k in m09.PHENOS if k.startswith("BMD_")]


def load_oa_eas(path):
    """解析 EAS OA GWAMA/tsv -> CHR,BP,A1,A2,beta,se（合并效应）。"""
    # 兼容 GO gwama 格式（含 Male./Female.）与通用格式
    d = pd.read_csv(path, sep="\t")
    d = d.rename(columns={"Chromosome": "CHR", "Position": "BP",
                          "EffectAllele": "A1", "AlternateAllele": "A2"})
    d = d.dropna(subset=["CHR", "BP", "A1", "A2"])
    d.CHR = d.CHR.astype(int); d.BP = d.BP.astype(int)
    d.A1 = d.A1.astype(str).str.upper(); d.A2 = d.A2.astype(str).str.upper()
    # 优先用合并 EffectSize/Pvalue；若无则取 Male.（同为 OR）
    es = d.get("EffectSize", d.get("Male.EffectSize"))
    pv = d.get("Pvalue", d.get("Male.Pvalue"))
    es = pd.to_numeric(es, errors="coerce").clip(lower=1e-6)
    pv = pd.to_numeric(pv, errors="coerce").clip(lower=1e-300, upper=1 - 1e-12)
    beta = np.log(es.values)
    z = norm.ppf(1 - pv.values / 2.0)
    z = np.where(np.abs(z) < 1e-6, 1e-6, z)
    se = np.abs(beta) / z
    se = np.where(se <= 0, 1e-6, se)
    return pd.DataFrame({"CHR": d.CHR, "BP": d.BP, "A1": d.A1, "A2": d.A2,
                         "beta": beta, "se": se})


def wmean(b, s):
    w = 1.0 / s**2
    return float((w * b).sum() / w.sum()), float(1.0 / np.sqrt(w.sum()))


def run(outcome_label, path):
    print(f"\n########## T2 祖先敏感性 结局={outcome_label} (EAS OA) ##########")
    d_oa = load_oa_eas(path)
    print(f"  EAS OA 位点总数: {len(d_oa):,}")
    rows = []
    per_site = {}
    for site in BMD_SITES:
        inst = m13.select_instruments(m13.load_munged(site))
        if len(inst) == 0:
            continue
        out = d_oa[["CHR", "BP", "A1", "A2", "beta", "se"]]
        h = m13.harmonise(inst, out)
        if len(h) < 3:
            continue
        r = m13.run_mr(h.beta_e.values, h.se_e.values, h.beta_o.values, h.se_o.values)
        if "error" in r:
            continue
        per_site[site] = dict(beta=r["ivw"], se=r["ivw_se"], re_p=r.get("ivw_re_p"),
                              n=r.get("n_snp"), egger_int_p=r.get("egger_intercept_p"))
        rows.append(dict(exposure=site, beta=r["ivw"], se=r["ivw_se"],
                         re_p=r.get("ivw_re_p"), n_snp=r.get("n_snp")))
        print(f"  {site:16s} β={r['ivw']:+.4f} (SE {r['ivw_se']:.4f}) "
              f"RE_P={r.get('ivw_re_p'):.2e} n={r.get('n_snp')} Egger_int_P={r.get('egger_intercept_p'):.3f}")
    df = pd.DataFrame(rows)
    if len(df) >= 2:
        bf, sf = wmean(df.beta.values, df.se.values)
        n_pos = int((df.beta > 0).sum())
        print(f"  [跨位点 meta] β={bf:+.4f} (SE {sf:.4f}) RE_P={2*norm.sf(abs(bf/sf)):.2e} "
              f"n_site={len(df)} 阳性位点={n_pos}/{len(df)}")
        res = dict(meta_beta=bf, meta_se=sf, meta_re_p=2*norm.sf(abs(bf/sf)),
                   n_sites=len(df), n_pos=n_pos,
                   conclusion=("EAS 方向一致为正(支持 SOS)" if (bf > 0 and 2*norm.sf(abs(bf/sf)) < 0.05)
                               else "EAS 未达显著/方向需复核"))
    else:
        res = dict(note="harmonised sites <2, 无法 meta")
    df.to_csv(os.path.join(OUT, f"tier1_ancestry_{outcome_label}.csv"), index=False)
    return dict(outcome=outcome_label, meta=res, per_site=per_site)


if __name__ == "__main__":
    p = os.path.join(DATA, "EAS_OA.tsv.gz")
    if not os.path.exists(p):
        print(f"!! 缺少 EAS OA 文件: {p}\n"
              f"   请先获取 Hatzikotoulas 2025 Nature 东亚(EAS)分层 OA 汇总统计"
              f"（EBI GWAS GCST，ancestry-stratified），放到该路径后运行。")
    else:
        out = run("EAS_OA", p)
        with open(os.path.join(OUT, "tier1_ancestry_summary.json"), "w") as f:
            json.dump(out, f, indent=2, default=str)
        print("\n✅ T2 输出:", os.path.join(OUT, "tier1_ancestry_*.csv / _summary.json"))
