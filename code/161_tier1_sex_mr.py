#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
code/161_tier1_sex_mr.py — Tier-1 F5 电池之 T1：性别分层 OA MR（女 vs 男）
(零新样本；暴露=10 区域 BMD，结局=GO Consortium 性别分层 OA，公开 broadinstitute)

数据：KP.Format.KneeOA / AllOA .gwama.plotting.gz.minus9sout.gz（含 Male.* / Female.* 列）
  EffectSize=OR(每等位基因) -> beta_o = ln(OR)；se_o = |ln(OR)| / Z, Z=norm.ppf(1-P/2)
伪证逻辑：BMD→OA 方向应在两性均为正（方向一致）；两性估计不应显著分歧（无性别交互推翻 SOS）。

尺度统一（G7，2026-10-20）：
  · 暴露侧：BMD munged 文件 N 为占位值 400,000（真实 N=31,873，Qian 2025），
    code/13 按 SE=1/sqrt(N·2p(1-p)) 重建会把效应低估 1/3.5428 倍。
    本脚本与主分析（code/21_mr_corrected_N.py）同尺度：暴露 beta/se ×F_BMD=3.5428。
  · 结局侧：Boer 性别特异 OA 的 beta=ln(OR)、se=|ln(OR)|/Z 由 P 值直接重建，
    已是真 log-OR 尺度（等价于 G6 对 OA_knee 的 K(1-K) 修正后尺度），无需再校正。
  · 注意：F8 用性别特异结局 GWAS + 跨 10 个 BMD 部位的随机效应 meta，
    其绝对量级与单一部位主分析不可直接比拟；推断依据为两性方向一致为正、
    且无显著性别交互（F vs M）。
"""
import os, sys, importlib, json
import numpy as np
import pandas as pd
from scipy.stats import norm, chi2

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
m09 = importlib.import_module("09_rg_matrix")
m13 = importlib.import_module("13_mr_analysis")

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "genetics")
os.makedirs(OUT, exist_ok=True)
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "gwas_stratified")

BMD_SITES = [k for k in m09.PHENOS if k.startswith("BMD_")]

# G7：BMD 暴露尺度校正因子，与 code/21_mr_corrected_N.py 的 F_BMD 完全一致
F_BMD = float(np.sqrt(400000.0 / 31873.0))     # = 3.5428


def load_exposure(site):
    """BMD 暴露工具变量：select_instruments 后施加 ×F_BMD 尺度校正（G7）。"""
    inst = m13.select_instruments(m13.load_munged(site))
    if len(inst):
        inst = inst.copy()
        inst["beta"] = inst["beta"] * F_BMD
        inst["se"] = inst["se"] * F_BMD
    return inst


def load_oa_sex_filtered(path, inst_pos, chunksize=500000):
    """分块读取性别分层 GWAMA，仅保留属于工具变量 (CHR,BP) 的位点 -> 内存安全。
    避免把整份 ~1.6GB gz GWAS 读入内存（原版在受限内存下会 OOM）。"""
    cols = ["Chromosome", "Position", "EffectAllele", "AlternateAllele",
            "Male.EffectSize", "Male.Pvalue", "Female.EffectSize", "Female.Pvalue"]
    out_cols = ["CHR", "BP", "A1", "A2", "beta_f", "se_f", "beta_m", "se_m", "p_f", "p_m"]
    parts = []
    for chunk in pd.read_csv(path, sep="\t", usecols=cols, chunksize=chunksize):
        chunk = chunk.rename(columns={"Chromosome": "CHR", "Position": "BP",
                                      "EffectAllele": "A1", "AlternateAllele": "A2"})
        chunk = chunk.dropna(subset=["CHR", "BP", "A1", "A2"])
        chunk.CHR = chunk.CHR.astype(int); chunk.BP = chunk.BP.astype(int)
        chunk["key"] = list(zip(chunk.CHR.tolist(), chunk.BP.tolist()))
        chunk = chunk[chunk["key"].isin(inst_pos)]
        if len(chunk) == 0:
            continue
        chunk.A1 = chunk.A1.astype(str).str.upper(); chunk.A2 = chunk.A2.astype(str).str.upper()
        for sex in ("Male", "Female"):
            es = chunk[f"{sex}.EffectSize"].clip(lower=1e-6)
            beta = np.log(es.values)
            p = chunk[f"{sex}.Pvalue"].clip(lower=1e-300, upper=1 - 1e-12).values
            z = norm.ppf(1 - p / 2.0)
            z = np.where(np.abs(z) < 1e-6, 1e-6, z)
            se = np.abs(beta) / z
            se = np.where(se <= 0, 1e-6, se)
            tag = "m" if sex == "Male" else "f"
            chunk[f"beta_{tag}"] = beta
            chunk[f"se_{tag}"] = se
            chunk[f"p_{tag}"] = p   # 保留结局 P 值：m13.harmonise() 需访问 out.P
        parts.append(chunk[out_cols])
    if not parts:
        return pd.DataFrame(columns=out_cols)
    return pd.concat(parts, ignore_index=True)


def wmean(b, s):
    w = 1.0 / s**2
    return float((w * b).sum() / w.sum()), float(1.0 / np.sqrt(w.sum()))


def dl_re_meta(b, s):
    """DerSimonian-Laird random-effects meta across the 10 BMD sites.
    Non-estimable (degenerate-SE) sites are excluded upstream before calling."""
    b = np.asarray(b, float); s = np.asarray(s, float)
    if len(b) == 0:
        return dict(meta_beta=float('nan'), meta_se=float('nan'),
                    meta_re_p=float('nan'), tau2=float('nan'), n_sites=0)
    if len(b) == 1:
        p = float(2 * norm.sf(abs(b[0] / s[0])))
        return dict(meta_beta=float(b[0]), meta_se=float(s[0]),
                    meta_re_p=p, tau2=0.0, n_sites=1)
    w = 1.0 / s**2
    beta_fe = float((w * b).sum() / w.sum())
    Q = float((w * (b - beta_fe) ** 2).sum())
    k = len(b); df_ = k - 1
    tau2 = max(0.0, (Q - df_) / (w.sum() - (w ** 2).sum() / w.sum()))
    w2 = 1.0 / (s ** 2 + tau2)
    beta_re = float((w2 * b).sum() / w2.sum())
    se_re = float(np.sqrt(1.0 / w2.sum()))
    p_re = float(2 * norm.sf(abs(beta_re / se_re)))
    return dict(meta_beta=beta_re, meta_se=se_re, meta_re_p=p_re,
                tau2=tau2, n_sites=k)


def run_sex(outcome_label, path):
    print(f"\n########## T1 结局={outcome_label} ##########")
    # 预收集全部 BMD 工具变量的 (CHR,BP)，仅加载这些位点到内存（避免 OOM）
    inst_pos = set()
    for site in BMD_SITES:
        inst = load_exposure(site)
        if len(inst):
            inst_pos.update(zip(inst.CHR.astype(int).tolist(), inst.BP.astype(int).tolist()))
    print(f"  工具变量唯一 (CHR,BP) 数: {len(inst_pos):,}")
    d_oa = load_oa_sex_filtered(path, inst_pos)
    print(f"  OA 匹配位点总数: {len(d_oa):,}")
    rows = []
    per_site = {sx: {} for sx in ("Female", "Male")}
    for site in BMD_SITES:
        inst = load_exposure(site)
        if len(inst) == 0:
            continue
        for sex, bcol, scol, pcol in (("Female", "beta_f", "se_f", "p_f"),
                                       ("Male", "beta_m", "se_m", "p_m")):
            out = d_oa[["CHR", "BP", "A1", "A2", bcol, scol, pcol]].rename(
                columns={bcol: "beta", scol: "se", pcol: "P"})
            h = m13.harmonise(inst, out)
            if len(h) < 3:
                continue
            r = m13.run_mr(h.beta_e.values, h.se_e.values, h.beta_o.values, h.se_o.values)
            if "error" in r:
                continue
            # 非可估计位点：结局 GWAS 在该位点 OR≈1 且 p≈1 -> beta=0, SE 被 floor 到 ~1e-6。
            # 逆方差权重会将其放大为无穷权，主导 meta 使其虚假归零。此类位点不提供方向信息，
            # 从跨位点 meta 中剔除（仍保留于 per_site 以备透明披露），符合 meta 分析常规。
            degenerate = (abs(r["ivw"]) < 1e-3) and (r["ivw_se"] <= 1e-3)
            per_site[sex][site] = dict(beta=r["ivw"], se=r["ivw_se"],
                                       re_p=r.get("ivw_re_p"), n=r.get("n_snp"),
                                       egger_int_p=r.get("egger_intercept_p"),
                                       non_informative=bool(degenerate))
            if degenerate:
                print(f"    [skip meta] {site} {sex}: non-informative (beta={r['ivw']:.2e}, se={r['ivw_se']:.2e})")
                continue
            rows.append(dict(outcome=outcome_label, exposure=site, sex=sex,
                             beta=r["ivw"], se=r["ivw_se"], re_p=r.get("ivw_re_p"),
                             n_snp=r.get("n_snp"), egger_int_p=r.get("egger_intercept_p")))

    df = pd.DataFrame(rows)
    # 每性别跨可估计位点 meta（随机效应 DerSimonian-Laird）
    res = {}
    for sex in ("Female", "Male"):
        g = df[df.sex == sex]
        if len(g) >= 2:
            m = dl_re_meta(g.beta.values, g.se.values)
            m["n_snp_total"] = int(g.n_snp.sum())
            res[sex] = m
            print(f"  [{sex}] RE meta β={m['meta_beta']:+.4f} (SE {m['meta_se']:.4f}) "
                  f"RE_P={m['meta_re_p']:.2e} (n_site={m['n_sites']}, tau²={m['tau2']:.4f})")
    # 性别异质性：meta_F vs meta_M
    if "Female" in res and "Male" in res:
        bf, sf = res["Female"]["meta_beta"], res["Female"]["meta_se"]
        bm, sm = res["Male"]["meta_beta"], res["Male"]["meta_se"]
        d_ = bf - bm; s_ = np.sqrt(sf**2 + sm**2)
        p_ = float(2*norm.sf(abs(d_/s_)))
        print(f"  [F vs M] diff={d_:+.4f} (SE {s_:.4f}) 双侧 P={p_:.3f} -> "
              f"{'两性显著分歧(性别交互)' if p_<0.05 else '两性一致(无显著性别交互)'}")
        res["F_vs_M"] = dict(diff=round(d_, 4), se=round(s_, 4), p_two_sided=round(p_, 4))
    # 每位点 F vs M 差异 + 位点级性别异质性 Cochran Q（H0: 所有位点 F-M 差=0）
    sites_both = sorted(set(per_site["Female"]) & set(per_site["Male"]))
    site_het = []
    for s in sites_both:
        bf, sf = per_site["Female"][s]["beta"], per_site["Female"][s]["se"]
        bm, sm = per_site["Male"][s]["beta"], per_site["Male"][s]["se"]
        dd = bf - bm; ss = np.sqrt(sf**2 + sm**2)
        site_het.append(dict(site=s, beta_F=bf, beta_M=bm, diff=dd,
                             se=ss, p_two=2*norm.sf(abs(dd/ss))))
    Q_het = None
    if len(site_het) >= 2:
        sh = pd.DataFrame(site_het)
        # 加权均值对比（meta 的 F-M 差）
        w = 1.0 / sh["se"]**2
        c_bar = float((w * sh["diff"]).sum() / w.sum())
        Q_het = float((w * (sh["diff"] - c_bar)**2).sum())
        Qp_het = float(chi2.sf(Q_het, len(sh) - 1))
        print(f"  [位点级性别异质性] Q={Q_het:.2f} (df={len(sh)-1}, P={Qp_het:.3f}) -> "
              f"{'存在显著位点性别交互' if Qp_het<0.05 else '无显著位点性别交互'}")
        res["site_het_Q"] = round(Q_het, 2); res["site_het_P"] = round(Qp_het, 4)
    df.to_csv(os.path.join(OUT, f"tier1_sex_mr_{outcome_label}.csv"), index=False)
    return dict(outcome=outcome_label, meta=res,
                per_site_Female=per_site["Female"], per_site_Male=per_site["Male"],
                site_both=len(sites_both),
                site_het=site_het if site_het else None)


results = {}
for label, fn in (("KneeOA", "KneeOA_sexspecific.tsv.gz"),
                  ("AllOA", "AllOA_sexspecific.tsv.gz")):
    p = os.path.join(DATA, fn)
    if os.path.exists(p):
        results[label] = run_sex(label, p)
    else:
        print(f"!! 缺少 {p}")

with open(os.path.join(OUT, "tier1_sex_mr_summary.json"), "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\n✅ T1 输出:", os.path.join(OUT, "tier1_sex_mr_*.csv / _summary.json"))
