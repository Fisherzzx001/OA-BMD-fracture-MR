#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
312_f10_eur_sensitivity.py  --  F10 跨祖先敏感性（欧裔限定结局）

目标：用与正文完全一致的 31 个 BMD 工具变量 + IVW 引擎（code/13），
将 BMD_total → knee-OA 的两样本 MR 在「欧裔限定结局 GCST90566811」上重跑，
与正文使用的「混合祖先结局 GCST90566800（ headlines OR/SD 1.0812）」比较。

方法学要点：
  · 31 个工具变量取自 bmd_instruments.json['inst']（chr/bp/beta_e/se_e/exp_allele）。
  · 结局效应取 GWAS Catalog 的 odds_ratio，log-OR = ln(OR)，SE 由 95% CI 反解
    se = (ln(ci_upper) - ln(ci_lower)) / 3.92 —— 与 legacy knee_matched.tsv 一致，
    绕开 GWAS Catalog 的 standard_error 列（OR 尺度歧义）。
  · 等位基因按 exp_allele（BMD 效应等位基因）定向：E==exp → 不翻；O==exp → 翻；都不等 → 丢弃。
  · 先对混合结局重跑并与 harmonised/BMD_total__OA_knee.csv 的 beta_o 逐 SNP 校验，
    确认解析+定向无误，再报告欧裔限定结果。
  · 不改任何已报 β/SE/p/OR/CI：本脚本产物是新增敏感性，正文 headline 数值保持不变。
"""
import os, json, gzip, math, importlib.util
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.abspath(__file__))
CACHE = "G:/OA_B_PROJECT_cache/l7"
INST = os.path.join(CACHE, "bmd_instruments.json")
MIXED = os.path.join(CACHE, "GCST90566800.tsv.gz")          # 混合祖先 knee OA（正文 headline）
EUR = os.path.join(CACHE, "joints", "GCST90566811.tsv.gz")   # 欧裔限定 knee OA（F10）
HARM = os.path.join(ROOT, "output", "genetics", "harmonised", "BMD_total__OA_knee.csv")

# ---- 载入 code/13 的 MR 引擎（IVW/Egger/加权中位数/Radial/留一法） ----
spec = importlib.util.spec_from_file_location("mr13", os.path.join(ROOT, "13_mr_analysis.py"))
mr13 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mr13)

COLS = ["chromosome", "base_pair_location", "effect_allele", "other_allele",
        "odds_ratio", "standard_error", "effect_allele_frequency", "p_value",
        "variant_id", "ci_upper", "ci_lower", "NCASES", "NCONTROLS"]


def load_inst():
    d = json.load(open(INST))
    return d["inst"]   # rsid -> {chr,bp,beta_e,se_e,beta_o,se_o,eaf,exp_allele}


def parse_outcome(path, targets):
    """流式读取 gz GWAS Catalog 文件，返回 { (chr,bp): row } 与样本量统计。"""
    out, tot_c, tot_ctrl, n_read = {}, 0.0, 0.0, 0
    with gzip.open(path, "rt") as f:
        for chunk in pd.read_csv(f, sep="\t", usecols=COLS, chunksize=400000):
            tot_c += pd.to_numeric(chunk["NCASES"], errors="coerce").fillna(0).sum()
            tot_ctrl += pd.to_numeric(chunk["NCONTROLS"], errors="coerce").fillna(0).sum()
            n_read += len(chunk)
            key = list(zip(chunk["chromosome"].astype(str), chunk["base_pair_location"].astype(str)))
            m = pd.Series(key).isin(targets).values
            if m.any():
                sub = chunk[m].copy()
                for _, r in sub.iterrows():
                    try:
                        out[(str(r["chromosome"]), str(r["base_pair_location"]))] = dict(
                            E=str(r["effect_allele"]).lower(),
                            O=str(r["other_allele"]).lower(),
                            OR=float(r["odds_ratio"]), cu=float(r["ci_upper"]),
                            cl=float(r["ci_lower"]), p=float(r["p_value"]))
                    except Exception:
                        continue
            if len(out) >= len(targets):
                break
    return out, tot_c, tot_ctrl, n_read


def orient_and_mr(inst, rows):
    be, se_e, bo, so, snps = [], [], [], [], []
    for rs, s in inst.items():
        k = (str(s["chr"]), str(s["bp"]))
        if k not in rows:
            continue
        r = rows[k]
        beta = math.log(r["OR"])
        if r["cu"] > 0 and r["cl"] > 0:
            se = (math.log(r["cu"]) - math.log(r["cl"])) / 3.92
        else:
            se = None
        if se is None or not np.isfinite(se) or se <= 0:
            continue
        if r["E"] == s["exp_allele"].lower():
            pass
        elif r["O"] == s["exp_allele"].lower():
            beta = -beta
        else:
            continue
        be.append(float(s["beta_e"])); se_e.append(float(s["se_e"]))
        bo.append(beta); so.append(se); snps.append(rs)
    res = mr13.run_mr(be, se_e, bo, so, snps=snps)
    return res, snps


def main():
    inst = load_inst()
    targets = {(str(s["chr"]), str(s["bp"])) for s in inst.values()}
    print(f"[inst] 31 工具变量；目标位点 = {len(targets)}")

    # ---------- 混合结局（校验） ----------
    print("\n=== 混合祖先结局 GCST90566800（校验解析+定向） ===")
    mix, mc, mctrl, mr = parse_outcome(MIXED, targets)
    print(f"  样本量 NCASES={mc:,.0f}  NCONTROLS={mctrl:,.0f}  合计={mc+mctrl:,.0f}  匹配 SNP={len(mix)}")
    res_mix, snps_mix = orient_and_mr(inst, mix)
    print(mr13.fmt(res_mix))
    if "ivw" in res_mix:
        print(f"  >>> OR/SD(混合) = exp({res_mix['ivw']:.4f}) = {math.exp(res_mix['ivw']):.4f}  "
              f"（正文 headline 1.0812，方向应为正）")

    # 逐 SNP 校验：与 harmonised CSV 的 beta_o 比对
    print("\n--- 逐 SNP 校验（复算混合 beta_o vs harmonised CSV） ---")
    h = pd.read_csv(HARM)
    hmap = {(str(r.SNP)): r for _, r in h.iterrows()}
    diffs = []
    for rs in snps_mix:
        if rs in hmap:
            # 复算该 SNP 的混合 log-OR（定向到 exp allele）
            k = (str(inst[rs]["chr"]), str(inst[rs]["bp"]))
            r = mix[k]
            b = math.log(r["OR"])
            if r["O"] == inst[rs]["exp_allele"].lower():
                b = -b
            diffs.append(b - float(hmap[rs].beta_o))
    if diffs:
        mad = max(abs(d) for d in diffs)
        print(f"  比对 {len(diffs)} 个 SNP；最大 |Δbeta_o| = {mad:.5f}  "
              f"（<1e-3 即解析/定向正确）")

    # ---------- 欧裔限定结局（F10） ----------
    print("\n=== 欧裔限定结局 GCST90566811（F10 执行） ===")
    eur, ec, ectrl, er = parse_outcome(EUR, targets)
    print(f"  样本量 NCASES={ec:,.0f}  NCONTROLS={ectrl:,.0f}  合计={ec+ectrl:,.0f}  匹配 SNP={len(eur)}")
    res_eur, snps_eur = orient_and_mr(inst, eur)
    print(mr13.fmt(res_eur))
    if "ivw" in res_eur:
        print(f"  >>> OR/SD(欧裔限定) = exp({res_eur['ivw']:.4f}) = {math.exp(res_eur['ivw']):.4f}")

    # ---------- 比较 ----------
    print("\n=== F10 敏感性结论 ===")
    if "ivw" in res_mix and "ivw" in res_eur:
        dm = math.exp(res_mix["ivw"]); de = math.exp(res_eur["ivw"])
        print(f"  混合 OR/SD = {dm:.4f}  vs  欧裔限定 OR/SD = {de:.4f}")
        print(f"  差异 = {abs(de-dm):.4f}（{100*abs(de-dm)/dm:.1f}%）→ "
              f"{'方向一致且量级相近，跨祖先稳健' if (res_mix['ivw']>0)==(res_eur['ivw']>0) and abs(de-dm)/dm<0.2 else '需复核'}")

    out = dict(
        mixed=dict(n_snp=len(snps_mix), **{k: v for k, v in res_mix.items() if k != "snps"}),
        eur=dict(n_snp=len(snps_eur), **{k: v for k, v in res_eur.items() if k != "snps"}),
        eur_sample=dict(NCASES=ec, NCONTROLS=ectrl, total=ec+ectrl),
        mixed_sample=dict(NCASES=mc, NCONTROLS=mctrl, total=mc+mctrl),
        per_snp_betao_max_abs_diff_vs_harmonised=float(mad) if diffs else None,
    )
    os.makedirs(os.path.join(ROOT, "output", "genetics"), exist_ok=True)
    with open(os.path.join(ROOT, "output", "genetics", "L7_F10_eur_sensitivity.json"), "w") as f:
        json.dump(out, f, indent=2, default=str)
    print("\n[saved] output/genetics/L7_F10_eur_sensitivity.json")


if __name__ == "__main__":
    main()
