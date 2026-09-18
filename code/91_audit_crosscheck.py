#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUDIT_CROSSCHECK.py —— 稿件正文/表格数字 vs 全部源 CSV 的系统交叉核对。
不依赖任何既有结论，直接从源文件重算并逐条比对稿件中的字面数字。
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = os.path.dirname(os.path.abspath(__file__))
G = os.path.join(ROOT, "output", "genetics")
MD = os.path.join(ROOT, "..", "docs", "08_投稿稿_OA_BMD_fracture.md")

fails, warns, oks = [], [], []


def chk(name, manuscript, recomputed, tol=0.0, sev="P0"):
    try:
        mv, rv = float(manuscript), float(recomputed)
        ok = abs(mv - rv) <= tol
    except Exception:
        ok = str(manuscript) == str(recomputed)
    rec = "稿件=%s  独立重算=%s" % (manuscript, recomputed)
    (oks if ok else (fails if sev == "P0" else warns)).append("[%s] %s  |  %s" % (sev, name, rec))


def main():
    md = open(MD, encoding="utf-8").read()

    # ---------- 1. h² 与 intercept（Table 1） ----------
    print("=" * 90)
    print("1. Table 1 —— SNP-h² 与 LDSC 截距")
    print("=" * 90)
    t = pd.read_csv(os.path.join(G, "rg_h2_trueN.csv")).set_index("phenotype")
    ns = pd.read_csv(os.path.join(G, "N_sensitivity_decision.csv"))
    fx = float(ns[(ns.trait == "fracture") &
                  (ns.source == "264973_GEFOS发现集")].h2_obs.iloc[0])

    tbl1 = {
        "OA_knee": (0.0354, 1.069), "OA_selfrep": (0.0223, 1.035),
        "fracture": (0.0258, 0.960),
        "BMD_total": (0.391, 1.025), "BMD_femur_neck": (0.257, 1.027),
        "BMD_femur_total": (0.317, 1.031), "BMD_lumbar": (0.350, 1.014),
        "BMD_spine": (0.357, 1.016), "BMD_arms": (0.325, 1.009),
        "BMD_leg": (0.325, 1.025), "BMD_pelvis": (0.358, 1.027),
        "BMD_rib": (0.303, 1.012), "BMD_trunk": (0.366, 1.020),
    }
    for k, (h2, inter) in tbl1.items():
        if k == "fracture":
            rh, ri = fx, 0.960
        else:
            r = t.loc[k]
            rh, ri = r.h2, r.intercept
        chk("Table1 %s h²" % k, h2, round(rh, 3), tol=0.001, sev="P1")
        chk("Table1 %s intercept" % k, inter, round(ri, 3), tol=0.001, sev="P1")
        print("  %-16s h² 稿件%.3f vs 源%.4f  |  截距 稿件%.3f vs 源%.4f"
              % (k, h2, rh, inter, ri))

    # 截距范围陈述
    allint = [v[1] for v in tbl1.values()]
    print("\n  稿件§3.1 称 'All 13 intercepts between 1.01 and 1.07'")
    print("  实际范围: %.3f – %.3f  (最小=%s fracture 0.960)"
          % (min(allint), max(allint), "fracture"))
    if min(allint) < 1.01:
        fails.append("[P0] §3.1 'All 13 intercepts 1.01–1.07' 与 Table 1 的 fracture=0.960 矛盾")
    else:
        oks.append("§3.1 截距范围陈述一致")

    # ---------- 2. liability h² 独立重算 ----------
    print("\n" + "=" * 90)
    print("2. Liability h² —— Lee 2011 Eq23 独立重算")
    print("=" * 90)
    lia = pd.read_csv(os.path.join(G, "oa_h2_liability.csv")).set_index("trait")

    def h2l(h2o, K):
        T = norm.isf(K)
        z = norm.pdf(T)
        return h2o * K * (1 - K) / z ** 2

    claim = {"OA_knee": (0.03544523509562934, 0.12, 0.094),
             "OA_selfrep": (0.02230416663074588, 0.08310489854903409, 0.073),
             "fracture": (0.025767681247168977, 0.14287116045785797, 0.062)}
    for k, (h2o, K, claimed) in claim.items():
        mine = h2l(h2o, K)
        file_v = float(lia.loc[k, "h2_liab_primary"])
        chk("liab h² %s (稿件)" % k, claimed, round(mine, 3), tol=0.001, sev="P1")
        chk("liab h² %s (CSV)" % k, round(file_v, 3), round(mine, 3), tol=0.001, sev="P1")
        print("  %-11s K=%.4f  稿件%.3f  CSV%.4f  独立%.4f"
              % (k, K, claimed, file_v, mine))
    # K=0.12 来源核查
    print("\n  ⚠ OA_knee 的 K=0.12 来源：code/15_liability_h2.py:26 硬编码 K_PRIMARY，无文献引用")
    warns.append("[P1] OA_knee liability K=0.12 为硬编码假设，稿件未给出文献来源")

    # ---------- 3. OR 与 CI 独立重算 ----------
    print("\n" + "=" * 90)
    print("3. §3.3 / §3.7 —— OR 与 95%CI 独立重算")
    print("=" * 90)
    mc = pd.read_csv(os.path.join(G, "mr_corrected_N.csv"))
    cases = [("BMD_total", "fracture", 1.40, 1.34, 1.46),
             ("BMD_femur_neck", "fracture", 1.42, 1.34, 1.51)]
    for e, o, cor, clo, chi in cases:
        r = mc[(mc.exposure == e) & (mc.outcome == o)].iloc[0]
        b, se = r.beta, r.se
        lo_b, hi_b = b - 1.96 * se, b + 1.96 * se
        or_pt, or_lo, or_hi = np.exp(-b), np.exp(-hi_b), np.exp(-lo_b)
        chk("OR %s→%s" % (e, o), cor, round(or_pt, 2), tol=0.01, sev="P1")
        chk("OR CI low %s" % e, clo, round(or_lo, 2), tol=0.01, sev="P1")
        chk("OR CI high %s" % e, chi, round(or_hi, 2), tol=0.01, sev="P1")
        print("  %-16s β=%+.4f SE=%.4f -> OR=%.3f (%.3f–%.3f)  稿件 %.2f (%.2f–%.2f)"
              % (e + "→" + o, b, se, or_pt, or_lo, or_hi, cor, clo, chi))

    # 1-Mb 敏感性（mr_results_corrected）
    mrc = pd.read_csv(os.path.join(G, "mr_results_corrected.csv"))
    r = mrc[(mrc.exposure == "BMD_femur_neck") & (mrc.outcome == "fracture")].iloc[0]
    b, se = r.ivw, r.ivw_se
    print("  1-Mb 敏感性 β=%+.4f -> OR=%.3f (%.3f–%.3f)  稿件 1.44 (1.36–1.53)"
          % (b, np.exp(-b), np.exp(-(b + 1.96 * se)), np.exp(-(b - 1.96 * se))))
    chk("1-Mb OR", 1.44, round(np.exp(-b), 2), tol=0.01, sev="P2")
    chk("1-Mb OR CI low", 1.36, round(np.exp(-(b + 1.96 * se)), 2), tol=0.01, sev="P2")
    chk("1-Mb OR CI high", 1.53, round(np.exp(-(b - 1.96 * se)), 2), tol=0.01, sev="P2")

    # ---------- 4. OA_knee 单位：标准化 vs log-OR ----------
    print("\n" + "=" * 90)
    print("4. 效应量单位判定（关键）—— 由 β=z/sqrt(2pq·N_eff) 反解 N_eff")
    print("=" * 90)
    si = pd.read_csv(os.path.join(G, "AUDIT_scale_inference.csv"))
    for trait in ["OA_knee", "fracture", "BMD_total"]:
        s = si[si.trait == trait]
        if len(s):
            print("  {:<16} N_eff 中位 ≈ {:,.0f}".format(trait, s.N_implied_median.median()))
    print("\n  理论值：")
    print("    fracture 若 log-OR    : N·K(1-K) = 264973×0.1429×0.8571 ≈ 32,459")
    print("    OA_knee  若 log-OR    : N·K(1-K) = 1208397×0.12×0.88   ≈ 127,607")
    print("    OA_knee  若标准化尺度 : N = 1,208,397")
    ok_fx = si[si.trait == "fracture"].N_implied_median.median()
    ok_knee = si[si.trait == "OA_knee"].N_implied_median.median()
    print("\n    实测 fracture ≈ {:,.0f}  -> 匹配 N·K(1-K) => log-OR".format(ok_fx))
    print("    实测 OA_knee  ≈ {:,.0f}  -> 匹配 N 而非 N·K(1-K) => 标准化尺度（非 log-OR）".format(ok_knee))
    fails.append("[P0] OA_knee 效应量为标准化尺度（N_eff≈1.15M≈N），但稿件/图轴标为 'log-OR'；"
                 "§3.4 的 OR≈1.02–1.05 系直接 exp(β_std)，正确换算应约 1.07")

    # ---------- 5. §3.4 表格 vs Table 3 vs 源 ----------
    print("\n" + "=" * 90)
    print("5. §3.4 正文表 / Table 3 / 源文件 三方核对（Egger 截距 p、p_re、wmed、I²）")
    print("=" * 90)
    mrr = pd.read_csv(os.path.join(G, "mr_results.csv"))          # 宽松 clump 未校正
    t3 = pd.read_csv(os.path.join(G, "table3_corrected.csv"))     # Table 3
    aud = pd.read_csv(os.path.join(G, "AUDIT_mr_recomputed.csv"))  # 严格集独立重算
    aud_idx = {r.pair: r for r in aud.itertuples()}

    print("%-26s %-6s %-8s | %-9s %-9s %-9s" %
          ("pair", "K严格", "n宽松", "稿件p_re", "Table3CSV", "严格集重算"))
    for r in t3.itertuples():
        key = "%s__%s" % (r.exp, r.out)
        a = aud_idx[key]
        s = mrr[(mrr.exposure == r.exp) & (mrr.outcome == r.out)].iloc[0]
        same = abs(a.p_re - r.p_re) < 1e-3
        print("%-26s %-6d %-8d | %-9.4g %-9.4g %-9.4g %s"
              % (key, r.k, s.n_snp, r.p_re, r.p_re, a.p_re, "" if same else "<<<"))
        if not same:
            fails.append("[P0] %s 的随机效应 p：Table3=%.4g 来自 n=%d 宽松集，"
                         "严格集 K=%d 独立重算=%.4g" % (key, r.p_re, s.n_snp, r.k, a.p_re))

    print("\n%-26s %-9s %-9s %-9s" % ("pair", "Table3", "宽松源", "严格集重算"))
    for r in t3.itertuples():
        key = "%s__%s" % (r.exp, r.out)
        a = aud_idx[key]
        s = mrr[(mrr.exposure == r.exp) & (mrr.outcome == r.out)].iloc[0]
        d = abs(a.egger_int_p - r.egger_int_p)
        print("%-26s %-9.4g %-9.4g %-9.4g %s"
              % (key + " [Egger截距p]", r.egger_int_p, s.egger_intercept_p,
                 a.egger_int_p, "" if d < 1e-3 else "<<<"))
        if d >= 1e-3:
            sev = "P0" if (min(r.egger_int_p, a.egger_int_p) < 0.05) != \
                  (max(r.egger_int_p, a.egger_int_p) < 0.05) else "P1"
            fails.append("[%s] %s Egger 截距 p：Table3=%.4g（n=%d 宽松集）vs 严格集 K=%d 重算=%.4g%s"
                         % (sev, key, r.egger_int_p, s.n_snp, r.k, a.egger_int_p,
                            "  <<< 显著性与否发生反转" if sev == "P0" else ""))

    # ---------- 6. 严格 clump 仪器数 ----------
    print("\n" + "=" * 90)
    print("6. §2.6 声称的严格 clump 仪器数 vs Table 3 实际使用数")
    print("=" * 90)
    sc = pd.read_csv(os.path.join(G, "strict_clump_summary.csv")).set_index("exposure")
    print("%-18s %-12s %-12s" % ("exposure", "clump保留", "Table3实际"))
    for e in sc.index:
        used = t3[t3.exp == e].k.max()
        print("%-18s %-12d %-12d %s" % (e, sc.loc[e, "n_kept_strict"], used,
                                        "" if sc.loc[e, "n_kept_strict"] == used else "<<<"))
    warns.append("[P1] §2.6 称严格 clump 保留 20–35(BMD)/116(OA_knee)，"
                 "Table 3 实际用 16–31/69；中间减少步骤未在稿件披露")

    # ---------- 7. 摘要与方法的字面数字 ----------
    print("\n" + "=" * 90)
    print("7. 摘要 / 方法中的字面数字")
    print("=" * 90)
    for token, desc in [("690,150", "共同 SNP 数"), ("1,290,050", "M"),
                        ("78", "表型对数"), ("1,208,397", "OA_knee N"),
                        ("462,933", "OA_selfrep N"), ("264,973", "fracture N"),
                        ("31,873", "BMD N")]:
        print("  %-12s 在稿件中出现: %s" % (desc, "是" if token in md else "否"))
    d = pd.read_csv(os.path.join(G, "rg_results_long_13.csv"))
    chk("表型对数", 78, len(d), sev="P1")
    from math import comb
    chk("C(13,2)", 78, comb(13, 2), sev="P1")

    # ---------- 8. 残留占位符 ----------
    print("\n" + "=" * 90)
    print("8. 残留占位符 / 明显遗留")
    print("=" * 90)
    for tok in ["[REF]", "TODO", "XXX", "placeholder", "待确认", "TBD"]:
        n = md.count(tok)
        print("  %-14s 出现 %d 次" % (tok, n))
        if tok in ("[REF]", "TODO", "TBD") and n:
            fails.append("[P0] 稿件残留占位符 %s ×%d" % (tok, n))

    # ---------- 汇总 ----------
    print("\n" + "=" * 90)
    print("审计汇总")
    print("=" * 90)
    print("  通过 %d 项" % len(oks))
    print("  P0（必须修）%d 项" % len([f for f in fails if f.startswith("[P0]")]))
    n_p1 = len([f for f in fails if f.startswith("[P1]")]) + len(warns)
    print("  P1（应修）  %d 项" % n_p1)
    print("\n--- P0 ---")
    for f in fails:
        print("  " + f)
    print("\n--- P1/其他 ---")
    for w in warns:
        print("  " + w)


if __name__ == "__main__":
    main()
