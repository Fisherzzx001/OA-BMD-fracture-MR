# -*- coding: utf-8 -*-
"""
129 · Table 3 稳健估计量的可复现性核验 + weighted median bootstrap 95% CI
---------------------------------------------------------------------------
背景（第十二轮整改 P1-3）：正文明确声明因果读数 "rests on the concordance of
pleiotropy-robust estimators"，但 Table 3 只给 IVW 的 SE/p，MR-Egger 斜率与
weighted median 只给点估计、无不确定性。

本脚本做三件事：
  1. 从沉积的 harmonised 工具表（code/output/genetics/harmonised/*.csv）独立
     复算 IVW(fixed/random)、MR-Egger(slope/intercept)、Cochran Q、I²，
     与 AUDIT_mr_recomputed.csv（即 Table 3 的源）逐项比对 —— 证明沉积数据与
     已报告数值一致；
  2. 用 Bowden 2016 加权中位数（权重 w_j = (beta_e/se_o)² = 1/Var(ratio)，
     与 code/13_mr_analysis.py::_wmed、code/100::weighted_median 同一定义）
     复算 WM，并与 Table 3 已列的 WM 比对；
  3. 对 WM 做 10,000 次按工具重抽样的 bootstrap，给出 95% 百分位 CI。

输出：code/output/genetics/wm_recomputed_vs_tabulated.csv
       文档：code/output/genetics/wm_recomputed_report.txt
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import chi2, norm

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HARM = os.path.join(PROJ, "code", "output", "genetics", "harmonised")
AUDIT = os.path.join(PROJ, "code", "output", "genetics", "AUDIT_mr_recomputed.csv")
OUT = os.path.join(PROJ, "code", "output", "genetics", "wm_recomputed_vs_tabulated.csv")
TXT = os.path.join(PROJ, "code", "output", "genetics", "wm_recomputed_report.txt")

NBOOT = 10000
SEED = 20261017

LINES = []


def say(s=""):
    print(s)
    LINES.append(str(s))


# ---------------------------------------------------------------- estimators
def wm_bowden(be, so, bo):
    """Bowden 2016 加权中位数：ratio = bo/be，权重 w = (be/so)^2，取 50%
    累积权处的次序统计量（落在两点之间时取均值）。"""
    ratio = bo / be
    w = (be / so) ** 2
    o = np.argsort(ratio)
    r, w = ratio[o], w[o]
    cw = np.cumsum(w)
    cw = cw / cw[-1]
    idx = int(np.searchsorted(cw, 0.5))
    idx = min(idx, len(r) - 1)
    if idx == 0:
        return float(r[0])
    if abs(cw[idx - 1] - 0.5) < 1e-12:
        return float(0.5 * (r[idx - 1] + r[idx]))
    return float(r[idx])


def ivw(be, so, bo):
    w = 1.0 / so ** 2
    b = float(np.sum(w * be * bo) / np.sum(w * be ** 2))
    se = float(1.0 / np.sqrt(np.sum(w * be ** 2)))
    z = b / se
    p = 2 * norm.sf(abs(z))
    return b, se, z, p


def egger(be, so, bo):
    X = np.column_stack([np.ones_like(be), be])
    W = 1.0 / so ** 2
    A = np.linalg.inv((X.T * W) @ X)
    bh = A @ ((X.T * W) @ bo)
    resid = bo - X @ bh
    k = len(bo)
    disp = float(np.sum(W * resid ** 2) / (k - 2))
    cov = A * max(disp, 1.0) if disp > 1 else A
    se = np.sqrt(np.diag(cov))
    from scipy.stats import t as tdist
    p_slope = 2 * tdist.sf(abs(bh[1] / se[1]), k - 2)
    p_int = 2 * tdist.sf(abs(bh[0] / se[0]), k - 2)
    return float(bh[1]), float(se[1]), float(p_slope), float(bh[0]), float(se[0]), float(p_int)


def het(be, so, bo, b):
    w = 1.0 / so ** 2
    resid = bo - b * be
    Q = float(np.sum(w * resid ** 2))
    df = len(be) - 1
    I2 = max(0.0, (Q - df) / Q) * 100 if Q > 0 else 0.0
    return Q, df, float(chi2.sf(Q, df)), I2


# ---------------------------------------------------------------- main
def main():
    audit = pd.read_csv(AUDIT)
    rows = []
    rng_master = np.random.default_rng(SEED)
    say("=" * 96)
    say("129 · Table 3 稳健估计量可复现性核验（沉积 harmonised 工具表 → 独立复算）")
    say("=" * 96)

    for _, a in audit.iterrows():
        pair = a["pair"]
        p = os.path.join(HARM, pair + ".csv")
        if not os.path.exists(p):
            say(f"  [MISSING] {pair}")
            continue
        d = pd.read_csv(p)
        be = d["beta_e"].to_numpy(float)
        so = d["se_o"].to_numpy(float)
        bo = d["beta_o"].to_numpy(float)

        b, se, z, pval = ivw(be, so, bo)
        eb, ese, ep, ei, eise, eip = egger(be, so, bo)
        Q, qdf, qp, I2 = het(be, so, bo, b)
        wm = wm_bowden(be, so, bo)

        # bootstrap CI for the weighted median
        seeds = rng_master.integers(0, 2 ** 31 - 1, NBOOT)
        boots = np.empty(NBOOT)
        n = len(be)
        for i in range(NBOOT):
            idx = np.random.default_rng(seeds[i]).integers(0, n, n)
            boots[i] = wm_bowden(be[idx], so[idx], bo[idx])
        lo, hi = np.percentile(boots, [2.5, 97.5])

        rows.append(dict(
            pair=pair, n=len(be),
            ivw_recomp=b, ivw_tab=float(a["beta"]), d_ivw=b - float(a["beta"]),
            se_recomp=se, se_tab=float(a["se"]),
            egger_recomp=eb, egger_tab=float(a["egger_beta"]), d_egger=eb - float(a["egger_beta"]),
            egger_se=ese, egger_p=ep, egger_int_p_tab=float(a["egger_int_p"]),
            egger_int_p_recomp=eip,
            Q_recomp=Q, Q_tab=float(a["Q"]), I2_recomp=I2, I2_tab=float(a["I2"]),
            wm_recomp=wm, wm_tab=float(a["wmed"]), d_wm=wm - float(a["wmed"]),
            wm_lo=float(lo), wm_hi=float(hi),
        ))
        say(f"  {pair:<28} n={len(be):>3}  IVW {b:+.6f} (tab {float(a['beta']):+.6f}, "
            f"Δ={b-float(a['beta']):+.2e})")
        say(f"  {'':<28}          Egger {eb:+.6f} (tab {float(a['egger_beta']):+.6f}, "
            f"Δ={eb-float(a['egger_beta']):+.2e})   "
            f"I² {I2:5.1f} (tab {float(a['I2']):5.1f})")
        say(f"  {'':<28}          WM  {wm:+.6f} (tab {float(a['wmed']):+.6f}, "
            f"Δ={wm-float(a['wmed']):+.6f})  boot95% [{lo:+.6f}, {hi:+.6f}]")

    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False, encoding="utf-8")

    say()
    say("---- 结论性判据 ----")
    ok_ivw = np.allclose(out.ivw_recomp, out.ivw_tab, atol=5e-6)
    ok_egger = np.allclose(out.egger_recomp, out.egger_tab, atol=5e-6)
    # AUDIT 表中的 I² 以小数存放（0.808 = 80.8%），正文以百分数报告
    ok_i2 = np.allclose(out.I2_recomp / 100.0, out.I2_tab, atol=6e-3)
    ok_wm = np.allclose(out.wm_recomp, out.wm_tab, atol=5e-6)
    say(f"  IVW  全部可复现: {ok_ivw}")
    say(f"  Egger 全部可复现: {ok_egger}")
    say(f"  I²   全部可复现: {ok_i2}")
    say(f"  WM   全部可复现: {ok_wm}   （不一致对数 = {int((out.d_wm.abs() > 5e-6).sum())}/{len(out)}）")
    if not ok_wm:
        say("  → 已列表的 weighted median 与沉积工具表不一致；Table 3 采用复算值。")
    assert ok_ivw and ok_egger and ok_i2, "IVW/Egger/I² 与沉积数据不一致，先查数据管线"

    with open(TXT, "w", encoding="utf-8") as fh:
        fh.write("\n".join(LINES))
    say(f"\nwrote {OUT}")
    say(f"wrote {TXT}")


if __name__ == "__main__":
    main()
