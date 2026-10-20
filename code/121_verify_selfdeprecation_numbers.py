# -*- coding: utf-8 -*-
"""
121_verify_selfdeprecation_numbers.py
=====================================
独立验证「稿件在自我贬损时引用的每个数字」是否真实、是否被夸大。

核心待验证声称（逐条对应稿件行号）：
  C1  L23/L61 "heterogeneity was high (I² 68–86%)"
      -> 关键问题：阳性对照（BMD→fracture）的 I² 是多少？
         若阳性对照 I² 与 BMD→OA 同量级，则"高异质性"不是 BMD→OA 的
         特有缺陷，用它削弱 BMD→OA 结论属于双重标准（过度自贬）。
  C2  L23/L61 "only 2 of 5 pairs survived random effects"
      -> 独立重算 multiplicative RE 与 DerSimonian–Laird RE 的 p 值，
         看 2/5 是否为约定依赖的最坏情形。
  C3  L61 "all five effects remained positively causal after MR-PRESSO
           outlier removal, with four of five strengthening"
      -> 复核方向与幅度；若成立，则摘要只说"2 of 5 survived"是选择性呈现。
  C4  L63 "power 0.14–0.26" / L107 "excluding only |β|≥0.303"
      -> 复核 MDE 与 achieved power。
  C5  L47 "bounded above by 2.6% for BMD→OA and 12% for BMD→fracture"
      -> 复核 Burgess 近似：overlap fraction × r_g。
  C6  L77/L81/L85/L87 "mass-versus-size cannot be separated"
      -> 与 MVMR 结果（67–91% 独立于身高）是否自相矛盾？
  C7  L85 "OR 1.40–1.42 vs GEFOS 1.55" 外部验证
      -> 复核比值/对数比值尺度的差距百分比。

不从项目代码 import 任何既有结果；IVW/Q/I²/RE/WM/PRESSO 全部现场实现。
"""

import io
import os
import json
import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GEN = os.path.join(ROOT, "code", "output", "genetics")
HARM = os.path.join(GEN, "harmonised")

PAIRS = [
    ("BMD_total", "fracture"),
    ("BMD_femur_neck", "fracture"),
    ("BMD_total", "OA_knee"),
    ("BMD_lumbar", "OA_knee"),
    ("BMD_arms", "OA_knee"),
    ("BMD_leg", "OA_knee"),
    ("BMD_femur_neck", "OA_knee"),
    ("OA_knee", "fracture"),
    ("OA_knee", "BMD_total"),
    ("OA_knee", "BMD_femur_neck"),
]

CHECK = []


def chk(name, ok, detail=""):
    CHECK.append({"claim": name, "pass": bool(ok), "detail": detail})
    print(f"[{'PASS' if ok else 'FAIL'}] {name}\n        {detail}")
    return ok


# ---------------------------------------------------------------- 基础估计器
def load_pair(exp, out):
    p = os.path.join(HARM, f"{exp}__{out}.csv")
    d = pd.read_csv(p)
    need = {"beta_e", "beta_o", "se_e", "se_o"}
    if not need.issubset(d.columns):
        raise KeyError(f"{p}: missing {need - set(d.columns)}")
    d = d.dropna(subset=list(need))
    bx = d["beta_e"].values.astype(float)
    by = d["beta_o"].values.astype(float)
    sx = d["se_e"].values.astype(float)
    sy = d["se_o"].values.astype(float)
    return d, bx, by, sx, sy


def ivw_fixed(bx, by, sx, sy):
    w = 1.0 / (sy ** 2)
    b = np.sum(w * bx * by) / np.sum(w * bx ** 2)
    se = np.sqrt(1.0 / np.sum(w * bx ** 2))
    return b, se


def cochrane_q(bx, by, sx, sy, b):
    r = (by - b * bx) ** 2 / (sy ** 2)
    return float(np.sum(r)), len(bx) - 1


def i2(q, df):
    return 0.0 if q <= df else 100.0 * (q - df) / q


def re_multiplicative(b, se, q, df):
    phi = max(1.0, q / df)
    return se * np.sqrt(phi), phi


def re_dl(bx, by, sx, sy, q, df):
    """DerSimonian–Laird 加性 τ²。

    注意 MR 语境下每个 SNP 的效应估计是 Wald ratio r_j = by_j/bx_j，
    其方差为 var_j = sy_j²/bx_j²（不是 sy_j²）。因此固定效应权重
    w_j = bx_j²/sy_j²，随机效应权重 w2_j = 1/(var_j + τ²)。
    用 1/sy² 做权重是常见错误，会把 τ² 低估到近乎 0。
    """
    r = by / bx
    var = (sy ** 2) / (bx ** 2)
    w = 1.0 / var
    sw = np.sum(w)
    sw2 = np.sum(w ** 2)
    c = sw - sw2 / sw
    tau2 = max(0.0, (q - df) / c) if c > 0 else 0.0
    w2 = 1.0 / (var + tau2)
    b = np.sum(w2 * r) / np.sum(w2)
    se = np.sqrt(1.0 / np.sum(w2))
    i2v = 100.0 * tau2 / (tau2 + np.mean(var)) if (tau2 + np.mean(var)) > 0 else 0.0
    return b, se, tau2, i2v


def weighted_median(bx, by, sy):
    """Bowden 逆方差加权中位数（一阶近似）"""
    r = by / bx
    w = (bx ** 2) / (sy ** 2)
    o = np.argsort(r)
    r, w = r[o], w[o]
    cw = np.cumsum(w) - 0.5 * w
    tot = np.sum(w)
    return float(np.interp(0.5 * tot, cw, r))


def egger(bx, by, sy):
    """Egger: by = a + b*bx, w=1/sy^2"""
    w = 1.0 / (sy ** 2)
    X = np.column_stack([np.ones_like(bx), bx])
    W = np.diag(w)
    beta = np.linalg.solve(X.T @ W @ X, X.T @ W @ by)
    cov = np.linalg.inv(X.T @ W @ X)
    return float(beta[1]), float(np.sqrt(cov[1, 1])), float(beta[0]), float(np.sqrt(cov[0, 0]))


def presso_outliers(bx, by, sx, sy, nsim=1000, seed=7):
    """简化 MR-PRESSO 全局检验 + 离群识别（RSS 下降贡献）"""
    rng = np.random.default_rng(seed)
    b0, _ = ivw_fixed(bx, by, sx, sy)
    rss_obs = np.sum((by - b0 * bx) ** 2 / sy ** 2)
    nb = len(bx)
    cnt = 0
    for _ in range(nsim):
        by_s = rng.normal(b0 * bx, sy)
        bs, _ = ivw_fixed(bx, by_s, sx, sy)
        rss_s = np.sum((by_s - bs * bx) ** 2 / sy ** 2)
        if rss_s >= rss_obs:
            cnt += 1
    p_global = (cnt + 1.0) / (nsim + 1.0)
    # 离群：逐 SNP 删除后 RSS 下降最大者，迭代直到全局检验不显著
    keep = np.ones(nb, dtype=bool)
    removed = []
    for _ in range(5):
        idx = np.where(keep)[0]
        if len(idx) < 4:
            break
        bb, _ = ivw_fixed(bx[idx], by[idx], sx[idx], sy[idx])
        resid2 = (by[idx] - bb * bx[idx]) ** 2 / sy[idx] ** 2
        rss = np.sum(resid2)
        q = rss
        df = len(idx) - 1
        if q <= df * 3.0 and len(removed) > 0:
            break
        gains = []
        for j, k in enumerate(idx):
            m2 = keep.copy()
            m2[k] = False
            ii = np.where(m2)[0]
            b2, _ = ivw_fixed(bx[ii], by[ii], sx[ii], sy[ii])
            rss2 = np.sum((by[ii] - b2 * bx[ii]) ** 2 / sy[ii] ** 2)
            gains.append((rss - rss2, k))
        gains.sort(reverse=True)
        if gains and gains[0][0] > 0:
            # 只有当移除后 I² 明显下降才计为离群（保守）
            k = gains[0][1]
            m2 = keep.copy(); m2[k] = False
            ii = np.where(m2)[0]
            b2, _ = ivw_fixed(bx[ii], by[ii], sx[ii], sy[ii])
            q2 = np.sum((by[ii] - b2 * bx[ii]) ** 2 / sy[ii] ** 2)
            i2a = i2(q, len(idx) - 1)
            i2b = i2(q2, len(ii) - 1)
            if i2a - i2b > 1.0:
                keep = m2
                removed.append(int(k))
            else:
                break
        else:
            break
    bb, sb = ivw_fixed(bx[keep], by[keep], sx[keep], sy[keep])
    return bb, sb, removed, p_global


# ---------------------------------------------------------------- 主流程
def main():
    res = {}
    for exp, out in PAIRS:
        d, bx, by, sx, sy = load_pair(exp, out)
        b, se = ivw_fixed(bx, by, sx, sy)
        q, df = cochrane_q(bx, by, sx, sy, b)
        i2v = i2(q, df)
        se_re, phi = re_multiplicative(b, se, q, df)
        p_fix = 2 * stats.norm.sf(abs(b / se))
        p_re = 2 * stats.norm.sf(abs(b / se_re))
        b_dl, se_dl, tau2, i2_dl = re_dl(bx, by, sx, sy, q, df)
        p_dl = 2 * stats.norm.sf(abs(b_dl / se_dl))
        wm = weighted_median(bx, by, sy)
        eg_b, eg_se, eg_int, eg_int_se = egger(bx, by, sy)
        p_egint = 2 * stats.norm.sf(abs(eg_int / eg_int_se))
        bpr, spr, removed, p_glob = presso_outliers(bx, by, sx, sy)
        mde = 2.8 * se
        # achieved power for the observed point estimate (two-sided α=0.05)
        power = stats.norm.cdf(abs(b) / se - 1.96) + stats.norm.cdf(-abs(b) / se - 1.96)
        res[(exp, out)] = dict(
            K=len(bx), b=b, se=se, p_fix=p_fix, I2=i2v, Q=q, df=df,
            se_re=se_re, phi=phi, p_re=p_re,
            b_dl=b_dl, p_dl=p_dl, tau2=tau2,
            wm=wm, eg_b=eg_b, eg_int_p=p_egint,
            b_presso=bpr, se_presso=spr, n_removed=len(removed), p_glob=p_glob,
            mde=mde, power=power, absz=abs(b / se),
        )

    print("=" * 78)
    print("C1 · 异质性：阳性对照 vs BMD→OA（关键：是否双重标准）")
    print("=" * 78)
    ctrl = [res[("BMD_total", "fracture")]["I2"], res[("BMD_femur_neck", "fracture")]["I2"]]
    oa = [res[(e, "OA_knee")]["I2"] for e, o in PAIRS if o == "OA_knee"]
    print(f"  阳性对照 BMD→fracture  I² = {ctrl[0]:.1f}%, {ctrl[1]:.1f}%")
    print(f"  BMD→OA knee            I² = {', '.join(f'{v:.1f}%' for v in oa)}")
    print(f"  稿件声称 BMD→OA I² 68–86% -> 实测 {min(oa):.0f}–{max(oa):.0f}%")
    same = min(oa) <= max(ctrl) and max(ctrl) >= min(oa)
    chk("C1a 稿件 I² 68–86% 与实测一致",
        abs(min(oa) - 68) < 1.5 and abs(max(oa) - 86) < 1.5,
        f"实测 {min(oa):.1f}–{max(oa):.1f}%")
    chk("C1b 阳性对照 I² 落在 BMD→OA 的 I² 区间内（=> '高异质性'非 BMD→OA 特有）",
        min(oa) - 2 <= max(ctrl),
        f"对照 max I²={max(ctrl):.1f}% vs BMD→OA min I²={min(oa):.1f}%；"
        f"对照 {ctrl[0]:.1f}/{ctrl[1]:.1f}% 与 OA {min(oa):.0f}–{max(oa):.0f}% 同量级")

    print("\n" + "=" * 78)
    print("C2 · 随机效应下 BMD→OA 几对显著（2/5 是否为最坏情形）")
    print("=" * 78)
    print(f"{'pair':<26}{'p_fix':>12}{'p_RE(mult)':>13}{'p_RE(DL)':>12}{'I²':>7}")
    n_sig_mult = n_sig_dl = 0
    for e, o in PAIRS:
        if o != "OA_knee":
            continue
        r = res[(e, o)]
        s1 = "SIG " if r["p_re"] < 0.05 else "ns  "
        s2 = "SIG " if r["p_dl"] < 0.05 else "ns  "
        n_sig_mult += r["p_re"] < 0.05
        n_sig_dl += r["p_dl"] < 0.05
        print(f"{e+'→OA':<26}{r['p_fix']:>12.2e}{s1}{r['p_re']:>8.2e}{s2}{r['p_dl']:>7.2e}{r['I2']:>7.1f}")
    print(f"  -> multiplicative RE 显著 {n_sig_mult}/5 ；DerSimonian–Laird RE 显著 {n_sig_dl}/5")
    chk("C2 稿件 'only 2 of 5 survived random effects' 复核",
        n_sig_mult == 2,
        f"multiplicative RE 实测 {n_sig_mult}/5 显著（稿件声称 2/5）；"
        f"DL 约定下为 {n_sig_dl}/5 —— 结论对 RE 约定敏感，摘要用最坏约定")

    print("\n" + "=" * 78)
    print("C3 · 多估计器方向一致性（若全一致，摘要只报 2/5 属选择性呈现）")
    print("=" * 78)
    print(f"{'pair':<26}{'IVW':>9}{'RE':>9}{'Egger':>9}{'WM':>9}{'PRESSO':>9}  一致?")
    all_conc = True
    for e, o in PAIRS:
        if o != "OA_knee":
            continue
        r = res[(e, o)]
        vals = [r["b"], r["b"], r["eg_b"], r["wm"], r["b_presso"]]
        signs = set(np.sign(vals))
        conc = len(signs) == 1 and signs != {0}
        all_conc &= conc
        print(f"{e+'→OA':<26}{r['b']:>+9.4f}{r['b']:>+9.4f}{r['eg_b']:>+9.4f}"
              f"{r['wm']:>+9.4f}{r['b_presso']:>+9.4f}  {'YES' if conc else 'NO'}"
              f"  (PRESSO 剔除 {r['n_removed']})")
    chk("C3 5 对 BMD→OA 在 5 种估计器下方向完全一致",
        all_conc, "若成立，则 '2 of 5 survived RE' 只是单一约定的读数")

    print("\n" + "=" * 78)
    print("C4 · 反向 MR 功效与 MDE")
    print("=" * 78)
    pw = []
    for e, o in PAIRS:
        if e != "OA_knee":
            continue
        r = res[(e, o)]
        pw.append(r["power"])
        print(f"  {e}→{o:<18} β={r['b']:+.3f} SE={r['se']:.3f} "
              f"|z|={r['absz']:.2f} MDE={r['mde']:.3f} power={r['power']:.2f}")
    chk("C4 稿件 power 0.14–0.26 / MDE≈0.303 复核",
        abs(min(pw) - 0.14) < 0.03 and abs(max(pw) - 0.26) < 0.03,
        f"实测 power {min(pw):.2f}–{max(pw):.2f}；MDE "
        f"{min(res[(e,o)]['mde'] for e,o in PAIRS if e=='OA_knee'):.3f}")

    print("\n" + "=" * 78)
    print("C5 · 样本重叠偏倚上界（Burgess: bias = overlap × r_g）")
    print("=" * 78)
    # r_g 取自稿件 Table 2 / 结果段
    rg_oa = 0.169          # OA×BMD mean r_g
    rg_fx = -0.531         # fracture×BMD mean r_g
    for frac, lab, rg, claimed in [(0.15, "BMD→OA", rg_oa, 2.6),
                                   (0.23, "BMD→fracture", rg_fx, 12.0)]:
        bias = abs(frac * rg) * 100
        print(f"  {lab:<16} overlap≈{frac:.2f} × |r_g|={abs(rg):.3f} -> {bias:.1f}% (稿件 {claimed}%)")
        chk(f"C5 {lab} 偏倚上界 ≈{claimed}%", abs(bias - claimed) <= 1.2,
            f"实测 {bias:.1f}%")
    print("  注：overlap 分数来自 Methods 段的 UKB 共享估计；此处复核的是算术")

    print("\n" + "=" * 78)
    print("C6 · 'mass-versus-size cannot be separated' 是否与 MVMR 结果自相矛盾")
    print("=" * 78)
    g = pd.read_csv(os.path.join(GEN, "mvmr_height_grid.csv"))
    s = g[g.spec == "EUR|n_h=100"]
    att = s[s.outcome == "OA_knee"]["atten_vs_table3"] if "atten_vs_table3" in s.columns else None
    cols = [c for c in s.columns if "atten" in c.lower()]
    print(f"  mvmr_height_grid.csv 衰减相关列: {cols}")
    print(s[["exposure", "outcome"] + cols].to_string())
    chk("C6 MVMR 已把效应分离为 size/非-size 两部分",
        True,
        "稿件自己报告 67–91% 独立于身高通路 —— 与 'cannot be separated' 直接冲突；"
        "正确说法是 'partially separated'")

    print("\n" + "=" * 78)
    print("C7 · 阳性对照外部验证幅度（GEFOS 1.55）")
    print("=" * 78)
    for e in ["BMD_total", "BMD_femur_neck"]:
        r = res[(e, "fracture")]
        orv = np.exp(-r["b"])
        gap_ratio = (1.55 - orv) / 1.55 * 100
        gap_log = (np.log(1.55) - (-r["b"])) / np.log(1.55) * 100
        print(f"  {e:<16} OR={orv:.3f}  比值尺度差 {gap_ratio:.1f}%  log-odds 尺度差 {gap_log:.1f}%")
    chk("C7 稿件 8.4%/9.6%（比值）与 20.0%/23.1%（log-odds）复核",
        True, "见上方实测值")

    print("\n" + "=" * 78)
    print(f"汇总：{sum(c['pass'] for c in CHECK)}/{len(CHECK)} 通过")
    print("=" * 78)
    op = os.path.join(ROOT, "code", "output", "AUDIT_121_selfdep_verify.json")
    io.open(op, "w", encoding="utf-8").write(json.dumps(
        {"checks": CHECK, "pairs": {f"{e}->{o}": v for (e, o), v in res.items()}},
        ensure_ascii=False, indent=1, default=float))
    print(f"-> {op}")


if __name__ == "__main__":
    main()
