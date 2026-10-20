#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
170_from_scratch_recompute.py — 从零独立重算（full-trust-but-verify）

设计原则（用户硬约束）：
  1. 不 import 项目内任何分析模块；不读取任何既有结果 CSV（mr_results*.csv, table3*.csv 等）。
  2. 输入唯一来源：code/output/genetics/harmonised/<exposure>__<outcome>.csv（SNP 级）。
  3. 所有统计量与闭合式/定义式逐一手写：IVW、SE、z、p、OR+CI、Cochran Q、I^2、tau^2、
     随机效应 IVW、MR-Egger（含截距）、加权中位数、留一法、PRESSO 式 outlier 检出。
  4. 输出 JSON，供与稿件/门禁逐条对账。

Date: 2026-10-20
"""
import csv, json, math, os, sys
from statistics import NormalDist, median

HERE = os.path.dirname(os.path.abspath(__file__))
HARM = os.path.join(HERE, "output", "genetics", "harmonised")
OUT = os.path.join(HERE, "output", "genetics", "AUDIT_170_from_scratch.json")

ND = NormalDist()


# ---------- 基础统计函数（手写，不依赖 scipy） ----------
def pnorm(z):
    return ND.cdf(z)


def p_two_sided(z):
    """two-sided p from z; 对 |z| 极大时用渐近式避免下溢"""
    a = abs(z)
    if a < 8.0:
        return 2.0 * (1.0 - pnorm(a))
    # Mills-ratio asymptotic: 1-Phi(a) ~ phi(a)/a * (1 - 1/a^2 + 3/a^4)
    logp = -0.5 * a * a - math.log(a) - 0.5 * math.log(2 * math.pi)
    corr = math.log(max(1.0 - 1.0 / a**2 + 3.0 / a**4, 1e-300))
    return math.exp(logp + corr) * 2.0


def chi2_sf(x, df):
    """survival function of chi-square. 用 Wilson-Hilferty / 级数；df 小且 x 大时用渐近。"""
    if df <= 0:
        return float("nan")
    # 对整数/半整数 df 用不完全 gamma 的级数/连分式（Numerical Recipes gammq）
    a = df / 2.0
    x2 = x / 2.0
    if x2 < a + 1.0:
        # 级数展开 P(a,x)
        s = 1.0 / a
        term = s
        for n in range(1, 1000):
            term *= x2 / (a + n)
            s += term
            if abs(term) < abs(s) * 1e-16:
                break
        return 1.0 - s * math.exp(-x2 + a * math.log(x2) - math.lgamma(a))
    else:
        # 连分式 Q(a,x) —— 自 Lanczos 版
        lbeta = math.lgamma(a)
        # Log of Q via continued fraction (Lentz)
        tiny = 1e-300
        b = x2 + 1.0 - a
        c = 1.0 / tiny
        d = 1.0 / b
        h = d
        for i in range(1, 1000):
            an = -i * (i - a)
            b += 2.0
            d = an * d + b
            if abs(d) < tiny:
                d = tiny
            c = b + an / c
            if abs(c) < tiny:
                c = tiny
            d = 1.0 / d
            de = d * c
            h *= de
            if abs(de - 1.0) < 1e-16:
                break
        return math.exp(-x2 + a * math.log(x2) - lbeta) * h


# ---------- 载入 ----------
def load(pair):
    p = os.path.join(HARM, f"{pair}.csv")
    with open(p, encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    out = []
    for r in rows:
        out.append(
            dict(
                snp=r["SNP"],
                be=float(r["beta_e"]),
                se=float(r["se_e"]),
                pe=float(r["p_e"]),
                bo=float(r["beta_o"]),
                so=float(r["se_o"]),
                po=float(r["p_o"]),
                eaf=float(r["eaf"]),
            )
        )
    return out


# ---------- IVW ----------
def ivw(d):
    """IVW fixed effect: beta = sum(be*bo/so^2)/sum(be^2/so^2); se = 1/sqrt(sum(be^2/so^2))"""
    num = 0.0
    den = 0.0
    for r in d:
        wgt = r["be"] * r["be"] / (r["so"] ** 2)
        num += r["be"] * r["bo"] / (r["so"] ** 2)
        den += wgt
    beta = num / den
    se = 1.0 / math.sqrt(den)
    return beta, se, den


def cochran_q(d, beta):
    """Cochran Q = sum( (bo - beta*be)^2 / so^2 )"""
    q = 0.0
    for r in d:
        resid = r["bo"] - beta * r["be"]
        q += resid * resid / (r["so"] ** 2)
    return q


def i2_and_tau2(q, k):
    """I^2 = max(0,(Q-(k-1))/Q); tau^2 (DerSimonian-Laird)"""
    df = k - 1
    i2 = max(0.0, (q - df) / q) * 100.0 if q > 0 else 0.0
    return i2, df


def random_effects(d, q, k):
    """DerSimonian-Laird random-effects IVW.
    tau2 = max(0, (Q-(k-1)) / (sum(w) - sum(w^2)/sum(w))), w_i = be_i^2/so_i^2
    then beta_RE = sum(w*_i * ratio_i)/sum(w*_i), w*_i = 1/(so_i^2/be_i^2 + tau2)
    """
    w = [r["be"] ** 2 / r["so"] ** 2 for r in d]
    sw = sum(w)
    sw2 = sum(x * x for x in w)
    denom = sw - sw2 / sw
    df = k - 1
    tau2 = max(0.0, (q - df) / denom) if denom > 0 else 0.0
    ws = [1.0 / (1.0 / wi + tau2) for wi in w]
    ratio = [r["bo"] / r["be"] for r in d]
    beta = sum(a * b for a, b in zip(ws, ratio)) / sum(ws)
    se = math.sqrt(1.0 / sum(ws))
    return beta, se, tau2


def egger(d):
    """MR-Egger regression: bo = b0 + b1 * be, weights 1/so^2 (weighted linear regression).
    Returns slope, se_slope, p_slope, intercept, se_int, p_int."""
    w = [1.0 / r["so"] ** 2 for r in d]
    xs = [r["be"] for r in d]
    ys = [r["bo"] for r in d]
    sw = sum(w)
    xbar = sum(a * b for a, b in zip(w, xs)) / sw
    ybar = sum(a * b for a, b in zip(w, ys)) / sw
    sxx = sum(a * (b - xbar) ** 2 for a, b in zip(w, xs))
    sxy = sum(a * (b - xbar) * (c - ybar) for a, b, c in zip(w, xs, ys))
    b1 = sxy / sxx
    b0 = ybar - b1 * xbar
    n = len(d)
    resid = [y - (b0 + b1 * x) for x, y in zip(xs, ys)]
    # weighted residual variance (chi2/df)
    chi2 = sum(wi * ri * ri for wi, ri in zip(w, resid))
    s2 = chi2 / (n - 2)
    se_b1 = math.sqrt(s2 / sxx)
    se_b0 = math.sqrt(s2 * (1.0 / sw + xbar * xbar / sxx))
    z1 = b1 / se_b1
    z0 = b0 / se_b0
    # t-test with n-2 df approximated by normal (n>=16 -> error small); use t via normal + note
    return b1, se_b1, p_two_sided(z1), b0, se_b0, p_two_sided(z0)


def weighted_median(d):
    """Bowden weighted median: weights (be/so)^2 as stated in manuscript footnote."""
    w = [ (r["be"] / r["so"]) ** 2 for r in d]
    ratio = [r["bo"] / r["be"] for r in d]
    order = sorted(range(len(d)), key=lambda i: ratio[i])
    rs = [ratio[i] for i in order]
    ws = [w[i] for i in order]
    tot = sum(ws)
    cum = 0.0
    for idx in range(len(rs)):
        cum += ws[idx]
        if cum >= tot / 2.0:
            return rs[idx]
    return rs[-1]


def loo(d):
    """leave-one-out: returns (min, max, n_sign_flips)"""
    full, _, _ = ivw(d)
    vals = []
    for i in range(len(d)):
        sub = d[:i] + d[i + 1:]
        b, _, _ = ivw(sub)
        vals.append((d[i]["snp"], b))
    bs = [v for _, v in vals]
    flips = sum(1 for b in bs if (b > 0) != (full > 0))
    return min(bs), max(bs), flips, vals


def presso_like(d):
    """MR-PRESSO-style global test: observed RSS vs expected; outlier = SNP whose removal
    reduces RSS most and whose residual|z|>3 (approximate, for cross-check only)."""
    beta, se, _ = ivw(d)
    rss = cochran_q(d, beta)
    # expected RSS ~ chi2(k-1)
    k = len(d)
    p_global = chi2_sf(rss, k - 1)
    outliers = []
    for i, r in enumerate(d):
        resid = r["bo"] - beta * r["be"]
        z = resid / r["so"]
        if abs(z) > 3.0:
            outliers.append((r["snp"], round(z, 2)))
    return rss, p_global, outliers


# ---------- 主流程 ----------
PAIRS = [
    ("BMD_total__fracture", "BMD total -> fracture", "positive_control", "lower"),
    ("BMD_femur_neck__fracture", "BMD femoral neck -> fracture", "positive_control", "lower"),
    ("BMD_total__OA_knee", "BMD total -> OA knee", "forward", "higher"),
    ("BMD_lumbar__OA_knee", "Lumbar spine -> OA knee", "forward", "higher"),
    ("BMD_arms__OA_knee", "Arms -> OA knee", "forward", "higher"),
    ("BMD_leg__OA_knee", "Legs -> OA knee", "forward", "higher"),
    ("BMD_femur_neck__OA_knee", "Femoral neck -> OA knee", "forward", "higher"),
    ("OA_knee__fracture", "OA knee -> fracture", "reverse", None),
    ("OA_knee__BMD_total", "OA knee -> BMD total", "reverse", None),
    ("OA_knee__BMD_femur_neck", "OA knee -> BMD femoral neck", "reverse", None),
]

results = {}
print(f"{'pair':34s} {'K':>3s} {'beta':>10s} {'se':>9s} {'z':>8s} {'p':>12s} {'OR':>7s} {'Q':>8s} {'I2%':>6s} {'p_RE':>10s}")
print("-" * 115)
for key, label, kind, ordir in PAIRS:
    d = load(key)
    k = len(d)
    beta, se, den = ivw(d)
    z = beta / se
    p = p_two_sided(z)
    q = cochran_q(d, beta)
    i2, df = i2_and_tau2(q, k)
    bre, sre, tau2 = random_effects(d, q, k)
    p_re = p_two_sided(bre / sre)
    b1, se1, p1, b0, se0, p0 = egger(d)
    wm = weighted_median(d)
    lmin, lmax, flips, _ = loo(d)
    rss, pg, outl = presso_like(d)
    if ordir is not None:
        ab = abs(beta)
        orv = math.exp(ab)
        lo, hi = math.exp(ab - 1.96 * se), math.exp(ab + 1.96 * se)
        orstr = f"{orv:.3f}({lo:.3f}-{hi:.3f})"
    else:
        orstr = "—"
    results[key] = dict(
        label=label, kind=kind, K=k,
        beta=beta, se=se, z=z, p=p,
        OR_str=orstr,
        Q=q, Q_df=df, Q_p=chi2_sf(q, df), I2=i2, tau2=tau2,
        beta_RE=bre, se_RE=sre, p_RE=p_re,
        egger_beta=b1, egger_se=se1, egger_p=p1,
        egger_intercept=b0, egger_intercept_se=se0, egger_intercept_p=p0,
        wmedian=wm,
        loo_min=lmin, loo_max=lmax, loo_sign_flips=flips,
        presso_global_p=pg, presso_outliers=outl,
        MDE=2.8 * se,
    )
    print(f"{label:34s} {k:3d} {beta:+10.5f} {se:9.5f} {z:8.3f} {p:12.3e} "
          f"{(orv if ordir else 0):7.3f} {q:8.2f} {i2:6.1f} {p_re:10.3e}")

with open(OUT, "w", encoding="utf-8") as f:
    json.dump(results, f, indent=2, ensure_ascii=False)
print(f"\nwritten -> {OUT}")
