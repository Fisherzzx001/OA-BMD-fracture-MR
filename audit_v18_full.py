# -*- coding: utf-8 -*-
"""
独立复盘审计 v18 —— 从零重算核心统计量
不依赖任何既有结论/中间结论，仅从源 CSV 重新实现算法后与稿件 Table 1-4 逐项比对。
"""
import os, json, math
import numpy as np
import pandas as pd
from scipy import stats

PROJ = os.path.dirname(os.path.abspath(__file__))
GEN = os.path.join(PROJ, "code", "output", "genetics")
HARM = os.path.join(GEN, "harmonised")
MD = os.path.join(PROJ, "docs", "08_投稿稿_OA_BMD_fracture.md")

md = open(MD, encoding="utf-8").read()

CHECKS = [0]
FAILS = []
WARNS = []


def ck(name, ok, detail=""):
    CHECKS[0] += 1
    if not ok:
        FAILS.append(f"{name}: {detail}")
    return ok


def close(name, a, b, tol=5e-3):
    """比对两个数值；相对容差"""
    CHECKS[0] += 1
    if a is None or b is None:
        FAILS.append(f"{name}: None value ({a} vs {b})")
        return False
    denom = max(abs(a), abs(b), 1e-12)
    ok = abs(a - b) / denom <= tol
    if not ok:
        FAILS.append(f"{name}: recomputed={a:.6g} reported={b:.6g} (reldiff={abs(a-b)/denom:.3%})")
    return ok


def fmt_ok(name, recomputed, reported_str, fmt):
    """稿件报告值是格式化字符串时的比对"""
    CHECKS[0] += 1
    s = format(recomputed, fmt)
    if s != reported_str:
        FAILS.append(f"{name}: manuscript='{reported_str}' recomputed={s}")
        return False
    return True


# ============================================================ 1. 独立重算 MR
PAIRS = [
    ("BMD_total", "fracture"), ("BMD_femur_neck", "fracture"),
    ("BMD_total", "OA_knee"), ("BMD_femur_neck", "OA_knee"),
    ("BMD_lumbar", "OA_knee"), ("BMD_arms", "OA_knee"), ("BMD_leg", "OA_knee"),
    ("OA_knee", "fracture"), ("OA_knee", "BMD_total"), ("OA_knee", "BMD_femur_neck"),
]


def ivw_fixed(x, y, se_y):
    w = 1.0 / se_y ** 2
    b = np.sum(w * x * y) / np.sum(w * x * x)
    se = math.sqrt(1.0 / np.sum(w * x * x))
    return float(b), float(se)


def flex_fmt(name, recomputed, reported_str, lo=2, hi=5, signed=False):
    """稿件报告值可能是任意小数位；任一精度匹配即通过（消除格式假阳性）"""
    CHECKS[0] += 1
    for nd in range(lo, hi + 1):
        f = f"+.{nd}f" if signed else f".{nd}f"
        if format(recomputed, f) == reported_str:
            return True
    # 无符号格式也试一下（稿件有时省略 +）
    for nd in range(lo, hi + 1):
        if format(recomputed, f".{nd}f") == reported_str.lstrip("+"):
            return True
    FAILS.append(f"{name}: manuscript='{reported_str}' recomputed={recomputed:.6f}")
    return False


def ivw_random(x, y, se_y, q, df):
    """随机效应：SE_re = SE_fe × sqrt(max(1, Q/df))（TwoSampleMR 标准做法）"""
    b, se_fe = ivw_fixed(x, y, se_y)
    mult = math.sqrt(max(1.0, q / df)) if df > 0 else 1.0
    se = se_fe * mult
    p = 2 * stats.norm.sf(abs(b / se))
    return float(b), float(se), p


def egger(x, y, se_y):
    """标准 MR-Egger (Bowden 2015)：WLS 权重 1/se_y²，SE 用 psi 残差方差估计"""
    w = 1.0 / se_y ** 2
    Sw = np.sum(w)
    xbar = np.sum(w * x) / Sw
    ybar = np.sum(w * y) / Sw
    Sxx = np.sum(w * (x - xbar) ** 2)
    b1 = np.sum(w * (x - xbar) * (y - ybar)) / Sxx
    b0 = ybar - b1 * xbar
    resid = y - (b0 + b1 * x)
    psi = np.sum(w * resid ** 2) / (len(x) - 2)
    se_b1 = math.sqrt(psi / Sxx)
    se_b0 = math.sqrt(psi / Sw + psi * xbar ** 2 / Sxx)
    p1 = 2 * stats.norm.sf(abs(b1 / se_b1))
    p0 = 2 * stats.norm.sf(abs(b0 / se_b0))
    return float(b0), float(b1), float(se_b1), float(se_b0), p1, p0


def wmedian(x, y, se_y):
    """Bowden 2016 weighted median — 中点百分位定义 p_j = 100*(cumsum(w_j) - 0.5*w_j)
    （MendelianRandomization::mr_median 第 3 步；本项目实现，无第 4 步外推）"""
    ratio = y / x
    se_r = se_y / np.abs(x)
    w = 1.0 / se_r ** 2
    o = np.argsort(ratio)
    r, ww = ratio[o], w[o]
    cum = np.cumsum(ww) - 0.5 * ww
    tot = np.sum(ww)
    idx = np.searchsorted(cum / tot, 0.5)
    idx = min(int(idx), len(r) - 1)
    return float(r[idx])


def wmedian_full(x, y, se_y):
    """Bowden 2016 完整实现 — 含第 4 步线性外推（标准 R 包行为）"""
    ratio = y / x
    se_r = se_y / np.abs(x)
    w = 1.0 / se_r ** 2
    o = np.argsort(ratio)
    r, ww = ratio[o], w[o]
    tot = np.sum(ww)
    s = np.cumsum(ww)
    k = int(np.searchsorted(s, 0.5 * tot)) - 1
    k = max(k, 0)
    k = min(k, len(r) - 2)
    sk, sk1 = s[k] / tot, s[k + 1] / tot
    return float(r[k] + (r[k + 1] - r[k]) * (0.5 - sk) / (sk1 - sk))


def loo(x, y, se_y):
    out = []
    for i in range(len(x)):
        m = np.ones(len(x), bool)
        m[i] = False
        b, _ = ivw_fixed(x[m], y[m], se_y[m])
        out.append(b)
    return np.array(out)


print("=" * 100)
print("PART 1 — 从 harmonised 原始数据独立重算 MR（10 对）")
print("=" * 100)

RE = {}
for e, o in PAIRS:
    d = pd.read_csv(os.path.join(HARM, f"{e}__{o}.csv"))
    x, y, se = d.beta_e.values, d.beta_o.values, d.se_o.values
    n = len(d)

    b, bse = ivw_fixed(x, y, se)
    z = b / bse
    p = 2 * stats.norm.sf(abs(z))
    resid = y - b * x
    Q = float(np.sum(resid ** 2 / se ** 2))
    qdf = n - 1
    Qp = float(stats.chi2.sf(Q, qdf))
    I2 = max(0.0, (Q - qdf) / Q) if Q > 0 else 0.0
    bre, bse_re, p_re = ivw_random(x, y, se, Q, qdf)
    b0, b1, se_b1, se_b0, p_egger, p_int = egger(x, y, se)
    wm = wmedian(x, y, se)
    wmi = wmedian_full(x, y, se)
    l = loo(x, y, se)

    # F 统计量
    F = float(np.mean((x / d.se_e.values) ** 2))
    Fmin = float(np.min((x / d.se_e.values) ** 2))

    RE[(e, o)] = dict(n=n, beta=b, se=bse, z=z, p=p, Q=Q, Qdf=qdf, Qp=Qp, I2=I2,
                      beta_re=bre, se_re=bse_re, p_re=p_re,
                      egger=b1, egger_se=se_b1, intercept=b0, intercept_se=se_b0,
                      p_int=p_int, wmed=wm, wmed_i=wmi,
                      loo_min=float(l.min()), loo_max=float(l.max()),
                      loo_flips=int(np.sum(np.sign(l) != np.sign(b))),
                      F=F, Fmin=Fmin,
                      n_presso=int(d.presso_outlier.sum()) if "presso_outlier" in d else -1)
    print(f"{e:16s} -> {o:14s} n={n:3d} beta={b:+.6f} se={bse:.6f} p={p:.3e} "
          f"I2={I2*100:5.1f}% wmed={wm:+.5f} F={F:.1f}")

# ---- 与 AUDIT_mr_recomputed.csv 比对（稿件 Table 3 的权威源）
A = pd.read_csv(os.path.join(GEN, "AUDIT_mr_recomputed.csv")).set_index(["exposure", "outcome"])
print("\n-- 与 AUDIT_mr_recomputed.csv（= 稿件 Table 3 源）逐项比对 --")
for (e, o), r in RE.items():
    a = A.loc[(e, o)]
    close(f"MR {e}->{o} n", float(r["n"]), float(a.K), 0)
    close(f"MR {e}->{o} beta", r["beta"], float(a.beta))
    close(f"MR {e}->{o} se", r["se"], float(a.se))
    close(f"MR {e}->{o} z", r["z"], float(a.z))
    close(f"MR {e}->{o} p", r["p"], float(a.p), 1e-6)
    close(f"MR {e}->{o} Q", r["Q"], float(a.Q))
    close(f"MR {e}->{o} I2", r["I2"], float(a.I2))
    close(f"MR {e}->{o} p_re", r["p_re"], float(a.p_re), 1e-5)
    close(f"MR {e}->{o} egger", r["egger"], float(a.egger_beta))
    close(f"MR {e}->{o} egger_int", r["intercept"], float(a.egger_int), 1e-4)
    close(f"MR {e}->{o} egger_int_p", r["p_int"], float(a.egger_int_p), 1e-3)
    close(f"MR {e}->{o} loo_min", r["loo_min"], float(a.loo_min))
    close(f"MR {e}->{o} loo_max", r["loo_max"], float(a.loo_max))
    close(f"MR {e}->{o} loo_flips", float(r["loo_flips"]), float(a.loo_flips), 0)
    # 加权中位数：两种定义至少一种匹配
    okw = (abs(r["wmed"] - float(a.wmed)) <= 1e-6 or abs(r["wmed_i"] - float(a.wmed)) <= 1e-6
           or abs(r["wmed"] - float(a.wmed)) / max(abs(float(a.wmed)), 1e-9) < 0.02)
    if not okw:
        WARNS.append(f"MR {e}->{o} wmed: recomputed {r['wmed']:+.6f} / interp {r['wmed_i']:+.6f} "
                     f"vs AUDIT {float(a.wmed):+.6f}")

# ---- 与稿件 Table 3 字符串逐一比对
print("\n-- 与稿件 Table 3 报告字符串比对 --")
T3 = [
    ("BMD_total", "fracture", 31, "-0.337", "0.022", "3.7×10⁻⁵¹", "2.4×10⁻¹²", "-0.337", "0.558", "-0.364", "78%", "0"),
    ("BMD_femur_neck", "fracture", 25, "-0.351", "0.030", "7.7×10⁻³¹", "2.3×10⁻⁶", "-0.343", "0.018", "-0.365", "83%", "0"),
    ("BMD_total", "OA_knee", 31, "+0.023", "0.004", "3.9×10⁻¹⁰", "6.2×10⁻³", "+0.023", "0.810", "+0.019", "81%", "0"),
    ("BMD_lumbar", "OA_knee", 21, "+0.028", "0.005", "7.7×10⁻⁹", "1.7×10⁻²", "+0.026", "0.295", "+0.037", "83%", "0"),
    ("BMD_arms", "OA_knee", 16, "+0.022", "0.005", "1.7×10⁻⁶", "7.0×10⁻²", "+0.022", "0.982", "+0.031", "86%", "0"),
    ("BMD_leg", "OA_knee", 24, "+0.018", "0.004", "2.3×10⁻⁵", "9.8×10⁻²", "+0.017", "0.815", "+0.025", "85%", "0"),
    ("BMD_femur_neck", "OA_knee", 25, "+0.017", "0.005", "6.4×10⁻⁴", "5.5×10⁻²", "+0.017", "0.875", "+0.015", "68%", "0"),
    ("OA_knee", "fracture", 69, "+0.092", "0.108", "0.395", "0.542", "+0.106", "0.476", "+0.118", "49%", "0"),
    ("OA_knee", "BMD_total", 69, "+0.106", "0.108", "0.326", "0.481", "+0.123", "0.376", "+0.112", "49%", "0"),
    ("OA_knee", "BMD_femur_neck", 69, "+0.142", "0.108", "0.190", "0.272", "+0.164", "0.169", "+0.113", "30%", "0"),
]


SUP = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")


def sci(p):
    """把 p 格式化成稿件风格 x×10⁻ⁿ（Unicode 上标）"""
    if p >= 0.1:
        return f"{p:.3f}"
    e = math.floor(math.log10(p))
    m = p / (10 ** e)
    return f"{m:.1f}×10⁻{str(abs(e)).translate(SUP)}"


for e, o, K, beta, se, p, pre, egg, pint, wm, i2, flips in T3:
    r = RE[(e, o)]
    flex_fmt(f"T3 {e}->{o} n", r["n"], str(K), 0, 2)
    flex_fmt(f"T3 {e}->{o} beta", r["beta"], beta, 2, 4, True)
    flex_fmt(f"T3 {e}->{o} se", r["se"], se, 2, 4)
    CHECKS[0] += 1
    if sci(r["p"]) != p:
        FAILS.append(f"T3 {e}->{o} p: manuscript='{p}' recomputed='{sci(r['p'])}' (raw {r['p']:.4e})")
    CHECKS[0] += 1
    if sci(r["p_re"]) != pre:
        FAILS.append(f"T3 {e}->{o} p_re: manuscript='{pre}' recomputed='{sci(r['p_re'])}' (raw {r['p_re']:.4e})")
    flex_fmt(f"T3 {e}->{o} egger", r["egger"], egg, 2, 4, True)
    flex_fmt(f"T3 {e}->{o} egger_int_p", r["p_int"], pint, 2, 4)
    flex_fmt(f"T3 {e}->{o} wmed", r["wmed"], wm, 2, 4, True)
    flex_fmt(f"T3 {e}->{o} I2", r["I2"] * 100, i2.rstrip("%") if i2.rstrip("%").isdigit() else i2[:-1], 0, 2)
    fmt_ok(f"T3 {e}->{o} loo_flips", r["loo_flips"], flips, "d")

# ============================================================ 2. Table 4 / power
print("\n" + "=" * 100)
print("PART 2 — 独立重算统计功效（Table 4）")
print("=" * 100)
PW = pd.read_csv(os.path.join(GEN, "mr_power.csv"))
for e, o in PAIRS:
    r = RE[(e, o)]
    mde = 2.80 * r["se"]
    power = float(stats.norm.sf(1.96 - abs(r["z"])) + stats.norm.cdf(-1.96 - abs(r["z"])))
    row = PW[(PW.exposure == e) & (PW.outcome == o)].iloc[0]
    close(f"PWR {e}->{o} MDE", mde, float(row.MDE_beta_80pct), 2e-2)
    close(f"PWR {e}->{o} power", power, float(row.power_observed), 2e-2)
    close(f"PWR {e}->{o} abs_z", abs(r["z"]), float(row.abs_z), 2e-2)
    RE[(e, o)]["mde"] = mde
    RE[(e, o)]["power"] = power

# ============================================================ 3. Table 1 / h2
print("\n" + "=" * 100)
print("PART 3 — 独立复算 h² 与责任尺度转换（Table 1）")
print("=" * 100)
H2 = pd.read_csv(os.path.join(GEN, "rg_h2_trueN.csv")).set_index("phenotype")
LIAB = pd.read_csv(os.path.join(GEN, "oa_h2_liability.csv")).set_index("trait")


def h2_liability(h2_obs, K, P=None):
    """Lee et al. 2011 Eq 23"""
    if P is None:
        P = K
    thr = stats.norm.ppf(1 - K)          # 阈值（上尾）
    z = stats.norm.pdf(thr)              # 阈值处正态密度
    return h2_obs * (K * (1 - K)) ** 2 / (z ** 2 * P * (1 - P))


T1 = {"OA_knee": ("0.0354", "0.0007", "1.069"), "OA_selfrep": ("0.0223", "0.0016", "1.035"),
      "fracture": ("0.0258", "0.0014", "0.960"), "BMD_total": ("0.391", "0.012", "1.025"),
      "BMD_femur_neck": ("0.257", "0.009", "1.027"), "BMD_femur_total": ("0.317", "0.011", "1.031"),
      "BMD_lumbar": ("0.350", "0.013", "1.014"), "BMD_spine": ("0.357", "0.014", "1.016"),
      "BMD_arms": ("0.325", "0.012", "1.009"), "BMD_leg": ("0.325", "0.012", "1.025"),
      "BMD_pelvis": ("0.358", "0.014", "1.027"), "BMD_rib": ("0.303", "0.011", "1.012"),
      "BMD_trunk": ("0.366", "0.013", "1.020")}
NS = pd.read_csv(os.path.join(GEN, "N_sensitivity_decision.csv"))
FR_H2 = float(NS[(NS.trait == "fracture") & (NS.source.str.startswith("264973"))].h2_obs.iloc[0])
FR_SE = float(NS[(NS.trait == "fracture") & (NS.source.str.startswith("264973"))].h2_se.iloc[0])
FR_INT = float(NS[(NS.trait == "fracture") & (NS.source.str.startswith("264973"))].intercept.iloc[0])

for t, (hs, ss, isx) in T1.items():
    if t == "fracture":
        flex_fmt(f"T1 {t} h2", FR_H2, hs, 3, 5)
        flex_fmt(f"T1 {t} se", FR_SE, ss, 3, 5)
        flex_fmt(f"T1 {t} intercept", FR_INT, isx, 3, 5)
        continue
    flex_fmt(f"T1 {t} h2", float(H2.loc[t, "h2"]), hs, 3, 5)
    flex_fmt(f"T1 {t} se", float(H2.loc[t, "h2_se"]), ss, 3, 5)
    flex_fmt(f"T1 {t} intercept", float(H2.loc[t, "intercept"]), isx, 3, 5)

# BMD h² 实际计算用的 N（rg_h2_trueN 的 N_used）与 Table 1 声称的 31,873 是否一致
N_USED = {p: float(H2.loc[p, "N_used"]) for p in H2.index}
print("rg_h2_trueN N_used:", {k: int(v) for k, v in N_USED.items()})
# 2026-10-16 第十一轮清理：该条长期悬挂为 WARN。实因为 h² 用的是 representative N=32,000
# （BMD 各位点真实 N=31,873），差 0.40%，远小于报告 SE。Table 1 的 "N used" 列已改标实际
# 使用的 32,000 并补 ‡ 脚注；门禁改为核验「列值 = 计算值」+「舍入已披露」，不再静默 WARN。
if abs(N_USED["BMD_total"] - 32000) > 1:
    WARNS.append(f"Table 1 标注 BMD N used=32,000，但 h² 实际按 N={int(N_USED['BMD_total'])} 计算")
else:
    ck("T1 BMD 'N used' 列 = 32,000（h² representative N，与 rg_h2_trueN 一致）",
       "| 32,000‡" in md, "Table 1 应把 h² 实际使用的 representative N=32,000 写入 N used 列")
    ck("T1 ‡ 脚注披露 31,873→32,000 的 0.40% 舍入（h² ∝ 1/N）",
       "0.40%" in md and "representative N of 32,000" in md,
       "应披露 per-site N=31,873 与 representative N=32,000 的 0.40% 差异")

for t, (rep, K) in {"OA_knee": ("0.094", 0.12), "OA_selfrep": ("0.073", None),
                    "fracture": ("0.062", None)}.items():
    row = LIAB.loc[t]
    k = row.K_primary
    calc = h2_liability(float(row.h2_obs), float(k))
    flex_fmt(f"T1 {t} h2_liab", calc, rep, 3, 5)
    close(f"T1 {t} h2_liab vs file", calc, float(row.h2_liab_primary), 1e-3)

# ============================================================ 4. Table 2 / r_g
print("\n" + "=" * 100)
print("PART 4 — 独立复算 r_g 矩阵与分块（Table 2）")
print("=" * 100)
RG = pd.read_csv(os.path.join(GEN, "rg_matrix_13.csv"), index_col=0)
RL = pd.read_csv(os.path.join(GEN, "rg_results_long_13.csv"))
print("matrix shape:", RG.shape, "| long rows:", len(RL))
ck("T2 78 pairs", len(RL) == 78, f"long={len(RL)}")
ck("T2 matrix symmetric", np.allclose(RG.values, RG.values.T, atol=1e-9),
   "asymmetry max=%.2g" % np.max(np.abs(RG.values - RG.values.T)))
ck("T2 diagonal = 1", np.allclose(np.diag(RG.values), 1.0, atol=1e-9),
   "diag range %.4f-%.4f" % (np.min(np.diag(RG.values)), np.max(np.diag(RG.values))))

# 矩阵与 long 表一致性
mism = 0
for _, r in RL.iterrows():
    a, b = r.trait1 if "trait1" in r else r.iloc[0], r.trait2 if "trait2" in r else r.iloc[1]
    v = r.rg
    if abs(RG.loc[a, b] - v) > 1e-9:
        mism += 1
ck("T2 matrix vs long consistency", mism == 0, f"{mism} mismatches")

BLK = pd.read_csv(os.path.join(GEN, "rg_blocks_summary.csv"))
print(BLK.to_string())
T2 = {"OA × OA": (1, "+0.838", "+0.838", "1/1"), "OA × BMD": (20, "+0.117", "+0.250", "20/20"),
      "Fracture × BMD": (10, "-0.596", "-0.446", "10/10"), "BMD × BMD": (45, "+0.689", "+0.960", "45/45"),
      "OA × Fracture": (2, "+0.148", "+0.211", "0/2")}
# blocks_summary 只含 3 块，用 long 表独立重算全部分块
OA = ["OA_knee", "OA_selfrep"]
BMD = [c for c in RG.columns if c.startswith("BMD")]
ck("T2 BMD site count = 10", len(BMD) == 10, f"found {len(BMD)}")


def block_rg(g1, g2, exclude_same=False):
    vals = []
    for i, a in enumerate(g1):
        for b in g2:
            if a == b:
                continue
            if exclude_same and b in g1 and a in g2 and RG.index.get_loc(a) > RG.index.get_loc(b):
                continue
            vals.append((a, b, RG.loc[a, b]))
    # 去重（对称）
    seen, out = set(), []
    for a, b, v in vals:
        key = tuple(sorted([a, b]))
        if key in seen:
            continue
        seen.add(key)
        out.append((a, b, v))
    return out


b_oo = block_rg(OA, OA)
b_ob = block_rg(OA, BMD)
b_fb = block_rg(["fracture"], BMD)
b_bb = block_rg(BMD, BMD)
b_of = block_rg(OA, ["fracture"])
print(f"blocks: OA×OA={len(b_oo)} OA×BMD={len(b_ob)} Frac×BMD={len(b_fb)} "
      f"BMD×BMD={len(b_bb)} OA×Frac={len(b_of)} total={len(b_oo)+len(b_ob)+len(b_fb)+len(b_bb)+len(b_of)}")
ck("T2 block counts sum=78",
   len(b_oo) + len(b_ob) + len(b_fb) + len(b_bb) + len(b_of) == 78,
   f"{len(b_oo)},{len(b_ob)},{len(b_fb)},{len(b_bb)},{len(b_of)}")
ck("T2 OA×OA n=1", len(b_oo) == 1, str(len(b_oo)))
ck("T2 OA×BMD n=20", len(b_ob) == 20, str(len(b_ob)))
ck("T2 Frac×BMD n=10", len(b_fb) == 10, str(len(b_fb)))
ck("T2 BMD×BMD n=45", len(b_bb) == 45, str(len(b_bb)))
ck("T2 OA×Frac n=2", len(b_of) == 2, str(len(b_of)))

for nm, b, lo, hi, mean in [("OA×OA", b_oo, "+0.838", "+0.838", "+0.838"),
                            ("OA×BMD", b_ob, "+0.117", "+0.250", "+0.169"),
                            ("Fracture×BMD", b_fb, "-0.596", "-0.446", "-0.531"),
                            ("BMD×BMD", b_bb, "+0.689", "+0.960", "+0.812"),
                            ("OA×Fracture", b_of, "+0.148", "+0.211", "+0.180")]:
    vs = [v for _, _, v in b]
    flex_fmt(f"T2 {nm} min", min(vs), lo, 2, 4, True)
    flex_fmt(f"T2 {nm} max", max(vs), hi, 2, 4, True)
    flex_fmt(f"T2 {nm} mean", float(np.mean(vs)), mean, 2, 4, True)

# ============================================================ 5. MRlap
print("\n" + "=" * 100)
print("PART 5 — 独立复算 MRlap 校正量（§2.8 表）")
print("=" * 100)
ML = pd.read_csv(os.path.join(PROJ, "code", "output", "mrlap_overlap_results.csv"))
T_ML = [("BMD_total", "fracture", "-0.122", "0.012", "-0.141", "0.014", "1.0×10⁻²⁴", "1.15", "1.7×10⁻¹³", 48, "-0.49"),
        ("BMD_femur_neck", "fracture", "-0.103", "0.023", "-0.123", "0.027", "5.2×10⁻⁶", "1.19", "6.7×10⁻⁴", 34, "-0.51"),
        ("BMD_total", "OA_knee", "+0.0241", "0.0063", "+0.0284", "0.0074", "1.4×10⁻⁴", "1.18", "2.0×10⁻⁴", 46, "+0.16"),
        ("BMD_femur_neck", "OA_knee", "+0.0222", "0.0090", "+0.0270", "0.0109", "0.013", "1.22", "0.013", 33, "+0.19"),
        ("BMD_lumbar", "OA_knee", "+0.0333", "0.0083", "+0.0406", "0.0103", "8.6×10⁻⁵", "1.22", "<10⁻¹⁰⁰", 33, "+0.15"),
        ("BMD_arms", "OA_knee", "+0.0281", "0.0090", "+0.0343", "0.0109", "1.6×10⁻³", "1.22", "7.2×10⁻⁴", 31, "+0.17"),
        ("BMD_leg", "OA_knee", "+0.0253", "0.0085", "+0.0302", "0.0101", "2.7×10⁻³", "1.20", "5.7×10⁻³", 35, "+0.20"),
        ("OA_knee", "fracture", "+0.0384", "0.0344", "+0.0398", "0.0418", "0.34", "1.03", "0.89", 156, "+0.12"),
        ("OA_knee", "BMD_total", "+0.2742", "0.1103", "+0.3435", "0.1359", "0.011", "1.25", "5.7×10⁻³", 156, "+0.16"),
        ("OA_knee", "BMD_femur_neck", "+0.2608", "0.1007", "+0.3210", "0.1237", "9.5×10⁻³", "1.23", "0.015", 156, "+0.19")]
for e, o, ob, os_, cb, cs_, cp, ratio, pd_, k, rg in T_ML:
    r = ML[(ML.exposure == e) & (ML.outcome == o)].iloc[0]
    flex_fmt(f"MRlap {e}->{o} obs", float(r.observed_effect), ob, 3, 5, True)
    flex_fmt(f"MRlap {e}->{o} obs_se", float(r.observed_effect_se), os_, 2, 5)
    flex_fmt(f"MRlap {e}->{o} corr", float(r.corrected_effect), cb, 3, 5, True)
    flex_fmt(f"MRlap {e}->{o} corr_se", float(r.corrected_effect_se), cs_, 2, 5)
    close(f"MRlap {e}->{o} ratio", abs(float(r.corrected_effect)) / abs(float(r.observed_effect)),
          float(ratio), 1e-2)
    flex_fmt(f"MRlap {e}->{o} IVs", int(r.m_IVs), str(k), 0, 2)
    flex_fmt(f"MRlap {e}->{o} rg", float(r.rg), rg, 2, 3, True)

pct = [(abs(float(r.corrected_effect)) / abs(float(r.observed_effect)) - 1) * 100 for _, r in ML.iterrows()]
print(f"MRlap pct change: min={min(pct):.2f}% median={np.median(pct):.2f}% max={max(pct):.2f}%")
close("MRlap median pct (legend 20.5%)", float(np.median(pct)), 20.5, 1e-2)
ck("MRlap range 3%-25%", min(pct) >= 3.0 and max(pct) <= 25.5, f"{min(pct):.2f}-{max(pct):.2f}")
ck("MRlap 9/10 difference sig", int((ML.p_difference < 0.05).sum()) == 9,
   f"{int((ML.p_difference<0.05).sum())}/10")
lam = max(abs(float(r.int_crosstrait)) for _, r in ML.iterrows())
print(f"MRlap |lambda| max = {lam:.6f} (legend claims <= 0.037)")
ck("MRlap |lambda| <= 0.037", lam <= 0.037, f"{lam:.6f}")
ck("MRlap IV range 31-156", ML.m_IVs.min() == 31 and ML.m_IVs.max() == 156,
   f"{ML.m_IVs.min()}-{ML.m_IVs.max()}")

# MRlap h2 vs 主 LDSC
MAP = {"BMD_total": "BMD_total", "BMD_femur_neck": "BMD_femur_neck", "BMD_lumbar": "BMD_lumbar",
       "BMD_arms": "BMD_arms", "BMD_leg": "BMD_leg", "OA_knee": "OA_knee", "fracture": "fracture"}
mlexp = {r.exposure: float(r.h2_exp) for _, r in ML.iterrows()}
mlout = {r.outcome: float(r.h2_out) for _, r in ML.iterrows()}
mlh2 = dict(mlexp)
mlh2.update(mlout)
prim = {}
for t in ["BMD_total", "BMD_femur_neck", "BMD_lumbar", "BMD_arms", "BMD_leg", "OA_knee", "fracture"]:
    if t == "fracture":
        ns = pd.read_csv(os.path.join(GEN, "N_sensitivity_decision.csv"))
        prim[t] = float(ns[(ns.trait == "fracture") & (ns.source.str.startswith("264973"))].h2_obs.iloc[0])
    else:
        prim[t] = float(H2.loc[t, "h2"])
xs = [prim[t] for t in mlh2 if t in prim]
ys = [mlh2[t] for t in mlh2 if t in prim]
rr = float(np.corrcoef(xs, ys)[0, 1])
ratios = [ys[i] / xs[i] for i in range(len(xs))]
print("MRlap h2 vs primary LDSC:")
for i, t in enumerate([k for k in mlh2 if k in prim]):
    print(f"   {t:16s} primary={xs[i]:.6f} mrlap={ys[i]:.6f} ratio={ratios[i]:.4f}")
print(f"   Pearson r = {rr:.6f}; mean ratio = {np.mean(ratios):.4f} "
      f"(systematic offset {100*(np.mean(ratios)-1):+.1f}%)")
ck("MRlap h2 pearson r = 0.999", round(rr, 3) == 0.999, f"r={rr:.6f}")

# ============================================================ 6. 稿件文本断言核验
print("\n" + "=" * 100)
print("PART 6 — 稿件正文/摘要数值断言核验（独立于表格）")
print("=" * 100)

# 6.1 摘要：BMD→OA "p_IVW < 1e-4"
pmax_oa = max(RE[(e, "OA_knee")]["p"] for e in ["BMD_total", "BMD_femur_neck", "BMD_lumbar", "BMD_arms", "BMD_leg"])
print(f"[摘要] BMD→OA 五对最大 p = {pmax_oa:.3e}  (摘要声称 p_IVW < 1×10⁻⁴)")
ck("摘要 p_IVW ≤ 6.4e-4 成立", pmax_oa <= 6.4e-4, f"实际最大 p={pmax_oa:.3e}")

# 6.2 §3.4：随机效应下 femoral-neck 是否仍显著
p_re_fn = RE[("BMD_femur_neck", "OA_knee")]["p_re"]
print(f"[§3.4] femoral-neck→OA 随机效应 p = {p_re_fn:.4f}  (正文称 'remained significant (all p<0.05)')")
n_re_sig = sum(1 for e in ["BMD_total","BMD_femur_neck","BMD_lumbar","BMD_arms","BMD_leg"] if RE[(e,"OA_knee")]["p_re"] < 0.05)
print(f"[§3.4] 随机效应下显著的对数 = {n_re_sig}/5  (修正后正文称 2/5)")
ck("§3.4 随机效应显著数=2", n_re_sig == 2, f"实际 {n_re_sig}/5")
for e in ["BMD_total", "BMD_lumbar", "BMD_arms", "BMD_leg"]:
    print(f"        {e:16s} p_re = {RE[(e,'OA_knee')]['p_re']:.4f}")

# 6.3 §3.6：移除离群后 β 报告值
PRESS = pd.read_csv(os.path.join(GEN, "mr_presso_results.csv"))
CLAIM36 = {("BMD_total", "fracture"): "-0.352", ("BMD_femur_neck", "fracture"): "-0.345",
           ("BMD_total", "OA_knee"): "+0.026", ("BMD_femur_neck", "OA_knee"): "+0.021",
           ("BMD_lumbar", "OA_knee"): "+0.033", ("BMD_arms", "OA_knee"): "+0.026",
           ("BMD_leg", "OA_knee"): "+0.015"}
print("\n[§3.6] 移除 PRESSO 离群后 β（按 presso 比值 × 严格集 β 重算）")
for (e, o), claim in CLAIM36.items():
    row = PRESS[(PRESS.exposure == e) & (PRESS.outcome == o)].iloc[0]
    after = RE[(e, o)]["beta"] * (row.b_ivw_out / row.b_ivw_full)
    s = format(after, "+.3f")
    ok = (s == claim)
    ck(f"§3.6 {e}->{o} post-outlier β", ok, f"稿件='{claim}' 重算={s} (raw {after:+.6f})")
    print(f"        {e:16s}->{o:12s} 稿件 {claim}  重算 {s}")
# null 对移除离群后变化（稿件未报告）
#
# FIX 2026-10-16 (trust-but-verify audit F3):
#   b0 came from the strict-set recomputation (RE[...]) while b1 came from
#   mr_presso_results.csv, and those two sources are NOT on the same numerical
#   scale for the OA_knee-exposure rows. The mismatched pair produced spurious
#   magnitudes AND reversed signs (e.g. -71.8% for a row from which zero
#   outliers had actually been removed).
#   Both quantities are now derived from the SAME harmonised SNP files on one
#   consistent scale, and the outlier count is cross-checked against the
#   MR-PRESSO result file.
for e, o in [("OA_knee", "fracture"), ("OA_knee", "BMD_total"), ("OA_knee", "BMD_femur_neck")]:
    _d = pd.read_csv(os.path.join(GEN, "harmonised", f"{e}__{o}.csv"))
    _d["_b"] = _d.beta_o / _d.beta_e
    _d["_s"] = (_d.se_o / _d.beta_e).abs()

    def _ivw(df):
        _w = 1.0 / df._s ** 2
        return float((_w * df._b).sum() / _w.sum())

    b0 = _ivw(_d)
    _n_out = int(_d.presso_outlier.astype(float).sum())
    b1 = _ivw(_d[_d.presso_outlier.astype(float) == 0])
    row = PRESS[(PRESS.exposure == e) & (PRESS.outcome == o)].iloc[0]
    ck(f"§3.6 {e}->{o} outlier count consistent",
       _n_out == int(row.n_outlier),
       f"harmonised={_n_out} presso_results={int(row.n_outlier)}")
    print(f"        {e:16s}->{o:12s} 全样本 β={b0:+.4f}  移除{_n_out}个离群后 β={b1:+.4f}"
          f"  (变化 {100*(b1/b0-1):+.1f}%)")
    if abs(b1 / b0 - 1) > 0.5:
        # 稿件 §3.4 已报告幅度（+0.171 / +0.200）；若不再含则恢复 WARN
        if ("total-body BMD to +0.171" not in md) or ("fracture to +0.200" not in md):
            WARNS.append(f"§3.6 {e}->{o}: 移除离群后 β 由 {b0:+.3f} 变为 {b1:+.3f} "
                         f"({100*(b1/b0-1):+.0f}%)，稿件未报告幅度")

# 6.4 §3.2：多重检验计数
RL2 = pd.read_csv(os.path.join(GEN, "rg_results_long_13.csv"))
n_fdr = int((RL2.q_fdr < 0.05).sum())
n_bonf = int((RL2.p_bonf < 0.05).sum())
print(f"\n[§3.2] 78 对中 FDR<0.05: {n_fdr}, Bonferroni<0.05: {n_bonf}  (正文称 76 / 75)")
ck("§3.2 FDR 76", n_fdr == 76, str(n_fdr))
ck("§3.2 Bonf 75", n_bonf == 75, str(n_bonf))
nb = RL2[RL2.p_bonf >= 0.05]
print("        Bonferroni 下不显著的对：")
for _, r in nb.iterrows():
    print(f"          {r.trait1} × {r.trait2}: rg={r.rg:+.3f} p={r.p:.2e} p_bonf={r.p_bonf:.4f}")
ck("§3.2 Bonferroni 不显著对数=3(修正后)", len(nb) == 3, f"实际 {len(nb)} 对")

# 6.5 §3.3 weighted median 与 Table 3
wm_strict = RE[("BMD_total", "fracture")]["wmed"]
print(f"\n[§3.3] total-body→fracture 加权中位数：严格集 {wm_strict:+.4f} "
      f"(稿件 §3.3 与 Table 3 均报 −0.364)")
ck("§3.3 weighted median 与 Table 3 一致", format(wm_strict, ".3f") == "-0.364",
   f"重算 {wm_strict:+.4f} 与 Table 3 −0.364 不一致")

# 6.6 over-claim 扫描
print("\n[全文] over-claim 扫描")
import re as _re
PATTERNS = [
    (r"does not raise", "确定性: does not raise"),
    (r"were null", "非显著当零效应: were null"),
    (r"no causal effect", "确定性: no causal effect"),
    (r"is not mediated", "确定性: is not mediated"),
]
for pat, label in PATTERNS:
    for m in _re.finditer(pat, md):
        ln = md[:m.start()].count("\n") + 1
        ctx = md[max(0, m.start() - 90):m.end() + 90].replace("\n", " ")
        print(f"        L{ln:4d} [{label}] ...{ctx}...")
        WARNS.append(f"L{ln} 残留断言 [{label}]")

# ============================================================ 6.7 附属文件 vs 稿件 一致性
# 本轮新增：此前门禁不校验附属文件，导致 cover letter / Table S11 落后稿件一整个版本而静默通过。
print("\n[\u00a76.7] \u9644\u5c5e\u6587\u4ef6 vs \u7a3f\u4ef6 \u6570\u503c\u4e00\u81f4\u6027")
import io as _io
_anc = {
    "cover letter": os.path.join(PROJ, "docs", "cover_letter_OAC.md"),
    "Table S11": os.path.join(PROJ, "docs", "Supplementary",
                              "01_predictions_and_falsification.md"),
}
_anc_txt = {}
for _k, _p in _anc.items():
    _anc_txt[_k] = _io.open(_p, encoding="utf-8").read() if os.path.exists(_p) else ""

_ck = lambda n, ok, d: ck(n, ok, d)

# --- G1/G2: MVMR \u8870\u51cf\u8303\u56f4\u5fc5\u987b\u4e0e\u7a3f\u4ef6\u540c\u6b65 ---
ck("\u00a76.7 G1 cover letter MVMR \u8303\u56f4\u4e0e\u7a3f\u4ef6\u4e00\u81f4",
   "9.5\u201333.2%" in _anc_txt.get("cover letter", "")
   and "67\u201391%" in _anc_txt.get("cover letter", ""),
   "cover letter \u5e94\u4e3a 9.5\u201333.2% / 67\u201391%")

ck("\u00a76.7 G2 Table S11 \u8870\u51cf\u8303\u56f4\u4e0e\u7a3f\u4ef6\u4e00\u81f4",
   "9.5\u201333.2%" in _anc_txt.get("Table S11", ""),
   "Table S11 P4 \u884c\u5e94\u4e3a 9.5\u201333.2%")

ck("\u00a76.7 G3 Table S11 \u9633\u6027\u5bf9\u7167\u4e0a\u754c\u4e0e\u7a3f\u4ef6\u4e00\u81f4",
   "\u22121.7% to \u22123.0%" in _anc_txt.get("Table S11", ""),
   "Table S11 \u5e94\u4e3a \u22121.7% to \u22123.0%")

# --- G4: \u65e7\u503c\u6b8b\u7559 ---
for _old, _lab in [("6\u201333%", "MVMR 6-33%"), ("67\u201394%", "MVMR 67-94%"),
                   ("22\u201337", "\u5de5\u5177\u6570 22-37"), ("118 for knee OA", "\u5de5\u5177\u6570 118")]:
    ck(f"\u00a76.7 G4 \u7a3f\u4ef6\u4e0d\u6b8b\u7559\u65e7\u503c {_lab}",
       _old not in md, f"\u627e\u5230\u6b8b\u7559: {_lab}")
    for _k, _t in _anc_txt.items():
        ck(f"\u00a76.7 G4 {_k} \u4e0d\u6b8b\u7559\u65e7\u503c {_lab}",
           _old not in _t, f"{_k} \u6b8b\u7559: {_lab}")

# --- G5: \u7a7a\u6d1e\u81ea\u8d2c\u8bcd ---
for _w in ["provisional", "cannot be fully excluded", "cannot decompose",
           "better described as absence of evidence"]:
    ck(f"\u00a76.7 G5 \u7a3f\u4ef6\u4e0d\u542b\u7a7a\u6d1e\u81ea\u8d2c \u201c{_w}\u201d",
       _w not in md, f"\u7a3f\u4ef6\u6b8b\u7559: {_w}")

# --- G6: Highlights \u4e0d\u5e94\u628a\u975e\u663e\u8457\u5f53\u96f6\u6548\u5e94 ---
_hl = md.split("## Highlights")[1].split("**Keywords**")[0] if "## Highlights" in md else ""
ck("\u00a76.7 G6 Highlights \u4e0d\u542b \u201cwere null\u201d\uff08\u975e\u663e\u8457\u5f53\u96f6\u6548\u5e94\uff09",
   "were null" not in _hl, "Highlights \u5e94\u7528 non-significant")

# --- G7: \u5b9e\u8d28\u62ab\u9732\u4e0d\u5f97\u88ab\u5220\uff08\u9632\u6b62\u4e0b\u8f6e\u6539\u52a8\u8d70\u5411\u53e6\u4e00\u6781\u7aef\uff09 ---
for _must in ["0.14\u20130.26", "I\u00b2 68\u201386%", "9.5\u201333.2%",
              "MR-Egger", "positive control"]:
    ck(f"\u00a76.7 G7 \u5b9e\u8d28\u62ab\u9732\u4fdd\u7559 \u201c{_must}\u201d",
       _must in md, f"\u7a3f\u4ef6\u7f3a\u5c11\u5fc5\u8981\u62ab\u9732: {_must}")

# --- G8: \u9633\u6027\u5bf9\u7167\u5f02\u8d28\u6027\u5bf9\u7167\u5df2\u8865\u5165\uff08\u672c\u8f6e\u6838\u5fc3\u53d1\u73b0\uff09 ---
ck("\u00a76.7 G8 \u5df2\u8865\u5165\u9633\u6027\u5bf9\u7167\u5f02\u8d28\u6027\u5bf9\u7167",
   "positive control was comparably heterogeneous" in md,
   "\u9633\u6027\u5bf9\u7167 I\u00b2 78/83% \u4e0e BMD\u2192OA 68\u201386% \u540c\u91cf\u7ea7\uff1b\u5fc5\u987b\u663e\u5f0f\u5199\u5165")


# --- G9: \u5168\u6587\u62ec\u53f7\u5e73\u8861\uff08\u6458\u8981\u7b2c\u4e00\u53e5\u66fe\u7f3a\u53f3\u62ec\u53f7\uff09 ---
_bad = []
for _ln, _line in enumerate(md.split("\n"), 1):
    if _line.startswith("|"):
        continue
    if _line.count("(") != _line.count(")"):
        _bad.append(f"L{_ln}")
ck("\u00a76.7 G9 \u5168\u6587\u62ec\u53f7\u5e73\u8861", not _bad,
   f"\u4e0d\u5e73\u8861\u884c: {_bad}")

# --- G10: \u8fc7\u5ea6\u65ad\u8a00\u7981\u7528\uff08\u7b2c\u5341\u8f6e\uff1aopposing-endpoints \u65b0\u610f\u951a\u5b9a + forming/size \u672a\u5206\u79bb\uff09---
for _w, _lab in [("mass-attributable", "\u6b8b\u5dee\u88ab\u5f53\u4f5c mass \u7684\u65ad\u8a00"),
                 ("explains the opposing", "set-point \u88ab\u8bf4\u6210\u89e3\u91ca\u5bf9\u7acb\u7aef\u70b9"),
                 ("independent of skeletal size", "height-independent \u88ab\u5077\u6362\u4e3a size-independent"),
                 ("explained by skeletal size", "\u8870\u51cf\u88ab\u8bf4\u6210 explained by skeletal size"),
                 ("mass effect intact", "\u9633\u6027\u5bf9\u7167\u88ab\u8bf4\u6210 mass effect intact"),
                 ("true mass effect survives", "\u5bf9\u7167\u88ab\u8bf4\u6210\u8bc1\u5b9e true mass effect")]:
    ck(f"\u00a76.7 G10 \u7a3f\u4ef6\u4e0d\u542b\u8fc7\u5ea6\u65ad\u8a00 \u201c{_w}\u201d\uff08{_lab}\uff09",
       _w not in md, f"\u7a3f\u4ef6\u6b8b\u7559: {_w}")
    for _k, _t in _anc_txt.items():
        ck(f"\u00a76.7 G10 {_k} \u4e0d\u542b\u8fc7\u5ea6\u65ad\u8a00 \u201c{_w}\u201d", _w not in _t, f"{_k} \u6b8b\u7559: {_w}")

# --- G11: \u4e25\u683c\u672f\u8bed\u4e0e forming/size \u672a\u5206\u79bb\u62ab\u9732\u5fc5\u987b\u5b58\u5728\uff08\u9632\u4e0b\u8f6e\u56de\u9000\u5230\u65e7\u53d9\u4e8b\uff09---
ck("\u00a76.7 G11 \u7a3f\u4ef6\u542b\u4e25\u683c\u672f\u8bed height-independent",
   "height-independent" in md, "\u5e94\u4f7f\u7528 height-independent \u800c\u975e size-independent")
ck("\u00a76.7 G11 \u7a3f\u4ef6\u62ab\u9732 height \u4e0d\u542b\u9aa8\u5bbd/\u51e0\u4f55\uff08bone width / geometry\uff09",
   "bone width" in md and "geometry" in md,
   "\u5e94\u62ab\u9732 height \u4e0d\u542b\u9aa8\u5bbd/\u51e0\u4f55\uff0c\u6b8b\u5dee\u672a\u9274\u522b\u4e3a mass")
ck("\u00a76.7 G11 \u7a3f\u4ef6\u542b size \u5f52\u56e0\u4e0b\u754c\u8868\u8ff0\uff08lower bound\uff09",
   "lower bound" in md, "9.5\u201333.2% \u5e94\u8868\u8ff0\u4e3a size \u5f52\u56e0\u4e0b\u754c")
ck("\u00a76.7 G11 Table S11 \u540c\u6b65 forming/size \u672a\u5206\u79bb\u8fb9\u754c",
   "bone width" in _anc_txt.get("Table S11", ""),
   "S11 \u6838\u5fc3\u547d\u9898/\u8fb9\u754c\u6bb5\u5e94\u540c\u6b65 height \u7c97\u7cd9\u4ee3\u7406\u62ab\u9732")

# --- G12–G18: 第十一轮整改项防回退（第九轮预审报告 P0/P1/P2 与逐条批评）---
# G12 首例声明必须带 hedge（R3#1）
ck("\u00a76.7 G12 \u9996\u4f8b\u58f0\u660e\u5df2\u8f6f\u5316\uff08among the first\uff09",
   "To our knowledge this is among the first" in md,
   "\u7a3d\u4ef6\u5e94\u4e3a 'To our knowledge this is among the first'")
ck("\u00a76.7 G12 \u4e0d\u6b8b\u7559\u7edd\u5bf9\u9996\u4f8b\u58f0\u660e",
   "This is the first direct MR test" not in md,
   "\u7edd\u5bf9\u9996\u4f8b\u58f0\u660e\u662f\u4f2a\u521b\u65b0\u89e6\u53d1\u70b9")

# G13 femoral-neck\u2192fracture \u7684 Egger intercept p=0.018 \u5df2\u5728\u6b63\u6587\u70b9\u660e\uff08R2#6\uff09
ck("\u00a76.7 G13 \u6b63\u6587\u70b9\u660e femoral-neck intercept p=0.018",
   "intercept was not null (p=0.018)" in md,
   "\u9633\u6027\u5bf9\u7167\u4e0d\u5e94\u53ea\u62a5 total-body \u7684 null intercept")

# G14 MR-Egger \u529f\u7387\u9650\u5b9a\uff08R2#5\uff09
ck("\u00a76.7 G14 MR-Egger \u622a\u8ddd\u68c0\u9a8c\u5df2\u9644\u529f\u7387\u9650\u5b9a",
   "limited power" in md and "16\u201331 instruments" in md,
   "intercept \u4e0d\u663e\u8457 \u2260 \u65e0\u591a\u6548\u6027\uff1b\u9700\u5199\u660e 16\u201331 \u5de5\u5177\u4e0b\u529f\u7387\u6709\u9650")

# G15 \u70b9\u4f4d\u7279\u5f02\u6027\u9650\u5b9a\u4e3a\u5df2\u5206\u6790\u4f4d\u70b9\uff08R2#4\uff09
ck("\u00a76.7 G15 \u70b9\u4f4d\u7ed3\u8bba\u9650\u5b9a\u4e3a analysed sites",
   "analysed sites" in md,
   "\u5e94\u4e3a 'all five analysed sites'\uff0c\u4e0d\u5f97\u88f8\u79f0 all five sites")
ck("\u00a76.7 G15 \u70b9\u4f4d\u9009\u62e9\u89c4\u5219\u5df2\u5199\u660e\u4e3a\u9884\u8bbe",
   "prespecified canonical DXA regions" in md and "a site enters MR only if" in md,
   "5/10 \u4f4d\u70b9\u4e0d\u5e94\u8bfb\u4e3a\u6570\u636e\u53ef\u7528\u6027\u9a71\u52a8")

# G16 P6/F1 \u524d\u5411\u5f15\u7528\u9884\u544a\uff08R4#1\uff09
ck("\u00a76.7 G16 P6/F1 \u9996\u6b21\u51fa\u73b0\u5904\u5df2\u52a0\u524d\u5411\u6307\u9488",
   "specified in Table S11 and adjudicated in" in md,
   "L84 \u9996\u6b21\u63d0 P6/F1 \u65f6\u5e94\u9884\u544a\u5b9a\u4e49\u4f4d\u7f6e")

# G17 [31] \u5df2\u66f4\u6362\u4e3a\u5173\u8282\u7279\u5f02\u6027\u5b9e\u8bc1\uff08R3#2/#3\uff09
ck("\u00a76.7 G17 [31] = Zhu 2025 \u5173\u8282\u7f6e\u6362 MR",
   "Clin Rheumatol" in md and "39738847" in md,
   "[31] \u5e94\u4e3a Zhu et al. 2025 Clin Rheumatol 44:823\u2013830")
ck("\u00a76.7 G17 [31] \u4e0d\u518d\u4e3a Tong 2025 \u8e1d\u8db3\u75be\u75c5",
   "Tong X" not in md, "Tong 2025 \u4e0e BMD\u2192OA \u5173\u8282\u7279\u5f02\u6027\u65e0\u5173")

# G18 \u6458\u8981\u56e0\u679c\u8bed\u6c14\u5bf9\u9f50 pleiotropy-robust \u6846\u67b6\uff08Editor #2\uff09
ck("\u00a76.7 G18 \u6458\u8981\u4e0d\u4ee5\u56fa\u5b9a\u6548\u5e94 p \u9886\u5934",
   "fixed-effect p \u2264" in md and "estimator concordance rather than" in md,
   "\u6458\u8981\u5e94\u4ee5\u4f30\u8ba1\u91cf\u4e00\u81f4\u9886\u5934\uff0c p_IVW \u964d\u4e3a\u9644\u8bc1")

# ============================================================ 7. 输出
print("\n" + "=" * 100)
print(f"CHECKS = {CHECKS[0]}   FAILS = {len(FAILS)}   WARNS = {len(WARNS)}")
print("=" * 100)
for f in FAILS:
    print("  FAIL:", f)
for w in WARNS:
    print("  WARN:", w)

json.dump({f"{e}->{o}": {k: (float(v) if isinstance(v, (int, float, np.floating)) else v)
                          for k, v in r.items()} for (e, o), r in RE.items()},
          open(os.path.join(PROJ, "audit_v18_recomputed.json"), "w", encoding="utf-8"), indent=1)
print("\nwrote audit_v18_recomputed.json")
