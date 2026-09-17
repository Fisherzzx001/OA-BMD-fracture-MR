# -*- coding: utf-8 -*-
"""
verify_figures.py
=========================================================================
独立核验：确认每一张投稿图上的每个数值都真实来自底层统计结果。

原则：
  1. 不 import build_figures.py / build_figures_mr.py 的任何函数；
     所有统计量在本文件内用独立实现重新推导。
  2. 每个被画到图上的量（点、误差棒、斜率、显著性标记、坐标截断）
     都与源 CSV 逐项比对。
  3. 顺带检查"图是否把数据画全了"（坐标轴截断导致误差棒不可见）。

用法： python verify_figures.py
"""
import os, csv, math, json
import numpy as np
import pandas as pd

GEN = "C:/Users/Administrator/WorkBuddy/2026-09-08-22-17-40/OA_B_PROJECT/code/output/genetics"
HARM = os.path.join(GEN, "harmonised")
MS = "C:/Users/Administrator/WorkBuddy/2026-09-08-22-17-40/OA_B_PROJECT/docs/08_投稿稿_OA_BMD_fracture.md"

FAIL, WARN, OK = [], [], []


def chk(cond, msg, sev="FAIL"):
    (OK if cond else (FAIL if sev == "FAIL" else WARN)).append(msg)
    return cond


def phi(x):
    return math.exp(-0.5 * x * x) / math.sqrt(2 * math.pi)


def norm_sf(z):
    return 0.5 * math.erfc(z / math.sqrt(2))


def norm_ppf(p):
    """Acklam 近似，用于独立核对（避免依赖 scipy）。"""
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00]
    pl, ph = 0.02425, 1 - 0.02425
    if p < pl:
        q = math.sqrt(-2 * math.log(p))
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    if p > ph:
        q = math.sqrt(-2 * math.log(1 - p))
        return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    q = p - 0.5
    r = q * q
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)


def read_mat(fname):
    with open(os.path.join(GEN, fname), encoding="utf-8") as f:
        r = csv.reader(f)
        header = next(r)
        rows = [x for x in r if x]
    labels = [h for h in header if h]
    M = np.array([[float(v) if v not in ("", "NA", "nan") else np.nan for v in row[1:]] for row in rows])
    return labels, M


# ---------------------------------------------------------------- BH-FDR（独立实现）
def bh_fdr(pvals):
    """Benjamini-Hochberg q 值。独立实现：显式循环 + 从后向前取最小。"""
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    q = [0.0] * m
    running = 1.0
    for rank in range(m - 1, -1, -1):
        i = order[rank]
        val = pvals[i] * m / (rank + 1)
        running = min(running, val)
        q[i] = min(running, 1.0)
    return q


# ---------------------------------------------------------------- 统计（独立实现）
def ivw_indep(be, bo, so):
    """逆方差加权（过原点）。用加权最小二乘闭式解独立实现。"""
    w = 1.0 / np.asarray(so, float) ** 2
    be = np.asarray(be, float)
    bo = np.asarray(bo, float)
    b = float(np.sum(w * be * bo) / np.sum(w * be * be))
    se = float(math.sqrt(1.0 / np.sum(w * be * be)))
    return b, se


def egger_indep(be, bo, so):
    w = 1.0 / np.asarray(so, float) ** 2
    be = np.asarray(be, float); bo = np.asarray(bo, float)
    mx = np.sum(w * be) / np.sum(w)
    my = np.sum(w * bo) / np.sum(w)
    sl = np.sum(w * (be - mx) * (bo - my)) / np.sum(w * (be - mx) ** 2)
    return float(my - sl * mx), float(sl)


print("=" * 78)
print("图 / 数据 一致性独立核验")
print("=" * 78)

# ================================================================ FIGURE 1
print("\n[Figure 1] LDSC 遗传相关热图")
labels, RG = read_mat("rg_matrix_13.csv")
_, SE = read_mat("rg_se_matrix_13.csv")
_, PV = read_mat("rg_p_matrix_13.csv")
n = len(labels)
li = {l: i for i, l in enumerate(labels)}   # 表型名 -> 矩阵下标
print(f"  表型数 n={n}  唯一对 = C(n,2) = {n*(n-1)//2}")

# 矩阵基本性质
chk(np.allclose(RG, RG.T, atol=1e-9, equal_nan=True), "r_g 矩阵对称")
chk(np.allclose(np.nan_to_num(SE), np.nan_to_num(SE.T), atol=1e-9), "SE 矩阵对称")
dg = np.diag(RG)
chk(np.allclose(dg[~np.isnan(dg)], 1.0, atol=1e-9), f"对角线 r_g = 1（{np.nanmin(dg):.4f}~{np.nanmax(dg):.4f}）")

# BH-FDR 独立重算
pu, idx = [], []
for i in range(n):
    for j in range(i + 1, n):
        pu.append(PV[i, j]); idx.append((i, j))
q = bh_fdr(pu)
Q = np.ones((n, n))
for k, (i, j) in enumerate(idx):
    Q[i, j] = Q[j, i] = q[k]
print(f"  BH-FDR: 78 对中 q<0.05 共 {sum(1 for v in q if v < 0.05)} 对，"
      f"q<0.01 {sum(1 for v in q if v < 0.01)} 对，q<0.001 {sum(1 for v in q if v < 0.001)} 对")

# 图注声称 "78 pairs"
ms_txt = open(MS, encoding="utf-8").read()
chk(n * (n - 1) // 2 == 78 and "across 78 pairs" in ms_txt,
    "图注 'across 78 pairs' 与 C(13,2)=78 一致")

# ---------------------------------------------------------- 溯源：矩阵本身来自真实 LDSC
print("  溯源：r_g 矩阵 vs LDSC 长表输出 rg_results_long_13.csv")
lg = pd.read_csv(os.path.join(GEN, "rg_results_long_13.csv"))
chk(len(lg) == 78, f"长表 {len(lg)} 行 = 78 个唯一表型对")
bad_rg, bad_se = [], []
for r in lg.itertuples():
    i, j = li[r.trait1], li[r.trait2]
    if abs(RG[i, j] - r.rg) > 1e-9:
        bad_rg.append(f"{r.trait1}-{r.trait2}: 矩阵{RG[i,j]:+.6f} vs 长表{r.rg:+.6f}")
    if abs(SE[i, j] - r.se) > 1e-9:
        bad_se.append(f"{r.trait1}-{r.trait2}")
chk(not bad_rg, f"矩阵中每个 r_g 都能在 LDSC 长表输出中逐位找到（不一致 {bad_rg[:3]}）")
chk(not bad_se, f"矩阵中每个 SE 都能在 LDSC 长表输出中逐位找到（不一致 {bad_se[:3]}）")
# 稿件自检用的 q_fdr 列与本次独立 BH 重算比较
if "q_fdr" in lg.columns:
    lgq = {frozenset((r.trait1, r.trait2)): r.q_fdr for r in lg.itertuples()}
    dq = max(abs(lgq[frozenset((labels[i], labels[j]))] - Q[i, j]) for i, j in idx)
    chk(dq < 1e-9, f"独立 BH-FDR 重算与 LDSC 自带的 q_fdr 列一致（最大差 {dq:.2e}）")

# 色阶截断：TwoSlopeNorm(vmin=-0.65, vcenter=0, vmax=1.0)
off = [RG[i, j] for i in range(n) for j in range(n) if i != j and not np.isnan(RG[i, j])]
oob = [v for v in off if v < -0.65 or v > 1.0]
chk(not oob, f"所有 r_g 落在色阶 [-0.65, 1.0] 内（越界 {len(oob)} 个；实测范围 "
             f"{min(off):+.3f} ~ {max(off):+.3f}）", sev="WARN" if oob else "OK")
# 图注声称色阶对称 —— 实际不对称
if "Colour scale is symmetric about zero" in ms_txt:
    chk(False, "图注称 'Colour scale is symmetric about zero'，但绘图用 "
               "TwoSlopeNorm(vmin=-0.65, vcenter=0, vmax=1.0)，数值范围并不对称（-0.65 ~ +1.00）", sev="WARN")

# 分块线位置应与 ORDER 的分组一致
OA = [i for i, l in enumerate(labels) if l.startswith("OA_")]
FR = [i for i, l in enumerate(labels) if l == "fracture"]
BM = [i for i, l in enumerate(labels) if l.startswith("BMD_")]
chk(OA == [0, 1] and FR == [2] and BM == list(range(3, 13)),
    f"分块线（画在 1.5 与 2.5）与真实分组一致：OA{OA} / fracture{FR} / BMD{BM[0]}~{BM[-1]}")

# ================================================================ FIGURE 2
print("\n[Figure 2] 遗传相关森林图")
bmd = [l for l in labels if l.startswith("BMD_")]
li = {l: i for i, l in enumerate(labels)}


def pair(a, b):
    i, j = li[a], li[b]
    return RG[i, j], SE[i, j], Q[i, j]


# 端到端：读取作图脚本实际落盘的坐标轴范围（而非复算），逐条检查误差棒是否完整可见
lim_path = os.path.join(GEN, "figure_axis_limits.json")
LIM = json.load(open(lim_path, encoding="utf-8")) if os.path.exists(lim_path) else {}
chk(bool(LIM), "读到作图脚本实际使用的坐标轴范围 figure_axis_limits.json")

panel_specs = [
    ("a", "fracture", bmd, "neg", "all negative"),
    ("b", "OA_knee", bmd, "pos", "all positive"),
    ("c", "OA_selfrep", bmd, "pos", "all positive"),
]
for lab, row, cols, expect, claim in panel_specs:
    vals, clipped = [], []
    xlo, xhi = LIM.get(f"fig2_{lab}", (None, None))
    for c in cols:
        rg, se, qv = pair(row, c)
        vals.append(rg)
        lo, hi = rg - 1.96 * se, rg + 1.96 * se
        if xlo is not None and (lo < xlo - 1e-9 or hi > xhi + 1e-9):
            clipped.append(f"{c}: [{lo:+.3f}, {hi:+.3f}]")
    if expect == "neg":
        chk(all(v < 0 for v in vals), f"面板 {lab} 标题声称 '{claim}'，实测全部 r_g<0 "
                                      f"（范围 {min(vals):+.3f} ~ {max(vals):+.3f}）")
    else:
        chk(all(v > 0 for v in vals), f"面板 {lab} 标题声称 '{claim}'，实测全部 r_g>0 "
                                      f"（范围 {min(vals):+.3f} ~ {max(vals):+.3f}）")
    chk(not clipped,
        f"面板 {lab} 全部 {len(cols)} 条 95% CI 完整落在实绘轴 [{xlo:+.3f}, {xhi:+.3f}] 内（截断 {clipped}）")

# 面板 d
d_pairs = [("OA_knee", "OA_selfrep"), ("OA_knee", "fracture"), ("OA_selfrep", "fracture")]
xlo, xhi = LIM.get("fig2_d", (None, None))
clipped = []
for a, b in d_pairs:
    rg, se, qv = pair(a, b)
    lo, hi = rg - 1.96 * se, rg + 1.96 * se
    if xlo is not None and (lo < xlo - 1e-9 or hi > xhi + 1e-9):
        clipped.append(f"{a}–{b}: [{lo:+.3f}, {hi:+.3f}]")
chk(not clipped, f"面板 d 全部 3 条 95% CI 完整落在实绘轴 [{xlo:+.3f}, {xhi:+.3f}] 内（截断 {clipped}）")
rg_kk, se_kk, _ = pair("OA_knee", "OA_selfrep")
print(f"  OA(knee) vs OA(self-rep): r_g = {rg_kk:+.4f} ± {se_kk:.4f}")
chk(abs(rg_kk - 0.84) < 0.01, f"图注声称 r_g=+0.84，实测 {rg_kk:+.4f}")
for a, b in [("OA_knee", "fracture"), ("OA_selfrep", "fracture")]:
    rg, se, qv = pair(a, b)
    chk(qv >= 0.05, f"图注声称 {a}–{b} 不显著，实测 r_g={rg:+.4f}, BH q={qv:.3f}")

# ================================================================ FIGURE 3
print("\n[Figure 3] MR 森林图")
mr = pd.read_csv(os.path.join(GEN, "mr_corrected_N.csv"))
mr_idx = {(r.exposure, r.outcome): r for r in mr.itertuples()}
print(f"  mr_corrected_N.csv 共 {len(mr)} 对")

panels3 = [
    ("a", [("BMD_total", "fracture"), ("BMD_femur_neck", "fracture")]),
    ("b", [("BMD_total", "OA_knee"), ("BMD_femur_neck", "OA_knee"),
           ("BMD_lumbar", "OA_knee"), ("BMD_arms", "OA_knee"), ("BMD_leg", "OA_knee")]),
    ("c", [("OA_knee", "fracture"), ("OA_knee", "BMD_total"), ("OA_knee", "BMD_femur_neck")]),
]
plotted3 = []
for lab, items in panels3:
    xlo, xhi = LIM.get(f"fig3_{lab}", (None, None))
    clipped = []
    for e, o in items:
        r = mr_idx[(e, o)]
        b, se = float(r.beta), float(r.se)
        lo, hi = b - 1.96 * se, b + 1.96 * se
        plotted3.append((e, o))
        if xlo is not None and (lo < xlo - 1e-9 or hi > xhi + 1e-9):
            clipped.append(f"{e}→{o}: [{lo:+.4f}, {hi:+.4f}]")
    chk(not clipped,
        f"面板 {lab} 全部 {len(items)} 条 95% CI 完整落在实绘轴 [{xlo:+.3f}, {xhi:+.3f}] 内（截断 {clipped}）")
# 覆盖性：MR 表里的每一对是否都画出来了
allpairs = set(map(tuple, mr[["exposure", "outcome"]].values.tolist()))
chk(set(plotted3) == allpairs,
    f"Figure 3 覆盖 mr_corrected_N.csv 全部 {len(allpairs)} 对（图中 {len(set(plotted3))} 对）")

# 图注/轴标签里硬编码的 OR 与 CI，必须从 β/SE 可复现
print("  硬编码 OR 标签复算：")
for e, o, claim in [("BMD_total", "fracture", "OR 1.40 (1.34–1.46)"),
                    ("BMD_femur_neck", "fracture", "OR 1.42 (1.34–1.51)")]:
    r = mr_idx[(e, o)]
    b, se = float(r.beta), float(r.se)
    orl = math.exp(-b)                      # 每 1 SD 更低 BMD 的 OR
    lo, hi = math.exp(-(b + 1.96 * se)), math.exp(-(b - 1.96 * se))
    print(f"    {e}→{o}: OR={orl:.4f} ({lo:.4f}–{hi:.4f})   稿件标签：{claim}")
    chk(abs(orl - float(r.OR_per_SD_lower)) < 1e-9,
        f"  {e} 复算 OR={orl:.4f} 与 CSV 列 OR_per_SD_lower={float(r.OR_per_SD_lower):.4f} 一致")
    chk(abs(lo - 1.34) < 0.006 and abs(hi - (1.46 if e == "BMD_total" else 1.51)) < 0.006,
        f"  {e} 复算 CI ({lo:.2f}–{hi:.2f}) 与稿件标签 {claim} 一致")
# 面板 b "all positively causal"
sig_pos = all(float(mr_idx[(e, o)].beta) > 0 and float(mr_idx[(e, o)].p) < 0.05 for e, o in panels3[1][1])
chk(sig_pos, "图注称面板 b 五个 BMD 暴露对膝 OA 'all positively causal'，实测全部 β>0 且 p<0.05")
# 面板 c 全部不显著
chk(all(float(mr_idx[(e, o)].p) > 0.05 for e, o in panels3[2][1]),
    "图注称面板 c 为无效对，实测全部 p>0.05")

# ================================================================ FIGURE 4
print("\n[Figure 4] SNP 遗传力条形图")
h2t = pd.read_csv(os.path.join(GEN, "rg_h2_trueN.csv"))
h2map = {r.phenotype: (float(r.h2), float(r.h2_se)) for r in h2t.itertuples()}
nsd = pd.read_csv(os.path.join(GEN, "N_sensitivity_decision.csv"))
fr_row = nsd[(nsd.trait == "fracture") & (nsd.source == "264973_GEFOS发现集")].iloc[0]
h2map["fracture"] = (float(fr_row.h2_obs), float(fr_row.h2_se))

order4 = ["OA_knee", "OA_selfrep", "fracture", "BMD_total", "BMD_femur_neck", "BMD_femur_total",
          "BMD_lumbar", "BMD_spine", "BMD_arms", "BMD_leg", "BMD_pelvis", "BMD_rib", "BMD_trunk"]
miss = [l for l in order4 if l not in h2map]
chk(not miss, f"13 个表型 h² 全部有源数据（缺失 {miss}）")
xlo4, xhi4 = LIM.get("fig4", (None, None))
clipped4 = []
for l in order4:
    v, s = h2map[l]
    if xhi4 is not None and v + 1.96 * s > xhi4 + 1e-9:
        clipped4.append(f"{l}={v:.4f}+{1.96*s:.4f}")
chk(not clipped4, f"所有 h² 及 95%CI 在实绘 x 轴 [{xlo4}, {xhi4}] 内（越界 {clipped4}）")
bmd_h2 = [h2map[l][0] for l in order4 if l.startswith("BMD_")]
chk(round(min(bmd_h2), 2) == 0.26 and round(max(bmd_h2), 2) == 0.39,
    f"稿件称 BMD h² 范围 0.26–0.39，实测 {min(bmd_h2):.4f}–{max(bmd_h2):.4f}"
    f"（两位小数 {min(bmd_h2):.2f}–{max(bmd_h2):.2f}）")
for l, tgt in [("OA_knee", 0.035), ("OA_selfrep", 0.022), ("fracture", 0.026)]:
    chk(abs(h2map[l][0] - tgt) < 0.002, f"{l} h²_obs={h2map[l][0]:.4f} 与稿件 {tgt} 一致")

# 易感性尺度 h²（Lee 2011 Eq23）独立复算 —— 图中硬编码的菱形标记
print("  易感性尺度 h² 独立复算（Lee 2011 Eq 23）：")
liab_hard = {"OA_knee": 0.0935, "OA_selfrep": 0.0726, "fracture": 0.0620}
liab_K = {"OA_knee": 0.12, "OA_selfrep": 0.083, "fracture": 0.143}
for t, K in liab_K.items():
    h2o = h2map[t][0]
    T = norm_ppf(1 - K)
    z = phi(T)
    h2l = h2o * K * (1 - K) / (z * z)
    print(f"    {t:12s} K={K:.3f}  h²_obs={h2o:.5f} → h²_liab={h2l:.5f}   "
          f"图中硬编码 {liab_hard[t]:.4f}")
    chk(abs(h2l - liab_hard[t]) < 0.0015,
        f"  {t} 复算 h²_liab={h2l:.4f} 与图中菱形标记 {liab_hard[t]} 一致")

# 交叉核对：Tables/ 里随包分发的 oa_h2_liability.csv 是否与稿件/图一致
oa = pd.read_csv(os.path.join(GEN, "oa_h2_liability.csv"))
print("  随包分发的 oa_h2_liability.csv 交叉核对：")
for r in oa.itertuples():
    print(f"    {r.trait:12s} h2_obs={r.h2_obs:.5f} K={r.K_primary} h2_liab={r.h2_liab_primary:.5f}")
r_self = oa[oa.trait == "OA_selfrep"].iloc[0]
chk(abs(float(r_self.h2_liab_primary) - 0.073) < 0.002,
    f"oa_h2_liability.csv 的 OA_selfrep h²_liab={float(r_self.h2_liab_primary):.4f} "
    f"应与稿件 0.073 一致（现为旧值，且 h²_obs={float(r_self.h2_obs):.5f} 用的是已废弃的 N=350,000，K=0.25）")
h2_selfrep_auth = float(nsd[(nsd.trait == "OA_selfrep") &
                            (nsd.source == "462933_ukb-b-14486权威")].iloc[0]["h2_obs"])
chk(abs(float(r_self.h2_obs) - h2_selfrep_auth) < 1e-12,
    f"oa_h2_liability.csv 的 OA_selfrep h²_obs={float(r_self.h2_obs):.6f} 与权威源 "
    f"N_sensitivity_decision(N=462,933) 的 {h2_selfrep_auth:.6f} 一致")
chk(abs(float(r_self.h2_obs) - h2map["OA_selfrep"][0]) < 1e-4,
    f"...且与 rg_h2_trueN.csv（Fig 4 绘图源）{h2map['OA_selfrep'][0]:.6f} 一致到绘图精度"
    f"（差 {abs(float(r_self.h2_obs) - h2map['OA_selfrep'][0]):.2e}，图上按 3 位小数显示为 0.022）")
chk("fracture" in set(oa.trait), "oa_h2_liability.csv 应含 fracture 行（现缺失，断裂的 h²_liab=0.062 无源）")

# ================================================================ FIGURE 5
print("\n[Figure 5] 因果网络概念图（图上文字标签核对）")
rg_oa_bmd = [(pair(a, b)[0], f"{a}–{b}") for a in ["OA_knee", "OA_selfrep"] for b in bmd]
vals_ob = [v for v, _ in rg_oa_bmd]
chk(all(v > 0 for v in vals_ob), f"图上标注 'Genetic correlation r_g +0.12 to +0.25'："
                                 f"实测 OA–BMD r_g 全部为正，范围 {min(vals_ob):+.3f} ~ {max(vals_ob):+.3f}")
inrange = all(0.12 - 0.005 <= v <= 0.25 + 0.005 for v in vals_ob)
chk(inrange, f"图上标注区间 +0.12~+0.25 覆盖全部 OA–BMD r_g（越界 "
             f"{[k for v, k in rg_oa_bmd if not (0.115 <= v <= 0.255)]}）")
or1 = math.exp(-float(mr_idx[("BMD_total", "fracture")].beta))
or2 = math.exp(-float(mr_idx[("BMD_femur_neck", "fracture")].beta))
chk(1.40 <= or1 <= 1.421 and 1.40 <= or2 <= 1.421,
    f"图上标注 'Protective OR 1.40–1.42'：复算 {or1:.4f} / {or2:.4f}")
chk(float(mr_idx[("OA_knee", "fracture")].p) > 0.05,
    "图上标注 OA→Fracture 'Null (ns)'：实测 p=%.3f > 0.05" % float(mr_idx[("OA_knee", "fracture")].p))
chk(all(float(mr_idx[(e, "OA_knee")].beta) > 0 for e in
        ["BMD_total", "BMD_femur_neck", "BMD_lumbar", "BMD_arms", "BMD_leg"]),
    "图上标注 BMD→OA 'Causal ↑ risk'：五个 BMD 暴露 β 全部为正")

# ================================================================ FIG 6 / 7 / S1 / S2
print("\n[Figure 6 / 7 / S1 / S2] 逐 SNP 散点、漏斗、留一法")
harm_files = sorted(f for f in os.listdir(HARM) if f.endswith(".csv") and not f.startswith("_"))
print(f"  harmonised 目录逐 SNP 文件 {len(harm_files)} 个")
rows = []
for fn in harm_files:
    exp, out = fn[:-4].split("__")
    d = pd.read_csv(os.path.join(HARM, fn))
    d = d[np.isfinite(d.beta_e) & np.isfinite(d.beta_o)]
    d = d[(np.abs(d.beta_e) > 0) & (d.se_o > 0)].reset_index(drop=True)
    be, bo, so = d.beta_e.values, d.beta_o.values, d.se_o.values
    b, se = ivw_indep(be, bo, so)
    z = b / se
    p = 2 * norm_sf(abs(z))
    eg_int, eg_sl = egger_indep(be, bo, so)
    # 逐 SNP 比值列是否可由 β 推出
    rr = bo / be
    sr = so / np.abs(be)
    ok_ratio = np.allclose(rr, d.ratio.values, atol=1e-8)
    ok_sr = np.allclose(sr, d.se_ratio.values, atol=1e-8)
    # 留一法
    loo = np.array([ivw_indep(np.delete(be, i), np.delete(bo, i), np.delete(so, i))[0]
                    for i in range(len(d))])
    flips = int(np.sum(np.sign(loo) != np.sign(b)))
    rows.append(dict(exp=exp, out=out, K=len(d), b=b, se=se, p=p, eg=eg_sl, eg_int=eg_int,
                     ok_ratio=ok_ratio, ok_sr=ok_sr,
                     n_out=int(d.presso_outlier.sum()),
                     loo_lo=loo.min(), loo_hi=loo.max(), flips=flips,
                     crosses0=bool((loo.min() < 0 < loo.max()))))
res = pd.DataFrame(rows)

# 与已发表 MR 主结果逐位比对
print("\n  逐对：重跑 harmonisation 的 IVW  vs  mr_corrected_N.csv（稿件 Table 3 口径）")
print(f"  {'pair':38s} {'K':>4s} {'β(重跑)':>10s} {'β(稿件)':>10s} {'Δβ':>10s} {'p(重跑)':>11s} {'p(稿件)':>11s}")
maxdb = 0.0
for r in res.itertuples():
    key = (r.exp, r.out)
    if key not in mr_idx:
        chk(False, f"  {key} 不在 mr_corrected_N.csv 中"); continue
    m = mr_idx[key]
    db = abs(r.b - float(m.beta))
    maxdb = max(maxdb, db)
    dp = abs(r.p - float(m.p)) / max(float(m.p), 1e-300)
    print(f"  {r.exp+'→'+r.out:38s} {r.K:4d} {r.b:+10.6f} {float(m.beta):+10.6f} {db:10.2e} "
          f"{r.p:11.3e} {float(m.p):11.3e}")
    chk(db < 1e-9, f"  {r.exp}→{r.out}: 重跑 β 与稿件 Δ={db:.2e}")
    chk(r.K == int(m.K), f"  {r.exp}→{r.out}: 工具变量数 K={r.K} 与稿件 {int(m.K)} 一致")
    chk(dp < 1e-6, f"  {r.exp}→{r.out}: 重跑 p 与稿件相对差 {dp:.2e}")
    chk(r.ok_ratio and r.ok_sr, f"  {r.exp}→{r.out}: 逐 SNP ratio / se_ratio 可由 β_o/β_e、se_o/|β_e| 复现")
print(f"  最大 Δβ = {maxdb:.3e}")

# MR-PRESSO 离群 SNP 数应与 presso 主结果一致
presso = pd.read_csv(os.path.join(GEN, "mr_presso_results.csv"))
pi = {(r.exposure, r.outcome): r for r in presso.itertuples()}
for r in res.itertuples():
    m = pi[(r.exp, r.out)]
    chk(r.n_out == int(m.n_outlier),
        f"  {r.exp}→{r.out}: 图中标红离群点 {r.n_out} 个 = MR-PRESSO n_outlier {int(m.n_outlier)}")

# 漏斗图纵轴 1/SE 与散点范围
print("\n  漏斗图坐标完整性（散点是否落在 y 轴范围内）：")
for r in res.itertuples():
    d = pd.read_csv(os.path.join(HARM, f"{r.exp}__{r.out}.csv"))
    d = d[np.isfinite(d.beta_e) & np.isfinite(d.beta_o)]
    d = d[(np.abs(d.beta_e) > 0) & (d.se_o > 0)]
    prec = 1.0 / d.se_ratio.values
    rr = d.ratio.values
    lo, hi = np.nanpercentile(rr, [2, 98])
    span = max(hi - lo, 4 * abs(r.b) + 0.05)
    c = 0.5 * (lo + hi)
    xlim = (min(c - 0.75 * span, r.b - 0.3 * span), max(c + 0.75 * span, r.b + 0.3 * span))
    # y 轴：prec.max()*1.08 为上界，prec.min()*0.85 为下界 → 全部点必然在内
    # x 轴：按 2–98 分位截断，统计落在轴外的点
    nclip = int(np.sum((rr < xlim[0]) | (rr > xlim[1])))
    chk(nclip <= 0.05 * len(rr),
        f"  {r.exp}→{r.out}: 漏斗图 x 轴按 2–98 分位截断，{nclip}/{len(rr)} 个点落在轴外"
        f"（{'可接受' if nclip else '无截断'}）", sev="WARN" if nclip else "OK")

# 留一法：图注/正文的断言
print("\n  留一法稳定性：")
for r in res.itertuples():
    print(f"    {r.exp+'→'+r.out:38s} β={r.b:+.4f}  LOO范围 [{r.loo_lo:+.4f}, {r.loo_hi:+.4f}]  "
          f"符号翻转={r.flips}  跨零={'是' if r.crosses0 else '否'}")
    chk(r.flips == 0, f"  {r.exp}→{r.out}: 留一法无符号翻转（flips={r.flips}）")
sig = res[res.p < 0.05]
chk(bool((~sig.crosses0).all()),
    f"显著对（n={len(sig)}）的留一法区间均不跨零 —— 支持图注 'No leave-one-out estimate crosses the null line'")
nscross = res[res.p >= 0.05].crosses0.sum()
if nscross:
    chk(False, f"但 {nscross} 个非显著对的留一法区间跨零（本应如此）；"
               f"图注 S2 的 'No leave-one-out estimate crosses the null line' 若泛指全部 10 对则不成立，"
               f"需限定为显著对", sev="WARN")

# ================================================================ 汇总
print("\n" + "=" * 78)
print(f"通过 {len(OK)} 项")
if WARN:
    print(f"\n⚠ 警告 {len(WARN)} 项：")
    for w in WARN:
        print("   -", w)
if FAIL:
    print(f"\n✗ 不通过 {len(FAIL)} 项：")
    for f in FAIL:
        print("   -", f)
else:
    print("\n✓ 无致命不一致")
print("=" * 78)
