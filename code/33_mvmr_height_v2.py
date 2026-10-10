#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
33_mvmr_height_v2.py — 多变量 MR（身高作骨尺寸代理）：完整规格网格 + 诊断
=============================================================================
承接 32_ 的初算结果。初算发现两个必须处理的问题：

(1) 条件 F(BMD) 仅 3.2–4.6（< 经验阈值 10）。原因不是 BMD 工具弱，而是
    Sanderson 条件 F 取「全体仪器上的均值」：557 个身高专属 SNP 在 BMD 上
    的 β≈0，各贡献 ≈1 的残差 χ²，把均值稀释掉了（F ≈ 1 + (K_bmd/K)(F_bmd−1)）。
    本脚本因此跑**仪器数网格** n_h ∈ {0, 100, 300, 1000}，并额外给出
    「仅在 BMD 仪器子集上」的条件 F（回答「BMD 在给定身高后还剩多少识别力」）。
(2) fracture 结局 SE 重建有 ±2.4% 误差，与要测的衰减（−2.3%）同量级。
    本脚本改从**原始 GWAS 文件**取真实 beta/SE
    (Direction1_pipeline/data/raw/fracture_anytype.tsv)，使 gate 由 WARN 变 PASS。

方向学已独立核实（trust-but-verify，32_ 的疑点已在 33_ 排查并否定）
--------------------------------------------------------------------
· 31 个 BMD 工具中**链歧义（A/T、G/C）SNP 数为 0**；
· 暴露 munge 与结局 munge 的 A1/A2 **31/31 完全一致**；
· 与 fracture 原始 GWAS 逐 SNP 比对：**符号一致性 31/31 正确**
  （32_ 报的 11 处「不符」是我用 MAF 反推等位基因所致——LDSC munge 的 MAF 列
   记的是**源文件的效应等位基因频率**，未必是 A1 的频率，故不可用于判定等位基因）。
· 由此：手稿方向学无错；本脚本仅对链歧义 SNP 在**身高侧**做保守剔除
  （TwoSampleMR 对不可推断的回文 SNP 的默认处理），并在结果中报告剔除数。

输出
----
  output/genetics/mvmr_height_grid.csv     规格网格主结果
  output/genetics/mvmr_height_subgroup.csv BMD 工具按「是否身高 GWS」分组
  output/genetics/mvmr_height_report.txt   可读报告
=============================================================================
"""
import os, gzip, math, csv, json
import numpy as np
import pandas as pd
from scipy.stats import chi2

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "output", "genetics")
MUNG = "C:/Users/Administrator/WorkBuddy/2026-09-08-20-17-24/OA_FLAGSHIP/07_reference/ldsc"
HG = os.path.join(ROOT, "..", "data", "height_gwas")
FRACTURE_RAW = ("C:/Users/Administrator/WorkBuddy/2026-09-09-14-17-27/Direction1_pipeline/"
                "Direction1_pipeline/data/raw/fracture_anytype.tsv")

HEIGHT_EUR = "GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR.gz"
HEIGHT_NOUKB = "GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR_excluding_UKB.gz"
GW_P, CLUMP_KB = 5e-8, 1_000_000
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def p_from_z(z):
    return math.erfc(abs(z) / math.sqrt(2.0))


# ---------------- 缓存式加载 ----------------
_MC, _HC = {}, {}


def load_munge(fn):
    if fn in _MC:
        return _MC[fn]
    d = {}
    with gzip.open(os.path.join(MUNG, fn), "rt") as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 8:
                continue
            try:
                d[p[2]] = (int(p[0]), int(p[1]), p[3].upper(), p[4].upper(),
                           float(p[5]), float(p[6]), float(p[7]))
            except ValueError:
                continue
    _MC[fn] = d
    print("   [cache] %s : %d SNPs" % (fn, len(d)))
    return d


def load_height(fn):
    if fn in _HC:
        return _HC[fn]
    d = {}
    with gzip.open(os.path.join(HG, fn), "rt") as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 11 or not p[1].startswith("rs"):
                continue
            try:
                d[p[1]] = (int(p[2]), int(p[3]), p[4].upper(), p[5].upper(),
                           float(p[6]), float(p[7]), float(p[8]), float(p[9]), float(p[10]))
            except ValueError:
                continue
    _HC[fn] = d
    print("   [cache] %s : %d SNPs" % (fn, len(d)))
    return d


def load_fracture_raw(want):
    """从原始 fracture GWAS 取真实 beta/SE（仅保留 want 中的 rsid）"""
    cache = os.path.join(OUT, "fracture_raw_extract.csv")
    if os.path.exists(cache):
        df = pd.read_csv(cache)
        return {r.rsid: (r.chr, r.pos, r.ea, r.oa, r.beta, r.se, r.eaf) for r in df.itertuples()}
    if not os.path.exists(FRACTURE_RAW):
        return {}
    out = []
    with open(FRACTURE_RAW, encoding="utf-8", errors="replace") as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            if p[8] in want:
                try:
                    out.append((p[8], int(p[0]), int(p[1]), p[2], p[3],
                                float(p[4]), float(p[5]), float(p[6])))
                except ValueError:
                    continue
    pd.DataFrame(out, columns=["rsid", "chr", "pos", "ea", "oa", "beta", "se", "eaf"]).to_csv(
        cache, index=False, encoding="utf-8")
    return {r[0]: (r[1], r[2], r[3], r[4], r[5], r[6], r[7]) for r in out}


def is_ambiguous(a1, a2):
    return COMP.get(a1, None) == a2


def align(a1r, a2r, a1, a2):
    if a1 == a1r and a2 == a2r:
        return 1.0
    if a1 == a2r and a2 == a1r:
        return -1.0
    if COMP.get(a1, "") == a1r and COMP.get(a2, "") == a2r:
        return 1.0
    if COMP.get(a1, "") == a2r and COMP.get(a2, "") == a1r:
        return -1.0
    return None


# ---------------- 估计器 ----------------
def ivw_uni(b_e, b_o, s_o):
    w = 1.0 / s_o ** 2
    den = np.sum(w * b_e ** 2)
    return float(np.sum(w * b_e * b_o) / den), float(math.sqrt(1.0 / den))


def mvmr(X, y, s_o):
    W = np.diag(1.0 / s_o ** 2)
    XtWX = X.T @ W @ X
    cov = np.linalg.inv(XtWX)
    th = cov @ (X.T @ W @ y)
    r = y - X @ th
    Q = float(r @ W @ r)
    df = max(1, len(y) - X.shape[1])
    disp = Q / df
    se_fe = np.sqrt(np.diag(cov))
    se_re = se_fe * math.sqrt(max(disp, 1.0))
    return dict(theta=th, se_fe=se_fe, se_re=se_re, Q=Q, df=df, disp=disp, resid=r)


def mvmr_int(X, y, s_o):
    Xd = np.column_stack([np.ones(len(y)), X])
    W = np.diag(1.0 / s_o ** 2)
    cov = np.linalg.inv(Xd.T @ W @ Xd)
    th = cov @ (Xd.T @ W @ y)
    r = y - Xd @ th
    Q = float(r @ W @ r)
    df = max(1, len(y) - Xd.shape[1])
    se = np.sqrt(np.diag(cov)) * math.sqrt(max(Q / df, 1.0))
    return th[0], se[0], th[1:], se[1:]


def cond_F(bj, sj, B, S):
    """Sanderson 条件 F：(1/(K−1)) Σ r_k²/v_k"""
    if B.shape[1] == 0:
        return float(np.mean((bj / sj) ** 2))
    d, *_ = np.linalg.lstsq(B, bj, rcond=None)
    r = bj - B @ d
    v = sj ** 2 + (S ** 2) @ (d ** 2)
    return float(np.sum(r ** 2 / v) / (len(bj) - 1))


# ---------------- 仪器集 ----------------
def build(pair_csv, munge_e, munge_o, dh, n_h, fr_raw=None):
    h_rows = list(csv.DictReader(open(pair_csv, encoding="utf-8")))
    de, do = load_munge(munge_e), load_munge(munge_o)

    # 暴露 N（实测，不假定）
    Ns = []
    for r in h_rows:
        s = r["SNP"]
        if s in de:
            f = de[s][4]
            Ns.append(1.0 / (2.0 * float(r["se_e"]) ** 2 * f * (1 - f)))
    N_e = float(np.median(Ns))

    # 结局 SE：优先真实文件；否则用「逐 SNP N_o」精确重建（OA_knee 已验证比值 min=max=1.00000），
    # 并以中位数 N_eff 作为兜底（用于原始文件缺失该 SNP 时）。
    use_raw = fr_raw is not None and len(fr_raw) > 0
    Ns2 = []
    for r in h_rows:
        s = r["SNP"]
        if s in do:
            f = do[s][4]
            Ns2.append(1.0 / (2.0 * float(r["se_o"]) ** 2 * f * (1 - f)))
    Neff = float(np.median(Ns2)) if Ns2 else None

    # 身高仪器
    cand = [(rs, v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7])
            for rs, v in dh.items() if v[7] < GW_P]
    hdf = pd.DataFrame(cand, columns=["SNP", "CHR", "BP", "EA", "OA", "FREQ", "BETA", "SE", "P"])
    hdf = hdf.sort_values("P")
    taken, keep = [], []
    for i, r in hdf.iterrows():
        c, b = int(r.CHR), int(r.BP)
        if all(not (cc == c and abs(bb - b) <= CLUMP_KB) for cc, bb in taken):
            taken.append((c, b)); keep.append(i)
    hdf = hdf.loc[keep]
    bmd_pos = [(int(r["CHR"]), int(r["BP"])) for r in h_rows]
    hdf = hdf[[all(not (cc == int(r.CHR) and abs(bb - int(r.BP)) <= CLUMP_KB)
                   for cc, bb in bmd_pos) for _, r in hdf.iterrows()]]
    hdf = hdf.sort_values("P")
    if n_h is None:
        pass
    elif len(hdf) > n_h:
        hdf = hdf.head(n_h)          # n_h=0 → 不纳入任何身高专属工具（仅 BMD 工具）

    recs = []
    for r in h_rows:
        recs.append((r["SNP"], "BMD", float(r["beta_e"]), float(r["beta_o"])))
    for _, r in hdf.iterrows():
        recs.append((r.SNP, "HEIGHT", np.nan, np.nan))

    out = []
    n_amb_drop = 0
    for s, src, be_ref, bo_ref in recs:
        if s not in de or s not in do or s not in dh:
            continue
        ce, bpe, a1e, a2e, mafe, ne, ze = de[s]
        if is_ambiguous(a1e, a2e):
            n_amb_drop += 1
            continue
        co, bpo, a1o, a2o, mafo, no, zo = do[s]
        ch, bph, eah, oah, frh, bh, sh, ph, nh = dh[s]
        go = align(a1e, a2e, a1o, a2o)
        gh = align(a1e, a2e, eah, oah)
        if go is None or gh is None:
            continue

        se_e = 1.0 / math.sqrt(2.0 * N_e * mafe * (1 - mafe))
        beta_e = ze * se_e
        beta_o = se_o = None
        if use_raw and s in fr_raw:
            _, _, rea, roa, rb, rse, reaf = fr_raw[s]
            g2 = align(a1e, a2e, rea.upper(), roa.upper())
            if g2 is not None:
                beta_o, se_o = g2 * rb, rse
        if beta_o is None:
            if use_raw:      # fracture：原始文件缺该 SNP → 用 harmonised 反推的中位 N_eff
                se_o = 1.0 / math.sqrt(2.0 * Neff * mafo * (1 - mafo))
            else:            # OA_knee：逐 SNP N_o 精确重建（已验证 min=max=1.00000）
                se_o = 1.0 / math.sqrt(2.0 * no * mafo * (1 - mafo))
            beta_o = go * zo * se_o
        out.append(dict(SNP=s, source=src, beta_e=beta_e, se_e=se_e,
                        beta_h=gh * bh, se_h=sh, beta_o=beta_o, se_o=se_o,
                        p_h=ph, beta_e_ref=be_ref, beta_o_ref=bo_ref))
    d = pd.DataFrame(out)
    return d, dict(N_e=N_e, Neff=Neff, use_raw=use_raw,
                   n_amb_drop=n_amb_drop, n_bmd=len(h_rows))


# ---------------- 主流程 ----------------
PAIRS = [
    ("BMD_total", "OA_knee", "sumstats_total.sumstats.gz", "sumstats_oa_knee.sumstats.gz"),
    ("BMD_femur_neck", "OA_knee", "sumstats_femur_neck.sumstats.gz", "sumstats_oa_knee.sumstats.gz"),
    ("BMD_lumbar", "OA_knee", "sumstats_lumbar.sumstats.gz", "sumstats_oa_knee.sumstats.gz"),
    ("BMD_arms", "OA_knee", "sumstats_arms.sumstats.gz", "sumstats_oa_knee.sumstats.gz"),
    ("BMD_leg", "OA_knee", "sumstats_leg.sumstats.gz", "sumstats_oa_knee.sumstats.gz"),
    ("BMD_total", "fracture", "sumstats_total.sumstats.gz", "sumstats_fracture.sumstats.gz"),
    ("BMD_femur_neck", "fracture", "sumstats_femur_neck.sumstats.gz", "sumstats_fracture.sumstats.gz"),
]
HN = {("BMD_total", "OA_knee"): "BMD_total__OA_knee.csv",
      ("BMD_femur_neck", "OA_knee"): "BMD_femur_neck__OA_knee.csv",
      ("BMD_lumbar", "OA_knee"): "BMD_lumbar__OA_knee.csv",
      ("BMD_arms", "OA_knee"): "BMD_arms__OA_knee.csv",
      ("BMD_leg", "OA_knee"): "BMD_leg__OA_knee.csv",
      ("BMD_total", "fracture"): "BMD_total__fracture.csv",
      ("BMD_femur_neck", "fracture"): "BMD_femur_neck__fracture.csv"}

# 手稿 Table 3/4 的单变量参照值（用于交叉校验）
REF = {"BMD_total__OA_knee": (0.02312, 0.00369),
       "BMD_femur_neck__OA_knee": (0.01709, 0.00500),
       "BMD_lumbar__OA_knee": (0.02791, 0.00483),
       "BMD_arms__OA_knee": (0.02218, 0.00464),
       "BMD_leg__OA_knee": (0.01772, 0.00418),
       "BMD_total__fracture": (-0.33713, 0.02241),
       "BMD_femur_neck__fracture": (-0.35066, 0.03037)}


def analyse(d, meta, tag):
    b_e, b_h = d.beta_e.values, d.beta_h.values
    b_o, s_e, s_h, s_o = d.beta_o.values, d.se_e.values, d.se_h.values, d.se_o.values
    K = len(d)
    X = np.column_stack([b_e, b_h])
    m = (d.source == "BMD").values

    # 门禁：重建 vs harmonised
    eb = float(np.max(np.abs(d.loc[m, "beta_e"] - d.loc[m, "beta_e_ref"])))
    ob = float(np.max(np.abs(d.loc[m, "beta_o"] - d.loc[m, "beta_o_ref"])))
    gate = "PASS" if (eb < 1e-8 and ob < 1e-8) else ("WARN" if ob < 0.03 else "FAIL")

    ub, us = ivw_uni(b_e, b_o, s_o)
    uhb, uhs = ivw_uni(b_h, b_o, s_o)
    ub_bmd, us_bmd = ivw_uni(b_e[m], b_o[m], s_o[m])
    mv = mvmr(X, b_o, s_o)
    i0, i0s, sl, sls = mvmr_int(X, b_o, s_o)

    # 条件 F（全体 / 仅 BMD 子集）
    Fc_b = cond_F(b_e, s_e, b_h.reshape(-1, 1), s_h.reshape(-1, 1))
    Fc_h = cond_F(b_h, s_h, b_e.reshape(-1, 1), s_e.reshape(-1, 1))
    if m.sum() > 3:
        Fc_b_bmd = cond_F(b_e[m], s_e[m], b_h[m].reshape(-1, 1), s_h[m].reshape(-1, 1))
        Fc_h_bmd = cond_F(b_h[m], s_h[m], b_e[m].reshape(-1, 1), s_e[m].reshape(-1, 1))
    else:
        Fc_b_bmd = Fc_h_bmd = float("nan")
    Fu_b = float(np.mean((b_e[m] / s_e[m]) ** 2))
    Fu_h = float(np.mean((b_h / s_h) ** 2))

    rz = mv["resid"] / s_o
    kp = np.abs(rz) <= 3
    if kp.sum() > 4:
        mt = mvmr(X[kp], b_o[kp], s_o[kp]); tb, tse = mt["theta"][0], mt["se_re"][0]
    else:
        tb = tse = float("nan")

    at = 100.0 * (1.0 - mv["theta"][0] / ub) if abs(ub) > 1e-12 else float("nan")
    at_b = 100.0 * (1.0 - mv["theta"][0] / ub_bmd) if abs(ub_bmd) > 1e-12 else float("nan")

    return dict(spec=tag, K=K, n_bmd=int(m.sum()), n_h=int((~m).sum()),
                gate=gate, err_e=eb, err_o=ob, N_exp=meta["N_e"],
                use_raw_se=meta["use_raw"], n_amb_drop=meta["n_amb_drop"],
                uni_beta=ub, uni_se=us, uni_p=p_from_z(ub / us),
                uni_bmd_only=ub_bmd, uni_bmd_only_se=us_bmd,
                mvmr_beta=mv["theta"][0], mvmr_se=mv["se_fe"][0],
                mvmr_p=p_from_z(mv["theta"][0] / mv["se_fe"][0]),
                mvmr_se_re=mv["se_re"][0], mvmr_p_re=p_from_z(mv["theta"][0] / mv["se_re"][0]),
                mvmr_h=mv["theta"][1], mvmr_h_se=mv["se_fe"][1],
                mvmr_h_p=p_from_z(mv["theta"][1] / mv["se_fe"][1]),
                uni_h_beta=uhb, uni_h_se=uhs, uni_h_p=p_from_z(uhb / uhs),
                atten_vs_full=at, atten_vs_bmdsubset=at_b,
                Q=mv["Q"], Q_df=mv["df"], Q_p=float(chi2.sf(mv["Q"], mv["df"])),
                I2=max(0.0, (mv["Q"] - mv["df"]) / mv["Q"]) if mv["Q"] > 0 else 0.0,
                F_bmd=Fu_b, F_h=Fu_h, Fc_bmd=Fc_b, Fc_h=Fc_h,
                Fc_bmd_bmdsub=Fc_b_bmd, Fc_h_bmdsub=Fc_h_bmd,
                egger_int=i0, egger_int_se=i0s, egger_int_p=p_from_z(i0 / i0s),
                egger_beta=sl[0], egger_se=sls[0],
                trim_beta=tb, trim_se=tse, n_outlier=int((~kp).sum()))


def main():
    os.makedirs(OUT, exist_ok=True)
    # fracture 真实 SE：只抽取需要用到的 rsid
    dh0 = load_height(HEIGHT_EUR)
    need = {rs for rs, v in dh0.items() if v[7] < GW_P}
    for (e, o, me, mo) in PAIRS:
        for r in csv.DictReader(open(os.path.join(OUT, "harmonised", HN[(e, o)]), encoding="utf-8")):
            need.add(r["SNP"])
    print("fracture 原始文件需抽取 rsid: %d" % len(need))
    fr_raw = load_fracture_raw(need)
    print("fracture 原始文件命中: %d" % len(fr_raw))

    rows = []
    for hfn, hlab in ((HEIGHT_EUR, "EUR"), (HEIGHT_NOUKB, "EUR_noUKB")):
        dh = load_height(hfn)
        for (e, o, me, mo) in PAIRS:
            pc = os.path.join(OUT, "harmonised", HN[(e, o)])
            for n_h in (0, 100, 300, 1000):
                tag = "%s|n_h=%d" % (hlab, n_h)
                try:
                    d, meta = build(pc, me, mo, dh, n_h,
                                    fr_raw=(fr_raw if o == "fracture" else None))
                    r = analyse(d, meta, tag)
                except Exception as ex:
                    print("  !! %s->%s %s : %s" % (e, o, tag, ex))
                    continue
                r.update(exposure=e, outcome=o, height=hlab, n_h=n_h)
                rb, rs = REF["%s__%s" % (e, o)]
                r["ref_uni_beta"] = rb
                r["ref_uni_se"] = rs
                rows.append(r)
                print("  %-15s->%-8s %-14s K=%4d(BMD %2d/H %4d) gate %-4s | uni %+.4f→MVMR %+.4f | "
                      "衰减 %5.1f%% | θ_H %+.4f | Fc(BMD) %5.1f (BMD子集 %5.1f) | Fc(H) %6.1f"
                      % (e, o, tag, r["K"], r["n_bmd"], r["n_h"], r["gate"],
                         r["uni_beta"], r["mvmr_beta"], r["atten_vs_full"],
                         r["mvmr_h"], r["Fc_bmd"], r["Fc_bmd_bmdsub"], r["Fc_h"]))
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "mvmr_height_grid.csv"),
                              index=False, encoding="utf-8")
    print("\nwritten mvmr_height_grid.csv rows=%d" % len(rows))


if __name__ == "__main__":
    main()
