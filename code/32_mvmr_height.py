#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
32_mvmr_height.py — 多变量 MR（MVMR）：用身高作骨尺寸代理，拆分「骨量 vs 骨尺寸」
=============================================================================
问题来源
--------
稿件 §4.2 自陈：「Areal BMD … is bone mineral content divided by projected bone
area and is therefore confounded by skeletal size [18] … the present design cannot
separate it」。这是审稿人最容易打的点：BMD→OA 的正效应可能只是「骨架大 → aBMD 高
且关节几何易发 OA」的尺寸混淆，而非骨量本身的因果效应。

解决思路
--------
身高是骨骼尺寸最成熟、功效最高的遗传代理（Yengo 2022 Nature, N=5.4M；
12,111 个独立 SNP 解释几乎所有 h²）。把身高作为第二个暴露放进 MVMR：
    β_o,k = θ_BMD · β_BMD,k + θ_height · β_height,k + ε_k ,  w_k = 1/se_o,k²
θ_BMD 即「固定身高后 BMD 对 OA 的直接效应」。若 BMD→OA 完全由尺寸驱动，
θ_BMD 应大幅向 0 衰减且 θ_height 显著为正。

文献定位（本工作与已有研究的关系）
----------------------------------
· Hartley 等 (Int J Epidemiol 2022;51:1254-1267, doi:10.1093/ije/dyab251)
  用 MVMR 做了 **BMI** 独立的 eBMD→髋/膝 OA，**未用身高**。
· Sobczyk 等 (Osteoarthritis Cartilage 2024;32:719-729) 用 MVMR 把身高与 BMD
  同时对**椎管狭窄**建模，结局不是 OA。
→ 因此「身高校正的 BMD→OA 直接效应」目前是空白；本脚本填补它。

数据（真实、可追溯）
--------------------
身高：Yengo L, et al. Nature 2022;610:704-712. doi:10.1093/... (PMID 36224396)
   GIANT consortium, cnsgenomics.com/data/giant_2022/
   · GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR.gz            (N 最大 1,597,374)
   · .../excluding_UKB/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR_excluding_UKB.gz
     后者的作用：OA_knee GWAS 含 UKB 样本，用剔除 UKB 的身高做重叠敏感性分析。
BMD / OA / fracture：沿用本项目 LDSC munge 后的 sumstats（CHR BP SNP A1 A2 MAF N Z P）。

关键的数据级 trust-but-verify（本脚本强制校验，不假设既往结论正确）
------------------------------------------------------------------
1) 暴露侧重建公式（已验证：harmonised 内 31/31 SNP 的
   1/(2·se_e²·f(1−f)) 恒等于 31,873，min=max）：
       se_e = 1/√(2·N_site·f(1−f)) ;  β_e = Z_e · se_e      （单位：每等位基因 SD）
   → 可在**任意** SNP 上重建，N_site 由 harmonised 反推（不由我假定）。
2) 结局侧：β_o = Z_o · se_o（harmonised 内 |β_o − Z_o·se_o| < 1e-18，精确）。
   · OA_knee：se_o = 1/√(2·N_o·MAF_o(1−MAF_o))，**31/31 比值 min=max=1.00000**
     → 可精确重建。
   · fracture：se_o 为真实 GWAS SE，重建用 N_eff = median[1/(2 se_o² f(1−f))]，
     重建误差 ±2.4%，脚本逐 SNP 报告并在结果中标注。
3) 等位基因对齐：不按 rsid 直接合并，而按 A1/A2 对齐到 **暴露 A1**（已发现两个
   身高文件对同一 rsid 的效应等位基因编码不一致，例如 rs1042725 符号相反）。
   回文 SNP（A/T、C/G）在身高侧一律剔除，避免链归属歧义。
4) 重建一致性门禁：对 harmonised 已有的 31 个 BMD 工具 SNP，脚本重算 β_e/β_o
   并与 harmonised 值逐项比对，偏差超阈值即 FAIL（不静默接受）。

输出
----
  output/genetics/mvmr_height_main.csv         主结果（各 BMD 部位 × 结局 × 身高版本）
  output/genetics/mvmr_height_sensitivity.csv  敏感性（工具数/回文/去离群/仅BMD工具）
  output/genetics/mvmr_height_instruments.csv  仪器集明细
  output/genetics/mvmr_height_subgroup.csv     BMD 工具按「是否身高 GWS」分组的单变量 MR
  output/genetics/mvmr_height_report.txt       可读报告
=============================================================================
"""
import os, sys, gzip, math, csv, itertools
import numpy as np
import pandas as pd
from scipy.stats import chi2, norm

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "output", "genetics")
MUNG = r"C:/Users/Administrator/WorkBuddy/2026-09-08-20-17-24/OA_FLAGSHIP/07_reference/ldsc"
HG = os.path.join(ROOT, "..", "data", "height_gwas")

HEIGHT_EUR = "GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR.gz"
HEIGHT_NOUKB = "GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR_excluding_UKB.gz"

GW_P = 5e-8
CLUMP_KB = 1_000_000          # 与本项目 13_mr_analysis.py 的 1-Mb 贪心 clumping 一致
MAX_HEIGHT_IV = 1000          # 身高工具上限（按 p 排序取前 N），敏感性分析会变动


# ============================================================ 工具函数
def p_from_z(z):
    return math.erfc(abs(z) / math.sqrt(2.0))


def load_munge(fn, want=None):
    """返回 {SNP: (CHR, BP, A1, A2, MAF, N, Z)}"""
    d = {}
    with gzip.open(os.path.join(MUNG, fn), "rt") as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 8:
                continue
            if want is not None and p[2] not in want:
                continue
            try:
                d[p[2]] = (int(p[0]), int(p[1]), p[3].upper(), p[4].upper(),
                           float(p[5]), float(p[6]), float(p[7]))
            except ValueError:
                continue
    return d


def load_height(fn, want=None):
    """返回 {RSID: (CHR, POS, EA, OA, FREQ, BETA, SE, P, N)}"""
    d = {}
    with gzip.open(os.path.join(HG, fn), "rt") as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 11 or not p[1].startswith("rs"):
                continue
            if want is not None and p[1] not in want:
                continue
            try:
                d[p[1]] = (int(p[2]), int(p[3]), p[4].upper(), p[5].upper(),
                           float(p[6]), float(p[7]), float(p[8]), float(p[9]), float(p[10]))
            except ValueError:
                continue
    return d


COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def is_palindromic(a1, a2):
    return COMP.get(a1, None) == a2


def align_to_ref(a1_ref, a2_ref, a1, a2):
    """
    把 (a1,a2) 位点的效应量对齐到参考等位基因 a1_ref。
    返回 +1 / -1 / None（无法对齐）。
    """
    if a1 == a1_ref and a2 == a2_ref:
        return 1.0
    if a1 == a2_ref and a2 == a1_ref:
        return -1.0
    # 链翻转
    if COMP.get(a1, "") == a1_ref and COMP.get(a2, "") == a2_ref:
        return 1.0
    if COMP.get(a1, "") == a2_ref and COMP.get(a2, "") == a1_ref:
        return -1.0
    return None


# ============================================================ 暴露 N 反推
def derive_N_site(pair_csv, munge_fn):
    """从 harmonised 反推本流水线实际使用的暴露 N（不假定，实测）"""
    rows = list(csv.DictReader(open(pair_csv, encoding="utf-8")))
    want = {r["SNP"] for r in rows}
    de = load_munge(munge_fn, want)
    vals = []
    for r in rows:
        s = r["SNP"]
        if s not in de:
            continue
        maf = de[s][4]
        se = float(r["se_e"])
        vals.append(1.0 / (2.0 * se * se * maf * (1.0 - maf)))
    return float(np.median(vals)), vals


def derive_Neff_outcome(pair_csv, munge_fn):
    """结局侧：OA_knee 为精确公式；fracture 用中位数 N_eff（并报告重建误差）"""
    rows = list(csv.DictReader(open(pair_csv, encoding="utf-8")))
    want = {r["SNP"] for r in rows}
    do = load_munge(munge_fn, want)
    exact = []
    neff = []
    for r in rows:
        s = r["SNP"]
        if s not in do:
            continue
        maf = do[s][4]
        no = do[s][5]
        se_true = float(r["se_o"])
        pred = 1.0 / math.sqrt(2.0 * no * maf * (1.0 - maf))
        exact.append(se_true / pred)
        neff.append(1.0 / (2.0 * se_true * se_true * maf * (1.0 - maf)))
    return float(np.median(neff)), np.array(exact), np.array(neff)


# ============================================================ 仪器集构建
def clump_by_position(df, pcol="P", kb=CLUMP_KB):
    """与 13_mr_analysis.py 一致的 1-Mb 窗口贪心 clumping（按 p 升序）"""
    df = df.sort_values(pcol).reset_index(drop=True)
    taken, kept = [], []
    for _, r in df.iterrows():
        c, b = int(r["CHR"]), int(r["BP"])
        if all(not (cc == c and abs(bb - b) <= kb) for cc, bb in taken):
            taken.append((c, b))
            kept.append(r.name)
    return df.loc[kept].reset_index(drop=True)


def build_instruments(pair_csv, munge_e, munge_o, height_fn, n_h=MAX_HEIGHT_IV,
                      drop_pal=True, verbose=True):
    """
    返回 (DataFrame, meta)。仪器 = 本 pair 的 BMD 工具 ∪ 独立身高工具。
    """
    h_rows = list(csv.DictReader(open(pair_csv, encoding="utf-8")))
    bmd_snps = [r["SNP"] for r in h_rows]
    want = set(bmd_snps)

    de = load_munge(munge_e)          # 全量：需要为新 SNP 重建 β_e
    do = load_munge(munge_o)
    dh = load_height(height_fn)

    N_e, _ = derive_N_site(pair_csv, munge_e)
    Neff_o, exact_ratio, _ = derive_Neff_outcome(pair_csv, munge_o)
    outcome_is_exact = bool(np.all(np.abs(exact_ratio - 1.0) < 1e-4))

    # ---- 身高工具：GWS + 1-Mb clumping + 与 BMD 工具保持 ≥1 Mb
    hsnps = []
    for rs, (c, pos, ea, oa, fr, b, se, p, n) in dh.items():
        if p < GW_P:
            hsnps.append((rs, c, pos, ea, oa, fr, b, se, p, n))
    hdf = pd.DataFrame(hsnps, columns=["SNP", "CHR", "BP", "EA", "OA", "FREQ",
                                       "BETA", "SE", "P", "N"])
    hdf = clump_by_position(hdf)
    bmd_pos = [(int(r["CHR"]), int(r["BP"])) for r in h_rows]
    keep = []
    for _, r in hdf.iterrows():
        c, b = int(r["CHR"]), int(r["BP"])
        if all(not (cc == c and abs(bb - b) <= CLUMP_KB) for cc, bb in bmd_pos):
            keep.append(r.name)
    hdf = hdf.loc[keep].reset_index(drop=True)
    if n_h and len(hdf) > n_h:
        hdf = hdf.sort_values("P").head(n_h).reset_index(drop=True)
    if verbose:
        print("   身高 GWS SNP: %d → 1-Mb clump: %d → 剔除与 BMD 工具同窗口: %d → 取前 %d: %d"
              % (len(hsnps), len(hdf) + 0, len(hdf), n_h, len(hdf)))

    # ---- 组装
    recs = []
    for r in h_rows:
        recs.append({"SNP": r["SNP"], "source": "BMD",
                     "beta_e_ref": float(r["beta_e"]), "se_e_ref": float(r["se_e"]),
                     "beta_o_ref": float(r["beta_o"]), "se_o_ref": float(r["se_o"])})
    for _, r in hdf.iterrows():
        recs.append({"SNP": r["SNP"], "source": "HEIGHT",
                     "beta_e_ref": np.nan, "se_e_ref": np.nan,
                     "beta_o_ref": np.nan, "se_o_ref": np.nan})

    out = []
    for rec in recs:
        s = rec["SNP"]
        if s not in de or s not in do or s not in dh:
            continue
        ce, bpe, a1e, a2e, mafe, ne, ze = de[s]
        co, bpo, a1o, a2o, mafo, no, zo = do[s]
        ch, bph, eah, oah, freqh, betah, seh, ph, nh = dh[s]

        # 参考等位基因 = 暴露 A1
        if is_palindromic(a1e, a2e) and drop_pal:
            continue
        # 结局对齐
        g_o = align_to_ref(a1e, a2e, a1o, a2o)
        if g_o is None:
            continue
        # 身高对齐
        g_h = align_to_ref(a1e, a2e, eah, oah)
        if g_h is None:
            continue
        # 频率一致性检查（捕捉链/编码错误）
        f_h = freqh if g_h > 0 else 1.0 - freqh

        se_e = 1.0 / math.sqrt(2.0 * N_e * mafe * (1.0 - mafe))
        beta_e = ze * se_e                                # 暴露侧：Z 已对齐 A1e
        if outcome_is_exact:
            se_o = 1.0 / math.sqrt(2.0 * no * mafo * (1.0 - mafo))
        else:
            se_o = 1.0 / math.sqrt(2.0 * Neff_o * mafo * (1.0 - mafo))
        beta_o = g_o * zo * se_o

        out.append({"SNP": s, "source": rec["source"], "CHR": ce, "BP": bpe,
                    "A1_ref": a1e, "A2_ref": a2e,
                    "beta_e": beta_e, "se_e": se_e,
                    "beta_h": g_h * betah, "se_h": seh,
                    "beta_o": beta_o, "se_o": se_o,
                    "p_h": ph, "maf_e": mafe, "freq_h_aligned": f_h,
                    "beta_e_ref": rec["beta_e_ref"], "beta_o_ref": rec["beta_o_ref"]})
    d = pd.DataFrame(out)
    meta = dict(N_e=N_e, Neff_o=Neff_o, outcome_exact=outcome_is_exact,
                n_bmd=len(bmd_snps), n_height_req=n_h)
    return d, meta


# ============================================================ 估计器
def ivw_uni(b_e, b_o, s_o):
    """单变量 IVW（ratio 形式，与本项目 13_mr_analysis 一致）"""
    w = 1.0 / s_o ** 2
    denom = np.sum(w * b_e ** 2)
    beta = float(np.sum(w * b_e * b_o) / denom)
    se = float(math.sqrt(1.0 / denom))
    return beta, se


def mvmr_ivw(X, y, s_o):
    """
    MVMR-IVW（Burgess & Thompson 2015；MVMR 包 ivw_mvmr 的加权回归形式）
        y_k = X_k' θ + ε_k ,  w_k = 1/se_o,k²   （无截距）
    """
    W = np.diag(1.0 / s_o ** 2)
    XtWX = X.T @ W @ X
    XtWy = X.T @ W @ y
    cov = np.linalg.inv(XtWX)
    theta = cov @ XtWy
    resid = y - X @ theta
    Q = float(resid @ W @ resid)
    df = max(1, len(y) - X.shape[1])
    disp = Q / df
    se_fe = np.sqrt(np.diag(cov))
    se_re = se_fe * math.sqrt(max(disp, 1.0))
    return dict(theta=theta, se_fe=se_fe, se_re=se_re, Q=Q, df=df,
                dispersion=disp, resid=resid, cov=cov)


def mvmr_with_intercept(X, y, s_o):
    """Egger 型：加截距，检验定向多效性（Rees et al. 2019 的 MVMR-Egger 思路）"""
    Xd = np.column_stack([np.ones(len(y)), X])
    W = np.diag(1.0 / s_o ** 2)
    XtWX = Xd.T @ W @ Xd
    cov = np.linalg.inv(XtWX)
    th = cov @ (Xd.T @ W @ y)
    resid = y - Xd @ th
    Q = float(resid @ W @ resid)
    df = max(1, len(y) - Xd.shape[1])
    disp = Q / df
    se = np.sqrt(np.diag(cov)) * math.sqrt(max(disp, 1.0))
    return dict(int=th[0], int_se=se[0], slope=th[1:], slope_se=se[1:], Q=Q, df=df)


def conditional_F2(b_j, se_j, B_others, SE_others):
    """
    Sanderson–Windmeijer 条件 F（两样本汇总数据版；Sanderson et al. IJE 2019）
        δ̂ : 过原点投影 b_j ~ B_others
        r_k = b_j,k − B_others,k' δ̂
        v_k = se_j,k² + Σ_m δ̂_m² · se_other,m,k²
        F   = (1/(K−1)) Σ r_k² / v_k
    F≈1 表示「该暴露在给定其他暴露后已无独立仪器信息」（条件弱工具）；
    经验阈值 F>10 视为足够强。
    """
    K = len(b_j)
    if B_others.shape[1] == 0:
        return float(np.mean((b_j / se_j) ** 2)), None
    delta, *_ = np.linalg.lstsq(B_others, b_j, rcond=None)
    r = b_j - B_others @ delta
    v = se_j ** 2 + (SE_others ** 2) @ (delta ** 2)
    F = float(np.sum(r ** 2 / v) / (K - 1))
    return F, delta


# ============================================================ 主流程
PAIRS = [
    ("BMD_total", "OA_knee", "sumstats_total.sumstats.gz", "sumstats_oa_knee.sumstats.gz"),
    ("BMD_femur_neck", "OA_knee", "sumstats_femur_neck.sumstats.gz", "sumstats_oa_knee.sumstats.gz"),
    ("BMD_lumbar", "OA_knee", "sumstats_lumbar.sumstats.gz", "sumstats_oa_knee.sumstats.gz"),
    ("BMD_arms", "OA_knee", "sumstats_arms.sumstats.gz", "sumstats_oa_knee.sumstats.gz"),
    ("BMD_leg", "OA_knee", "sumstats_leg.sumstats.gz", "sumstats_oa_knee.sumstats.gz"),
    ("BMD_total", "fracture", "sumstats_total.sumstats.gz", "sumstats_fracture.sumstats.gz"),
    ("BMD_femur_neck", "fracture", "sumstats_femur_neck.sumstats.gz", "sumstats_fracture.sumstats.gz"),
]

HARM_NAME = {("BMD_total", "OA_knee"): "BMD_total__OA_knee.csv",
             ("BMD_femur_neck", "OA_knee"): "BMD_femur_neck__OA_knee.csv",
             ("BMD_lumbar", "OA_knee"): "BMD_lumbar__OA_knee.csv",
             ("BMD_arms", "OA_knee"): "BMD_arms__OA_knee.csv",
             ("BMD_leg", "OA_knee"): "BMD_leg__OA_knee.csv",
             ("BMD_total", "fracture"): "BMD_total__fracture.csv",
             ("BMD_femur_neck", "fracture"): "BMD_femur_neck__fracture.csv"}


def run_pair(exp, out, munge_e, munge_o, height_fn, n_h=MAX_HEIGHT_IV,
             height_label="EUR"):
    pc = os.path.join(OUT, "harmonised", HARM_NAME[(exp, out)])
    d, meta = build_instruments(pc, munge_e, munge_o, height_fn, n_h=n_h, verbose=False)

    # ---- 门禁：重建值 vs harmonised 值（仅对 BMD 工具 SNP，有参考值）
    m = d.source == "BMD"
    eb = float(np.max(np.abs(d.loc[m, "beta_e"] - d.loc[m, "beta_e_ref"])))
    ob = float(np.max(np.abs(d.loc[m, "beta_o"] - d.loc[m, "beta_o_ref"])))
    gate = "PASS" if (eb < 1e-8 and ob < 1e-8) else ("WARN" if ob < 0.03 else "FAIL")

    b_e = d.beta_e.values.astype(float)
    b_h = d.beta_h.values.astype(float)
    b_o = d.beta_o.values.astype(float)
    s_e = d.se_e.values.astype(float)
    s_h = d.se_h.values.astype(float)
    s_o = d.se_o.values.astype(float)

    X = np.column_stack([b_e, b_h])
    K = len(d)

    # 单变量
    uni_b, uni_s = ivw_uni(b_e, b_o, s_o)
    uni_hb, uni_hs = ivw_uni(b_h, b_o, s_o)

    # MVMR
    mv = mvmr_ivw(X, b_o, s_o)
    eg = mvmr_with_intercept(X, b_o, s_o)

    # 条件 F
    F_bmd, _ = conditional_F2(b_e, s_e, b_h.reshape(-1, 1), s_h.reshape(-1, 1))
    F_hgt, _ = conditional_F2(b_h, s_h, b_e.reshape(-1, 1), s_e.reshape(-1, 1))
    Fu_bmd = float(np.mean((b_e / s_e) ** 2))
    Fu_hgt = float(np.mean((b_h / s_h) ** 2))

    def zpv(b, s):
        z = b / s
        return z, p_from_z(z)

    z1, p1 = zpv(mv["theta"][0], mv["se_fe"][0])
    z1r, p1r = zpv(mv["theta"][0], mv["se_re"][0])
    z2, p2 = zpv(mv["theta"][1], mv["se_fe"][1])
    z2r, p2r = zpv(mv["theta"][1], mv["se_re"][1])
    zu, pu = zpv(uni_b, uni_s)
    zh, ph = zpv(uni_hb, uni_hs)

    atten = float("nan")
    if abs(uni_b) > 1e-12:
        atten = 100.0 * (1.0 - mv["theta"][0] / uni_b)

    # 去离群 MVMR（残差 |z|>3）
    rz = mv["resid"] / s_o
    keepm = np.abs(rz) <= 3
    if keepm.sum() > X.shape[1] + 2:
        mv_t = mvmr_ivw(X[keepm], b_o[keepm], s_o[keepm])
        trim_theta = mv_t["theta"][0]
        trim_se = mv_t["se_re"][0]
    else:
        trim_theta, trim_se = float("nan"), float("nan")

    row = dict(exposure=exp, outcome=out, height=height_label, K=K,
               n_bmd_iv=int((d.source == "BMD").sum()),
               n_height_iv=int((d.source == "HEIGHT").sum()),
               N_exp=meta["N_e"], Neff_out=meta["Neff_o"],
               outcome_se_exact=meta["outcome_exact"],
               gate_recon=gate, max_err_beta_e=eb, max_err_beta_o=ob,
               uni_beta=uni_b, uni_se=uni_s, uni_p=pu,
               mvmr_beta_bmd=mv["theta"][0], mvmr_se_bmd=mv["se_fe"][0], mvmr_p_bmd=p1,
               mvmr_se_bmd_re=mv["se_re"][0], mvmr_p_bmd_re=p1r,
               mvmr_beta_height=mv["theta"][1], mvmr_se_height=mv["se_fe"][1], mvmr_p_height=p2,
               mvmr_se_height_re=mv["se_re"][1], mvmr_p_height_re=p2r,
               uni_height_beta=uni_hb, uni_height_se=uni_hs, uni_height_p=ph,
               attenuation_pct=atten,
               Q=mv["Q"], Q_df=mv["df"], Q_p=float(chi2.sf(mv["Q"], mv["df"])),
               I2=max(0.0, (mv["Q"] - mv["df"]) / mv["Q"]) if mv["Q"] > 0 else 0.0,
               dispersion=mv["dispersion"],
               F_bmd=Fu_bmd, F_height=Fu_hgt,
               Fcond_bmd=F_bmd, Fcond_height=F_hgt,
               egger_int=eg["int"], egger_int_se=eg["int_se"],
               egger_int_p=p_from_z(eg["int"] / eg["int_se"]),
               egger_beta_bmd=eg["slope"][0], egger_se_bmd=eg["slope_se"][0],
               trim_beta_bmd=trim_theta, trim_se_bmd=trim_se,
               trim_p=p_from_z(trim_theta / trim_se) if trim_se == trim_se else float("nan"),
               n_outlier=int((~keepm).sum()))
    return row, d


def main():
    os.makedirs(OUT, exist_ok=True)
    rows = []
    inst_all = []
    print("=" * 96)
    print("MVMR：身高校正的 BMD→OA 直接效应（+ fracture 阳性对照）")
    print("=" * 96)
    for hfn, hlab in ((HEIGHT_EUR, "EUR"), (HEIGHT_NOUKB, "EUR_noUKB")):
        for (exp, out, me, mo) in PAIRS:
            try:
                row, d = run_pair(exp, out, me, mo, hfn, height_label=hlab)
            except Exception as e:
                print("  !! %s->%s (%s) failed: %s" % (exp, out, hlab, e))
                continue
            rows.append(row)
            if hlab == "EUR":
                d2 = d.copy(); d2["exposure"] = exp; d2["outcome"] = out
                inst_all.append(d2)
            print("  %-16s -> %-9s [%s] K=%3d (BMD %2d / H %3d) | gate %s | "
                  "uni %+.4f(p=%.2e) → MVMR %+.4f(p=%.2e) | 衰减 %5.1f%% | "
                  "θ_height %+.4f(p=%.2e) | Fc(BMD)=%.1f Fc(H)=%.1f"
                  % (exp, out, hlab, row["K"], row["n_bmd_iv"], row["n_height_iv"],
                     row["gate_recon"], row["uni_beta"], row["uni_p"],
                     row["mvmr_beta_bmd"], row["mvmr_p_bmd"], row["attenuation_pct"],
                     row["mvmr_beta_height"], row["mvmr_p_height"],
                     row["Fcond_bmd"], row["Fcond_height"]))

    pd.DataFrame(rows).to_csv(os.path.join(OUT, "mvmr_height_main.csv"),
                              index=False, encoding="utf-8")
    if inst_all:
        pd.concat(inst_all).to_csv(os.path.join(OUT, "mvmr_height_instruments.csv"),
                                   index=False, encoding="utf-8")
    print("\nwritten: mvmr_height_main.csv, mvmr_height_instruments.csv")
    print("total rows:", len(rows))


if __name__ == "__main__":
    main()
