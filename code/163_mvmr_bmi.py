#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
163_mvmr_bmi.py — 多变量 MR（MVMR）：用 BMI 作第二暴露，拆分「BMD→OA 残差是否 BMI 中介」
====================================================================================
来源：docs/117 §4.2 K-b（核心理论结 = 骨量 vs 几何，height-independent 残差 67–94%
      在定义上就是 mass/geometry 混合；前沿是"残差中骨量占多少、几何占多少"）。
      docs/117 §6.3 第 3 项："补 BMI 条件 MVMR（K-b 残留）"——直接回应"残差仍可能
      BMI 中介"的局限（BMI 本身是 OA 风险因子，故稿件刻意只条件 height、避开 BMI）。

数据（真实、已缓存，零下载）
----------------------------
· BMI：sumstats_bmi.sumstats.gz（与本项目 BMD/OA 同一 munge 格式
       CHR BP SNP A1 A2 MAF N Z，位于 FLAGSHIP_LDSC）。
· 身高：data/height_gwas/（仅作敏感性：BMI 与身高高度相关，可做 BMI+身高 三暴露，
       但主分析为 BMD+BMI 双暴露，与稿件 height-MVMR 并列）。
· BMD / OA_knee / fracture：沿用本项目 LDSC munge 后的 sumstats。

方法学（与 code/32_mvmr_height.py 完全对齐，仅把"身高"换成"BMI"，且 BMI 走
munge 格式而非 GIANT 专用格式）
-------------------------------------------------------------------------
    β_o,k = θ_BMD · β_BMD,k + θ_BMI · β_BMI,k + ε_k ,  w_k = 1/se_o,k²
θ_BMD = 固定 BMI 后 BMD 对 OA 的直接效应（K-b 目标量）。
若 BMD→OA 完全由 BMI 中介 → θ_BMD 大幅向 0 衰减且 θ_BMI 显著为正。
若 θ_BMD 仍显著 >0 且衰减有限 → 残差非 BMI 中介，支持"非尺寸的骨生物学"（K-c）。

信任-但-验证（强制门禁，不假设既往结论）
----------------------------------------
1) 暴露侧重建：se_e = 1/√(2·N_site·f(1−f))，β_e = Z_e·se_e（N_site 由 harmonised 反推）。
2) 结局侧：OA_knee 精确公式；fracture 用中位数 N_eff（报告重建误差）。
3) 等位基因按 A1/A2 对齐到 BMD 暴露 A1；回文 SNP 剔除。
4) 重建一致性门禁：对 harmonised 31 个 BMD 工具 SNP 重算 β_e/β_o 与 harmonised 比对。

输出
----
  output/genetics/mvmr_bmi_main.csv          主结果（各 BMD 部位 × 结局）
  output/genetics/mvmr_bmi_instruments.csv   仪器集明细
  output/genetics/mvmr_bmi_report.txt        可读报告
"""
import os, sys, math, csv
import importlib.util
import numpy as np
import pandas as pd
from scipy.stats import chi2

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "output", "genetics")
MUNG = r"C:/Users/Administrator/WorkBuddy/2026-09-08-20-17-24/OA_FLAGSHIP/07_reference/ldsc"
BMI_FILE = "sumstats_bmi.sumstats.gz"

# 复用 code/32 的成熟函数
_spec = importlib.util.spec_from_file_location("m32", os.path.join(ROOT, "32_mvmr_height.py"))
m32 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m32)

GW_Z = 5.45                 # |Z| >= 5.45  ≈ p < 5e-8（munge 格式无 P 列，由 Z 反推）
CLUMP_KB = 1_000_000


def p_from_z(z):
    return math.erfc(abs(z) / math.sqrt(2.0))


def load_munge(fn, want=None):
    return m32.load_munge(fn, want)


def align_to_ref(a1_ref, a2_ref, a1, a2):
    return m32.align_to_ref(a1_ref, a2_ref, a1, a2)


def is_palindromic(a1, a2):
    return m32.is_palindromic(a1, a2)


def clump_by_position(df, pcol="P", kb=CLUMP_KB):
    return m32.clump_by_position(df, pcol=pcol, kb=kb)


COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


# ============================================================ 仪器集构建（BMI 版）
def build_instruments(pair_csv, munge_bmd, munge_bmi, munge_o,
                       n_bmi=800, verbose=True):
    """
    返回 (DataFrame, meta)。仪器 = 本 pair 的 BMD 工具 ∪ 独立 BMI GWS 工具。
    """
    h_rows = list(csv.DictReader(open(pair_csv, encoding="utf-8")))
    bmd_snps = [r["SNP"] for r in h_rows]
    want = set(bmd_snps)

    de = load_munge(munge_bmd)          # BMD 全量
    do = load_munge(munge_o)            # 结局
    db = load_munge(munge_bmi)          # BMI 全量

    N_e, _ = m32.derive_N_site(pair_csv, munge_bmd)
    Neff_o, exact_ratio, _ = m32.derive_Neff_outcome(pair_csv, munge_o)
    outcome_is_exact = bool(np.all(np.abs(exact_ratio - 1.0) < 1e-4))

    # ---- BMI 工具：GWS(|Z|>=5.45) + 1-Mb clumping + 与 BMD 工具保持 ≥1 Mb ----
    bsnps = []
    for s, (c, bp, a1, a2, maf, n, z) in db.items():
        if abs(z) >= GW_Z:
            bsnps.append((s, c, bp, a1, a2, maf, n, z))
    bdf = pd.DataFrame(bsnps, columns=["SNP", "CHR", "BP", "A1", "A2", "MAF", "N", "Z"])
    bdf = clump_by_position(bdf, pcol="Z")
    bmd_pos = [(int(r["CHR"]), int(r["BP"])) for r in h_rows]
    keep = []
    for _, r in bdf.iterrows():
        c, b = int(r["CHR"]), int(r["BP"])
        if all(not (cc == c and abs(bb - b) <= CLUMP_KB) for cc, bb in bmd_pos):
            keep.append(r.name)
    bdf = bdf.loc[keep].reset_index(drop=True)
    if n_bmi and len(bdf) > n_bmi:
        bdf = bdf.sort_values("Z", key=lambda s: s.abs(), ascending=False).head(n_bmi).reset_index(drop=True)
    if verbose:
        print("   BMI GWS SNP: %d → 1-Mb clump: %d → 剔除与 BMD 工具同窗口: %d → 取前 %d: %d"
              % (len(bsnps), len(bdf) + 0, len(bdf), n_bmi, len(bdf)))

    # ---- 组装：每个仪器 SNP 同时有 β_BMD 与 β_BMI（参考等位基因 = BMD A1）----
    recs = []
    for r in h_rows:
        recs.append({"SNP": r["SNP"], "source": "BMD",
                     "beta_e_ref": float(r["beta_e"]), "se_e_ref": float(r["se_e"]),
                     "beta_o_ref": float(r["beta_o"]), "se_o_ref": float(r["se_o"])})
    for _, r in bdf.iterrows():
        recs.append({"SNP": r["SNP"], "source": "BMI",
                     "beta_e_ref": np.nan, "se_e_ref": np.nan,
                     "beta_o_ref": np.nan, "se_o_ref": np.nan})

    out = []
    for rec in recs:
        s = rec["SNP"]
        if s not in de or s not in do or s not in db:
            continue
        ce, bpe, a1e, a2e, mafe, ne, ze = de[s]
        co, bpo, a1o, a2o, mafo, no, zo = do[s]
        cb, bpb, a1b, a2b, mafb, nb, zb = db[s]

        if is_palindromic(a1e, a2e):
            continue
        g_o = align_to_ref(a1e, a2e, a1o, a2o)
        if g_o is None:
            continue
        g_b = align_to_ref(a1e, a2e, a1b, a2b)
        if g_b is None:
            continue

        se_e = 1.0 / math.sqrt(2.0 * N_e * mafe * (1.0 - mafe))
        beta_e = ze * se_e
        if outcome_is_exact:
            se_o = 1.0 / math.sqrt(2.0 * no * mafo * (1.0 - mafo))
        else:
            se_o = 1.0 / math.sqrt(2.0 * Neff_o * mafo * (1.0 - mafo))
        beta_o = g_o * zo * se_o

        freq_b = mafb if g_b > 0 else 1.0 - mafb
        se_bmi = 1.0 / math.sqrt(2.0 * nb * freq_b * (1.0 - freq_b))
        beta_bmi = g_b * zb * se_bmi

        out.append({"SNP": s, "source": rec["source"], "CHR": ce, "BP": bpe,
                    "A1_ref": a1e, "A2_ref": a2e,
                    "beta_e": beta_e, "se_e": se_e,
                    "beta_bmi": beta_bmi, "se_bmi": se_bmi,
                    "beta_o": beta_o, "se_o": se_o,
                    "beta_e_ref": rec["beta_e_ref"], "beta_o_ref": rec["beta_o_ref"]})
    d = pd.DataFrame(out)
    meta = dict(N_e=N_e, Neff_o=Neff_o, outcome_exact=outcome_is_exact,
                n_bmd=len(bmd_snps), n_bmi_req=n_bmi)
    return d, meta


# ============================================================ 估计器（复用 m32）
ivw_uni = m32.ivw_uni
mvmr_ivw = m32.mvmr_ivw
mvmr_with_intercept = m32.mvmr_with_intercept
conditional_F2 = m32.conditional_F2


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

HARM_NAME = m32.HARM_NAME

# Table 2（稿件主表）BMD→结局单变量基线（BMD 仪器专属估计）——与身高校正表 S9 同口径
TABLE2_BASELINE = {
    ("BMD_total", "OA_knee"): 0.0231,
    ("BMD_femur_neck", "OA_knee"): 0.0171,
    ("BMD_lumbar", "OA_knee"): 0.0279,
    ("BMD_arms", "OA_knee"): 0.0222,
    ("BMD_leg", "OA_knee"): 0.0177,
    ("BMD_total", "fracture"): -0.3371,
}


def run_pair(exp, out, munge_e, munge_o, munge_bmi, n_bmi=800):
    pc = os.path.join(OUT, "harmonised", HARM_NAME[(exp, out)])
    d, meta = build_instruments(pc, munge_e, munge_bmi, munge_o, n_bmi=n_bmi, verbose=False)

    m = d.source == "BMD"
    eb = float(np.max(np.abs(d.loc[m, "beta_e"] - d.loc[m, "beta_e_ref"]))) if m.any() else 0.0
    ob = float(np.max(np.abs(d.loc[m, "beta_o"] - d.loc[m, "beta_o_ref"]))) if m.any() else 0.0
    gate = "PASS" if (eb < 1e-8 and ob < 1e-8) else ("WARN" if ob < 0.03 else "FAIL")

    b_e = d.beta_e.values.astype(float)
    b_b = d.beta_bmi.values.astype(float)
    b_o = d.beta_o.values.astype(float)
    s_e = d.se_e.values.astype(float)
    s_b = d.se_bmi.values.astype(float)
    s_o = d.se_o.values.astype(float)

    X = np.column_stack([b_e, b_b])
    K = len(d)

    uni_b, uni_s = ivw_uni(b_e, b_o, s_o)
    uni_bb, uni_bs = ivw_uni(b_b, b_o, s_o)

    mv = mvmr_ivw(X, b_o, s_o)
    eg = mvmr_with_intercept(X, b_o, s_o)

    F_bmd, _ = conditional_F2(b_e, s_e, b_b.reshape(-1, 1), s_b.reshape(-1, 1))
    F_bmi, _ = conditional_F2(b_b, s_b, b_e.reshape(-1, 1), s_e.reshape(-1, 1))
    Fu_bmd = float(np.mean((b_e / s_e) ** 2))
    Fu_bmi = float(np.mean((b_b / s_b) ** 2))

    def zpv(b, s):
        z = b / s
        return z, p_from_z(z)

    z1, p1 = zpv(mv["theta"][0], mv["se_fe"][0])
    z1r, p1r = zpv(mv["theta"][0], mv["se_re"][0])
    z2, p2 = zpv(mv["theta"][1], mv["se_fe"][1])
    z2r, p2r = zpv(mv["theta"][1], mv["se_re"][1])
    zu, pu = zpv(uni_b, uni_s)
    zb_, pb = zpv(uni_bb, uni_bs)

    atten = float("nan")
    if abs(uni_b) > 1e-12:
        atten = 100.0 * (1.0 - mv["theta"][0] / uni_b)

    rz = mv["resid"] / s_o
    keepm = np.abs(rz) <= 3
    if keepm.sum() > X.shape[1] + 2:
        mv_t = mvmr_ivw(X[keepm], b_o[keepm], s_o[keepm])
        trim_theta = mv_t["theta"][0]
        trim_se = mv_t["se_re"][0]
    else:
        trim_theta, trim_se = float("nan"), float("nan")

    row = dict(exposure=exp, outcome=out, K=K,
               n_bmd_iv=int((d.source == "BMD").sum()),
               n_bmi_iv=int((d.source == "BMI").sum()),
               N_exp=meta["N_e"], Neff_out=meta["Neff_o"],
               outcome_se_exact=meta["outcome_exact"],
               gate_recon=gate, max_err_beta_e=eb, max_err_beta_o=ob,
               uni_beta_bmd=uni_b, uni_se_bmd=uni_s, uni_p_bmd=pu,
               uni_beta_bmi=uni_bb, uni_se_bmi=uni_bs, uni_p_bmi=pb,
               mvmr_beta_bmd=mv["theta"][0], mvmr_se_bmd=mv["se_fe"][0], mvmr_p_bmd=p1,
               mvmr_se_bmd_re=mv["se_re"][0], mvmr_p_bmd_re=p1r,
               mvmr_beta_bmi=mv["theta"][1], mvmr_se_bmi=mv["se_fe"][1], mvmr_p_bmi=p2,
               mvmr_se_bmi_re=mv["se_re"][1], mvmr_p_bmi_re=p2r,
               attenuation_pct=atten,
               Q=mv["Q"], Q_df=mv["df"], Q_p=float(chi2.sf(mv["Q"], mv["df"])),
               I2=max(0.0, (mv["Q"] - mv["df"]) / mv["Q"]) if mv["Q"] > 0 else 0.0,
               dispersion=mv["dispersion"],
               F_bmd=Fu_bmd, F_bmi=Fu_bmi,
               Fcond_bmd=F_bmd, Fcond_bmi=F_bmi,
               egger_int=eg["int"], egger_int_se=eg["int_se"],
               egger_int_p=p_from_z(eg["int"] / eg["int_se"]),
               egger_beta_bmd=eg["slope"][0], egger_se_bmd=eg["slope_se"][0],
               trim_beta_bmd=trim_theta, trim_se_bmd=trim_se,
               trim_p=p_from_z(trim_theta / trim_se) if trim_se == trim_se else float("nan"),
               n_outlier=int((~keepm).sum()))
    return row, d


def main():
    os.makedirs(OUT, exist_ok=True)
    GRID = [100, 300, 800]          # 主配置 n_bmi=100（与身高校正表 S9 的 n_h=100 对齐）
    rows = []
    inst_all = []
    print("=" * 104)
    print("MVMR：BMI 条件的 BMD→OA 直接效应（K-b 残留；+ fracture 阳性对照）")
    print("仪器网格 n_bmi ∈ {100,300,800}（主=100，与身高校正表 S9 的 n_h=100 对齐）")
    print("atten_T2 = 对稿件 Table 2 单变量基线（与 S9 同口径）；atten_mvset = 对 BMD+BMI 联合集")
    print("=" * 104)
    for (exp, out, me, mo) in PAIRS:
        for n_bmi in GRID:
            try:
                row, d = run_pair(exp, out, me, mo, BMI_FILE, n_bmi=n_bmi)
            except Exception as e:
                print("  !! %s->%s (n_bmi=%d) failed: %s" % (exp, out, n_bmi, e))
                continue
            base2 = TABLE2_BASELINE.get((exp, out))
            if base2 is not None and abs(base2) > 1e-12:
                row["atten_vs_table2"] = 100.0 * (1.0 - row["mvmr_beta_bmd"] / base2)
            else:
                row["atten_vs_table2"] = float("nan")
            row["n_bmi_config"] = n_bmi
            rows.append(row)
            if n_bmi == 100:
                d2 = d.copy(); d2["exposure"] = exp; d2["outcome"] = out
                inst_all.append(d2)
        prim = [r for r in rows if r["exposure"] == exp and r["outcome"] == out
                and r["n_bmi_config"] == 100]
        if prim:
            row = prim[-1]
            print("  %-16s -> %-9s [nBMI=100] K=%3d (BMD %2d / BMI %3d) | gate %s | "
                  "uni_BMD %+.4f(p=%.2e) → MVMR_BMD %+.4f(p=%.2e) | atten_T2 %5.1f%% | "
                  "θ_BMI %+.4f(p=%.2e) | Fc(BMD)=%.1f Fc(BMI)=%.1f"
                  % (exp, out, row["K"], row["n_bmd_iv"], row["n_bmi_iv"],
                     row["gate_recon"], row["uni_beta_bmd"], row["uni_p_bmd"],
                     row["mvmr_beta_bmd"], row["mvmr_p_bmd"], row["atten_vs_table2"],
                     row["mvmr_beta_bmi"], row["mvmr_p_bmi"],
                     row["Fcond_bmd"], row["Fcond_bmi"]))

    pd.DataFrame(rows).to_csv(os.path.join(OUT, "mvmr_bmi_main.csv"),
                              index=False, encoding="utf-8")
    if inst_all:
        pd.concat(inst_all).to_csv(os.path.join(OUT, "mvmr_bmi_instruments.csv"),
                                   index=False, encoding="utf-8")
    print("\nwritten: mvmr_bmi_main.csv (grid n_bmi∈{100,300,800}), "
          "mvmr_bmi_instruments.csv (primary n_bmi=100)")
    print("total rows:", len(rows))


if __name__ == "__main__":
    main()
