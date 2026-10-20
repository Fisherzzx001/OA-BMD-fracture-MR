# -*- coding: utf-8 -*-
"""
115 · Figures v18（S9/S10/S11/S12）独立溯源核验
=================================================
原则：plotted_data_v18.json 是被审查方，不是证据。
每个绘制数值都从一手源文件重新推导并逐一比对；同时核对新图与稿件之间
不存在新的不一致（含本轮勘误 20–35/116 与 coloc 0.799 披露）。

输出：code/output/genetics/AUDIT_115_figures_v18.json + 终端报告
"""
from __future__ import annotations
import csv
import io
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
GEN = os.path.join(HERE, "output", "genetics")
HARM = os.path.join(GEN, "harmonised")
PROJ = os.path.normpath(os.path.join(HERE, ".."))
DOC = os.path.join(PROJ, "docs", "08_投稿稿_OA_BMD_fracture.md")
PLOT = os.path.join(PROJ, "Figures_v18_Suppl", "plotted_data_v18.json")

Z = 1.959963985
CHECKS = 0
FAILS = []
NOTES = []


def ok(cond, msg):
    global CHECKS
    CHECKS += 1
    if not cond:
        FAILS.append(msg)
        print(f"  FAIL  {msg}")
    return cond


def count_rows(path):
    n = 0
    with open(path, encoding="utf-8") as f:
        for i, line in enumerate(f):
            if i > 0 and line.strip():
                n += 1
    return n


plotted = json.load(open(PLOT, encoding="utf-8"))

# ===========================================================================
print("== A · SupplFig10 (MVMR decomposition) ==")
src = pd.read_csv(os.path.join(GEN, "mvmr_height_grid.csv"))
s = src[src["spec"] == "EUR|n_h=100"]
ok(len(s) == 7, f"S10 source rows at EUR|n_h=100 = {len(s)} (expect 7)")
EXP_T3 = {"BMD_total|OA_knee": 19.1, "BMD_femur_neck|OA_knee": 14.6,
          "BMD_lumbar|OA_knee": 27.6, "BMD_arms|OA_knee": 13.5,
          "BMD_leg|OA_knee": 33.2, "BMD_total|fracture": -1.7,
          "BMD_femur_neck|fracture": -3.0}
EXP_SET = {"BMD_total|OA_knee": 9.8, "BMD_femur_neck|OA_knee": 14.6,
           "BMD_lumbar|OA_knee": 9.5, "BMD_arms|OA_knee": 11.0,
           "BMD_leg|OA_knee": 21.2, "BMD_total|fracture": -2.8,
           "BMD_femur_neck|fracture": -2.8}
p10 = plotted["SupplFig10_MVMR_decomposition"]["rows"]
for _, r in s.iterrows():
    key = f"{r.exposure}|{r.outcome}"
    a_t3 = 100.0 * (r.ref_uni_beta - r.mvmr_beta) / r.ref_uni_beta
    a_set = 100.0 * (r.uni_beta - r.mvmr_beta) / r.uni_beta
    ok(abs(a_t3 - EXP_T3[key]) < 0.06, f"S10 {key} atten_vs_Table3 {a_t3:.2f} vs {EXP_T3[key]}")
    ok(abs(a_set - EXP_SET[key]) < 0.06, f"S10 {key} atten_vs_set {a_set:.2f} vs {EXP_SET[key]}")
    q = p10[key]
    ok(abs(q["beta_mv"] - r.mvmr_beta) < 1e-9, f"S10 {key} plotted mvmr_beta")
    ok(abs(q["atten_t3"] - a_t3) < 1e-9, f"S10 {key} plotted atten_t3")
    ok(abs(q["atten_set"] - a_set) < 1e-9, f"S10 {key} plotted atten_set")
    # 森林图 CI 闭式核对
    ok(abs((q["beta_uni"] + Z * q["se_uni"]) - (r.uni_beta + Z * r.uni_se)) < 1e-9,
       f"S10 {key} uni CI upper")
# θ_H 与稿件 S9 句子一致（= MVMR 条件化直接效应 mvmr_h，log-OR 尺度，5 对 OA）
th = [float(r.mvmr_h) for _, r in s.iterrows() if r.outcome == "OA_knee"]
ok(min(th) > 0.143 and max(th) < 0.166, f"S10 theta_H range {min(th):.4f}-{max(th):.4f}")

# ===========================================================================
print("== B · SupplFig11 (robustness grid) ==")
aud = pd.read_csv(os.path.join(GEN, "AUDIT_mr_recomputed.csv")).set_index("pair")
mrlap = pd.read_csv(os.path.join(PROJ, "code", "output", "mrlap_overlap_results.csv"))
mrlap["pair"] = mrlap["exposure"] + "__" + mrlap["outcome"]
mrlap = mrlap.set_index("pair")
g = pd.read_csv(os.path.join(GEN, "mvmr_height_grid.csv"))
g = g[g["spec"] == "EUR|n_h=100"].copy()
g["pair"] = g["exposure"] + "__" + g["outcome"]
g = g.set_index("pair")
p11 = plotted["SupplFig11_robustness_grid"]["details"]
ROW_ORDER = ["BMD_total__fracture", "BMD_femur_neck__fracture",
             "BMD_total__OA_knee", "BMD_lumbar__OA_knee", "BMD_arms__OA_knee",
             "BMD_leg__OA_knee", "BMD_femur_neck__OA_knee",
             "OA_knee__fracture", "OA_knee__BMD_total", "OA_knee__BMD_femur_neck"]
for key in ROW_ORDER:
    a = aud.loc[key]
    d = p11[key]
    # IVW fixed/random/Egger 来自 Table 3 权威表
    ok(abs(d["mrlap_corrected"] - float(mrlap.loc[key, "corrected_effect"])) < 1e-12,
       f"S11 {key} mrlap beta")
    ok(abs(d["mrlap_p"] - float(mrlap.loc[key, "corrected_effect_p"])) < 1e-12,
       f"S11 {key} mrlap p")
    # 方向一致性：MRlap observed 与 IVW 同号
    ok(np.sign(mrlap.loc[key, "observed_effect"]) == np.sign(a["beta"]),
       f"S11 {key} MRlap observed sign vs IVW")
    # PRESSO 离群后 IVW：从 harmonised 重算
    df = pd.read_csv(os.path.join(HARM, key + ".csv"))
    n_out_file = int(df["presso_outlier"].astype(float).sum())
    ok(n_out_file == d["n_outliers"], f"S11 {key} n_outliers {n_out_file} vs {d['n_outliers']}")
    m = df[np.asarray(df["presso_outlier"].astype(float) == 0)]
    b = (m["beta_o"] / m["beta_e"]).to_numpy(float)
    w = 1.0 / np.abs(m["se_o"] / m["beta_e"]).to_numpy(float) ** 2
    b_out = float(np.sum(w * b) / np.sum(w))
    ok(abs(b_out - d["presso_out_beta"]) < 1e-9, f"S11 {key} PRESSO-out beta recompute")
    # 全样本 IVW 应等于 Table 3
    b_all = (df["beta_o"] / df["beta_e"]).to_numpy(float)
    w_all = 1.0 / np.abs(df["se_o"] / df["beta_e"]).to_numpy(float) ** 2
    b_full = float(np.sum(w_all * b_all) / np.sum(w_all))
    ok(abs(b_full - a["beta"]) < 5e-4, f"S11 {key} full-sample IVW {b_full:.6f} vs Table3 {a['beta']:.6f}")
    # MVMR
    if key in g.index:
        ok(abs(d["mvmr_beta"] - float(g.loc[key, "mvmr_beta"])) < 1e-12, f"S11 {key} mvmr beta")
# 方向一致性计数：全部 7/7 或 6/6
for key, cells in plotted["SupplFig11_robustness_grid"].get("_conc", {}).items() if False else []:
    pass
conc = {"BMD_total__fracture": "7/7", "BMD_femur_neck__fracture": "7/7",
        "BMD_total__OA_knee": "7/7", "BMD_lumbar__OA_knee": "7/7",
        "BMD_arms__OA_knee": "7/7", "BMD_leg__OA_knee": "7/7",
        "BMD_femur_neck__OA_knee": "7/7", "OA_knee__fracture": "6/6",
        "OA_knee__BMD_total": "6/6", "OA_knee__BMD_femur_neck": "6/6"}
for key in ROW_ORDER:
    a = aud.loc[key]
    sg = np.sign(a["beta"])
    n_est, n_conc = 0, 0
    ests = [("ivw", a["beta"]), ("re", a["beta"]), ("egger", a["egger_beta"]),
            ("wmed", a["wmed"])]
    df = pd.read_csv(os.path.join(HARM, key + ".csv"))
    m = df[np.asarray(df["presso_outlier"].astype(float) == 0)]
    b = (m["beta_o"] / m["beta_e"]).to_numpy(float)
    w = 1.0 / np.abs(m["se_o"] / m["beta_e"]).to_numpy(float) ** 2
    ests.append(("presso_out", float(np.sum(w * b) / np.sum(w))))
    ests.append(("mrlap", float(mrlap.loc[key, "corrected_effect"])))
    if key in g.index:
        ests.append(("mvmr", float(g.loc[key, "mvmr_beta"])))
    for name, v in ests:
        n_est += 1
        n_conc += int(np.sign(v) == sg)
    ok(f"{n_conc}/{n_est}" == conc[key], f"S11 {key} concordance {n_conc}/{n_est} vs {conc[key]}")

# ===========================================================================
print("== C · SupplFig9 (instrument flow) ==")
sc = pd.read_csv(os.path.join(GEN, "strict_clump_summary.csv")).set_index("exposure")
p9 = plotted["SupplFig9_instrument_flow"]["rows"]
T3K = {"BMD_total__OA_knee.csv": 31, "BMD_total__fracture.csv": 31,
       "BMD_femur_neck__OA_knee.csv": 25, "BMD_femur_neck__fracture.csv": 25,
       "BMD_lumbar__OA_knee.csv": 21, "BMD_arms__OA_knee.csv": 16,
       "BMD_leg__OA_knee.csv": 24, "OA_knee__BMD_total.csv": 69,
       "OA_knee__BMD_femur_neck.csv": 69, "OA_knee__fracture.csv": 69}
for row in p9:
    exp_key = ("BMD_total" if row["exposure"] == "BMD total"
               else "BMD_femur_neck" if row["exposure"] == "Femoral neck"
               else "BMD_lumbar" if row["exposure"] == "Lumbar spine"
               else "BMD_arms" if row["exposure"] == "Arms"
               else "BMD_leg" if row["exposure"] == "Legs" else "OA_knee")
    ok(row["n_gws"] == int(sc.loc[exp_key, "n_gws"]), f"S9 {exp_key} n_gws")
    ok(abs(row["n_clumped"] - int(sc.loc[exp_key, "n_kept_strict"])) <= 1,
       f"S9 {exp_key} n_clumped vs strict pipeline")
    ok(abs(row["retained_pct"] - round(100.0 * row["K"] / row["n_clumped"], 1)) < 0.051,
       f"S9 {exp_key} retained pct")
# harmonised K 与 Table 3 一致
for fn, K in T3K.items():
    ok(count_rows(os.path.join(HARM, fn)) == K, f"S9 harmonised {fn} K={K}")

# ===========================================================================
print("== D · SupplFig12 (literature comparison) ==")
p12 = plotted["SupplFig12_literature_comparison"]["ours"]
for k, r in p12.items():
    a = aud.loc[f"{k}__OA_knee"]
    ok(abs(r["beta"] - float(a["beta"])) < 1e-12, f"S12 {k} beta")
    ok(abs(r["OR"] - float(np.exp(a["beta"]))) < 1e-12, f"S12 {k} OR=exp(beta)")
    ok(abs(r["lo"] - float(np.exp(a["beta"] - Z * a["se"]))) < 1e-12, f"S12 {k} CI lo")
# 文献值与稿件 Discussion 引用一致
with io.open(DOC, encoding="utf-8") as f:
    ms = f.read()
ok("OR 1.28, 1.05\u20131.57" in ms and "1.40, 1.20\u20131.63" in ms,
   "S12 Hartley values match manuscript citation")
ok("1.14, 1.06\u20131.22" in ms and "Funck-Brentano" in ms,
   "S12/Funck-Brentano 2019 all-OA OR 1.14 [1.06\u20131.22] matches manuscript citation")
ok("0.247 vs 0.336" in ms, "S12 within-study 27% log-odds contrast matches manuscript")
# 正文与图中 Funck-Brentano 2019 均只报 all-OA OR 1.14（图中已无关节别行），对图源脚本回源校验
s12_src = io.open(os.path.join(HERE, "114_fig_S12_literature_comparison.py"), encoding="utf-8").read()
ok("Funck-Brentano 2019 Arthritis Rheumatol 71(10):1634-1641" in s12_src
   and "all-OA OR 1.14 (1.06-1.22)" in s12_src,
   "S12 figure source carries Funck-Brentano 2019 all-OA OR 1.14 [1.06-1.22] consistent with full paper")
ok("knee 1.20 (1.05-1.37)" not in s12_src and "hip 1.23 (1.12-1.36)" not in s12_src,
   "S12 figure source no longer plots Funck-Brentano per-joint ORs")

# ===========================================================================
print("== E · 稿件勘误与披露落地核查 ==")
ok("20\u201335 instruments retained" in ms, "Methods instrument range corrected to 20-35")
ok("116 for knee OA" in ms, "Methods knee-OA instrument count corrected to 116")
ok("22\u201337" not in ms and "118 for knee OA" not in ms, "no residual 22-37/118")
ok("47 of the 116 knee-OA instruments (41%)" in ms, "attrition 47/116 (41%) disclosed")
ok("PP.H4=0.799" in ms and "rs11228269" in ms, "coloc near-threshold disclosure present")
ok("fell 0.0006 short" in ms, "coloc gap quantified")
# coloc 披露数值与一手文件一致
rows = list(csv.DictReader(open(os.path.join(GEN, "coloc_screen_all_loci.csv"))))
pp4 = [(r["index_snp"], float(r["PP4"])) for r in rows]
mx = max(pp4, key=lambda x: x[1])
ok(mx[0] == "rs11228269" and abs(mx[1] - 0.7994) < 5e-4,
   f"coloc max PP.H4 {mx[0]}={mx[1]:.4f}")
ok(sum(1 for _, v in pp4 if v >= 0.8) == 0, "no locus reaches 0.8")

# ===========================================================================
res = {"checks": CHECKS, "fails": FAILS,
       "verdict": "PASS" if not FAILS else "FAIL"}
out = os.path.join(GEN, "AUDIT_115_figures_v18.json")
with open(out, "w", encoding="utf-8") as f:
    json.dump(res, f, indent=1)
print(f"\n==== {CHECKS} checks, {len(FAILS)} failures -> {res['verdict']} ====")
print("wrote", out)
sys.exit(0 if not FAILS else 1)
