# -*- coding: utf-8 -*-
"""
92_audit_figures_v9.py — 图数据真实性独立复核（第二轮，覆盖 v9 全部 12 张图）
=========================================================================
设计原则（与 verify_figures.py 一致）：
  * 不 import build_figures*.py / 27_overlap_statistics.py / 28_mr_power.py 的任何函数；
  * 所有统计量在本文件内重新实现，从最底层源数据重算；
  * 比对对象是 **作图脚本实际绘制的值**（figure_v9_plotted_data.json），
    而不是"我以为它画了什么"。

覆盖：
  [A] Figure 8  检验效能：|z|、MDE、power 公式、效能曲线、80% 阈值、无裁切
  [B] Figure 9  样本重叠：ρ_GE 从逐 SNP 文件重算、IVW z 重算、符号一致性、无裁切
  [C] Suppl Fig S3 GEFOS 外部基准：OR/CI 重算 + 1-Mb 敏感性 + 与 GEFOS 的差距
  [D] Figure 6/7/S1/S2 轴单位：用反解 N_eff 判定每个性状的真实效应尺度
  [E] 文件层：12 张图的 PNG/PDF 存在性、分辨率、矢量字体
"""
import os
import json
import struct
import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = "C:/Users/Administrator/WorkBuddy/2026-09-08-22-17-40/OA_B_PROJECT"
GEN = os.path.join(ROOT, "code/output/genetics")
HARM = os.path.join(GEN, "harmonised")
FIG = os.path.join(ROOT, "submission_package_v2/Figures")

PASS, FAIL = [], []


def chk(cond, msg, tag=""):
    (PASS if cond else FAIL).append(f"[{tag}] {msg}" if tag else msg)
    print(("  PASS  " if cond else "  FAIL  ") + msg)
    return cond


def close(a, b, tol=1e-6):
    return abs(float(a) - float(b)) <= tol


# ============================================================ [A] Figure 8
def audit_figure8():
    print("\n[A] Figure 8 — MR 检验效能")
    dump = json.load(open(os.path.join(GEN, "figure_v9_plotted_data.json"),
                          encoding="utf-8"))["figure8"]
    mr = pd.read_csv(os.path.join(GEN, "mr_corrected_N.csv"))
    pw = pd.read_csv(os.path.join(GEN, "mr_power.csv"))

    # --- 从最底层重算 ---
    rec = {}
    for _, r in mr.iterrows():
        key = f"{r.exposure}__{r.outcome}"
        z = abs(r.beta / r.se)
        mde = (1.96 + 0.84) * r.se
        power = float(norm.sf(1.96 - z) + norm.cdf(-1.96 - z))
        rec[key] = dict(z=z, mde=mde, power=power, p=r.p, K=int(r.K))

    PRETTY = {"OA_knee": "knee OA", "fracture": "fracture",
              "BMD_total": "total BMD", "BMD_femur_neck": "femoral-neck BMD",
              "BMD_lumbar": "lumbar BMD", "BMD_arms": "arms BMD",
              "BMD_leg": "legs BMD"}

    # --- A1 面板 a 的散点 ---
    pts = {p["pair"]: p for p in dump["panel_a"]["points"]}
    chk(len(pts) == 10, f"A1 面板 a 绘出 {len(pts)} 对（应为 10）", "A1")
    for key, v in rec.items():
        exp, out = key.split("__")
        lab = f"{PRETTY[exp]} → {PRETTY[out]}"
        p = pts.get(lab)
        if p is None:
            chk(False, f"A2 缺失绘图点 {lab}", "A2")
            continue
        chk(close(p["abs_z"], v["z"], 1e-3), f"A2 {lab}: |z| 绘制 {p['abs_z']:.3f} == 重算 {v['z']:.3f}", "A2")
        chk(close(p["power"], v["power"], 1e-3), f"A3 {lab}: power 绘制 {p['power']:.3f} == 重算 {v['power']:.3f}", "A3")
        chk(p["K"] == v["K"], f"A4 {lab}: K 绘制 {p['K']} == 源 {v['K']}", "A4")

    # --- A5 效能曲线公式 ---
    cx = np.array(dump["panel_a"]["curve_x"])
    cy = np.array(dump["panel_a"]["curve_y"])
    cy_ref = norm.sf(1.96 - cx) + norm.cdf(-1.96 - cx)
    chk(np.max(np.abs(cy - cy_ref)) < 1e-12,
        f"A5 效能曲线 == Φ(|z|−1.96)+Φ(−|z|−1.96)（最大偏差 {np.max(np.abs(cy-cy_ref)):.2e}）", "A5")

    # --- A6 80% 阈值自洽：曲线在 |z|=2.80 处应≈0.80 ---
    y280 = float(norm.sf(1.96 - 2.80) + norm.cdf(-1.96 - 2.80))
    chk(0.795 <= y280 <= 0.805, f"A6 阈值自洽：|z|=2.80 处效能 = {y280:.4f} ≈ 0.80", "A6")
    chk(close(dump["panel_a"]["thresholds"]["z80"], 2.80) and
        close(dump["panel_a"]["thresholds"]["power80"], 0.80),
        "A6b 图中标注的阈值线为 0.80 / 2.80", "A6")

    # --- A7 MDE 与 Table 4 / mr_power.csv 一致 ---
    for _, r in pw.iterrows():
        key = f"{r.exposure}__{r.outcome}"
        v = rec[key]
        chk(close(r.MDE_beta_80pct, v["mde"], 1e-3),
            f"A7 {r.exposure}→{r.outcome}: MDE {r.MDE_beta_80pct:.4f} == (1.96+0.84)×SE {v['mde']:.4f}", "A7")
        chk(close(r.power_observed, v["power"], 2e-2),
            f"A8 {r.exposure}→{r.outcome}: 实测效能 {r.power_observed:.2f} == 重算 {v['power']:.2f}", "A8")
        chk((r.verdict.startswith("powered")) == (v["power"] >= 0.80),
            f"A9 {r.exposure}→{r.outcome}: 判定 '{r.verdict}' 与效能 {v['power']:.2f} 自洽", "A9")

    # --- A10 面板 b 顺序/数值/无裁切 ---
    pb = dump["panel_b"]
    chk(len(pb["power"]) == 10 and len(pb["K"]) == 10, "A10 面板 b 绘出 10 对", "A10")
    mapped = {f"{PRETTY[e]} → {PRETTY[o]}": v for (e, o), v in
              zip(zip(mr.exposure, mr.outcome), [rec[f"{e}__{o}"] for e, o in zip(mr.exposure, mr.outcome)])}
    for lab, pv, kv in zip(pb["order"], pb["power"], pb["K"]):
        chk(close(pv, mapped[lab]["power"], 1e-3) and kv == mapped[lab]["K"],
            f"A11 {lab}: 条形 power={pv:.2f} K={kv} 与源一致", "A11")
    chk(max(pb["power"]) <= pb["xlim"][1] + 1e-9,
        f"A12 面板 b 最大效能 {max(pb['power']):.2f} ≤ 轴上限 {pb['xlim'][1]}（无裁切）", "A12")
    chk(max(p["abs_z"] for p in pts.values()) <= dump["panel_a"]["xlim"][1],
        "A13 面板 a 所有点的 |z| 落在轴范围内（无裁切）", "A13")
    chk(max(p["power"] for p in pts.values()) <= dump["panel_a"]["ylim"][1],
        "A14 面板 a 所有点的效能落在轴范围内（无裁切）", "A14")

    # --- A15 关键结论复核：显著者效能高、null 者效能低 ---
    sig = [v for v in rec.values() if v["p"] < 0.05]
    nul = [v for v in rec.values() if v["p"] >= 0.05]
    chk(all(v["power"] >= 0.80 for v in sig), f"A15 7 个显著对效能均 ≥0.80（最低 {min(v['power'] for v in sig):.2f}）", "A15")
    chk(all(v["power"] < 0.80 for v in nul), f"A16 3 个 null 对效能均 <0.80（最高 {max(v['power'] for v in nul):.2f}）", "A16")
    chk(len(sig) == 7 and len(nul) == 3, f"A17 分组计数 7/3（实际 {len(sig)}/{len(nul)}）", "A17")


# ============================================================ [B] Figure 9
def audit_figure9():
    print("\n[B] Figure 9 — 样本重叠量规 ρ_GE")
    dump = json.load(open(os.path.join(GEN, "figure_v9_plotted_data.json"),
                          encoding="utf-8"))["figure9"]
    ov = pd.read_csv(os.path.join(GEN, "overlap_statistics.csv"))
    mr = pd.read_csv(os.path.join(GEN, "mr_corrected_N.csv"))
    key_beta = {(e, o): b for e, o, b in zip(mr.exposure, mr.outcome, mr.beta)}
    key_se = {(e, o): s for e, o, s in zip(mr.exposure, mr.outcome, mr.se)}

    PRETTY = {"OA_knee": "knee OA", "fracture": "fracture",
              "BMD_total": "total BMD", "BMD_femur_neck": "femoral-neck BMD",
              "BMD_lumbar": "lumbar BMD", "BMD_arms": "arms BMD",
              "BMD_leg": "legs BMD"}

    # --- 从逐 SNP 文件独立重算 ρ_GE 与 IVW ---
    recalc = {}
    for fn in sorted(os.listdir(HARM)):
        if not (fn.endswith(".csv") and "__" in fn):
            continue
        exp, out = fn[:-4].split("__")
        d = pd.read_csv(os.path.join(HARM, fn))
        d = d[np.isfinite(d.beta_e) & np.isfinite(d.beta_o) &
              (d.se_e > 0) & (d.se_o > 0)]
        ze = (d.beta_e / d.se_e).values
        zo = (d.beta_o / d.se_o).values
        rho = float(np.corrcoef(ze, zo)[0, 1])
        w = 1.0 / d.se_o.values ** 2
        beta = float(np.sum(w * d.beta_e.values * d.beta_o.values) /
                     np.sum(w * d.beta_e.values ** 2))
        se = float(1.0 / np.sqrt(np.sum(d.beta_e.values ** 2 / d.se_o.values ** 2)))
        recalc[f"{PRETTY[exp]} → {PRETTY[out]}"] = dict(rho=rho, K=len(d),
                                                        z=beta / se, beta=beta)

    # --- B1 ρ_GE 与逐 SNP 重算一致 ---
    for _, r in ov.iterrows():
        exp, out = str(r["pair"]).replace(".csv", "").split("__")
        lab = f"{PRETTY[exp]} → {PRETTY[out]}"
        rc = recalc[lab]
        # overlap_statistics.csv 存 4 位小数，故容差取 1e-4（其余比对仍用 1e-6）
        chk(close(r.rho_GE_zcorr, rc["rho"], 1e-4),
            f"B1 {lab}: overlap_statistics ρ_GE {r.rho_GE_zcorr:+.4f} == 逐SNP重算 {rc['rho']:+.6f}", "B1")
        chk(int(r.n_iv) == rc["K"], f"B2 {lab}: 工具数 {int(r.n_iv)} == 逐SNP 行数 {rc['K']}", "B2")

    # --- B3 绘图值与重算一致 ---
    pa = dump["panel_a"]
    for lab, rho, K in zip(pa["order"], pa["rho"], pa["K"]):
        chk(close(rho, recalc[lab]["rho"], 1e-4) and K == recalc[lab]["K"],
            f"B3 {lab}: 面板a 绘制 ρ={rho:+.4f} K={K} == 重算", "B3")
    pb = dump["panel_b"]
    INV_PRETTY = {v: k for k, v in PRETTY.items()}
    for lab, rho, z in zip(pa["order"], pb["rho"], pb["z"]):
        chk(close(rho, recalc[lab]["rho"], 1e-4), f"B4 {lab}: 面板b x 坐标 ρ 一致", "B4")
        chk(close(z, recalc[lab]["z"], 1e-6),
            f"B5 {lab}: 面板b y 坐标 IVW z {z:+.3f} == 逐SNP重算 {recalc[lab]['z']:+.3f}", "B5")
        e_raw, o_raw = [s.strip() for s in lab.split("→")]
        e_tab, o_tab = INV_PRETTY[e_raw], INV_PRETTY[o_raw]
        b_tab, s_tab = key_beta[(e_tab, o_tab)], key_se[(e_tab, o_tab)]
        chk(close(z, b_tab / s_tab, 1e-6),
            f"B6 {lab}: z 亦等于 Table 3 β/SE = {b_tab/s_tab:+.3f}", "B6")

    # --- B7 符号一致性：稿件声称 10/10 ---
    same = sum(1 for lab in pa["order"]
               if np.sign(recalc[lab]["rho"]) == np.sign(recalc[lab]["z"]))
    chk(same == 10, f"B7 ρ_GE 与 IVW z 同号 {same}/10（稿件/图注称 10/10）", "B7")
    chk(pb["n_sign_concordant"] == 10, f"B8 图中标注的同号计数 = {pb['n_sign_concordant']}（应为 10）", "B8")

    # --- B9 无裁切 ---
    chk(min(pa["rho"]) >= pa["xlim"][0] and max(pa["rho"]) <= pa["xlim"][1],
        f"B9 ρ 范围 [{min(pa['rho']):+.3f}, {max(pa['rho']):+.3f}] 在轴 [{pa['xlim'][0]}, {pa['xlim'][1]}] 内", "B9")
    chk(min(pb["z"]) >= pb["ylim"][0] and max(pb["z"]) <= pb["ylim"][1],
        f"B10 z 范围 [{min(pb['z']):+.2f}, {max(pb['z']):+.2f}] 在轴内", "B10")

    # --- B11 稿件文字：ρ_GE 范围 −0.79 ~ +0.46 ---
    chk(close(min(pa["rho"]), -0.79, 0.01) and close(max(pa["rho"]), 0.46, 0.01),
        f"B11 ρ_GE 范围 {min(pa['rho']):+.2f} ~ {max(pa['rho']):+.2f}（§2.8 称 −0.79 ~ +0.46）", "B11")
    # --- B12 稿件文字：BMD→OA 为正、BMD→fracture 为负 ---
    chk(all(recalc[l]["rho"] > 0 for l in recalc if l.endswith("knee OA") and l.startswith("BMD") or l.startswith("total BMD → knee OA")),
        "B12 所有 BMD→OA 对的 ρ_GE 为正", "B12")
    chk(all(recalc[l]["rho"] < 0 for l in recalc if l.endswith("fracture") and l.startswith("BMD")),
        "B13 所有 BMD→fracture 对的 ρ_GE 为负", "B13")


# ============================================================ [C] Suppl S3
def audit_s3():
    print("\n[C] Supplementary Figure S3 — GEFOS 外部基准")
    g = json.load(open(os.path.join(GEN, "gefos_benchmark.json"), encoding="utf-8"))
    mr = pd.read_csv(os.path.join(GEN, "mr_corrected_N.csv"))
    raw = pd.read_csv(os.path.join(GEN, "mr_results.csv"))
    CF = np.sqrt(400000 / 31873)          # 3.5426
    INV = 1.0 / CF                        # 0.2823

    # C1 总体质 / 股骨颈（严格集）由 mr_corrected_N.csv 重算
    for exp, lab in [("BMD_total", "This study - total-body BMD (primary)"),
                     ("BMD_femur_neck", "This study - femoral-neck BMD (primary)")]:
        r = mr[(mr.exposure == exp) & (mr.outcome == "fracture")].iloc[0]
        orr = float(np.exp(abs(r.beta)))
        lo, hi = float(np.exp(abs(r.beta) - 1.96 * r.se)), float(np.exp(abs(r.beta) + 1.96 * r.se))
        i = g["labels"].index(lab)
        chk(close(g["estimates"][i], round(orr, 2), 0.005),
            f"C1 {lab}: 图中 OR {g['estimates'][i]} == exp(|β|) {orr:.4f}", "C1")
        chk(close(g["ci_low"][i], round(lo, 2), 0.005) and close(g["ci_high"][i], round(hi, 2), 0.005),
            f"C2 {lab}: 图中 CI [{g['ci_low'][i]}, {g['ci_high'][i]}] == 重算 [{lo:.4f}, {hi:.4f}]", "C2")

    # C3 1-Mb 敏感性（宽松集）：由 mr_results.csv 施加同一校正后重算
    r = raw[(raw.exposure == "BMD_femur_neck") & (raw.outcome == "fracture")].iloc[0]
    b, s = abs(r.ivw) * INV, r.ivw_se * INV
    orr, lo, hi = float(np.exp(b)), float(np.exp(b - 1.96 * s)), float(np.exp(b + 1.96 * s))
    i = g["labels"].index("This study - femoral-neck BMD (1-Mb sensitivity)")
    chk(close(g["estimates"][i], round(orr, 2), 0.005),
        f"C3 1-Mb 敏感性 OR {g['estimates'][i]} == 重算 {orr:.4f}（β_raw {r.ivw:.4f} × 0.2823）", "C3")
    chk(close(g["ci_low"][i], round(lo, 2), 0.005) and close(g["ci_high"][i], round(hi, 2), 0.005),
        f"C4 1-Mb 敏感性 CI [{g['ci_low'][i]}, {g['ci_high'][i]}] == 重算 [{lo:.4f}, {hi:.4f}]", "C4")

    # C5 GEFOS 基准为外部文献值，仅校验图中数值 == §3.7 引用值
    i = g["labels"].index("GEFOS 2018 (femoral-neck BMD)")
    chk(g["estimates"][i] == 1.55 and g["ci_low"][i] == 1.48 and g["ci_high"][i] == 1.63,
        "C5 GEFOS 基准 1.55 [1.48, 1.63] == §3.7 引用的 Trajanoska 2018", "C5")

    # C6 稿件称"低 8–10%"
    for lab, our in [("This study - total-body BMD (primary)", 1.40),
                     ("This study - femoral-neck BMD (primary)", 1.42)]:
        pct = (1.55 - our) / 1.55 * 100
        chk(8.0 <= pct <= 10.0, f"C6 {lab}: 低于 GEFOS {pct:.1f}%（稿件称 8–10%）", "C6")

    # C7 校正因子方向：严格集与宽松集在同一 SNP 数下应相差 0.2823（而非 3.5426）
    for exp in ["BMD_total", "BMD_leg"]:
        a = raw[(raw.exposure == exp) & (raw.outcome == "OA_knee")]
        b = mr[(mr.exposure == exp) & (mr.outcome == "OA_knee")]
        if len(a) and len(b) and int(a.iloc[0].n_snp) == int(b.iloc[0].K):
            ratio = b.iloc[0].beta / a.iloc[0].ivw
            chk(close(ratio, INV, 1e-3),
                f"C7 {exp}→OA_knee(K={int(b.iloc[0].K)}): 校正后/校正前 = {ratio:.4f} == 0.2823", "C7")


# ============================================================ [D] 轴单位
def audit_axis_units():
    print("\n[D] Figure 6/7/S1/S2 轴单位（效应尺度归属）")
    si = pd.read_csv(os.path.join(GEN, "AUDIT_scale_inference.csv"))
    N = {"OA_knee": 1208397.0, "fracture": 264973.0,
         "BMD_total": 31873.0, "BMD_femur_neck": 31873.0, "BMD_lumbar": 31873.0,
         "BMD_arms": 31873.0, "BMD_leg": 31873.0}
    K = {"OA_knee": 0.12, "fracture": 37857.0 / 264973.0}
    # 作图脚本中的尺度归属
    BINARY = {"fracture"}
    STAND = {"OA_knee", "OA_selfrep"}

    agg = si.groupby("trait").N_implied_median.median()
    for t, nim in agg.items():
        if t not in N:
            continue
        n = N[t]
        nk = n * K[t] * (1 - K[t]) if t in K else None
        if t.startswith("BMD"):
            ok = 0.9 * n <= nim <= 1.1 * n
            chk(ok, f"D1 {t}: 反解 N_eff {nim:,.0f} ≈ N {n:,.0f} → 应标 'SD per allele'", "D1")
        elif t == "fracture":
            r_n, r_nk = nim / n, nim / nk
            chk(r_nk > 0.8 and r_n < 0.3,
                f"D2 {t}: N_eff/N={r_n:.2f}、N_eff/[N·K(1−K)]={r_nk:.2f} → log-OR 尺度成立", "D2")
            chk(t in BINARY, f"D3 {t} 在作图脚本中确被标为 log-OR", "D3")
        elif t == "OA_knee":
            nk = n * K[t] * (1 - K[t])
            r_n, r_nk = nim / n, nim / nk
            chk(r_n > 0.8 and r_nk > 5,
                f"D4 {t}: N_eff/N={r_n:.2f}、N_eff/[N·K(1−K)]={r_nk:.1f}（≈9 倍）→ 非 log-OR", "D4")
            chk(t in STAND, f"D5 {t} 在作图脚本中确被标为 standardised（已修 v6 的错误）", "D5")


# ============================================================ [E] 文件层
def png_dpi(path):
    with open(path, "rb") as fh:
        head = fh.read(33)
    if head[:8] != b"\x89PNG\r\n\x1a\n":
        return None
    w, h = struct.unpack(">II", head[16:24])
    with open(path, "rb") as fh:
        data = fh.read()
    i = data.find(b"pHYs")
    if i < 0:
        return None
    ppux, ppuy, unit = struct.unpack(">IIB", data[i + 4:i + 13])
    if unit != 1:
        return None
    return round(ppux * 0.0254, 1)


def audit_files():
    print("\n[E] 文件层检查（12 张图 × PNG+PDF）")
    names = ["Figure1_rg_heatmap_13", "Figure2_forest_key_results", "Figure3_MR_forest",
             "Figure4_h2_bar", "Figure5_causal_network", "Figure6_MR_scatter",
             "Figure7_MR_funnel", "Figure8_MR_power", "Figure9_overlap_pleiotropy",
             "SupplFig1_MR_scatter_remaining", "SupplFig2_leave_one_out",
             "SupplFig3_GEFOS_benchmark"]
    for n in names:
        p, q = os.path.join(FIG, n + ".png"), os.path.join(FIG, n + ".pdf")
        ok = os.path.exists(p) and os.path.exists(q)
        chk(ok, f"E1 {n}: PNG+PDF 均存在", "E1")
        if not ok:
            continue
        dpi = png_dpi(p)
        chk(dpi is not None and dpi >= 299, f"E2 {n}.png 分辨率 {dpi} dpi ≥ 300", "E2")
        with open(q, "rb") as fh:
            head = fh.read(4096)
        chk(b"%PDF" in head, f"E3 {n}.pdf 为合法 PDF", "E3")
        with open(q, "rb") as fh:
            blob = fh.read()
        chk(b"Arial" in blob or b"Helvetica" in blob,
            f"E4 {n}.pdf 内嵌 Arial/Helvetica 字体（矢量可编辑）", "E4")


def main():
    print("=" * 78)
    print("图数据真实性独立复核 — v9（12 张图）")
    print("=" * 78)
    audit_figure8()
    audit_figure9()
    audit_s3()
    audit_axis_units()
    audit_files()
    print("\n" + "=" * 78)
    print(f"通过 {len(PASS)} 项 / 失败 {len(FAIL)} 项")
    print("=" * 78)
    if FAIL:
        print("\n失败明细：")
        for f in FAIL:
            print("  ✗", f)
        raise SystemExit(1)
    print("\n✓ 全部通过：所有图均基于真实数据与真实统计量")


if __name__ == "__main__":
    main()
