#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G6 未传播修复 · 补丁 3/13（F3）
=============================
从修正后的 mvmr_height_grid.csv 整体重生 Supplementary Table S10
（docs/Supplementary/02_mvmr_height_grid.md）。

背景：该表生成于 G6 修正之前，OA_knee 行的全部 β 列停留在旧尺度
（×3.3882 之前），且 θ_H 列在 2026-10-20 已核实为 mvmr_h（MVMR 直接
效应），须与稿件 S9 表（θ_H=+0.145~+0.165）收敛。衰减两列是比值，
尺度不变，重算后应与旧表一致（OA_knee 19.1/14.6/27.6/13.5/33.2）。

规则：
  - β_uni (Table 2)   <- ref_uni_beta（F1 已修：OA_knee 新尺度，fracture 原值）
  - β_uni (MVMR set)  <- uni_beta
  - β_MVMR            <- mvmr_beta
  - Atten% vs Table 2 = 100*(ref_uni_beta - mvmr_beta)/ref_uni_beta
  - Atten% vs BMD-subset = 100*(uni_beta - mvmr_beta)/uni_beta
  - θ_H               <- mvmr_h（n_h=0 行不可估计，记 —）
  - Fc(BMD)           <- Fc_bmd；Fc(H) <- Fc_h（n_h=0 记 —）
"""
import csv
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GRID = os.path.join(ROOT, "code", "output", "genetics", "mvmr_height_grid.csv")
OUT = os.path.join(ROOT, "docs", "Supplementary", "02_mvmr_height_grid.md")

EXP_NAME = {
    "BMD_total": "BMD total",
    "BMD_femur_neck": "BMD femur neck",
    "BMD_lumbar": "BMD lumbar",
    "BMD_arms": "BMD arms",
    "BMD_leg": "BMD leg",
}
OUT_NAME = {"OA_knee": "OA knee", "fracture": "fracture"}
GROUP_ORDER = [
    ("BMD_total", "OA_knee"), ("BMD_femur_neck", "OA_knee"),
    ("BMD_lumbar", "OA_knee"), ("BMD_arms", "OA_knee"),
    ("BMD_leg", "OA_knee"), ("BMD_total", "fracture"),
    ("BMD_femur_neck", "fracture"),
]
NH_ORDER = {"0": 0, "100": 1, "300": 2, "1000": 3}
HS_ORDER = {"EUR": 0, "EUR_noUKB": 1}

HEADER = ("# Supplementary Table S10. Multivariable MR grid — BMD→OA knee and "
          "BMD→fracture, conditioning on adult height\n\n"
          "Full sensitivity grid underlying the primary n_h=100 / EUR configuration reported in "
          "Table S9. Instruments were selected at p<5×10⁻⁸ and clumped at r²<0.001/10 Mb; the MVMR "
          "set overlaps the primary strict-clumped set except that six BMD instruments are absent "
          "from the height GWAS (see Table S9 note). `β_uni (Table 2)` is the univariable IVW "
          "estimate from the main Table 2; `β_uni (MVMR set)` is the univariable estimate "
          "recomputed in the MVMR pipeline; `Atten% vs Table 2` = 100×(β_uni(Table 2) − "
          "β_MVMR)/β_uni(Table 2); `Atten% vs BMD-subset` = 100×(β_uni(MVMR set) − "
          "β_MVMR)/β_uni(MVMR set). θ_H is the direct effect of height (MVMR); Fc = "
          "Sanderson–Smith conditional F (≥10 adequate). β per 1 SD on the N-corrected scale.\n\n")

COLS = ("| Exposure→Outcome | Height src | n_h (req) | n_bmd | β_uni (Table 2) | "
        "β_uni (MVMR set) | β_MVMR | Atten% vs Table 2 | Atten% vs BMD-subset | θ_H | "
        "Fc(BMD) | Fc(H) |\n"
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |\n")

NOTES = (
    "\n\n---\n\n"
    "**Note on the attenuation columns (corrected 2026-10-16; β columns re-scaled 2026-10-20 "
    "with the G6 N·K(1−K) correction).** `Atten% vs Table 2` = 100×(β_uni(Table 2) − "
    "β_MVMR)/β_uni(Table 2); `Atten% vs BMD-subset` = 100×(β_uni(MVMR set) − "
    "β_MVMR)/β_uni(MVMR set). Both are ratios and therefore invariant to the N-scale "
    "correction; they agree with Table S9 for the primary EUR / n_h=100 configuration "
    "(19.1 / 14.6 / 27.6 / 13.5 / 33.2% vs Table 2 and 9.8 / 14.6 / 9.5 / 11.0 / 21.2% vs "
    "the MVMR-set baseline). **Negative values mean the MVMR estimate moved _away_ from the "
    "null relative to that baseline** (i.e. the effect strengthened under conditioning), "
    "which occurs in the degenerate n_h=0 rows where no height instrument is conditioned "
    "on but the restricted instrument subset changes the univariable reference; the primary "
    "n_h=100 rows are positive throughout for the OA-knee outcomes.\n\n"
    "> **Note (2026-10-20 audit).** For the n_h=0 rows the MVMR model contains no height "
    "instruments (K = n_bmd), so no direct height effect θ_H or conditional F for height "
    "is estimable; those cells are shown as '—'. Values previously printed in these cells "
    "were rank-deficiency artefacts of the regression and are withdrawn. The n_h=0 rows "
    "are retained only as the BMD-only univariable baseline.\n"
)


def main():
    with open(GRID, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    by_key = {}
    for r in rows:
        by_key[(r["exposure"], r["outcome"], r["height"], r["n_h"])] = r

    out = [HEADER, COLS]
    n_rows = 0
    for exp, out_ in GROUP_ORDER:
        for hsrc in ("EUR", "EUR_noUKB"):
            for nh in ("0", "100", "300", "1000"):
                r = by_key[(exp, out_, hsrc, nh)]
                ref = float(r["ref_uni_beta"])
                uni = float(r["uni_beta"])
                mv = float(r["mvmr_beta"])
                a_t2 = 100.0 * (ref - mv) / ref
                a_set = 100.0 * (uni - mv) / uni
                if nh == "0":
                    theta, fch = "—", "—"
                else:
                    theta = f"{float(r['mvmr_h']):.4f}"
                    fch = f"{float(r['Fc_h']):.1f}"
                cells = [
                    f"{EXP_NAME[exp]} → {OUT_NAME[out_]}",
                    hsrc, nh, r["n_bmd"],
                    f"{ref:.5f}", f"{uni:.5f}", f"{mv:.5f}",
                    f"{a_t2:.1f}", f"{a_set:.1f}",
                    theta, f"{float(r['Fc_bmd']):.1f}", fch,
                ]
                out.append("| " + " | ".join(cells) + " |\n")
                n_rows += 1
    out.append(NOTES)

    with open(OUT, "w", encoding="utf-8") as f:
        f.write("".join(out))
    print(f"S10 重生完成：{n_rows} 行 -> {OUT}")

    # ---- 自验：主配置 vs 稿件 S9 ----
    print("\n自验（EUR n_h=100，OA_knee）：")
    exp_t2 = {"BMD_total": 19.1, "BMD_femur_neck": 14.6, "BMD_lumbar": 27.6,
              "BMD_arms": 13.5, "BMD_leg": 33.2}
    exp_set = {"BMD_total": 9.8, "BMD_femur_neck": 14.6, "BMD_lumbar": 9.5,
               "BMD_arms": 11.0, "BMD_leg": 21.2}
    exp_theta = {"BMD_total": 0.1445, "BMD_femur_neck": 0.1646, "BMD_lumbar": 0.1529,
                 "BMD_arms": 0.1585, "BMD_leg": 0.1522}
    ok = True
    for exp in exp_t2:
        r = by_key[(exp, "OA_knee", "EUR", "100")]
        a2 = 100.0 * (float(r["ref_uni_beta"]) - float(r["mvmr_beta"])) / float(r["ref_uni_beta"])
        as_ = 100.0 * (float(r["uni_beta"]) - float(r["mvmr_beta"])) / float(r["uni_beta"])
        th = float(r["mvmr_h"])
        good = (abs(a2 - exp_t2[exp]) < 0.06 and abs(as_ - exp_set[exp]) < 0.06
                and abs(th - exp_theta[exp]) < 0.0006)
        ok &= good
        print(f"  {exp:18s} a_t2={a2:5.1f}% (期望 {exp_t2[exp]})  a_set={as_:5.1f}% "
              f"(期望 {exp_set[exp]})  θ_H={th:.4f} (期望 {exp_theta[exp]})  {'✓' if good else '✗'}")
    #  fracture 主配置
    r = by_key[("BMD_total", "fracture", "EUR", "100")]
    a2 = 100.0 * (float(r["ref_uni_beta"]) - float(r["mvmr_beta"])) / float(r["ref_uni_beta"])
    print(f"  fracture         a_t2={a2:5.1f}% (期望 -1.7)  θ_H={float(r['mvmr_h']):.4f} (期望 0.1663)")
    ok &= abs(a2 - (-1.7)) < 0.06 and abs(float(r["mvmr_h"]) - 0.1663) < 0.0006
    assert ok, "S10 主配置自验失败"
    print("\n✅ F3 完成：S10 与稿件 S9 收敛")


if __name__ == "__main__":
    main()
