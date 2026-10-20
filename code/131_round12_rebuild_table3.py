# -*- coding: utf-8 -*-
"""
131 · 第十二轮整改：重建 Table 3
------------------------------------------------------------------
P1-1 加 OR (95% CI) 列（BMD→OA 每 SD 升高；BMD→fracture 每 SD 降低），
     并统一效应尺度表述（删 "standardised outcome scale"）；
P1-3 加 MR-Egger 斜率的 SE 与 p、weighted median 的 bootstrap 95% CI；
     同时把 weighted median 列换为可从沉积工具表复算的值（见 129 脚本结论：
     IVW 与 MR-Egger 逐位可复现，而已列表的 WM 有 7/10 对不可复现）。

数值来源：code/output/genetics/AUDIT_mr_recomputed.csv（IVW/Egger/Q/I²/p，
与沉积 harmonised 表逐位一致）+ code/output/genetics/wm_recomputed_vs_tabulated.csv
（WM 复算值 + bootstrap CI）+ 本脚本按 exp(β ± 1.96·SE) 计算的 OR/CI。
"""
import io
import math
import os
import re
import pandas as pd
from scipy.stats import norm

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MD = os.path.join(ROOT, "docs", "08_投稿稿_OA_BMD_fracture.md")
AUDIT = os.path.join(ROOT, "code", "output", "genetics", "AUDIT_mr_recomputed.csv")
WM = os.path.join(ROOT, "code", "output", "genetics", "wm_recomputed_vs_tabulated.csv")

SUP = {"-": "⁻", "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴",
       "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹"}

# 展示顺序与标签（与现 Table 3 完全一致）
ORDER = [
    ("BMD_total__fracture", "**BMD total → fracture** (positive control)"),
    ("BMD_femur_neck__fracture", "BMD femoral neck → fracture"),
    ("BMD_total__OA_knee", "**BMD total → OA knee**"),
    ("BMD_lumbar__OA_knee", "Lumbar spine → OA knee"),
    ("BMD_arms__OA_knee", "Arms → OA knee"),
    ("BMD_leg__OA_knee", "Legs → OA knee"),
    ("BMD_femur_neck__OA_knee", "Femoral neck → OA knee"),
    ("OA_knee__fracture", "OA knee → fracture"),
    ("OA_knee__BMD_total", "OA knee → BMD total"),
    ("OA_knee__BMD_femur_neck", "OA knee → BMD femoral neck"),
]


def sup(i: int) -> str:
    return "".join(SUP[c] for c in str(i))


def fmt_p(p: float) -> str:
    if p < 0.01:
        e = math.floor(math.log10(p))
        m = p / 10.0 ** e
        return f"{m:.1f}×10{sup(e)}"
    return f"{p:.3f}"


def fmt_b(b: float, d: int = 3) -> str:
    s = f"{b:+.{d}f}"
    return s.replace("+", "+").replace("-", "−") if False else s


def signed(b: float, d: int = 3) -> str:
    """按稿件风格：负号用 U+2212，正号保留 +；舍入为 0 时不带符号。"""
    v = f"{abs(b):.{d}f}"
    if float(v) == 0.0:
        return v
    return ("−" if b < 0 else "+") + v


def main():
    a = pd.read_csv(AUDIT).set_index("pair")
    w = pd.read_csv(WM).set_index("pair")

    rows = []
    for pair, label in ORDER:
        r = a.loc[pair]
        wm = w.loc[pair]
        beta, se = float(r["beta"]), float(r["se"])
        # OR：BMD→OA 每 SD 升高；BMD→fracture 每 SD 降低（高 BMD 为保护）
        if pair.endswith("OA_knee") and pair.startswith("BMD"):
            orv, lo, hi = math.exp(beta), math.exp(beta - 1.96 * se), math.exp(beta + 1.96 * se)
            or_txt = f"{orv:.3f} ({lo:.3f}–{hi:.3f})"
        elif pair.endswith("fracture") and pair.startswith("BMD"):
            orv, lo, hi = math.exp(-beta), math.exp(-beta - 1.96 * se), math.exp(-beta + 1.96 * se)
            or_txt = f"{orv:.2f} ({lo:.2f}–{hi:.2f})"
        else:
            or_txt = "—"
        eg_se = float(r["egger_se"])
        eg_p = float(2 * norm.sf(abs(float(r["egger_beta"]) / eg_se)))
        rows.append(
            f"| {label} | {int(r['K'])} | {signed(beta)} | {se:.3f} | {fmt_p(float(r['p']))} "
            f"| {fmt_p(float(r['p_re']))} | {or_txt} "
            f"| {signed(float(r['egger_beta']))} ({eg_se:.3f}) | {fmt_p(eg_p)} "
            f"| {fmt_p(float(r['egger_int_p']))} "
            f"| {signed(float(wm['wm_recomp']))} ({signed(float(wm['wm_lo']))} to "
            f"{signed(float(wm['wm_hi']))}) | {float(r['I2'])*100:.0f}% "
            f"| {int(r['loo_flips'])} |"
        )

    header = ("| Exposure → Outcome | n_snp | IVW β (fixed) | SE | p | IVW p (random) "
              "| OR (95% CI) | MR-Egger β (SE) | Egger p | Egger intercept p "
              "| Weighted median (95% CI) | I² | LOO sign flips |")
    sep = ("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    body = "\n".join(rows)

    footnote = (
        "Univariate instrument strength (F, R²_exposure) and the Steiger directionality test for "
        "every pair are in Table S13. β is the IVW estimate per 1 SD higher exposure; for the "
        "binary outcomes β is on the log-odds scale and OR = exp(|β|) with CI exp(β ± 1.96 SE), "
        "quoted per 1 SD *higher* BMD for the BMD→OA pairs and per 1 SD *lower* BMD for the "
        "BMD→fracture pairs, where higher BMD is protective; no OR is given for the reverse "
        "pairs, whose exposure is a binary liability rather than a 1-SD continuous exposure. "
        "MR-Egger p is a two-sided z-test on the slope. Weighted median is Bowden's estimator "
        "with weights (β_exposure/SE_outcome)², recomputed from the deposited harmonised "
        "instrument tables with 10,000-replicate bootstrap CIs (code/129); the IVW, MR-Egger, "
        "Cochran Q and I² columns reproduce from the same deposit to machine precision, and "
        "seven weighted-median values previously tabulated did not, so that column has been "
        "corrected (audit trail: "
        "code/output/genetics/wm_recomputed_vs_tabulated.csv). Forward BMD-exposure betas are "
        "per-SD on the verified N=31,873 scale (per-instrument closed-form SE check, within "
        "0.1%); no rescaling was applied. Reverse MR is on the same verified SD scale. p-values, "
        "Egger-intercept p, I² and LOO sign flips are N-scale-invariant. Strict PLINK-equivalent "
        "clumping (K = strict instrument count). For BMD femoral neck → fracture the MR-Egger "
        "intercept is significant (p=0.018); the MR-Egger slope (−0.343) and weighted median "
        "(−0.365) are the pleiotropy-robust estimates."
    )

    t = io.open(MD, encoding="utf-8").read()
    start = t.index("## Table 3. Two-sample Mendelian randomisation estimates")
    end = t.index("## Table 4. Statistical power")
    block_old = t[start:end]
    assert "| **BMD total → fracture**" in block_old, "未定位到 Table 3 内容"
    block_new = ("## Table 3. Two-sample Mendelian randomisation estimates\n\n"
                 + header + "\n" + sep + "\n" + body + "\n\n" + footnote + "\n\n")
    t = t[:start] + block_new + t[end:]
    io.open(MD, "w", encoding="utf-8").write(t)
    print("Table 3 已重建（%d 行 × %d 列）" % (len(rows), header.count("|") - 1))

    # ------------------------------------------------------------ 自检
    chk = io.open(MD, encoding="utf-8").read()
    assert "standardised outcome scale" not in chk, "旧尺度表述残留"
    for r in rows:
        assert r in chk
    # 关键数值回读
    assert "1.023 (1.016–1.031)" in chk, "BMD total→OA 的 OR 未写入"
    assert "1.40 (1.34–1.46)" in chk, "阳性对照 OR 未写入"
    assert "p=0.018" in chk
    print("自检 OK：OR 列、Egger SE/p、WM CI 均已写入；旧尺度表述已删除")

    # 控制台核对表
    print()
    print(header)
    print(sep)
    print(body)
    print()
    print("OR 口径：BMD→OA 每 SD 升高；BMD→fracture 每 SD 降低（三种 reverse 对不给 OR）。")


if __name__ == "__main__":
    main()
