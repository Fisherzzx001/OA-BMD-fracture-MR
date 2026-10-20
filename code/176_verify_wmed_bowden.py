# -*- coding: utf-8 -*-
"""独立第三实现：Bowden 加权中位数，验证 arms/legs/lumbar/total->fracture/OA->* 的 wmed。
不依赖 m13 / code/129，纯手工实现，交叉核对 mr_results_corrected(wm_recomp) 与
table3_corrected(wm_tab) 谁对。"""
import csv
import os

HARM = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                    "code", "output", "genetics", "harmonised")


def bowden_wm(be, so, bo):
    """Bowden 2016: ratio=bo/be, w=(be/so)^2, 加权 50% 分位（线性插值版）。"""
    ratios = sorted(zip([b / e for b, e in zip(bo, be)],
                        [(e / s) ** 2 for e, s in zip(be, so)]))
    total = sum(w for _, w in ratios)
    cum = 0.0
    for i, (r, w) in enumerate(ratios):
        cum += w
        if cum >= 0.5 * total:
            if i == 0:
                return r
            prev_cum = cum - w
            # 在 r[i-1] 与 r[i] 之间线性插值
            r0 = ratios[i - 1][0]
            if cum == prev_cum:
                return r
            frac = (0.5 * total - prev_cum) / (cum - prev_cum)
            return r0 + frac * (r - r0)
    return ratios[-1][0]


def bowden_wm_step(be, so, bo):
    """Bowden 2016 次序统计量版（code/129 的定义：searchsorted 取该处比值）。"""
    pairs = sorted(zip([b / e for b, e in zip(bo, be)],
                       [(e / s) ** 2 for e, s in zip(be, so)]))
    total = sum(w for _, w in pairs)
    cw, s = [], 0.0
    for _, w in pairs:
        s += w
        cw.append(s / total)
    import bisect
    idx = bisect.bisect_left(cw, 0.5)
    idx = min(idx, len(pairs) - 1)
    return pairs[idx][0]


PAIRS = [("BMD_total", "OA_knee"), ("BMD_femur_neck", "OA_knee"),
         ("BMD_lumbar", "OA_knee"), ("BMD_arms", "OA_knee"),
         ("BMD_leg", "OA_knee"), ("BMD_total", "fracture"),
         ("BMD_femur_neck", "fracture"), ("OA_knee", "fracture"),
         ("OA_knee", "BMD_total"), ("OA_knee", "BMD_femur_neck")]

# 参照值
REF = {  # (mr_results_corrected=m13/Bowden, table3_corrected=旧表)
    ("BMD_total", "OA_knee"): (0.0653, 0.0653),
    ("BMD_femur_neck", "OA_knee"): (0.0520, 0.0520),
    ("BMD_lumbar", "OA_knee"): (0.1146, 0.1255),
    ("BMD_arms", "OA_knee"): (0.0257, 0.1039),
    ("BMD_leg", "OA_knee"): (0.0261, 0.0840),
    ("BMD_total", "fracture"): (-0.3926, -0.3642),
    ("BMD_femur_neck", "fracture"): (-0.3646, -0.3646),
    ("OA_knee", "fracture"): (0.0185, 0.0349),
    ("OA_knee", "BMD_total"): (0.0256, 0.0330),
    ("OA_knee", "BMD_femur_neck"): (0.0231, 0.0334),
}

print(f"{'pair':<34} {'插值版':>10} {'次序版':>10} {'m13/ Bowden':>12} {'table3旧值':>11}  判定")
for exp, out in PAIRS:
    fp = os.path.join(HARM, f"{exp}__{out}.csv")
    with open(fp, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    be = [float(r["beta_e"]) for r in rows]
    se = [float(r["se_e"]) for r in rows]
    bo = [float(r["beta_o"]) for r in rows]
    so = [float(r["se_o"]) for r in rows]
    w_interp = bowden_wm(be, so, bo)
    w_step = bowden_wm_step(be, so, bo)
    ref_bow, ref_tab = REF[(exp, out)]
    verdict = "复算值正确" if abs(w_step - ref_bow) < 5e-4 else (
        "旧表值正确" if abs(w_step - ref_tab) < 5e-4 else "？？")
    print(f"{exp+' -> '+out:<34} {w_interp:>10.4f} {w_step:>10.4f} {ref_bow:>12.4f} {ref_tab:>11.4f}  {verdict}")
