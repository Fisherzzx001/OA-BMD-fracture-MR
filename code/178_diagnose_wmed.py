# -*- coding: utf-8 -*-
"""权威裁决：从【当前】harmonised 文件直接计算 BMD→OA_knee / OA_knee→* 的加权中位数
（Bowden 2016），分别用 步近约定(cum>=0.5) 与 中点约定(cum-0.5w)。
对照：稿件 Table 2、mr_results_corrected.csv、table3_corrected.csv、wm_recomputed_vs_tabulated.csv。
目的：判定哪种约定 + 哪个文件是 FINAL 真源。"""
import csv
import bisect
import os

H = "code/output/genetics/harmonised"

PAIRS = [
    ("BMD_total", "fracture"),
    ("BMD_femur_neck", "fracture"),
    ("BMD_total", "OA_knee"),
    ("BMD_lumbar", "OA_knee"),
    ("BMD_arms", "OA_knee"),
    ("BMD_leg", "OA_knee"),
    ("BMD_femur_neck", "OA_knee"),
    ("OA_knee", "fracture"),
    ("OA_knee", "BMD_total"),
    ("OA_knee", "BMD_femur_neck"),
]


def load_pairs(e, o):
    p = os.path.join(H, f"{e}__{o}.csv")
    if not os.path.exists(p):
        return None
    with open(p, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    be = [float(r["beta_e"]) for r in rows]
    se = [float(r["se_e"]) for r in rows]
    bo = [float(r["beta_o"]) for r in rows]
    so = [float(r["se_o"]) for r in rows]
    return be, se, bo, so


def wm_step(be, se, bo, so):
    pairs = sorted(zip([b / e for b, e in zip(bo, be)],
                       [(e / s) ** 2 for e, s in zip(be, so)]))
    w = [x[1] for x in pairs]
    cum = [sum(w[:i + 1]) for i in range(len(w))]
    total = cum[-1]
    idx = bisect.bisect_left([c / total for c in cum], 0.5)
    idx = min(idx, len(pairs) - 1)
    return pairs[idx][0]


def wm_mid(be, se, bo, so):
    pairs = sorted(zip([b / e for b, e in zip(bo, be)],
                       [(e / s) ** 2 for e, s in zip(be, so)]))
    w = [x[1] for x in pairs]
    cum = [sum(w[:i + 1]) for i in range(len(w))]
    total = cum[-1]
    for i in range(len(pairs)):
        if (cum[i] - 0.5 * w[i]) / total <= 0.5 <= cum[i] / total:
            if i + 1 < len(pairs):
                return 0.5 * (pairs[i][0] + pairs[i + 1][0])
            return pairs[i][0]
    return pairs[-1][0]


# 加载对照源
def load_csv(p):
    with open(p, newline="", encoding="utf-8") as f:
        return {(r.get("exposure") or r.get("exp"), r.get("outcome") or r.get("out")): r
                for r in csv.DictReader(f)}

mr = load_csv("code/output/genetics/mr_results_corrected.csv")
t3 = load_csv("code/output/genetics/table3_corrected.csv")
wm = {}
with open("code/output/genetics/wm_recomputed_vs_tabulated.csv", newline="", encoding="utf-8") as f:
    for r in csv.DictReader(f):
        wm[r["pair"]] = r

# 稿件 Table 2 wmed（手工转录自当前 docs/08）
MAN = {
    ("BMD_total", "fracture"): -0.364,
    ("BMD_femur_neck", "fracture"): -0.365,
    ("BMD_total", "OA_knee"): 0.065,
    ("BMD_lumbar", "OA_knee"): 0.126,
    ("BMD_arms", "OA_knee"): 0.104,
    ("BMD_leg", "OA_knee"): 0.084,
    ("BMD_femur_neck", "OA_knee"): 0.052,
    ("OA_knee", "fracture"): 0.035,
    ("OA_knee", "BMD_total"): 0.033,
    ("OA_knee", "BMD_femur_neck"): 0.033,
}

print(f"{'pair':32s} {'|b_o|max':>8s} {'step':>9s} {'mid':>9s} | {'man':>7s} {'mr':>8s} {'t3':>8s} {'wm_rec':>8s} {'wm_tab':>8s}")
for (e, o) in PAIRS:
    d = load_pairs(e, o)
    if d is None:
        print(f"{e}__{o:24s}  (harmonised 缺失)")
        continue
    be, se, bo, so = d
    s = wm_step(be, se, bo, so)
    m = wm_mid(be, se, bo, so)
    maxbo = max(abs(x) for x in bo)
    man = MAN.get((e, o), float("nan"))
    mrw = float(mr[(e, o)]["wmed"]) if (e, o) in mr else float("nan")
    t3w = float(t3[(e, o)]["wmed"]) if (e, o) in t3 else float("nan")
    wr = float(wm[f"{e}__{o}"]["wm_recomp"])
    wt = float(wm[f"{e}__{o}"]["wm_tab"])
    print(f"{e+'__'+o:32s} {maxbo:8.3f} {s:9.4f} {m:9.4f} | {man:7.3f} {mrw:8.4f} {t3w:8.4f} {wr:8.4f} {wt:8.4f}")
