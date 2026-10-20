# -*- coding: utf-8 -*-
"""
105 · coloc 阈值临界位点独立复核
目的：核验稿件"无位点达到 PP.H4>=0.8"是否字面成立，并定位最高 PP.H4 位点，
      为主动披露提供精确数字（证据来自一手 coloc_screen_all_loci.csv）。
"""
import csv
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "output", "genetics", "coloc_screen_all_loci.csv")
OUT = os.path.join(HERE, "output", "genetics", "AUDIT_105_coloc_threshold.json")

rows = list(csv.DictReader(open(SRC)))
assert len(rows) == 31, f"expected 31 loci, got {len(rows)}"

pp4 = [(r["index_snp"], float(r["PP4"]), float(r["PP3"])) for r in rows]
pp4_sorted = sorted(pp4, key=lambda x: -x[1])

max4 = pp4_sorted[0]
n_ge_08 = sum(1 for _, v, _ in pp4 if v >= 0.8)
n_pp3_ge_08 = sum(1 for _, _, v in pp4 if v >= 0.8)

vals = sorted(v for _, v, _ in pp4)
mid = vals[len(vals) // 2] if len(vals) % 2 == 1 else (vals[len(vals) // 2 - 1] + vals[len(vals) // 2]) / 2

res = {
    "n_loci": len(rows),
    "n_PP4_ge_0.8": n_ge_08,
    "n_PP3_ge_0.8": n_pp3_ge_08,
    "max_PP4": round(max4[1], 4),
    "max_PP4_locus": max4[0],
    "max_PP4_locus_PP3": round(max4[2], 4),
    "gap_to_threshold": round(0.8 - max4[1], 4),
    "median_PP4": round(mid, 4),
    "top5_PP4": [(s, round(v, 4)) for s, v, _ in pp4_sorted[:5]],
}

print(json.dumps(res, indent=2))
with open(OUT, "w", encoding="utf-8") as f:
    json.dump(res, f, indent=2)

# --- 断言：稿件现有声称是否字面成立 ---
assert n_ge_08 == 0, "FAIL: 存在 PP.H4>=0.8 的位点，稿件声称不成立"
assert n_pp3_ge_08 == 6, f"FAIL: PP.H3>=0.8 计数={n_pp3_ge_08}，稿件写 6"
assert abs(mid - 0.0052) < 0.001, f"FAIL: 中位数 PP.H4={mid}，稿件写 0.005"
print("\nOK: 稿件共定位三项声称全部字面成立；最高 PP.H4 =", round(max4[1], 4), "at", max4[0])
