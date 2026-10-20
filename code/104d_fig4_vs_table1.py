# -*- coding: utf-8 -*-
"""104d — Figure 4 h2 (all 13 traits) vs manuscript Table 1."""
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
f = json.load(open(os.path.join(ROOT, "Figures_v17_Nature", "plotted_data.json"), encoding="utf-8"))["Figure4"]

T1 = {
    "OA (knee)": 0.0354, "OA (self-reported)": 0.0223, "Any fracture": 0.0258,
    "BMD total body": 0.391, "BMD femoral neck": 0.257, "BMD femoral total": 0.317,
    "BMD lumbar spine": 0.350, "BMD spine": 0.357, "BMD arms": 0.325,
    "BMD legs": 0.325, "BMD pelvis": 0.358, "BMD ribs": 0.303, "BMD trunk": 0.366,
}

print("trait".ljust(22), "Fig4 h2".rjust(10), "Table1".rjust(10), "match".rjust(9))
bad = 0
for t, h in zip(f["traits"], f["h2_obs"]):
    exp = T1.get(t)
    if exp is None:
        print(t.ljust(22), f"{h:10.4f}", "??".rjust(10), "UNKNOWN".rjust(9))
        bad += 1
        continue
    ok = abs(h - exp) <= 0.0006
    bad += (not ok)
    print(t.ljust(22), f"{h:10.4f}", f"{exp:10.3f}", ("OK" if ok else "MISMATCH").rjust(9))
print()
print("Fig4 vs Table 1 mismatches:", bad)
