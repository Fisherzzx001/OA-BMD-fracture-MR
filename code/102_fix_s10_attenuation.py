# -*- coding: utf-8 -*-
"""
102_fix_s10_attenuation.py
==========================
Recompute BOTH attenuation columns of Supplementary Table S10 from each row's
own beta columns, using exactly the definitions printed in the table header:

    Atten% vs Table 3      = 100 * ( beta_uni(Table 3)      - beta_MVMR ) / beta_uni(Table 3)
    Atten% vs BMD-subset   = 100 * ( beta_uni(MVMR set)     - beta_MVMR ) / beta_uni(MVMR set)

The beta columns themselves are NOT modified (they were verified against
Table S9: beta_uni 0.02312 vs S9 0.0231; beta_MVMR 0.01871 vs S9 0.0187).

Rewrites the file in place and reports the resulting ranges so that the
manuscript headline ranges ("6-33%", "13.5-33.2%") can be re-checked.
"""
from __future__ import annotations

import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
TARGET = os.path.join(HERE, "..", "docs", "Supplementary", "02_mvmr_height_grid.md")
TARGET = os.path.normpath(TARGET)


def cell(row: str):
    return [c.strip() for c in row.strip().strip("|").split("|")]


def to_f(s: str) -> float:
    return float(s.replace("−", "-").replace("\u2212", "-"))


def main() -> int:
    with open(TARGET, "r", encoding="utf-8") as fh:
        lines = fh.read().split("\n")

    out = []
    n_fixed = 0
    recs = {}          # (height_src, n_h) -> list of dict
    for ln in lines:
        if not ln.strip().startswith("|"):
            out.append(ln)
            continue
        parts = ln.strip().strip("|").split("|")
        if len(parts) < 12:
            out.append(ln)
            continue
        head = parts[0].strip()
        if head.startswith("---") or head.lower().startswith("exposure") or set(head) <= set("-: "):
            out.append(ln)
            continue

        exposure = head
        hsrc = parts[1].strip()
        try:
            n_h = int(parts[2].strip())
        except ValueError:
            out.append(ln)
            continue
        b_t3 = to_f(parts[4])
        b_ms = to_f(parts[5])
        b_mv = to_f(parts[6])
        old1, old2 = to_f(parts[7]), to_f(parts[8])

        new1 = 100.0 * (b_t3 - b_mv) / b_t3
        new2 = 100.0 * (b_ms - b_mv) / b_ms

        if abs(new1 - old1) > 0.05 or abs(new2 - old2) > 0.05:
            n_fixed += 1
        parts[7] = f"{new1:.1f}"
        parts[8] = f"{new2:.1f}"
        out.append("| " + " | ".join(parts) + " |")

        recs.setdefault((hsrc, n_h), []).append(
            dict(exposure=exposure, old1=old1, new1=new1, old2=old2, new2=new2)
        )

    with open(TARGET, "w", encoding="utf-8") as fh:
        fh.write("\n".join(out))

    print(f"rows rewritten: {n_fixed} (of {sum(len(v) for v in recs.values())})")
    print()
    print("=== corrected 'Atten% vs Table 3' for the OA-knee rows, primary n_h=100 ===")
    for hsrc in ("EUR", "EUR_noUKB"):
        for r in recs.get((hsrc, 100), []):
            if "fracture" in r["exposure"].lower():
                continue
            print(f"  {hsrc:<10} {r['exposure']:<24} old={r['old1']:>7.1f}%   new={r['new1']:>7.1f}%")
    print()
    print("=== positive-control rows (fracture), primary n_h=100 ===")
    for hsrc in ("EUR", "EUR_noUKB"):
        for r in recs.get((hsrc, 100), []):
            if "fracture" not in r["exposure"].lower():
                continue
            print(f"  {hsrc:<10} {r['exposure']:<24} old={r['old1']:>7.1f}%   new={r['new1']:>7.1f}%")
    print()

    oa1 = [r["new1"] for k, v in recs.items() for r in v
           if "fracture" not in r["exposure"].lower()]
    oa2 = [r["new2"] for k, v in recs.items() for r in v
           if "fracture" not in r["exposure"].lower()]
    print(f"OA-knee rows, corrected 'vs Table 3'      : min={min(oa1):.1f}%  max={max(oa1):.1f}%")
    print(f"OA-knee rows, corrected 'vs BMD-subset'   : min={min(oa2):.1f}%  max={max(oa2):.1f}%")
    pos = [x for x in oa1 + oa2 if x > 0]
    print(f"positive subset across both baselines     : min={min(pos):.1f}%  max={max(pos):.1f}%")
    print(f"negative (opposite-direction) values present: "
          f"{sorted(set(round(x,1) for x in oa1+oa2 if x < 0))}")

    oa1_100 = [r["new1"] for (h, n), v in recs.items() if n == 100 for r in v
               if "fracture" not in r["exposure"].lower() and h == "EUR"]
    print(f"\nEUR n_h=100 'vs Table 3' range            : "
          f"min={min(oa1_100):.1f}%  max={max(oa1_100):.1f}%   <- compare to S9 claim 13.5-33.2%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
