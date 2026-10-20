# -*- coding: utf-8 -*-
"""
106 · 把 coloc 临界位点（PP.H4=0.7994 @ rs11228269）主动披露进稿件
两处：Results 共定位段 + Table S8 脚注。每处 text.count(old)==1 断言。
"""
import io
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DOC = os.path.join(HERE, "..", "docs", "08_投稿稿_OA_BMD_fracture.md")
DOC = os.path.normpath(DOC)

with io.open(DOC, encoding="utf-8") as f:
    t = f.read()

# ---------- (1) Results 共定位段 ----------
old1 = (
    "Across all 31 BMD instrument loci, none reached PP.H4\u22650.8 (median 0.005); "
    "six favoured distinct variants (PP.H3\u22650.8)."
)
new1 = (
    "Across all 31 BMD instrument loci, none reached PP.H4\u22650.8 (median 0.005); "
    "six favoured distinct variants (PP.H3\u22650.8). One locus fell marginally short: "
    "rs11228269 (chr11:68.29 Mb, in the LRP5 region; BMD z=\u22128.65, knee-OA z=\u22124.06, "
    "direction-concordant) gave PP.H4=0.799, i.e. 0.0006 below the pre-specified threshold, "
    "with the residual posterior mass on a BMD-only signal (PP.H1=0.156). We retain the "
    "pre-specified threshold and do not count this as colocalisation, but report it explicitly "
    "because the posterior is continuous and the 0.8 cut-off is a convention; a reader applying "
    "a 0.75 threshold would read this locus as shared."
)

# ---------- (2) Table S8 脚注 ----------
old2 = (
    "No locus reaches PP.H4\u22650.8; six favour distinct causal variants (PP.H3\u22650.8)."
)
new2 = (
    "No locus reaches PP.H4\u22650.8; six favour distinct causal variants (PP.H3\u22650.8). "
    "The highest single-locus PP.H4 was 0.799 (rs11228269, chr11:68.29 Mb; BMD z=\u22128.65, "
    "knee-OA z=\u22124.06, direction-concordant; PP.H1=0.156), i.e. one locus fell 0.0006 short "
    "of the threshold; the pre-specified 0.8 cut-off is retained."
)

for i, (old, new) in enumerate([(old1, new1), (old2, new2)], 1):
    c = t.count(old)
    assert c == 1, f"replacement {i}: expected 1 occurrence, found {c}\n---\n{old}\n---"
    t = t.replace(old, new)
    print(f"[ok] replacement {i} applied (occurrences={c})")

with io.open(DOC, "w", encoding="utf-8") as f:
    f.write(t)

print("\nWrote:", DOC)
print("Total chars:", len(t))
