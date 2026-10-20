# -*- coding: utf-8 -*-
"""
113 · 勘误：Methods 工具数 22–37/118 → 20–35/116
原因：上一轮审计用 `wc -l - 1` 统计 PLINK .clumped 文件，但文件尾部有 2 个空行，
导致把数据行数 20–35/116 误报为 22–37/118 并"改正"了原本正确的稿件。
一手证据：`awk 'NF && NR>1' *_plink.clumped | wc -l` = 20/28/24/26/35 与 116，
且与 strict_clump_summary.csv 的 n_kept_strict 逐项一致。
"""
import io
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DOC = os.path.normpath(os.path.join(HERE, "..", "docs", "08_投稿稿_OA_BMD_fracture.md"))

with io.open(DOC, encoding="utf-8") as f:
    t = f.read()

old = ("22\u201337 instruments retained per analysed BMD site and 118 for knee OA). "
       "Harmonisation against each outcome \u2014 exclusion of palindromic and "
       "strand-ambiguous variants and of instruments absent from the outcome GWAS \u2014 "
       "reduced the analysed sets to 16\u201331 instruments per forward pair and 69 per "
       "knee-OA-exposure pair, an attrition we report explicitly.")
new = ("20\u201335 instruments retained per analysed BMD site and 116 for knee OA). "
       "Harmonisation against each outcome \u2014 exclusion of palindromic and "
       "strand-ambiguous variants and of instruments absent from the outcome GWAS \u2014 "
       "reduced the analysed sets to 16\u201331 instruments per forward pair and 69 per "
       "knee-OA-exposure pair; 47 of the 116 knee-OA instruments (41%) were lost at this "
       "step, an attrition we report explicitly.")

c = t.count(old)
assert c == 1, f"expected 1 occurrence, found {c}"
t = t.replace(old, new)

with io.open(DOC, "w", encoding="utf-8") as f:
    f.write(t)
print("[ok] Methods instrument counts corrected: 20-35 / 116; knee-OA attrition 47/116 (41%)")
