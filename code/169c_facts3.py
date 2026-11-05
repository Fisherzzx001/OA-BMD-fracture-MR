# -*- coding: utf-8 -*-
"""取证3: (1) 正文中的作者年式引用精确定位; (2) 跨段重复句检测(压缩产物?)."""
import re
from collections import OrderedDict

P = "docs/08_投稿稿_OA_BMD_fracture.md"
t = open(P, encoding="utf-8").read()
main, rest = t.split("## References", 1)
refs_block = rest.split("\n## ", 1)[0]

print("=" * 78)
print("E. 作者年式引用 (需转为编号或补入参考表)")
print("=" * 78)
# 形如: Name et al., Journal 2021  或  Name et al. 2014
pat = re.compile(r"((?:[A-Z][A-Za-z\u00C0-\u017F\-]+\s+)?(?:[A-Z][A-Za-z\u00C0-\u017F\-]+\s+)?et\s+al\.?[^.]{0,40}?((?:19|20)\d{2})[a-z]?)")
seen = OrderedDict()
for m in pat.finditer(main):
    key = m.group(1).strip()
    if key in seen:
        continue
    seen[key] = m.start()
for k, pos in seen.items():
    ctx = main[max(0, pos - 60):pos + 45].replace("\n", " ")
    print("  %-46s | ...%s..." % (k[:46], ctx[-100:]))

print()
print("  参考表未登记的关键方法来源: Skrivankova/STROBE-MR=%s, Giambartolomei/coloc=%s"
      % ("Skrivankova" in refs_block, "Giambartolomei" in refs_block))

print()
print("=" * 78)
print("I. 重复句检测 (>=10 词, 归一化后完全相同)")
print("=" * 78)
lines = main.split("\n")
sent_sources = []
for i, ln in enumerate(lines, 1):
    if ln.startswith(("#", "|", "-", "**Figure", "**Supplementary", "**Table", "**Graphical")):
        continue
    for s in re.split(r"(?<=[.!?])\s+", ln):
        words = re.findall(r"[A-Za-z][A-Za-z\-']*", s)
        if len(words) >= 10:
            norm = " ".join(w.lower() for w in words)
            sent_sources.append((norm, i, s.strip()[:80]))

dup = {}
for norm, ln_no, orig in sent_sources:
    dup.setdefault(norm, []).append((ln_no, orig))

ndup = 0
for norm, occ in dup.items():
    if len(occ) > 1:
        lines_nums = [o[0] for o in occ]
        # 同一行内的重复不计(表格/脚注除外); 跨行重复才是问题
        if len(set(lines_nums)) > 1 or len(occ) > 1 and max(lines_nums) - min(lines_nums) > 0:
            ndup += 1
            print("  ×%d  行 %s" % (len(occ), sorted(set(lines_nums))))
            print("       %s" % occ[0][1])
print("  重复句组数:", ndup)

print()
print("=" * 78)
print("J. 个别句子是否在相邻位置被粘贴两次")
print("=" * 78)
# 直接检查关键句出现次数
for probe in ["the near-1.0 intercept does not exclude",
              "The 13\u00d713 r_g matrix showed three blocks",
              "Because LDSC uses a European-only reference"]:
    hits = [m.start() for m in re.finditer(re.escape(probe), main)]
    print("  '%s' -> %d 次, 位置 %s" % (probe[:44], len(hits), hits))
    for h in hits:
        print("       ...%s..." % main[h - 40:h + 110].replace("\n", " "))
