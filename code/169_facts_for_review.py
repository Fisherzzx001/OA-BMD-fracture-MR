# -*- coding: utf-8 -*-
"""本轮预审的事实取证: 只报可复现的具体缺陷, 不做印象式评价."""
import re

P = "docs/08_投稿稿_OA_BMD_fracture.md"
t = open(P, encoding="utf-8").read()
main, rest = t.split("## References", 1)
refs_block = rest.split("\n## ", 1)[0]

print("=" * 78)
print("A. 参考文献编号顺序 (Elsevier/OAC = Vancouver, 必须按首次出现顺序编号)")
print("=" * 78)
seq = []
for m in re.finditer(r"\[([0-9,\u2013\-\s]+)\]", main):
    s = m.group(1)
    # expand ranges written with hyphen or en-dash
    for part in re.findall(r"(\d+)\s*[\u2013-]\s*(\d+)|(\d+)", s):
        if part[2]:
            seq.append(int(part[2]))
        else:
            lo, hi = int(part[0]), int(part[1])
            if hi >= lo and hi - lo < 30:
                seq.extend(range(lo, hi + 1))
first = []
seen = set()
for c in seq:
    if c not in seen:
        seen.add(c)
        first.append(c)
print("  正文首次出现顺序:", first)
print("  是否为升序 1..N:", first == sorted(first))
if first != sorted(first):
    print("  ✗ 违反 Vancouver 按出现顺序编号 → 编辑初筛可退回")

reflist = []
for ln in refs_block.splitlines():
    m = re.match(r"^(\d+)\.\s+", ln)
    if m:
        reflist.append(int(m.group(1)))
print("  参考文献表条目数:", len(reflist))
uncited = [r for r in reflist if r not in seen]
print("  全文(含区间)从未编号引用的条目:", uncited)
print("  被引但表中无此号:", sorted(n for n in seen if n not in reflist))

print()
print("=" * 78)
print("B. 可疑/畸形数字与排版瑕疵")
print("=" * 78)
checks = [
    ("N≈1.873 (疑似压缩时把 31,873 压坏)", r"N\u22481\.873"),
    (": 后缺空格", r":[A-Za-z]{2,}"),
    ("句点后缺空格", r"[a-z,;\)]\.[A-Z][a-z]"),
    ("重复标点", r",,|;;|\.\.(?!\.)"),
    ("括号不匹配的行", None),
]
for label, pat in checks[:-1]:
    for m in re.finditer(pat, main):
        s = max(0, m.start() - 70)
        ctx = main[s:m.end() + 40].replace("\n", " ")
        print("  [%s] ...%s..." % (label, ctx[-95:]))

print()
print("=" * 78)
print("C. 关键数值一致性")
print("=" * 78)
# beta range claim in Discussion
m = re.search(r"against (0\.\d{3}\u20130\.\d{3}) here", main)
print("  L91 声称的 beta 区间:", m.group(1) if m else "未找到")
betas = ["total-body β=+0.078", "lumbar-spine +0.095", "arms +0.075",
         "legs +0.060", "femoral-neck +0.058"]
vals = [float(re.search(r"=([+-]?\d\.\d+)", b).group(1).replace('+','')) if '=' in b
        else float(re.findall(r"([+-]?\d\.\d+)", b)[-1]) for b in betas]
print("  Results/Table 2 实际 beta:", dict(zip([b.split()[0] for b in betas], vals)))
print("  实际区间: %.3f–%.3f" % (min(vals), max(vals)))

# main absolute numeric claims
for pat in [r"OR 1\.060\u20131\.099", r"OR 1\.017\u20131\.028", r"3\.6\u20135\.8-fold",
            r"9\.5\u201333\.2", r"67\u201391%", r"\-3\.7% to 46\.0%", r"11\.9"]:
    hits = len(re.findall(pat, main))
    print("  '%s' 出现 %d 次" % (pat, hits))

print()
print("=" * 78)
print("D. 声称的工作量与 MR 覆盖面的匹配")
print("=" * 78)
for kw in ["We addressed all three", "site specificity is unknown", "only if strict clumping completes",
           "reflects compute availability"]:
    for m in re.finditer(re.escape(kw), main):
        print("  [%s] ...%s..." % (kw, main[max(0, m.start()-90):m.end()+90].replace("\n"," ")[-150:]))
