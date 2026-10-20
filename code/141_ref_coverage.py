# -*- coding: utf-8 -*-
"""参考文献覆盖核查（修正版）：抓复合引用 [18,32] / [27–30]，并按作者姓氏二次确认"""
import io, os, re

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
md = io.open(os.path.join(PROJ, "Manuscript_Ready", "01_Manuscript.md"), encoding="utf-8").read()
ref_head = md.index("## References")
body = md[:ref_head]
refs = md[ref_head:]

# 1) 所有方括号内的数字（含复合、含 27–30 这类范围，须展开中间数字）
cited = set()
for m in re.finditer(r"\[([\d\s,\u2013\u2014\-]+)\]", body):
    inner = m.group(1)
    # 先展开范围 27–30
    def _expand(mo):
        a, b = int(mo.group(1)), int(mo.group(2))
        return ",".join(str(i) for i in range(a, b + 1))
    inner = re.sub(r"(\d+)\s*[\u2013\u2014-]\s*(\d+)", _expand, inner)
    for n in re.findall(r"\d+", inner):
        cited.add(int(n))
# 2) "Refs 1–3,10"
for m in re.finditer(r"Refs?\s+([\d,\u2013\u2014\-and ]+)", body):
    inner = re.sub(r"(\d+)\s*[\u2013\u2014-]\s*(\d+)", _expand, m.group(1))
    for n in re.findall(r"\d+", inner):
        cited.add(int(n))

listed = {}
for m in re.finditer(r"^(\d+)\.\s+(.+?)(?=^\d+\.\s|\Z)", refs, re.M | re.S):
    listed[int(m.group(1))] = " ".join(m.group(2).split())

truly_uncited = []
for n, txt in sorted(listed.items()):
    if n in cited:
        continue
    # 作者姓氏是否出现在正文（文本式引用）
    first = re.match(r"([A-Z][A-Za-z\u00C0-\u024F'\-]+)\s", txt)
    sur = first.group(1) if first else None
    # 词边界匹配，避免 "Zhen" 命中作者名 "Zhenxin"
    txt_hit = bool(sur) and re.search(rf"\b{re.escape(sur)}\b", body) is not None
    truly_uncited.append((n, sur, txt_hit, txt[:95]))
    flag = "文本式引用(作者名)" if txt_hit else ">>> 完全未引用 <<<"
    print(f"[{n:>2}] {sur or '?':<12} {flag}")
    print(f"        {txt[:110]}")

print()
print(f"参考文献总数 {len(listed)}；正文编号引用 {len(cited & set(listed))}；"
      f"未以编号引用 {len(truly_uncited)}")
hard = [n for n, s, t, _ in truly_uncited if not t]
print(f"其中正文既无编号也无作者名的（硬挂名）: {hard}")
