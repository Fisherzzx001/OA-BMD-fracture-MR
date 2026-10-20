# -*- coding: utf-8 -*-
import re
t = open("docs/08_投稿稿_OA_BMD_fracture.md", encoding="utf-8").read()
def words(s):
    return len(re.findall(r"[A-Za-z][A-Za-z0-9'\-]*", s))

ai = t.find("**Objective.**")
ci = t.find("**Conclusions.**", ai)
tail = t[ci:]
end_marker = re.search(r"\n---", tail)
abstract_end = ci + (end_marker.start() if end_marker else len(tail))
leg = t.find("## Figure legends")
body = t[abstract_end:leg]
# split into paragraphs (blank-line separated), drop table rows
paras = [p for p in body.split("\n\n") if p.strip() and not p.strip().startswith("|")]
ranked = []
for p in paras:
    # strip markdown bold markers for counting but keep text
    txt = p
    w = words(txt)
    if w >= 60:
        ranked.append((w, p[:90].replace("\n"," ")))
ranked.sort(reverse=True)
print("Paragraphs >=60 words:", len(ranked), " total words in those:", sum(w for w,_ in ranked))
tot = sum(words(p) for p in paras)
print("TOTAL body words:", tot)
print("---- top 20 longest ----")
for w, head in ranked[:20]:
    print(f"{w:4d}  {head}")
