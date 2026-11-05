# Sentence/paragraph-level word audit using the SAME counting rule as code/148,
# so I can price each candidate edit before applying it.
# rule: words(s) = len(re.findall(r"[A-Za-z][A-Za-z0-9'\-]*", s))
# NOTE: numeric citation markers such as [35-38] cost ZERO words under this rule.
import re, sys

F = "docs/08_投稿稿_OA_BMD_fracture.md"
t = open(F, encoding="utf-8").read()

def words(s):
    return len(re.findall(r"[A-Za-z][A-Za-z0-9'\-]*", s))

leg = t.find("## Figure legends")
body_lines = [l for l in t.splitlines() if not l.strip().startswith("|")]

# split body into paragraphs (non-empty, non-table)
start_line = None
end_line = None
lines = t.splitlines()
for i, l in enumerate(lines):
    if l.startswith("## Figure legends"):
        end_line = i
        break
body = "\n".join(lines[:end_line])

paras = [p for p in body.split("\n\n") if p.strip() and not p.strip().startswith("|")]

print("=" * 78)
print("PARAGRAPH-LEVEL WORD COST (body, tables excluded)")
print("=" * 78)
idx = 0
for p in paras:
    # skip headings and everything before Discussion start marker
    first = p.strip().split("\n")[0][:70]
    w = words(p)
    if w >= 40:  # only show substantial paragraphs
        idx += 1
        print(f"[P{idx:02d}] {w:4d}w | {first}...")

print()
print("=" * 78)
print("TARGET: Discussion eta paragraphs (L85-L101) - sentence breakdown")
print("=" * 78)
targets = ["The architecture is consistent across 13",
           "The opposing-endpoint pattern is an old",
           "Zhou et al. [10] reported",
           "The directional agreement with Hartley",
           "OA is a disorder of the osteochondral unit",
           "**Pre-registered SOS predictions",
           "A direct MR test of OA",
           "**Strengths and limitations.**",
           "**Conclusion.**"]
for p in paras:
    ps = p.strip()
    for tg in targets:
        if ps.startswith(tg):
            print(f"\n### TOTAL {words(p):4d}w :: {tg}")
            for s in re.split(r"(?<=[.!?])\s+", ps):
                s = s.strip()
                if s:
                    print(f"    {words(s):4d}w | {s[:150]}")
            break
