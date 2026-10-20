# -*- coding: utf-8 -*-
import re
t = open("docs/08_投稿稿_OA_BMD_fracture.md", encoding="utf-8").read()

def words(s):
    return len(re.findall(r"[A-Za-z][A-Za-z0-9'\-]*", s))

# --- abstract block ---
ai = t.find("**Objective.**")
# abstract ends at the end of the Conclusions paragraph; find "**Conclusions.**" then take to next blank-line section
ci = t.find("**Conclusions.**", ai)
# grab from ci to the next "\n\n---\n" or "\n\n## "
tail = t[ci:]
end_marker = re.search(r"\n---", tail)
abstract_end = ci + (end_marker.start() if end_marker else len(tail))
abstract_txt = t[ai:abstract_end]
print("ABSTRACT words:", words(abstract_txt))

# --- main body = after abstract_end, before '## Figure legends' ---
leg = t.find("## Figure legends")
body = t[abstract_end:leg]
# drop table rows
body_lines = [l for l in body.splitlines() if not l.strip().startswith("|")]
body_nt = "\n".join(body_lines)
print("MAIN BODY words (excl abstract, tables, legends):", words(body_nt))

# --- references count ---
ri = t.find("## References")
refblock = t[ri:]
# reference entries start with optional number then '.' at line start
refs = re.findall(r"^\n?(\d+)\. ", refblock, re.M)
# count max index
nums = [int(x) for x in refs]
print("REFERENCE count (max index):", max(nums) if nums else 0)

# --- display items in MAIN text (figures/tables actually shown, not S-prefixed) ---
# NOTE: in this manuscript the tables are placed after '## Figure legends', so the
# count must span the WHOLE file and exclude supplementary (S-prefixed) items.
figs = sorted(set(int(x) for x in re.findall(r"Figure (\d+)", t)), key=int)
tabs = sorted(set(int(x) for x in re.findall(r"## Table (\d+)\.", t)), key=int)
print("MAIN figures (Figure N):", ["Figure " + str(i) for i in figs])
print("MAIN tables (## Table N):", ["Table " + str(i) for i in tabs])
print("DISPLAY ITEMS TOTAL:", len(figs) + len(tabs), "(OAC limit 8)")
# supplementary (should not count toward 8)
sfigs = len(re.findall(r"Supplementary Figure S\d+", t))
stabs = len(re.findall(r"## Supplementary Table S\d+", t))
print("SUPPL figures:", sfigs, "SUPPL tables:", stabs)
