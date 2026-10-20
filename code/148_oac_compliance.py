# -*- coding: utf-8 -*-
# Full OAC compliance check for docs/08.
# Body = prose between end of Abstract(Conclusions) and '## Figure legends', table rows excluded.
# Mirrors the repository's live checker code/144 for word counts, and adds the
# OAC Highlights rule (3-5 bullets, each <=85 characters incl. spaces).
import re, sys

F = "docs/08_投稿稿_OA_BMD_fracture.md"
t = open(F, encoding="utf-8").read()

def words(s):
    return len(re.findall(r"[A-Za-z][A-Za-z0-9'\-]*", s))

# --- abstract ---
ai = t.find("**Objective.**")
ci = t.find("**Conclusions.**", ai)
tail = t[ci:]
em = re.search(r"\n---", tail)
ae = ci + (em.start() if em else len(tail))
abstract_txt = t[ai:ae]
aw = words(abstract_txt)

# --- body ---
leg = t.find("## Figure legends")
body = t[ae:leg]
body_lines = [l for l in body.splitlines() if not l.strip().startswith("|")]
body_nt = "\n".join(body_lines)
bw = words(body_nt)

# --- references ---
ri = t.find("## References")
refblock = t[ri:]
nums = [int(x) for x in re.findall(r"^\n?(\d+)\. ", refblock, re.M)]
rcount = max(nums) if nums else 0

# --- display items ---
figs = sorted(set(int(x) for x in re.findall(r"Figure (\d+)", t)), key=int)
tabs = sorted(set(int(x) for x in re.findall(r"## Table (\d+)\.", t)), key=int)
display = len(figs) + len(tabs)

# --- highlights ---
hi = t.find("## Highlights")
hl = t.find("**Keywords.**")
hl_block = t[hi:hl]
bullets = [m.group(1) for m in re.finditer(r"^- (.+)$", hl_block, re.M)]
hl_lens = [len(b) for b in bullets]
hl_ok = all(L <= 85 for L in hl_lens) and (3 <= len(bullets) <= 5)

# --- keywords ---
kw = t.find("**Keywords.**")
ke = t.find("\n", kw)
kwt = t[kw:ke]
kw_body = kwt.split("**Keywords.**", 1)[1] if "**Keywords.**" in kwt else ""
n_kw = len([x for x in kw_body.split(";") if x.strip()])

print("=== OAC COMPLIANCE (docs/08) ===")
print(f"Abstract words      : {aw:4d}  (limit 350)  {'OK' if aw<=350 else 'FAIL'}")
print(f"Body words          : {bw:4d}  (limit 4000) {'OK' if bw<=4000 else 'FAIL'}")
print(f"References          : {rcount:4d}  (limit 50)   {'OK' if rcount<=50 else 'FAIL'}")
print(f"Display items       : {display:4d}  (limit 8)    {'OK' if display<=8 else 'FAIL'}")
print(f"Highlights          : {len(bullets)} bullets; lengths {hl_lens} (limit 85 each) {'OK' if hl_ok else 'FAIL'}")
print(f"Keywords            : {n_kw} (limit 4-6) {'OK' if 4<=n_kw<=6 else 'FAIL'}")
overall = (aw<=350 and bw<=4000 and rcount<=50 and display<=8 and hl_ok and 4<=n_kw<=6)
print("OVERALL OAC COMPLIANT:", "YES" if overall else "NO")
sys.exit(0 if overall else 1)
