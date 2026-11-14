# -*- coding: utf-8 -*-
# Code/352d — small margin trim so the post independent-test scale-validation
# update (which adds words) keeps body <=4000. Removes only redundant prose.
import re, os, shutil, sys
F = "docs/08_投稿稿_OA_BMD_fracture.md"
BAK = F + ".bak_before_352d"
if not os.path.exists(BAK):
    shutil.copyfile(F, BAK)
t = open(F, encoding="utf-8").read()
def words(s): return len([w for w in s.split() if re.search(r"[^\W_]", w, re.UNICODE)])
def body_words(txt):
    ai=txt.find("**Objective.**"); ci=txt.find("**Conclusions.**",ai)
    tail=txt[ci:]; em=re.search(r"\n---",tail); ae=ci+(em.start() if em else len(tail))
    leg=txt.find("## Figure legends"); body=txt[ae:leg]
    return words("\n".join(l for l in body.splitlines() if not l.strip().startswith("|")))
before=body_words(t)
EDITS=[
 ("Height and BMI were conditioned separately, each adjusted one at a time.",
  "", 1, "M1"),
 ("six loci favoured distinct variants (PP.H3≥0.8) (Fig. S18, S19).",
  "six loci favoured distinct variants (PP.H3≥0.8).", 1, "M2"),
 ("F1 demands a positive, knee-like hip effect.",
  "F1 demands a knee-like hip effect.", 1, "M3"),
]
fail=[]; rem=0
for old,new,exp,desc in EDITS:
    c=t.count(old)
    if c!=exp: fail.append((desc,exp,c,old[:60])); continue
    rem+=words(old)-words(new); t=t.replace(old,new)
if fail:
    print("ABORT:",fail); sys.exit(1)
open(F,"w",encoding="utf-8").write(t)
after=body_words(t)
print(f"body {before} -> {after} (removed {rem}); <=4000: {after<=4000}")
