# -*- coding: utf-8 -*-
# Definitive numeric-integrity check: compare every digit-bearing token between the
# last committed docs/08 and the working tree. Any difference = a changed number.
import re, subprocess, sys
from collections import Counter

def load(path):
    if path.startswith("HEAD:"):
        return subprocess.check_output(["git", "show", path]).decode("utf-8")
    return open(path, encoding="utf-8").read()

tok_re = re.compile(r"[\w%×−–\-+./()]*\d[\w%×−–\-+./()]*")

def toks(t):
    return Counter(tok_re.findall(t))

head = load("HEAD:docs/08_投稿稿_OA_BMD_fracture.md")
work = load("docs/08_投稿稿_OA_BMD_fracture.md")
ch, cw = toks(head), toks(work)

only_head = ch - cw   # present in committed but not working
only_work = cw - ch   # present in working but not committed
changed = {k: (ch[k], cw[k]) for k in (only_head | only_work) if k in ch and k in cw and ch[k] != cw[k]}

print("committed digit-tokens:", sum(ch.values()), "| working digit-tokens:", sum(cw.values()))
print("tokens only in committed:", dict(only_head))
print("tokens only in working  :", dict(only_work))
print("tokens with changed count:", changed)
if only_head or only_work or changed:
    print("RESULT: NUMERIC DIFFERENCE DETECTED -> DO NOT COMMIT")
    sys.exit(1)
else:
    print("RESULT: numeric multiset identical (only non-numeric prose differs) -> SAFE")
    sys.exit(0)
