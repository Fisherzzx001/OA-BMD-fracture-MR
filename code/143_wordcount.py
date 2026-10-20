# -*- coding: utf-8 -*-
import re
F = "docs/08_投稿稿_OA_BMD_fracture.md"
t = open(F, encoding="utf-8").read()

# Split off the supplementary legend / supplementary tables if any (they are part of same md).
# OAC word limit applies to the full article body; we count the whole main md.
# Count words = whitespace-separated tokens containing at least one letter.
tokens = re.findall(r"[A-Za-z][A-Za-z0-9'\-]*", t)
n_words = len(tokens)
print("main md word count (alpha tokens):", n_words)

# Also count excluding the figure-legend block (## Figure legends ... through Graphical Abstract)
# and excluding supplementary table blocks, to estimate the "prose" body.
idx_leg = t.find("## Figure legends")
body_prose = t[:idx_leg] if idx_leg != -1 else t
ptokens = re.findall(r"[A-Za-z][A-Za-z0-9'\-]*", body_prose)
print("prose before figure legends (words):", len(ptokens))

# References block
idx_ref = t.find("## References")
if idx_ref != -1:
    refblock = t[idx_ref:]
    rtok = re.findall(r"[A-Za-z][A-Za-z0-9'\-]*", refblock)
    print("references block words:", len(rtok))
