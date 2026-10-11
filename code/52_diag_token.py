#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Diagnostic: locate token-like material in the .doc WITHOUT printing it."""
import re
PATH = r"C:\Users\Administrator\github_token.doc"
raw = open(PATH, "rb").read()
print("file size:", len(raw))
print("magic:", raw[:8].hex())

def redact(s):
    s = s.strip()
    if len(s) < 8:
        return "***(" + str(len(s)) + "ch)***"
    return s[:4] + "…" + s[-4:]

# 1) raw latin-1 scan for github-ish keywords
for kw in [b"ghp_", b"github", b"gho_", b"ghu_", b"ghs_", b"ghr_", b"pat_", b"token", b"PAT", b"TOKEN"]:
    idx = 0
    hits = 0
    while True:
        i = raw.find(kw, idx)
        if i < 0:
            break
        seg = raw[i:i+60]
        # take printable ascii run
        run = re.sub(rb"[^\x20-\x7e]", b" ", seg).decode("ascii", "ignore").strip()
        print(f"  kw={kw.decode()} @ {i}: {redact(run)}")
        hits += 1
        idx = i + 1
        if hits >= 3:
            break

# 2) UTF-16LE scan (Word sometimes stores unicode)
try:
    u = raw.decode("utf-16-le", "ignore")
    for kw in ["ghp_", "github", "github_pat_", "token"]:
        i = u.find(kw)
        if i >= 0:
            seg = re.sub(r"[^ -~]", " ", u[i:i+60]).strip()
            print(f"  utf16 kw={kw} @ {i}: {redact(seg)}")
except Exception as e:
    print("utf16 decode err:", e)

# 3) longest printable ascii runs overall (top 10) redacted
runs = re.findall(rb"[\x20-\x7e]{12,}", raw)
runs = sorted(runs, key=len, reverse=True)[:10]
print("--- top ascii runs (redacted) ---")
for r in runs:
    s = r.decode("ascii", "ignore")
    print("  ", redact(s), f"({len(s)}ch)")
