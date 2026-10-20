#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""152 · Diagnose GitHub token type + OAuth scopes (never prints the token)."""
import re, urllib.request, urllib.error

PATH = r"C:\Users\Administrator\github_token.doc"
raw = open(PATH, "rb").read()
cands = []
try:
    cands.append(raw.decode("utf-16-le", "ignore"))
except Exception:
    pass
cands.append(raw.decode("latin-1", "ignore"))
pat = re.compile(r"(github_pat_[A-Za-z0-9_]{20,}|gh[pousr]_[A-Za-z0-9]{20,})")
tok = None
for t in cands:
    m = pat.search(t)
    if m:
        tok = m.group(1)
        break

if not tok:
    print("no token found")
    raise SystemExit(1)

kind = "fine-grained (github_pat_)" if tok.startswith("github_pat_") else "classic (ghp_)"
print("token type :", kind, "| prefix:", tok[:4])
print("token len  :", len(tok))

req = urllib.request.Request("https://api.github.com/user")
req.add_header("Authorization", "Bearer " + tok)
req.add_header("User-Agent", "scope-diag")
try:
    with urllib.request.urlopen(req, timeout=30) as r:
        import json
        info = json.loads(r.read().decode())
        print("login      :", info.get("login"))
        print("type       :", info.get("type"))
        print("X-OAuth-Scopes        :", r.headers.get("X-OAuth-Scopes", "(empty)"))
        print("X-Accepted-OAuth-Scopes:", r.headers.get("X-Accepted-OAuth-Scopes", "(n/a)"))
except urllib.error.HTTPError as e:
    print("HTTPError", e.code, e.read().decode()[:300])
