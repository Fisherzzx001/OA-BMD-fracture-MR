#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Diagnostic: report token scopes + whether target repo already exists. Read-only."""
import re, json, urllib.request, urllib.error
PATH = r"C:\Users\Administrator\github_token.doc"
API = "https://api.github.com"
REPO = "OA-BMD-fracture-MR"

raw = open(PATH, "rb").read()
text = raw.decode("utf-16-le", "ignore")
m = re.search(r"(github_pat_[A-Za-z0-9_]{20,}|gh[pousr]_[A-Za-z0-9]{20,})", text)
tok = m.group(1)

def get(url):
    req = urllib.request.Request(url)
    req.add_header("Authorization", "Bearer " + tok)
    req.add_header("User-Agent", "oa-check")
    req.add_header("Accept", "application/vnd.github+json")
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            return r.status, dict(r.getheaders()), json.loads(r.read().decode() or "{}")
    except urllib.error.HTTPError as e:
        return e.code, dict(e.headers), {}

st, hdrs, body = get(API + "/user")
print("GET /user status:", st)
print("  X-OAuth-Scopes:", hdrs.get("X-OAuth-Scopes", "(none)"))
print("  X-Accepted-OAuth-Scopes:", hdrs.get("X-Accepted-OAuth-Scopes", "(none)"))
login = body.get("login")
print("  login:", login)

st2, _, body2 = get(f"{API}/repos/{login}/{REPO}")
print(f"GET /repos/{login}/{REPO} status:", st2)
if st2 == 200:
    print("  repo EXISTS -> can push directly")
    print("  private:", body2.get("private"), "| default_branch:", body2.get("default_branch"))
elif st2 == 404:
    print("  repo does NOT exist -> need 'repo'/'public_repo' scope to create it")
else:
    print("  unexpected:", st2)
