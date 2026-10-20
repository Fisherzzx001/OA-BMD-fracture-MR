#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""153 · Diagnose write-permission failure on push.

`git push` returning "Permission to <repo> denied to <user>" (403) means the
credential authenticated but is not allowed to write. Causes differ:
  (a) the token lacks Contents: write  -> check repo.permissions.push
  (b) branch protection / ruleset blocks force-push -> check protection state
This script reports both using read-only API calls (no side effects).
"""
import json
import re
import urllib.error
import urllib.request

TOKEN_DOC = r"C:\Users\Administrator\github_token.doc"
API = "https://api.github.com"
OWNER = "Fisherzzx001"
REPO = "OA-BMD-fracture-MR"


def extract_token(path):
    raw = open(path, "rb").read()
    cands = []
    if raw[:4] == b"PK\x03\x04":
        import zipfile
        z = zipfile.ZipFile(path)
        cands.append(re.sub(r"<[^>]+>", " ", z.read("word/document.xml").decode("utf-8", "ignore")))
    try:
        cands.append(raw.decode("utf-16-le", "ignore"))
    except Exception:
        pass
    cands.append(raw.decode("latin-1", "ignore"))
    pat = re.compile(r"(github_pat_[A-Za-z0-9_]{20,}|gh[pousr]_[A-Za-z0-9]{20,})")
    for t in cands:
        m = pat.search(t)
        if m:
            return m.group(1)
    raise SystemExit("no token found")


def call(method, url, tok, data=None):
    req = urllib.request.Request(url, method=method)
    req.add_header("Accept", "application/vnd.github+json")
    req.add_header("User-Agent", "perm-diag")
    req.add_header("Authorization", "Bearer " + tok)
    if data is not None:
        req.add_header("Content-Type", "application/json")
        req.data = json.dumps(data).encode()
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            return r.status, json.loads(r.read().decode() or "{}")
    except urllib.error.HTTPError as e:
        body = e.read().decode("utf-8", "ignore")
        try:
            msg = json.loads(body).get("message", body)
        except Exception:
            msg = body
        return e.code, {"message": msg}


def main():
    tok = extract_token(TOKEN_DOC)
    print("token prefix:", tok[:11], "| suffix:", tok[-4:],
          "| type:", "fine-grained" if tok.startswith("github_pat_") else "classic")

    st, me = call("GET", API + "/user", tok)
    if st != 200:
        print("auth FAILED:", st, me.get("message"))
        return
    print("authenticated as:", me["login"])

    # (a) can this credential push?  permissions.push is authoritative
    st, r = call("GET", API + f"/repos/{OWNER}/{REPO}", tok)
    perms = r.get("permissions") or {}
    print("\n--- repo permissions for this token ---")
    for k in ("admin", "maintain", "push", "triage", "pull"):
        print(f"  {k:9}: {perms.get(k)}")
    if st != 200:
        print("repo lookup failed:", st, r.get("message"))
        return

    if not perms.get("push"):
        print("\n>>> ROOT CAUSE (a): token has NO Contents/Write permission on this repo.")
        print("    Fix: regenerate the fine-grained token with 'Contents: Read and write'")
        print("    (owner must be Fisherzzx001, repository access must include this repo),")
        print("    or use a classic token with the 'repo' scope.")
        return

    print("\npush permission present -> failure is NOT a token-permission issue")

    # (b) is force-push blocked by protection / rulesets?
    st, br = call("GET", API + f"/repos/{OWNER}/{REPO}/branches/main", tok)
    print("\n--- branch 'main' protection state ---")
    if st == 200:
        print("  protected:", br.get("protected"))
        pr = br.get("protection") or {}
        for key in ("allow_force_pushes", "required_status_checks", "restrictions",
                    "required_pull_request_reviews", "enforce_admins", "required_linear_history"):
            v = pr.get(key)
            if isinstance(v, dict):
                v = {k2: v[k2].get("enabled") for k2 in v if isinstance(v[k2], dict)}
            if v not in (None, {}, False):
                print(f"  {key}: {v}")
    else:
        print("  lookup failed:", st, br.get("message"))

    st, rs = call("GET", API + f"/repos/{OWNER}/{REPO}/rulesets?includes_parents=true", tok)
    if st == 200 and rs:
        print("  rulesets:")
        for rule in rs:
            print("   -", rule.get("name"), "| enforcement:", rule.get("enforcement"),
                  "| id:", rule.get("id"))


if __name__ == "__main__":
    main()
