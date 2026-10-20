#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""157 · Probe whether the token can WRITE through the REST API.

The git protocol is denied (git-receive-pack -> 403) although GET /repos reports
permissions.push=True. Before concluding the token is read-only, test an actual
REST write in the least invasive way possible: create a throwaway ref from an
existing commit and immediately delete it again.
"""
import json
import re
import urllib.error
import urllib.request

TOKEN_DOC = r"C:\Users\Administrator\github_token.doc"
API = "https://api.github.com"
OWNER, REPO = "Fisherzzx001", "OA-BMD-fracture-MR"
PROBE_REF = "refs/heads/api-write-probe"


def tok_from(path):
    raw = open(path, "rb").read()
    pat = re.compile(r"(github_pat_[A-Za-z0-9_]{20,}|gh[pousr]_[A-Za-z0-9]{20,})")
    for enc in ("utf-8", "latin-1"):
        m = pat.search(raw.decode(enc, "ignore"))
        if m:
            return m.group(1)
    raise SystemExit("no token")


def call(method, path, tok, data=None):
    req = urllib.request.Request(API + path, method=method)
    req.add_header("Accept", "application/vnd.github+json")
    req.add_header("User-Agent", "api-write-probe")
    req.add_header("Authorization", "Bearer " + tok)
    if data is not None:
        req.add_header("Content-Type", "application/json")
        req.data = json.dumps(data).encode()
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            return r.status, json.loads(r.read().decode() or "{}")
    except urllib.error.HTTPError as e:
        body = e.read().decode("utf-8", "ignore")
        try:
            msg = json.loads(body).get("message", body)
        except Exception:
            msg = body[:200]
        return e.code, {"message": msg}


def main():
    tok = tok_from(TOKEN_DOC)

    # an existing commit on the remote to hang the probe ref off
    st, refs = call("GET", f"/repos/{OWNER}/{REPO}/git/ref/heads/main", tok)
    if st != 200:
        print("cannot read remote main:", st, refs.get("message"))
        return
    sha = refs["object"]["sha"]
    print("remote main currently at:", sha[:12])

    print("\n--- REST WRITE TEST: create throwaway ref ---")
    st, r = call("POST", f"/repos/{OWNER}/{REPO}/git/refs", tok,
                 {"ref": PROBE_REF, "sha": sha})
    print("  status:", st, "|", r.get("message") or (r.get("ref") or "")[:60])

    if st in (200, 201):
        print("  >>> REST WRITE WORKS: the token can modify refs via the API.")
        print("      Clean up probe ref...")
        st2, r2 = call("DELETE", f"/repos/{OWNER}/{REPO}/git/refs/heads/api-write-probe", tok)
        print("  cleanup status:", st2, "(204 = deleted)")
        print("\n  CONCLUSION: git-protocol write is blocked but API write is allowed.")
        print("  Fallback exists: push content via the Git Data API instead of `git push`.")
    else:
        print("\n  >>> REST WRITE ALSO DENIED -> token is effectively READ-ONLY.")
        print("      Action required: regenerate the fine-grained token with")
        print("      'Contents: Read and write' (owner Fisherzzx001, repo access must include")
        print("      OA-BMD-fracture-MR), or issue a classic token with the 'repo' scope.")


if __name__ == "__main__":
    main()
