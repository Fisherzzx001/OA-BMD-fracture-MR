#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""154 · Probe whether ANY git push works, then perform the force push.

Diagnostic context (2026-10-21): API reports permissions.push = True and branch
`main` is unprotected, yet `git push` returns 403. This probes a harmless
throwaway branch first to separate "credentials/transport" from "force-push
policy", disables the credential helper (may hold a stale cached token), and
redacts the token from every line of output.
"""
import os
import re
import subprocess
import sys
import urllib.error
import urllib.request

REPO_DIR = r"C:\Users\Administrator\WorkBuddy\2026-09-08-22-17-40\OA_B_PROJECT"
TOKEN_DOC = r"C:\Users\Administrator\github_token.doc"
OWNER, REPO = "Fisherzzx001", "OA-BMD-fracture-MR"
PROBE = "refs/heads/probe-write-check"


def extract_token(path):
    raw = open(path, "rb").read()
    cands = []
    if raw[:4] == b"PK\x03\x04":
        import zipfile
        cands.append(re.sub(r"<[^>]+>", " ",
                            zipfile.ZipFile(path).read("word/document.xml").decode("utf-8", "ignore")))
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


def git(args, tok, url=None):
    """Run git with credential helpers disabled so the URL token is authoritative."""
    env = dict(os.environ)
    env["GIT_TERMINAL_PROMPT"] = "0"
    env["GIT_ASKPASS"] = "true"
    cmd = ["git", "-c", "credential.helper=", "-c", "credential.autoCache=false"] + args
    p = subprocess.run(cmd, cwd=REPO_DIR, capture_output=True, text=True, env=env)
    out = (p.stdout or "") + (p.stderr or "")
    out = out.replace(tok, "***TOKEN***")
    return p.returncode, out.strip()


def main():
    tok = extract_token(TOKEN_DOC)
    auth_url = f"https://{OWNER}:{tok}@github.com/{OWNER}/{REPO}.git"

    print("=== PROBE 1: plain (non-force) push to throwaway branch ===")
    rc, out = git(["push", auth_url, f"master:{PROBE}"], tok)
    print("exit:", rc)
    print(out[:1200])
    if rc == 0:
        print(">>> plain push WORKS -> credentials fine; issue was specific to main/force-push")
        # clean up probe branch
        rc2, out2 = git(["push", auth_url, "--delete", PROBE], tok)
        print("probe cleanup:", "OK" if rc2 == 0 else out2[:300])
    else:
        print(">>> plain push ALSO fails -> credential/transport problem, not branch policy")
        return 1

    print("\n=== PROBE 2: FORCE push master -> main ===")
    rc, out = git(["push", "--force", auth_url, "master:refs/heads/main"], tok)
    print("exit:", rc)
    print(out[:2000])
    return rc


if __name__ == "__main__":
    sys.exit(main())
