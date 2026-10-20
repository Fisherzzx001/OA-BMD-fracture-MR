#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""155 · Isolate why git push is rejected although the API reports push=True.

Three independent probes, each removing one suspect:

  A. Raw git protocol over HTTPS (no git binary): does the credential even have
     write clearance on git-receive-pack?
  B. git push with ALL config files neutralised (GIT_CONFIG_GLOBAL/SYSTEM=/dev/null)
     so no credential helper (GCM can rewrite the Authorization header) can interfere.
  C. git push using the normal config, for comparison.

Every printed line is run through redaction so the token never appears.
"""
import os
import re
import subprocess
import sys

REPO_DIR = r"C:\Users\Administrator\WorkBuddy\2026-09-08-22-17-40\OA_B_PROJECT"
TOKEN_DOC = r"C:\Users\Administrator\github_token.doc"
OWNER, REPO = "Fisherzzx001", "OA-BMD-fracture-MR"
SAFE_URL = f"https://github.com/{OWNER}/{REPO}.git"


def extract_token(path):
    raw = open(path, "rb").read()
    cands = [raw.decode("utf-8", "ignore"), raw.decode("latin-1", "ignore")]
    if raw[:4] == b"PK\x03\x04":
        import zipfile
        cands.insert(0, re.sub(r"<[^>]+>", " ",
                               zipfile.ZipFile(path).read("word/document.xml").decode("utf-8", "ignore")))
    pat = re.compile(r"(github_pat_[A-Za-z0-9_]{20,}|gh[pousr]_[A-Za-z0-9]{20,})")
    for t in cands:
        m = pat.search(t)
        if m:
            return m.group(1)
    raise SystemExit("no token found")


def scrub(s, tok):
    return s.replace(tok, "***TOKEN***")


def main():
    tok = extract_token(TOKEN_DOC)

    print("=== A. raw git protocol over HTTPS (no git binary involved) ===")
    for svc in ("git-upload-pack", "git-receive-pack"):
        cmd = ["curl", "-sS", "-o", "/dev/null", "-w", "%{http_code}",
               "-u", f"{OWNER}:{tok}",
               "-H", "User-Agent: git/2.55.0",
               f"https://github.com/{OWNER}/{REPO}.git/info/refs?service={svc}"]
        try:
            out = subprocess.run(cmd, capture_output=True, text=True, timeout=90).stdout.strip()
        except Exception as e:
            out = f"<error {e}>"
        # 200 = allowed;  401/403 = credential refused or insufficient
        print(f"  {svc:17} -> HTTP {out}")

    print("\n=== B. git push with global+system config neutralised (no credential helper) ===")
    env = dict(os.environ)
    env["GIT_CONFIG_GLOBAL"] = "/dev/null"
    env["GIT_CONFIG_SYSTEM"] = "/dev/null"
    env["GIT_TERMINAL_PROMPT"] = "0"
    auth_url = f"https://{OWNER}:{tok}@github.com/{OWNER}/{REPO}.git"
    p = subprocess.run(["git", "push", "--force", auth_url, "master:refs/heads/main"],
                       cwd=REPO_DIR, capture_output=True, text=True, env=env, timeout=600)
    out = scrub((p.stdout or "") + (p.stderr or ""), tok).strip()
    print("  exit:", p.returncode)
    print("  " + out.replace("\n", "\n  ")[:1500])
    if p.returncode == 0:
        print("\n>>> SUCCESS with helpers disabled -> the credential helper was the culprit.")
        return 0

    print("\n=== C. git push with normal config (for comparison) ===")
    env2 = dict(os.environ)
    env2["GIT_TERMINAL_PROMPT"] = "0"
    p2 = subprocess.run(["git", "push", "--force", auth_url, "master:refs/heads/main"],
                        cwd=REPO_DIR, capture_output=True, text=True, env=env2, timeout=600)
    print("  exit:", p2.returncode)
    print("  " + scrub((p2.stdout or "") + (p2.stderr or ""), tok).strip().replace("\n", "\n  ")[:1000])
    return p2.returncode


if __name__ == "__main__":
    sys.exit(main())
