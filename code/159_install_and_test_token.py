#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""159 · Install a (classic) token and verify WRITE clearance before publishing.

Order matters: prove the credential can write FIRST, then let code/151 publish.
This avoids discovering a 403 halfway through. Nothing secret is printed.
"""
import os
import re
import shutil
import subprocess
import sys
import time

TOKEN_DOC = r"C:\Users\Administrator\github_token.doc"
OWNER, REPO = "Fisherzzx001", "OA-BMD-fracture-MR"
UA = "git/2.55.0"

TOKEN = os.environ.get("NEW_TOKEN") or sys.argv[1] if len(sys.argv) > 1 else None


def mask(s):
    return (s[:4] + "..." + s[-4:]) if s and len(s) > 10 else "(none)"


def retry_curl(url, userpwd=None, tries=10, wait=2.0):
    last = "000"
    for _ in range(tries):
        cmd = ["curl", "-sS", "-o", "/dev/null", "-w", "%{http_code}", "-H", f"User-Agent: {UA}"]
        if userpwd:
            cmd += ["-u", userpwd]
        cmd.append(url)
        try:
            out = subprocess.run(cmd, capture_output=True, text=True, timeout=60).stdout.strip()
        except Exception:
            out = "000"
        if out and out != "000":
            return out
        last = out
        time.sleep(wait)
    return last


def main():
    tok = TOKEN
    if not tok:
        raise SystemExit("no token supplied")
    print("token length:", len(tok), "| masked:", mask(tok),
          "| type:", "classic" if tok.startswith("ghp_") else "fine-grained")

    if os.path.exists(TOKEN_DOC):
        bak = TOKEN_DOC + ".prev.bak"
        shutil.copy2(TOKEN_DOC, bak)
        print("previous token backed up ->", bak)
    open(TOKEN_DOC, "w", encoding="utf-8").write(tok)
    print("installed to", TOKEN_DOC)

    base = f"https://github.com/{OWNER}/{REPO}.git/info/refs"
    creds = f"{OWNER}:{tok}"

    print("\n=== write-clearance probe (with retry past proxy 502 flakes) ===")
    up = retry_curl(base + "?service=git-upload-pack", creds)
    rp = retry_curl(base + "?service=git-receive-pack", creds)
    print(f"  git-upload-pack  -> HTTP {up}   (read)")
    print(f"  git-receive-pack -> HTTP {rp}   (write)")

    if rp == "200":
        print("\n>>> CLEARED: this token can write over the git protocol.")
        print("    Now run: code/151_publish_github_release.py")
        return 0
    print("\n>>> BLOCKED (expected 200). Do NOT attempt the push yet.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
