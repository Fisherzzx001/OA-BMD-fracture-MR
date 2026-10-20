#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""156 · Authoritative test of *write* clearance on the git protocol.

Context: the sandbox proxy intermittently returns 502 on CONNECT to github.com,
so every probe must retry until the host is reachable before drawing conclusions.
Uses HTTP Basic auth exactly as `git push` does.
"""
import re
import subprocess
import time

TOKEN_DOC = r"C:\Users\Administrator\github_token.doc"
OWNER, REPO = "Fisherzzx001", "OA-BMD-fracture-MR"
UA = "git/2.55.0"


def tok_from(path, prefer_utf8=True):
    raw = open(path, "rb").read()
    cands = [raw.decode("utf-8", "ignore"), raw.decode("latin-1", "ignore")]
    pat = re.compile(r"(github_pat_[A-Za-z0-9_]{20,}|gh[pousr]_[A-Za-z0-9]{20,})")
    for t in cands:
        m = pat.search(t)
        if m:
            return m.group(1)
    raise SystemExit("no token in " + path)


def curl_status(url, userpwd=None, tries=8, wait=2.0):
    """Retry past proxy 502 flakes; return the first meaningful HTTP status."""
    last = "000"
    for i in range(tries):
        cmd = ["curl", "-sS", "-o", "/dev/null", "-w", "%{http_code}", "-H", f"User-Agent: {UA}"]
        if userpwd:
            cmd += ["-u", userpwd]
        cmd.append(url)
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=60).stdout.strip()
        if out and out != "000":
            return out, i + 1
        last = out
        time.sleep(wait)
    return last, tries


def main():
    tok = tok_from(TOKEN_DOC)
    base = f"https://github.com/{OWNER}/{REPO}.git/info/refs"
    creds = f"{OWNER}:{tok}"

    print("--- anonymous read (sanity: is the host reachable at all?) ---")
    st, n = curl_status(base + "?service=git-upload-pack")
    print(f"  upload-pack  anon   -> HTTP {st}  (after {n} attempt(s))")

    print("\n--- WITH this token (Basic auth) ---")
    st, n = curl_status(base + "?service=git-upload-pack", creds)
    print(f"  upload-pack  auth   -> HTTP {st}  (after {n} attempt(s))")
    st2, n2 = curl_status(base + "?service=git-receive-pack", creds)
    print(f"  receive-pack auth   -> HTTP {st2}  (after {n2} attempt(s))")

    print("\n--- interpretation ---")
    if st2 == "200":
        print("  receive-pack 200 => this credential MAY WRITE over git; any 403 from")
        print("  `git push` therefore comes from the push itself (branch policy), not auth.")
    elif st2 in ("401", "403"):
        print("  receive-pack denied => the token lacks git write clearance even though the")
        print("  REST API reported permissions.push=True. Regenerate with 'Contents: Read and write'.")
    else:
        print("  receive-pack inconclusive (network flake). Re-run.")


if __name__ == "__main__":
    main()
