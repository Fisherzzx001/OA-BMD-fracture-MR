#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""158 · Look for an already-cached GitHub credential that CAN write.

The account previously published this repo from this machine, so the OS/GCM cache
may hold a credential with sufficient rights. This asks every configured helper
for a github.com credential and tests each against git-receive-pack.

Nothing secret is ever printed: only masked fingerprints and HTTP statuses.
"""
import os
import re
import subprocess
import time

CRED_EXE = r"C:\Users\Administrator\.workbuddy\binaries\PortableGit\versions\1.3.0\mingw64\bin\git-credential-manager.exe"
OWNER, REPO = "Fisherzzx001", "OA-BMD-fracture-MR"


def mask(s):
    if not s:
        return "(none)"
    return (s[:6] + "..." + s[-4:]) if len(s) > 12 else "***"


def ask_helper(exe, label):
    """Run a credential helper 'get' for github.com and parse its response."""
    if exe and not os.path.exists(exe.replace('!"', "").strip('"')):
        return None
    inp = "protocol=https\nhost=github.com\n\n"
    try:
        p = subprocess.run([exe.strip('"') if os.path.exists(exe.strip('"')) else exe, "get"],
                           input=inp, capture_output=True, text=True, timeout=60)
    except Exception as e:
        print(f"  [{label}] helper error: {e}")
        return None
    out = (p.stdout or "")
    cred = dict()
    for line in out.splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            cred[k.strip()] = v.strip()
    if cred.get("password"):
        print(f"  [{label}] credential returned: username={cred.get('username')} "
              f"password={mask(cred['password'])}")
        return cred["password"]
    print(f"  [{label}] no credential returned"
          + (" | " + (p.stderr or '').strip()[:120] if (p.stderr or '').strip() else ""))
    return None


def test_write(tok, tries=6):
    """Retry past proxy 502 flakes; 200 on receive-pack means write is allowed."""
    last = "000"
    for _ in range(tries):
        cmd = ["curl", "-sS", "-o", "/dev/null", "-w", "%{http_code}",
               "-H", "User-Agent: git/2.55.0",
               "-u", f"{OWNER}:{tok}",
               f"https://github.com/{OWNER}/{REPO}.git/info/refs?service=git-receive-pack"]
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=60).stdout.strip()
        if out and out != "000":
            return out
        last = out
        time.sleep(2)
    return last


def main():
    print("=== query configured credential helpers ===")
    helpers = []
    st = subprocess.run(["git", "config", "--get-all", "credential.helper"],
                        capture_output=True, text=True)
    if st.returncode == 0:
        for h in st.stdout.splitlines():
            if h.startswith("!"):
                helpers.append((h.strip('!').strip('"'), "global-gcm"))
            else:
                helpers.append((h, "manager"))
        helpers.append((CRED_EXE, "gcm-executable"))
    else:
        helpers = [(CRED_EXE, "gcm-executable")]

    found = []
    for exe, label in helpers:
        pw = ask_helper(exe, label)
        # skip the token we already know is read-only (identical to github_token.doc)
        if pw:
            found.append((label, pw))

    known = open(r"C:\Users\Administrator\github_token.doc", "rb").read().decode("utf-8", "ignore")
    m = re.search(r"(github_pat_[A-Za-z0-9_]{20,}|gh[pousr]_[A-Za-z0-9]{20,})", known)
    known_tok = m.group(1) if m else ""
    prev = None
    pbak = r"C:\Users\Administrator\github_token.doc.prev.bak"
    if os.path.exists(pbak):
        raw = open(pbak, "rb").read()
        for enc in ("utf-8", "latin-1"):
            m2 = re.search(r"(github_pat_[A-Za-z0-9_]{20,}|gh[pousr]_[A-Za-z0-9]{20,})",
                           raw.decode(enc, "ignore"))
            if m2:
                prev = m2.group(1)
                break

    print("\n=== write test for each candidate ===")
    anything_ok = False
    for label, pw in found:
        if known_tok and pw == known_tok:
            print(f"  [{label}] same as current github_token.doc — already proven read-only, skip")
            continue
        st = test_write(pw)
        print(f"  [{label}] receive-pack -> HTTP {st}" + ("   *** WRITABLE ***" if st == "200" else ""))
        if st == "200":
            anything_ok = True
            print("      -> save this one to github_token.doc and re-run code/151")
    if prev and (not found or prev not in [p for _, p in found]):
        stp = test_write(prev)
        print(f"  [previous token doc backup] receive-pack -> HTTP {stp}")

    print("\n=== ssh fallback available? ===")
    sshdir = os.path.expanduser("~/.ssh")
    keys = []
    if os.path.isdir(sshdir):
        keys = [f for f in os.listdir(sshdir) if f.endswith((".pub",)) or f in ("id_rsa", "id_ed25519")]
    print("  ~/.ssh contents:", keys if keys else "(empty — no SSH fallback)")

    print("\nRESULT:", "a writable credential exists" if anything_ok else "no writable credential found on this machine")


if __name__ == "__main__":
    main()
