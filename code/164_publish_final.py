#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""164 · Final publish: push the corrected commit and mint a clean Zenodo version.

 Why a second release (v1.2) after v1.1?
   v1.1 was minted from the pre-correction snapshot, so its archive still carries
   the old (and, after the git cleanup, inaccurate) Data-availability wording.
   v1.2 archives the corrected text; the manuscript now cites the CONCEPT DOI,
   which stays valid across versions, so no further DOI churn is needed.

Steps:
   1. plain push master -> main (histories share a base now, so no force needed)
   2. create Release v1.2
   3. attach Submission_Package_2026-10-17.zip as a release asset
The token is injected into the remote only during step 1 and always scrubbed in
a finally block; it is also redacted from every printed line.
"""
import json
import os
import re
import subprocess
import sys
import time
import urllib.error
import urllib.request

REPO_DIR = r"C:\Users\Administrator\WorkBuddy\2026-09-08-22-17-40\OA_B_PROJECT"
TOKEN_DOC = r"C:\Users\Administrator\github_token.doc"
ZIP = os.path.join(REPO_DIR, "Submission_Package_2026-10-17.zip")
CLEAN_REMOTE = "https://github.com/Fisherzzx001/OA-BMD-fracture-MR.git"
OWNER, REPO = "Fisherzzx001", "OA-BMD-fracture-MR"
API = "https://api.github.com"


def extract(path):
    raw = open(path, "rb").read()
    pat = re.compile(r"(github_pat_[A-Za-z0-9_]{20,}|gh[pousr]_[A-Za-z0-9]{20,})")
    for enc in ("utf-8", "utf-16-le", "latin-1"):
        try:
            m = pat.search(raw.decode(enc, "ignore"))
            if m:
                return m.group(1)
        except Exception:
            pass
    raise SystemExit("no token")


def api(method, url, tok, data=None, raw=None, ctype=None):
    req = urllib.request.Request(url, method=method)
    req.add_header("Accept", "application/vnd.github+json")
    req.add_header("User-Agent", "final-publish")
    req.add_header("Authorization", "Bearer " + tok)
    body = None
    if data is not None:
        req.add_header("Content-Type", "application/json")
        body = json.dumps(data).encode()
    elif raw is not None:
        req.add_header("Content-Type", ctype or "application/octet-stream")
        body = raw
    if body is not None:
        req.data = body
    try:
        with urllib.request.urlopen(req, timeout=180) as r:
            return r.status, json.loads(r.read().decode() or "{}")
    except urllib.error.HTTPError as e:
        b = e.read().decode("utf-8", "ignore")
        try:
            return e.code, {"message": json.loads(b).get("message", b)}
        except Exception:
            return e.code, {"message": b[:300]}


def git(args, tok):
    env = dict(os.environ)
    env["GIT_TERMINAL_PROMPT"] = "0"
    # retry past the sandbox proxy's intermittent 502 on CONNECT
    for attempt in range(6):
        p = subprocess.run(["git"] + args, cwd=REPO_DIR, capture_output=True,
                           text=True, env=env, timeout=600)
        out = ((p.stdout or "") + (p.stderr or "")).replace(tok, "***TOKEN***").strip()
        if p.returncode == 0:
            return 0, out, attempt + 1
        if "CONNECT tunnel" in out or "000" in out:
            time.sleep(5)
            continue
        return p.returncode, out, attempt + 1
    return p.returncode, out, 6


def main():
    tok = extract(TOKEN_DOC)
    auth_url = f"https://{OWNER}:{tok}@github.com/{OWNER}/{REPO}.git"

    print("=== 1. push corrected commit ===")
    try:
        subprocess.run(["git", "remote", "set-url", "origin", auth_url],
                       cwd=REPO_DIR, check=True, capture_output=True,
                       env={**os.environ, "GIT_TERMINAL_PROMPT": "0"})
        rc, out, attempts = git(["push", "origin", "master:refs/heads/main"], tok)
        print(f"exit={rc} (after {attempts} attempt(s))")
        print(out[:800])
        if rc != 0:
            return 1
    finally:
        subprocess.run(["git", "remote", "set-url", "origin", CLEAN_REMOTE],
                       cwd=REPO_DIR, capture_output=True)
        print("token scrubbed from remote: OK")

    rc, out = subprocess.run(["git", "ls-remote", "--heads", "origin", "main"],
                             cwd=REPO_DIR, capture_output=True, text=True).returncode, None
    head = subprocess.run(["git", "ls-remote", "--heads", "origin", "main"],
                          cwd=REPO_DIR, capture_output=True, text=True).stdout.strip()
    print("remote main now:", head.split()[0] if head else "?")

    print("\n=== 2. create Release v1.2 ===")
    body = (
        "Corrected, OAC-compliant version.\n\n"
        "- Fixes three independent root causes: rg-matrix reproducibility, effect-scale "
        "propagation (log-OR factor 1/sqrt(K(1-K)) = 3.389), and the weighted-median "
        "step-near convention.\n"
        "- Corrected figures and Table 2 replace the earlier OR 1.017-1.028 values.\n"
        "- Manuscript body 3973 words (OAC limit 4000); all reported beta/SE/p/OR/CI "
        "verified unchanged.\n"
        "- Repository history cleaned: regenerable bulky intermediates removed "
        "(221 MB -> 26 MB); content verified byte-identical.\n\n"
        "- Attached asset `Submission_Package_2026-10-17.zip` contains the full "
        "submission package (manuscript docx, 19 main + 9 supplementary figures, "
        "17 supplementary files and 19 table-source CSVs).\n"
    )
    st, rel = api("POST", API + f"/repos/{OWNER}/{REPO}/releases", tok, {
        "tag_name": "v1.2",
        "name": "v1.2 - corrected figures, Table 2 and OAC-compliant manuscript",
        "body": body, "draft": False, "prerelease": False,
    })
    if st not in (200, 201):
        print("release creation failed:", st, rel.get("message"))
        return 1
    print("release created:", rel.get("html_url"))

    print("\n=== 3. attach submission package asset ===")
    if not os.path.exists(ZIP):
        print("zip not found:", ZIP)
        return 1
    upload_url = (rel.get("upload_url") or "").split("{")[0]
    with open(ZIP, "rb") as f:
        data = f.read()
    st, up = api("POST", upload_url + "?name=" + os.path.basename(ZIP), tok,
                 raw=data, ctype="application/zip")
    if st in (200, 201):
        print(f"asset uploaded: {up.get('name')} ({up.get('size')} bytes)")
        print("download url:", up.get("browser_download_url"))
    else:
        print("asset upload failed:", st, up.get("message"))
        return 1

    print("\nZenodo will mint a new version shortly; the manuscript cites the concept")
    print("DOI https://doi.org/10.5281/zenodo.23400860 which resolves to the latest version.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
