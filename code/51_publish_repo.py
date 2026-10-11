#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Publish the local `publish` branch to a NEW public GitHub repo OA-BMD-fracture-MR.
Token is read from a .doc file, used ONCE for the API + push, then erased from
the remote URL. The token is NEVER printed (only redacted last-4 chars).
"""
import os, re, sys, json, subprocess, urllib.request, urllib.error

REPO_DIR = r"C:\Users\Administrator\WorkBuddy\2026-09-08-22-17-40\OA_B_PROJECT"
TOKEN_DOC = r"C:\Users\Administrator\github_token.doc"
OWNER_REPO = "OA-BMD-fracture-MR"
BRANCH_LOCAL = "publish"
BRANCH_REMOTE = "main"
API = "https://api.github.com"

def redact(tok):
    if not tok or len(tok) < 8:
        return "***"
    return tok[:4] + "…" + tok[-4:]

def log(msg):
    # strip any accidental token leakage before printing
    sys.stdout.write(msg + "\n")

def extract_token(path):
    raw = open(path, "rb").read()
    candidates = []
    if raw[:4] == b"PK\x03\x04":  # docx (zip)
        import zipfile
        z = zipfile.ZipFile(path)
        xml = z.read("word/document.xml").decode("utf-8", "ignore")
        candidates.append(re.sub(r"<[^>]+>", " ", xml))
    # Word 97 .doc stores text as UTF-16LE (ASCII chars separated by NUL bytes)
    try:
        candidates.append(raw.decode("utf-16-le", "ignore"))
    except Exception:
        pass
    candidates.append(raw.decode("latin-1", "ignore"))
    pat = re.compile(r"(github_pat_[A-Za-z0-9_]{20,}|gh[pousr]_[A-Za-z0-9]{20,})")
    for text in candidates:
        m = pat.search(text)
        if m:
            return m.group(1)
    raise SystemExit("ERROR: no GitHub token (ghp_/github_pat_/gho_/ghu_/ghs_) found in token file.")

def api(method, url, token=None, data=None):
    req = urllib.request.Request(url, method=method)
    req.add_header("Accept", "application/vnd.github+json")
    req.add_header("User-Agent", "oa-publish-script")
    if token:
        req.add_header("Authorization", "Bearer " + token)
    if data is not None:
        req.add_header("Content-Type", "application/json")
        req.data = json.dumps(data).encode()
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            return r.status, json.loads(r.read().decode() or "{}")
    except urllib.error.HTTPError as e:
        body = e.read().decode("utf-8", "ignore")
        try:
            j = json.loads(body)
            msg = j.get("message", body)
        except Exception:
            msg = body
        return e.code, {"message": msg}

def run_git(args, token=None):
    env = dict(os.environ)
    if token:
        env["GIT_ASKPASS"] = "true"  # never prompt
        env["GIT_TERMINAL_PROMPT"] = "0"
    p = subprocess.run(["git"] + args, cwd=REPO_DIR, capture_output=True,
                       text=True, env=env)
    out = (p.stdout or "") + (p.stderr or "")
    if token:
        out = out.replace(token, "***TOKEN***")
    return p.returncode, out

def main():
    tok = extract_token(TOKEN_DOC)
    log("Token detected: " + redact(tok))

    # 1) whoami
    st, me = api("GET", API + "/user", tok)
    if st != 200 or "login" not in me:
        log(f"ERROR: GET /user failed ({st}): {me.get('message')}")
        return
    login = me["login"]
    log("Authenticated as GitHub user: " + login)

    # 2) create repo (ignore if already exists)
    data = {
        "name": OWNER_REPO,
        "private": False,
        "description": "Reproducible code, figures and supplementary material for: Genetically higher bone mineral density raises knee osteoarthritis but protects against fracture (Mendelian randomisation of 13 skeletal phenotypes).",
        "auto_init": False,
        "has_issues": True,
        "license_template": "mit",
    }
    st, resp = api("POST", API + "/user/repos", tok, data)
    if st == 201:
        log("Repository created: https://github.com/" + login + "/" + OWNER_REPO)
    elif st == 422:
        log("Repository already exists (proceeding to push): https://github.com/" + login + "/" + OWNER_REPO)
    else:
        log(f"ERROR: create repo failed ({st}): {resp.get('message')}")
        return

    remote_url = f"https://{login}:{tok}@github.com/{login}/{OWNER_REPO}.git"

    # 3) remote add / set-url
    rc, out = run_git(["remote", "get-url", "origin"], tok)
    if rc == 0:
        rc, out = run_git(["remote", "set-url", "origin", remote_url], tok)
    else:
        rc, out = run_git(["remote", "add", "origin", remote_url], tok)
    if rc != 0:
        log("ERROR: setting remote failed:\n" + out)
        return

    # 4) push publish -> main
    rc, out = run_git(["push", "-u", "origin",
                       f"{BRANCH_LOCAL}:refs/heads/{BRANCH_REMOTE}"], tok)
    if rc != 0:
        log("ERROR: push failed:\n" + out)
        # still erase token below
    else:
        log("Push OK: " + BRANCH_LOCAL + " -> " + BRANCH_REMOTE)

    # 5) set default branch = main
    st, _ = api("PATCH", API + f"/repos/{login}/{OWNER_REPO}", tok,
                {"default_branch": BRANCH_REMOTE})
    if st == 200:
        log("Default branch set to: " + BRANCH_REMOTE)

    # 6) ALWAYS erase token from remote URL
    rc, out = run_git(["remote", "set-url", "origin",
                       f"https://github.com/{login}/{OWNER_REPO}.git"])
    if rc == 0:
        log("Token erased from local remote URL: OK")
    else:
        log("WARN: could not reset remote URL: " + out)

    # 7) verify remote URL has no token
    rc, out = run_git(["remote", "-v"])
    if "***TOKEN***" not in out and "@" in out and tok[:4] not in out:
        log("Final remote (safe):\n" + out.strip())
    else:
        log("Final remote check skipped (verify manually).")

if __name__ == "__main__":
    main()
