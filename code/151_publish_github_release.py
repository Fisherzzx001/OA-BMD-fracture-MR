#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""151 · Re-publish corrected content to GitHub and trigger Zenodo auto-mint.

Steps (user-confirmed 2026-10-20):
 1. Extract GitHub token from github_token.doc (never printed; only redacted form).
 2. Verify auth (GET /user) and repo access.
 3. Force-push local master -> refs/heads/main  (histories are disjoint: no merge base,
    so no fast-forward path exists; the remote "Initial public release" content is
    superseded and master already carries the real URL + Zenodo DOI).
 4. Create a GitHub Release (new tag) -> triggers the GitHub-Zenodo integration.
 5. ALWAYS scrub the token from the local remote URL (try/finally).

Token redaction: only first4...last4 is ever printed.
"""
import os, re, sys, json, subprocess, urllib.request, urllib.error

REPO_DIR = r"C:\Users\Administrator\WorkBuddy\2026-09-08-22-17-40\OA_B_PROJECT"
TOKEN_DOC = r"C:\Users\Administrator\github_token.doc"
CLEAN_REMOTE = "https://github.com/Fisherzzx001/OA-BMD-fracture-MR.git"
API = "https://api.github.com"


def redact(tok):
    if not tok or len(tok) < 8:
        return "***"
    return tok[:4] + "..." + tok[-4:]


def extract_token(path):
    raw = open(path, "rb").read()
    candidates = []
    if raw[:4] == b"PK\x03\x04":
        import zipfile
        z = zipfile.ZipFile(path)
        xml = z.read("word/document.xml").decode("utf-8", "ignore")
        candidates.append(re.sub(r"<[^>]+>", " ", xml))
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
    raise SystemExit("ERROR: no GitHub token found in token file.")


def api(method, url, token, data=None):
    req = urllib.request.Request(url, method=method)
    req.add_header("Accept", "application/vnd.github+json")
    req.add_header("User-Agent", "oa-publish-script")
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
            msg = json.loads(body).get("message", body)
        except Exception:
            msg = body
        return e.code, {"message": msg}


def run_git(args, token=None):
    env = dict(os.environ)
    env["GIT_ASKPASS"] = "true"
    env["GIT_TERMINAL_PROMPT"] = "0"
    p = subprocess.run(["git"] + args, cwd=REPO_DIR, capture_output=True, text=True, env=env)
    out = (p.stdout or "") + (p.stderr or "")
    if token:
        out = out.replace(token, "***TOKEN***")
    return p.returncode, out


def next_tag(tok, login, repo):
    st, tags = api("GET", API + f"/repos/{login}/{repo}/tags?per_page=100", tok)
    have = [t["name"] for t in tags if isinstance(t, dict) and "name" in t] if st == 200 else []
    versions = []
    for nm in have:
        m = re.fullmatch(r"v?(\d+)\.(\d+)", nm or "")
        if m:
            versions.append((int(m.group(1)), int(m.group(2))))
    print("  existing tags:", have if have else "(none)")
    if versions:
        maj, minr = max(versions)
        return f"v{maj}.{minr + 1}"
    return "v1.1"


def main():
    tok = extract_token(TOKEN_DOC)
    print("Token detected:", redact(tok))

    st, me = api("GET", API + "/user", tok)
    if st != 200 or "login" not in me:
        print(f"ERROR: GET /user failed ({st}): {me.get('message')}")
        return
    login = me["login"]
    repo = "OA-BMD-fracture-MR"
    print("Authenticated as:", login)

    st, r = api("GET", API + f"/repos/{login}/{repo}", tok)
    if st != 200:
        print(f"ERROR: cannot access repo ({st}): {r.get('message')}")
        return
    print("Repo OK:", r.get("full_name"), "| default branch:", r.get("default_branch"))

    # inject token into remote, but guarantee scrub afterwards
    try:
        token_remote = f"https://{login}:***TOKEN***@github.com/{login}/{repo}.git"
        real_remote = f"https://{login}:{tok}@github.com/{login}/{repo}.git"
        rc, out = run_git(["remote", "set-url", "origin", real_remote], tok)
        if rc != 0:
            print("ERROR: set-url failed:\n", out)
            return
        print("Remote set to:", token_remote)

        rc, out = run_git(["push", "--force", "origin", "master:refs/heads/main"], tok)
        print("--- push ---")
        print(out.strip()[:2000])
        if rc != 0:
            print("ERROR: force push failed")
            return
        print("Push OK: master -> main (force)")
    finally:
        rc, out = run_git(["remote", "set-url", "origin", CLEAN_REMOTE])
        print("Token scrubbed from remote:", "OK" if rc == 0 else "FAILED")

    # verify remote state
    rc, out = run_git(["ls-remote", "--heads", "origin", "main"])
    print("--- ls-remote origin main ---")
    print(out.strip()[:300])

    tag = next_tag(tok, login, repo)
    body = (
        "Corrected and OAC-compliant version.\n\n"
        "- Fixes three independent root causes: rg-matrix reproducibility (G8), "
        "G6 effect-scale propagation (log-OR factor 3.3887), and the weighted-median "
        "step-near convention.\n"
        "- Corrected figures and Table 2 now replace the earlier OR 1.017-1.028 / wmed values.\n"
        "- Manuscript body compressed to 3973 words (OAC limit 4000); all reported "
        "beta/SE/p/OR/CI values verified unchanged.\n"
    )
    st, rel = api("POST", API + f"/repos/{login}/{repo}/releases", tok, {
        "tag_name": tag,
        "name": f"{tag} — corrected figures, Table 2 and OAC-compliant manuscript",
        "body": body,
        "draft": False,
        "prerelease": False,
    })
    if st in (200, 201):
        print(f"Release created: {tag}")
        print("Release URL:", rel.get("html_url"))
        print("The GitHub-Zenodo integration should mint a new Zenodo version shortly.")
    else:
        print(f"ERROR: release creation failed ({st}): {rel.get('message')}")


if __name__ == "__main__":
    main()
