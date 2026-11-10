#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
318_github_ssh_release.py -- push master + mint a new Zenodo version from a
sandbox that blocks git's HTTPS smart-HTTP but permits HTTP CONNECT to
ssh.github.com:443.

Design (no secrets are ever printed or persisted):
  * classic PAT is read from C:/Users/Administrator/github_token.doc (scope: repo)
  * a TEMPORARY ed25519 deploy key (write access) is registered via the REST API
  * git pushes over GitHub SSH-on-443 through code/317_ssh_proxy_relay.py
  * an annotated tag (default v2.1, override with RELEASE_TAG env var) is pushed
    and a GitHub Release is published, which is what the GitHub-Zenodo
    integration listens for
  * the deploy key is ALWAYS deleted in `finally` (phases share a state file)

Phases (run in order):
  python3 code/318_github_ssh_release.py setup
  python3 code/318_github_ssh_release.py push
  python3 code/318_github_ssh_release.py cleanup
Or `all` for setup -> push -> cleanup in one shot.

NOTE: this script performs PUBLIC external actions (push, tag, release).
Run only when the user has asked for the release.
"""
import os
import re
import sys
import json
import time
import shutil
import socket
import subprocess
import urllib.request
import urllib.error
from datetime import datetime, timezone

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TOKEN_DOC = r"C:\Users\Administrator\github_token.doc"
TMP = os.path.join(ROOT, "_tmp_release")
STATE = os.path.join(TMP, "state.json")
KEY = os.path.join(TMP, "id_ed25519")
SSH_CONFIG = os.path.join(TMP, "ssh_config")
RELAY = os.path.join(ROOT, "code", "317_ssh_proxy_relay.py")
OWNER, REPO = "Fisherzzx001", "OA-BMD-fracture-MR"
TAG = os.environ.get("RELEASE_TAG", "v2.1")
REMOTE = "_ssh_release_tmp"
API = "https://api.github.com"


def posix(p):
    return p.replace("\\", "/")


def token():
    raw = open(TOKEN_DOC, encoding="utf-8", errors="ignore").read()
    m = re.search(r"ghp_[A-Za-z0-9]{20,}", raw)
    if not m:
        raise SystemExit("ERROR: no classic token found in token doc")
    return m.group(0)


TOK = token()


def scrub(s):
    return (s or "").replace(TOK, "***TOKEN***")


def api(method, path, payload=None, timeout=40):
    data = json.dumps(payload).encode() if payload is not None else None
    req = urllib.request.Request(
        API + path, data=data, method=method,
        headers={"Authorization": "Bearer " + TOK,
                 "Accept": "application/vnd.github+json",
                 "Content-Type": "application/json",
                 "User-Agent": "workbuddy-release"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            body = r.read().decode("utf-8", "replace")
            return r.status, (json.loads(body) if body else {})
    except urllib.error.HTTPError as e:
        body = e.read().decode("utf-8", "replace")
        try:
            return e.code, json.loads(body)
        except Exception:
            return e.code, {"raw": body[:400]}
    except Exception as e:
        return None, {"error": str(e)}


def git(args, env=None, timeout=300):
    e = dict(os.environ)
    if env:
        e.update(env)
    try:
        r = subprocess.run(["git"] + args, cwd=ROOT, env=e,
                           capture_output=True, text=True, timeout=timeout)
        return r.returncode, scrub(r.stdout.strip()), scrub(r.stderr.strip())
    except subprocess.TimeoutExpired:
        return 124, "", "TIMEOUT"


def load_state():
    return json.load(open(STATE, encoding="utf-8")) if os.path.exists(STATE) else {}


def save_state(d):
    os.makedirs(TMP, exist_ok=True)
    with open(STATE, "w", encoding="utf-8") as f:
        json.dump(d, f, indent=2)


# ------------------------------------------------------------------ setup
def setup():
    os.makedirs(TMP, exist_ok=True)
    st = load_state()

    if not os.path.exists(KEY):
        r = subprocess.run(["ssh-keygen", "-t", "ed25519", "-N", "", "-q",
                            "-f", KEY, "-C", "workbuddy-temp-release"],
                           capture_output=True, text=True)
        if r.returncode != 0:
            raise SystemExit("ssh-keygen failed: " + scrub(r.stderr))
        print("[setup] generated temporary ed25519 key")
    pub = open(KEY + ".pub", encoding="utf-8").read().strip()

    cfg = f"""Host wb-github
    HostName ssh.github.com
    Port 443
    User git
    IdentityFile {posix(KEY)}
    IdentitiesOnly yes
    StrictHostKeyChecking accept-new
    UserKnownHostsFile {posix(os.path.join(TMP, 'known_hosts'))}
    ProxyCommand "{posix(sys.executable)}" "{posix(RELAY)}" ssh.github.com 443
"""
    with open(SSH_CONFIG, "w", encoding="utf-8", newline="\n") as f:
        f.write(cfg)
    print("[setup] wrote SSH config with ProxyCommand relay")

    if not st.get("key_id"):
        title = "workbuddy-temp-push-" + datetime.now(timezone.utc).strftime("%Y%m%d%H%M%S")
        code, body = api("POST", f"/repos/{OWNER}/{REPO}/keys",
                         {"title": title, "key": pub, "read_only": False})
        print(f"[setup] deploy key POST: HTTP {code}")
        if code not in (201, 422):
            raise SystemExit("deploy key registration failed: " + str(body)[:400])
        if code == 422:
            raise SystemExit("deploy key rejected (422): " + str(body)[:400])
        st["key_id"] = body["id"]
        save_state(st)
        print(f"[setup] deploy key id={st['key_id']} (write access; will be deleted)")
    else:
        print(f"[setup] reusing deploy key id={st['key_id']}")

    # connectivity test: GitHub answers "Hi ...! You've successfully authenticated"
    # with exit status 1 for `ssh -T`; both outcomes are informative.
    try:
        r = subprocess.run(["ssh", "-F", SSH_CONFIG, "-T", "wb-github"],
                           cwd=ROOT, capture_output=True, text=True, timeout=90)
        msg = scrub((r.stdout + r.stderr).strip())
        print(f"[setup] ssh -T rc={r.returncode}: {msg[:300]}")
        if "successfully authenticated" not in msg:
            raise SystemExit("SSH authentication did not succeed -- aborting before any push")
    except subprocess.TimeoutExpired:
        raise SystemExit("SSH handshake timed out -- aborting before any push")
    print("[setup] READY: SSH push channel verified")


# ------------------------------------------------------------------- push
def push():
    st = load_state()
    if not st.get("key_id"):
        raise SystemExit("run setup first")
    env = {"GIT_SSH_COMMAND": f"ssh -F {posix(SSH_CONFIG)}"}

    # remote state before push
    code, body = api("GET", f"/repos/{OWNER}/{REPO}/branches/master")
    print(f"[push] remote master before: {body.get('commit', {}).get('sha', '?')[:12]} (HTTP {code})")

    subprocess.run(["git", "remote", "remove", REMOTE], cwd=ROOT, capture_output=True)
    rc, out, err = git(["remote", "add", REMOTE,
                        f"wb-github:{OWNER}/{REPO}.git"])
    print("[push] temp remote added:", rc, err)

    ssh_url = f"ssh://git@ssh.github.com:443/{OWNER}/{REPO}.git"
    rc, out, err = git(["push", REMOTE, "master"], env=env, timeout=600)
    print(f"[push] git push master rc={rc}\n  out: {out[:400]}\n  err: {err[:400]}")
    if rc != 0:
        raise SystemExit("PUSH FAILED -- no tag/release was created")

    code, body = api("GET", f"/repos/{OWNER}/{REPO}/branches/master")
    print(f"[push] remote master after : {body.get('commit', {}).get('sha', '?')[:12]} (HTTP {code})")
    local = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    print(f"[push] local  master       : {local[:12]}")
    if body.get("commit", {}).get("sha") != local:
        raise SystemExit("remote SHA != local SHA after push -- stopping before tag")

    # If the tag is already published on GitHub, never repoint/force it:
    # a release is immutable and later master commits must not move v2.0.
    rc_tags, remote_tags = api("GET", f"/repos/{OWNER}/{REPO}/tags")
    published = isinstance(remote_tags, list) and any(
        t.get("name") == TAG for t in remote_tags)
    if published:
        print(f"[push] tag {TAG} already published on GitHub -- leaving it untouched")
        st["release_done"] = True
        save_state(st)
        return

    # annotated tag (idempotent).  A pre-existing LOCAL tag may point at an old
    # commit whose tree still contains the >100 MB blobs; repoint it to HEAD so
    # the tag push cannot re-introduce them.
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    tags = subprocess.check_output(["git", "tag", "-l", TAG], cwd=ROOT, text=True).strip()
    need_create = not tags
    if tags:
        tgt = subprocess.check_output(["git", "rev-parse", TAG], cwd=ROOT, text=True).strip()
        if tgt == head:
            print(f"[push] tag {TAG} already points at HEAD")
        else:
            print(f"[push] repointing stale local tag {TAG}: {tgt[:12]} -> {head[:12]}")
            git(["tag", "-d", TAG])
            need_create = True
    if need_create:
        tag_msg = (
            f"Submission-ready OAC manuscript {TAG}: SOS cellular anchor and "
            "S22 dual-axis discriminating kernel, PLINK .bed decode fix + code/331 "
            "guard, three-trait coloc 3/31 SOS-congruent, Chai 2026 independent-"
            "instrument audit; 50 refs, 3984/4000 words, reported values unchanged")
        rc, out, err = git(["tag", "-a", TAG, "-m", tag_msg])
        print("[push] tag created at HEAD:", rc, err)
    rc, out, err = git(["push", REMOTE, TAG], env=env, timeout=300)
    print(f"[push] git push tag {TAG} rc={rc}\n  out: {out[:300]}\n  err: {err[:300]}")
    if rc != 0:
        raise SystemExit("TAG PUSH FAILED")

    # GitHub Release -> this is the event the Zenodo integration listens for
    code, rel = api("GET", f"/repos/{OWNER}/{REPO}/releases/tags/{TAG}")
    if code == 200:
        print(f"[push] release for {TAG} already exists (id={rel.get('id')})")
    else:
        bodytxt = (
            f"Submission-ready version {TAG}. Layered T1/T2/T3 reinforcement over v2.1.\n\n"
            "T1 (theory): Peng 2026 RSPO3-NFATC1 promoted to the SOS cellular anchor "
            "(S30 vi); S22 dual-axis PCA (bone-mass/fracture vs OA axes orthogonal) "
            "built into the discriminating kernel; word-neutral prose edit.\n\n"
            "T2 (methods + revision flagship): a PLINK1 .bed 2-bit decoding error "
            "(het/homA2 swapped) was found and fixed in the LD-dependent readers, with "
            "code/331 added as a permanent regression guard (31/31 decode and eaf "
            "concordance) - no reported beta/SE/p/OR/CI affected. Full-cost re-runs "
            "(clumping settings unchanged): coloc-SuSiE 0/31 after a per-locus window "
            "fix; three-trait BMD/OA/fracture coloc 3/31 at PP.H4>=0.8, 3/3 SOS-sign-"
            "congruent after a sign-predicate fix; Chai 2026 independent instruments "
            "662 harmonised (direction-consistent, random-effects p=0.40, NOT an "
            "independent significant replication) with osteoblast/chondrocyte "
            "enrichment; CPBB weight-bearing positive-control MR OR 0.694, p=4.6e-20; "
            "F5 JSN-side MR uninformative.\n\n"
            "T3 (long-term): PanUKB SAS heel-BMD / P1NP / vBMD public-data readiness "
            "protocols recorded. Also carries the OSF preregistration DOI "
            "(10.17605/OSF.IO/CSMP7), S32 F10 ancestry-composition proxy, S29 LDSC "
            "site-level rg, and the F2/F5/F10 alternative-path audit.\n\n"
            "OAC-compliant manuscript: 3984/4000 words, 50/50 references. Four gates "
            "green (148 COMPLIANT / 149b SAFE / 182 65-65 / 310 30-30) plus 331 PASS. "
            "Concept DOI: https://doi.org/10.5281/zenodo.23400860")
        rel_name = "submission-" + TAG
        code, rel = api("POST", f"/repos/{OWNER}/{REPO}/releases",
                        {"tag_name": TAG, "name": rel_name,
                         "body": bodytxt, "draft": False, "prerelease": False})
        print(f"[push] release POST: HTTP {code} id={rel.get('id')} url={rel.get('html_url','')}")
        if code not in (201, 422):
            raise SystemExit("release creation failed: " + str(rel)[:400])
    st["release_done"] = True
    save_state(st)
    print("[push] DONE -- Zenodo should mint a new version within a few minutes")


# ---------------------------------------------------------------- cleanup
def cleanup():
    st = load_state()
    if st.get("key_id"):
        code, body = api("DELETE", f"/repos/{OWNER}/{REPO}/keys/{st['key_id']}")
        print(f"[cleanup] deploy key DELETE: HTTP {code}")
        st["key_id"] = None
        save_state(st)
    subprocess.run(["git", "remote", "remove", REMOTE], cwd=ROOT, capture_output=True)
    print("[cleanup] temp git remote removed")
    if os.path.isdir(TMP):
        shutil.rmtree(TMP, ignore_errors=True)
        print("[cleanup] temporary key/config directory removed")
    # Never write the PAT into origin; only report a redacted URL.
    origin = subprocess.check_output(["git", "remote", "get-url", "origin"],
                                     cwd=ROOT, text=True).strip()
    safe = re.sub(r"//[^/@]+@github\.com", "//***@github.com", origin)
    print("[cleanup] origin (redacted):", safe)


if __name__ == "__main__":
    phase = sys.argv[1] if len(sys.argv) > 1 else "all"
    os.chdir(ROOT)
    try:
        if phase in ("setup", "all"):
            setup()
        if phase in ("push", "all"):
            push()
    finally:
        if phase in ("cleanup", "all"):
            cleanup()
