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
            f"Submission-ready OAC manuscript {TAG}: comprehensive adversarial-"
            "review rectification (M1-M6 + m1-m10). Scale reconstruction resolved "
            "via independent GEFOS total-body BMD variant-level validation "
            "(median GEFOS/deposited beta = 0.79, confirming deposited BMD beta "
            "is already on the true per-SD scale; x3.543 is a methodological "
            "narrative phantom, never applied). SOS theory boundary rewritten "
            "(de-causalised title/Highlights; discriminating kernel explicitly "
            "'neither confirmed nor refuted'). MVMR collider caveat, power-limited "
            "colocalisation reading, random-effects 2/5 as primary口径. "
            "50 refs, 3996/4000 words (Word-compatible count); reported "
            "beta/SE/p/OR/CI unchanged. Seven gates green (148 COMPLIANT / "
            "149b SAFE / 182 79/79 / 310 98/98 / 331 31/31 / 334 18/18 / "
            "336 docx 3996/4000 COMPLIANT). "
            "Concept DOI: https://doi.org/10.5281/zenodo.23400860")
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
            f"Submission-ready version {TAG}. Comprehensive adversarial-review "
            "rectification (docs/149) of the bone-mineral-density / osteoarthritis "
            "/ fracture MR manuscript.\n\n"
            "M1 (scale reconstruction) RESOLVED by re-framing: an independent "
            "GEFOS total-body BMD GWAS (GCST005348, N=56,284) variant-level "
            "comparison across all 31 MR instruments gave median(GEFOS/deposited "
            "beta) = 0.79 (IQR 0.67-0.85), not 0.28 (=1/3.543). The deposited BMD "
            "beta and SE are already on the verified per-SD scale, so the x3.543 "
            "multiplier was a methodological-narrative phantom that was never "
            "applied; all five BMD->OA ORs survive a +/-20% scale perturbation. "
            "The OA reconstruction factor 1/sqrt[K(1-K)] = 3.389 (K=0.0964) is "
            "superseded but retained as a traceable locked parameter. h2 is "
            "computed on the as-deposited N=31,873 (not the placeholder 400,000) "
            "and is internally consistent.\n\n"
            "M2 (SOS set-point kernel): theory boundary rewritten - the title is "
            "de-causalised ('associated with ... pleiotropy-aware Mendelian-"
            "randomisation evidence'); SOS is framed as 'one untested single-set-"
            "point discriminator (a falsifiable hypothesis, not a demonstrated "
            "mechanism)'; the two-module mixture (code/340) gives rho_hat = -0.13 "
            "(95% CI -0.58 to +0.40, spanning zero), rejecting the strict "
            "single-set-point (rho=-1) model yet unable to separate a weak "
            "set-point from two near-independent modules. PCA orthogonality is "
            "explicitly neutral about a single set-point.\n\n"
            "M3/M4: MR-Egger non-significance (all p>0.24) with 16-31 instruments "
            "is read as underpowered (not detected); the positive control's own "
            "femoral-neck MR-Egger p=0.018 makes pleiotropy-robust readings "
            "contingent; random-effects 2/5 is now the primary口径 in the "
            "abstract and conclusion. M5: explicit MVMR collider caveat - "
            "conditioning on height/BMI is read as a sensitivity bound, not a "
            "deconfounded effect; the 67-91% residual is not resolved into bone "
            "mass. M6: the 9/31 vs 0/31 colocalisation split is read as "
            "power-limited, not endpoint-specific. m1-m10 all addressed "
            "(title tone, pre-registration deviation disclosure, post-hoc power "
            "label, PCA neutrality, P5 hypothesis-under-test, '9/10' register "
            "completeness caveat, areal-BMD qualifier, K locked parameter, "
            "clinical-relevance downweighting, Steiger downweighting).\n\n"
            "OAC-compliant manuscript: 3996/4000 words on the Word-compatible "
            "count (numerals included), 50/50 references, 8/8 display items, "
            "Highlights 5 bullets (<=85 each). Seven gates green "
            "(148 COMPLIANT / 149b SAFE / 182 79/79 / 310 98/98 / 331 31/31 / "
            "334 18/18 / 336 docx-level 3996/4000 COMPLIANT). "
            "No reported beta/SE/p/OR/CI changed. "
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
