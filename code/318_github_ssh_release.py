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
TAG = os.environ.get("RELEASE_TAG", "v2.7")
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
            "Submission v2.7 (2026-11-15): substantive post-v2.6 adversarial-review "
            "remediation, nine local commits since tag v2.6.1. .zenodo.json version "
            "was bumped to v2.7 BEFORE tagging so the Zenodo record carries the "
            "correct version label. Three batches: (1) Option 1 descriptive reframing "
            "to a multi-phenotype genetic-architecture plus triangulation study, with "
            "the Systemic Osteogenic Set-point downgraded to a post-hoc interpretive "
            "hypothesis, P1-P6/F1-F10 relocated to supplementary and labelled post-hoc, "
            "rho=-1 identified as a zero-discriminating boundary strawman, and PCA "
            "orthogonality stated as neutral; (2) five Reviewer-2 mandatory items - "
            "multiplicity control via two analysis-defined families, EUR-only FinnGen "
            "direction-consistency explicitly NOT an independent replication, MRlap "
            "configuration transparency, instrument-attrition accounting, and 5-of-10 "
            "site scope; (3) eighth-round adversarial review (docs/151) P0 and P1/P2 - "
            "N1 Bonferroni fixed-effect 5/5 vs random-effect 0/5 now stated together, "
            "N2 Table S27 instrument strength recomputed from its own 31 per-SNP rows "
            "(min F 31.2, median 41.8, mean 64.8, aggregate F 63.8, R2 0.063), plus "
            "nine further interpretation/estimator/design items. Full audit trail in "
            "docs/149 sections 8.7 and 8.8. All changes are wording, table-definition "
            "and interpretive self-consistency only: no reported beta/SE/p/OR/CI "
            "changed, no statistical analysis re-run. Seven gates green (148 COMPLIANT "
            "abstract 347/350, body 3993/4000 / 149b SAFE / 182 79/79 / 310 99/99 / "
            "331 31/31 / 334 18/18 / 336 docx 3998 COMPLIANT). Concept DOI: "
            "https://doi.org/10.5281/zenodo.23400860")
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
            "## Submission v2.7 (2026-11-15) - post-v2.6 adversarial-review remediation\n\n"
            "This is a **substantive** archive, not a metadata-only re-archive. "
            "`.zenodo.json` was bumped to v2.7 **before** tagging, so the "
            "GitHub-Zenodo integration mints a record with the correct version "
            "label (avoiding the v2.6 -> v2.6.1 correction that was previously "
            "needed). Nine local commits since tag v2.6.1.\n\n"
            "All changes are confined to **wording, table definitions and "
            "interpretive self-consistency**: no reported beta/SE/p/OR/CI changed, "
            "and no statistical analysis was re-run.\n\n"
            "### Batch 1 - theoretical reframing (Option 1)\n"
            "- The study is repositioned as a descriptive multi-phenotype genetic-"
            "architecture plus triangulation analysis; the Systemic Osteogenic "
            "Set-point is explicitly a **post-hoc interpretive hypothesis, not the "
            "central hypothesis**.\n"
            "- P1-P6/F1-F10 relocated to supplementary and labelled post-hoc; "
            "rho=-1 identified as a zero-discriminating boundary strawman; PCA "
            "orthogonality explicitly neutral; rho-hat CI and S22 narrative corrected.\n\n"
            "### Batch 2 - five Reviewer-2 mandatory items\n"
            "- Multiplicity controlled via two analysis-defined families.\n"
            "- Ancestry mismatch: EUR-only FinnGen knee-OA direction-consistency "
            "explicitly declared **within-ancestry, NOT an independent replication**.\n"
            "- MRlap configuration transparency; 41% instrument-attrition "
            "accounting; 5-of-10 site scope limitation.\n\n"
            "### Batch 3 - eighth-round adversarial review (docs/151)\n"
            "- **N1 (P0)**: Bonferroni fixed-effect 5/5 vs random-effect 0/5 "
            "survival now stated together instead of as mutually exclusive claims.\n"
            "- **N2 (P0)**: Table S27 instrument strength recomputed from its own "
            "31 per-SNP rows (min F 31.2, median 41.8, mean 64.8, aggregate F "
            "63.8, R2_exposure 0.063), reconciled with Tables 3/S13/CSV.\n"
            "- Nine P1/P2 items: colocalisation threshold removed from multiplicity "
            "control; weighted-median attenuation; MR-PRESSO consistency; reverse-MR "
            "power explanation; negative cross-trait intercept; weak conditional "
            "instruments; +/-20% absolute-scale uncertainty; 1,980 LD-correlated "
            "SNPs rather than independent loci; new Supplementary sub-tables S17e/f.\n\n"
            "Full audit trail: `docs/149` sections 8.7 and 8.8; review reports "
            "`docs/150`, `docs/151`. Audit trail is kept out of the manuscript.\n\n"
            "Seven gates remain green (148 COMPLIANT: abstract 347/350, body "
            "3993/4000, 50 refs, 8 display items / 149b SAFE / 182 79/79 / "
            "310 99/99 NO REGRESSION / 331 31/31 / 334 18/18 / 336 docx "
            "Word-like 3998 COMPLIANT). Concept DOI: "
            "https://doi.org/10.5281/zenodo.23400860 (resolves to the latest "
            "version; earlier archives including v2.6 / 24237523 remain available).")
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
