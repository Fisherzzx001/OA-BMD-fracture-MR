#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""160 · Verify the published GitHub content matches what we intend to ship.

Reads the remote tree straight from the API (not from anywhere local) and checks:
  * key deliverables present (docs/08, figures, Manuscript_Ready, table-source CSVs)
  * none of the stripped bulky artefacts leaked back in
  * total remote tree size, and the largest remaining file
"""
import json
import os
import re
import subprocess
import urllib.request

API = "https://api.github.com"
OWNER, REPO = "Fisherzzx001", "OA-BMD-fracture-MR"
TOKEN_DOC = r"C:\Users\Administrator\github_token.doc"


def tok():
    raw = open(TOKEN_DOC, "rb").read().decode("utf-8", "ignore")
    m = re.search(r"(github_pat_[A-Za-z0-9_]{20,}|gh[pousr]_[A-Za-z0-9]{20,})", raw)
    return m.group(1)


def get(path, token):
    req = urllib.request.Request(API + path)
    req.add_header("Accept", "application/vnd.github+json")
    req.add_header("User-Agent", "verify-remote")
    req.add_header("Authorization", "Bearer " + token)
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode())


BANNED = (".tif", ".zip", "_clump_input.txt", "build/")


def main():
    t = tok()
    tree = get(f"/repos/{OWNER}/{REPO}/git/trees/main?recursive=1", t)
    if tree.get("truncated"):
        print("WARNING: tree listing truncated by API")
    blobs = [e for e in tree["tree"] if e["type"] == "blob"]
    total = sum(e.get("size", 0) for e in blobs)
    print(f"remote 'main': {len(blobs)} files, {total/1048576:.1f} MB total")

    names = [e["path"] for e in blobs]
    print("\n--- key deliverables present? ---")
    must = ["docs/08_投稿稿_OA_BMD_fracture.md",
            "docs/123_根茎数据取证.md",
            "Manuscript_Ready/01_Manuscript.md",
            "submission_package/Cover_Letter.docx",
            "code/148_oac_compliance.py"]
    must = [m for m in must if True]
    for i, m in enumerate(must):
        if i == 1:  # docs/123 may have a different name; report rather than assert
            hit = [n for n in names if n.startswith("docs/123")]
            print(f"  docs/123*: {'YES (' + hit[0] + ')' if hit else 'NOT FOUND'}")
            continue
        print(f"  {m:52} {'YES' if m in names else 'NOT FOUND'}")

    counts = {
        "docs/": sum(1 for n in names if n.startswith("docs/")),
        "code/*.py": sum(1 for n in names if n.startswith("code/") and n.endswith(".py")),
        "Figures_v17_Nature/": sum(1 for n in names if n.startswith("Figures_v17_Nature/")),
        "Manuscript_Ready/": sum(1 for n in names if n.startswith("Manuscript_Ready/")),
        "table CSVs": sum(1 for n in names if n.startswith("code/output/genetics/") and n.endswith(".csv")),
    }
    print("\n--- counts ---")
    for k, v in counts.items():
        print(f"  {k:22} {v}")

    leaked = sorted([(e["path"], e.get("size", 0)) for e in blobs
                     if any(b in e["path"] for b in BANNED)])
    print("\n--- stripped artefacts leaked back? ---")
    if leaked:
        for p, s in leaked[:10]:
            print(f"  LEAKED: {p} ({s/1048576:.1f} MB)")
    else:
        print("  none — clean")

    top = sorted(blobs, key=lambda e: -e.get("size", 0))[:5]
    print("\n--- largest files on remote ---")
    for e in top:
        print(f"  {e.get('size',0)/1048576:6.2f} MB  {e['path']}")


if __name__ == "__main__":
    main()
