#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""161 · Check Zenodo for the new version minted by the GitHub integration.

The GitHub-Zenodo webhook usually creates a new VERSION under the same concept
DOI shortly after a Release. This reports every known version with its DOI and
publication date, so we can tell whether v1.1 (2026-10-21) has landed yet and
which DOI to write back into docs/08.
"""
import json
import sys
import time
import urllib.error
import urllib.request

RECORD_ID = "23400861"          # currently cited version DOI suffix
BASE = "https://zenodo.org/api/records"


def fetch(url, tries=6, wait=3.0):
    last = None
    for _ in range(tries):
        try:
            req = urllib.request.Request(url)
            req.add_header("User-Agent", "zenodo-check")
            req.add_header("Accept", "application/json")
            with urllib.request.urlopen(req, timeout=45) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            last = e
            time.sleep(wait)
    print("  fetch failed:", last)
    return None


def show(rec):
    if not rec:
        return
    md = rec.get("metadata", {})
    print(f"  id={rec.get('id')}  doi={rec.get('doi')}  version={rec.get('revision')} "
          f"  created={md.get('publication_date', '?')}")
    print(f"     title: {md.get('title', '?')}")


def main():
    print("=== currently cited record ===")
    cur = fetch(f"{BASE}/{RECORD_ID}")
    show(cur)

    concept = None
    if cur:
        concept = cur.get("conceptrecid")
        print("  conceptrecid:", concept, "| conceptdoi:", cur.get("conceptdoi"))

    if concept:
        print("\n=== all versions of this concept ===")
        vers = fetch(f"{BASE}/{concept}/versions?size=20&sort=version")
        hits = (vers or {}).get("hits", {}).get("hits", [])
        if not hits:
            print("  (none listed)")
        for h in hits:
            md = h.get("metadata", {})
            print(f"  v{h.get('revision')}: doi={h.get('doi')}  created={md.get('publication_date','?')}"
                  f"  {'  <-- NEW' if '2026-10-2' in str(md.get('publication_date','')) else ''}")

    print("\n=== note ===")
    print("  Zenodo can take minutes to hours after a Release to register the new version.")
    print("  Re-run this script until a version dated 2026-10-2x appears, then report its DOI.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
