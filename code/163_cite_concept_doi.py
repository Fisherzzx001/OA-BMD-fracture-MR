#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""163 · Cite the Zenodo CONCEPT DOI in Data availability (future-proof).

Why: the GitHub-Zenodo integration mints a new version DOI on every Release, so
hard-coding a version DOI inside the manuscript guarantees the text is one step
behind the archive (and the archived copy would carry a stale citation). The
concept DOI 10.5281/zenodo.23400860 always resolves to the latest version and
lists every version, which is what a reader needs.

Assertion-guarded: refuses to write unless the target string occurs exactly once.
"""
import io

PATH = r"C:\Users\Administrator\WorkBuddy\2026-09-08-22-17-40\OA_B_PROJECT\docs\08_投稿稿_OA_BMD_fracture.md"

OLD = ("https://doi.org/10.5281/zenodo.23505928 (all versions are listed under the "
       "concept DOI https://doi.org/10.5281/zenodo.23400860)")
NEW = ("https://doi.org/10.5281/zenodo.23400860; this concept DOI always resolves to the "
       "most recent version and lists every version of the deposit")


def main():
    text = io.open(PATH, encoding="utf-8").read()
    n = text.count(OLD)
    print("occurrences of target string:", n)
    if n != 1:
        raise SystemExit("ABORT: expected exactly 1 occurrence. File untouched.")

    text = text.replace(OLD, NEW)
    io.open(PATH, "w", encoding="utf-8", newline="\n").write(text)
    print("replaced with concept-DOI wording.")

    # show the resulting sentence fragment for eyeballing
    i = text.find("- **Data availability.**")
    print("\n--- Data availability (now) ---")
    print(text[i:i + 420])
    print("\nRebuild artefacts and re-run code/148 next.")


if __name__ == "__main__":
    main()
