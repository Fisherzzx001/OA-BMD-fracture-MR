#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""162 · Update the Data-availability statement after republishing + history cleanup.

Two precise, assertion-guarded prose edits in docs/08 line 294:

  1. Zenodo DOI -> the NEW version DOI minted by the v1.1 release
     (10.5281/zenodo.23505928), while also exposing the concept DOI (all versions).
  2. The claim "every output under code/output/ is deposited" stopped being true
     once the regenerable bulky intermediates (PLINK clumping inputs, harmonised
     files) were stripped from git history. Restate it accurately: those are
     regenerated deterministically by the deposited scripts.

Every replacement asserts count==1 first and refuses to write otherwise, so no
silent partial edit can occur. No reported statistic is touched.
"""
import io

PATH = r"C:\Users\Administrator\WorkBuddy\2026-09-08-22-17-40\OA_B_PROJECT\docs\08_投稿稿_OA_BMD_fracture.md"

OLD_DOI = "https://doi.org/10.5281/zenodo.23400861"
NEW_DOI = ("https://doi.org/10.5281/zenodo.23505928 (all versions are listed under the "
           "concept DOI https://doi.org/10.5281/zenodo.23400860)")

OLD_CLAIM = (
    "The deposit includes every script under `code/` and every output under "
    "`code/output/` needed to reproduce the tables and figures (including the MRlap "
    "sample-overlap correction), together with the supplementary files enumerated in "
    "the Supplementary Material section (STROBE-MR checklist, Table S10 MVMR grid, "
    "Table S11 SOS predictions and falsifying tests)."
)
NEW_CLAIM = (
    "The deposit includes every script under `code/` together with the result tables "
    "and source data under `code/output/` required to reproduce every table and figure "
    "(including the MRlap sample-overlap correction); bulk intermediate files that are "
    "regenerated deterministically from the source summary statistics by the deposited "
    "scripts (for example PLINK clumping inputs and harmonised files) are not archived, "
    "together with the supplementary files enumerated in the Supplementary Material "
    "section (STROBE-MR checklist, Table S10 MVMR grid, Table S11 SOS predictions and "
    "falsifying tests)."
)

PAIRS = [("zenodo DOI", OLD_DOI, NEW_DOI),
         ("deposit claim", OLD_CLAIM, NEW_CLAIM)]


def main():
    text = io.open(PATH, encoding="utf-8").read()

    # verify every pair first — nothing is written unless ALL assertions pass
    print("=== assertions ===")
    for label, old, new in PAIRS:
        n = text.count(old)
        print(f"  [{label}] occurrences: {n}")
        if n != 1:
            raise SystemExit(f"ABORT: expected exactly 1 occurrence of {label}, found {n}. File untouched.")
        if old == new:
            raise SystemExit(f"ABORT: {label} unchanged (no-op). File untouched.")

    # also confirm the two reported numbers we must never disturb are still there
    for guard, expect in ((("(" + "β=+0.078"), False), ("OR 1.060", True)):
        pass
    sentinels = ["OR 1.060", "3.3887", "PP.H4=0.799"]
    print("\n=== sentinel statistics still present? ===")
    for s in sentinels:
        print(f"  {s:14} -> {text.count(s)} occurrence(s)")

    print("\n=== applying ===")
    for label, old, new in PAIRS:
        text = text.replace(old, new)
        print("  replaced:", label)

    io.open(PATH, "w", encoding="utf-8", newline="\n").write(text)
    print("\nwritten OK — re-run code/148 and code/149b to confirm compliance + data integrity.")


if __name__ == "__main__":
    main()
