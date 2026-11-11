#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""332 · Export inline supplementary tables from docs/08 as standalone files.

Background (2026-11-11 full-text compliance audit)
--------------------------------------------------
docs/08_投稿稿_OA_BMD_fracture.md is the repository's single auditable source
of truth and carries twelve full supplementary tables INLINE.  The two
journal-facing docx files were built from the whole file, so those tables were
embedded in the submitted main text -- a packaging defect.  build_submission.py
/ build_maintext.py now strip that block before writing the docx.

If the block is stripped, every table must still be deliverable as a separate
supplementary file.  This script extracts EVERY inline table that needs a
standalone file, using the SAME numbering as the inline manuscript.  Two
historical defects are corrected here (both surfaced by the 2026-11-11 audit):

* Numbering drift: legacy files named the hip-replication S26, the geometry
  arm S27 and the LDSC site-rg S29, while docs/08 had already renumbered these
  to S27 / S28 / S30 (after inserting S25 tissue-enrichment and S26 eQTL).
  One table therefore carried two numbers across the delivery set.
* Stale/incorrect copies: standalone files generated before an earlier
  placeholder fix still contained "Table XX" cross-references, and the S32
  standalone file wrote a Cochran-Q p-value as "p=0.000" (statistically
  inadmissible; the source JSON gives Q_p=4.4e-19, correctly rendered "p<0.001").

Design rules
------------
* docs/08 is the ONLY source.  Each exported body is the verbatim inline slice
  (no statistic, heading or table cell is rewritten).
* The exported filename number == the inline heading number (asserted).
* No exported body may contain "Table XX" or "p=0.000" (asserted).
* Idempotent: identical output on re-run.
* The three legacy mis-numbered TRACKED files are removed only AFTER their
  correctly-numbered replacements have been written and verified.
* S32 keeps its richer hand-maintained standalone file; only its incorrect
  "p=0.000" cells are patched there (documented in the audit ledger).
"""
from __future__ import annotations

import hashlib
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "docs", "08_投稿稿_OA_BMD_fracture.md")
OUTDIR = os.path.join(ROOT, "docs", "Supplementary")

# Inline table number -> standalone filename.  These are the tables that either
# had no standalone file, a stale copy, or a mis-numbered legacy file.
TARGETS = {
    "S25": "Table_S25_tissue_enrichment.md",
    "S26": "Table_S26_eqtl_coloc.md",
    "S27": "Table_S27_hip_oa_replication.md",
    "S28": "Table_S28_geometry_arm.md",
    "S29": "Table_S29_osteophyte.md",
    "S30": "Table_S30_ldsc_site_rg.md",
    "S31": "Table_S31_mvmr_oa_fracture.md",
}

# Legacy mis-numbered files, removed once replacements are verified.  Content
# coverage was independently confirmed (digit-sequence identity for the hip and
# geometry tables; the LDSC legacy file's extra "Reproduction" section is
# re-attached to the S30 export below).
LEGACY = [
    "Table_S26_f1_hip_replication.md",   # content is inline S27 (hip)
    "Table_S27_f2_geometry_arm.md",      # content is inline S28 (geometry)
    "Table_S29_ldsc_site_rg.md",         # content is inline S30 (ldsc)
]

# Reproducibility block that existed only in the legacy LDSC file; re-attached
# to the S30 export so no provenance is lost.  It sits in the supplementary-file
# layer (not the manuscript body), consistent with the other standalone files.
S30_REPRO = (
    "\n## Reproduction\n\n"
    "`code/204_prep_oa_gwas_for_munge.py` → `code/205_run_ldsc_munge_rg.py "
    "--weights standard` (exploratory baseline, placeholder BMD N) → "
    "`code/206_summarise_ldsc_rg.py` → `code/207_bmd_verified_ldsc.py` "
    "(verified N + `--merge-alleles`, this table). Machine-readable values: "
    "`code/output/genetics/ldsc_rg_verified/verified_results.csv`; full "
    "comparison and boundary notes: `COMPARISON.md` in the same directory.\n"
)

HEADER = (
    "<!-- Standalone supplementary table extracted verbatim from\n"
    "     docs/08_投稿稿_OA_BMD_fracture.md by code/332_extract_inline_supplementary.py.\n"
    "     Source slice sha256: {sha} . No statistic was rewritten. -->\n\n"
)


def slices_from_source(text: str) -> dict[str, tuple[str, str]]:
    lines = text.splitlines(keepends=True)
    heads: list[tuple[int, str, str]] = []
    for i, ln in enumerate(lines):
        m = re.match(r"^## Supplementary Table (S\d+)[a-c]?\.\s", ln)
        if m:
            heads.append((i, m.group(1), ln.rstrip("\n")))
        elif ln.startswith("## References"):
            heads.append((i, "REFERENCES", ln.rstrip("\n")))
    out: dict[str, tuple[str, str]] = {}
    for j, (i, num, heading) in enumerate(heads):
        if num == "REFERENCES":
            continue
        end = heads[j + 1][0] if j + 1 < len(heads) else len(lines)
        body = "".join(lines[i:end]).rstrip() + "\n"
        out.setdefault(num, (heading, body))
    return out


def main() -> int:
    text = open(SRC, encoding="utf-8").read()
    slices = slices_from_source(text)

    print("=== inline supplementary sections detected ===")
    for num, (heading, body) in slices.items():
        print(f"  {num:4s} | {len(body.splitlines()):3d} lines | {heading[:74]}")

    missing = [n for n in TARGETS if n not in slices]
    if missing:
        print(f"ERROR: inline sections not found: {missing}", file=sys.stderr)
        return 1

    os.makedirs(OUTDIR, exist_ok=True)
    written: dict[str, str] = {}
    print("\n=== exporting correctly-numbered standalone files ===")
    for num, fname in TARGETS.items():
        heading, body = slices[num]
        if num == "S30":
            body = body.rstrip() + "\n" + S30_REPRO
        # Guard: never propagate a placeholder or an inadmissible p-value.
        for bad in ("Table XX", "p=0.000", "p = 0.000"):
            if bad in body:
                print(f"ERROR: {num} body contains forbidden token {bad!r}", file=sys.stderr)
                return 1
        sha = hashlib.sha256(body.encode("utf-8")).hexdigest()
        out = HEADER.format(sha=sha) + body
        path = os.path.join(OUTDIR, fname)
        prev = open(path, encoding="utf-8").read() if os.path.exists(path) else None
        if prev != out:
            open(path, "w", encoding="utf-8", newline="").write(out)
            state = "WRITTEN"
        else:
            state = "UNCHANGED (idempotent)"
        written[num] = fname
        print(f"  {fname:38s} {state:24s} sha256={sha[:16]}")

    # Remove legacy mis-numbered files only after replacements exist & verified.
    print("\n=== removing legacy mis-numbered files (replacements verified) ===")
    for fname in LEGACY:
        path = os.path.join(OUTDIR, fname)
        if os.path.exists(path):
            os.remove(path)
            print(f"  REMOVED {fname}")
        else:
            print(f"  (already absent) {fname}")

    # Verification: filename number == inline heading number; no forbidden tokens.
    print("\n=== verification ===")
    all_needed = ["S8", "S9", "S12", "S13", "S14", "S15", "S16", "S17",
                  "S25", "S26", "S27", "S28", "S29", "S30", "S31", "S32"]
    files = os.listdir(OUTDIR)
    ok = True
    # Count only Markdown deliverables; .csv files are machine-readable
    # companions to their .md table, not a second table with the same number.
    def md_files(n):
        return [f for f in files if f.startswith(f"Table_{n}_") and f.endswith(".md")]
    for n in all_needed:
        hit = md_files(n)
        status = ", ".join(hit) if hit else "*** MISSING ***"
        if len(hit) != 1:
            ok = False
        print(f"  {n:4s} -> {status}")
    # Duplicate-number check: no inline table may have two standalone .md files.
    for n in TARGETS:
        hit = md_files(n)
        if len(hit) != 1:
            print(f"ERROR: table {n} has {len(hit)} .md files: {hit}", file=sys.stderr)
            ok = False
    print("\nRESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
