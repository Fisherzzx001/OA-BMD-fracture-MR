#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""333 · Apply the 2026-11-11 full-text compliance edits to docs/08.

Scope (all verified by independent reads before editing):
  A. L49  "N≈1.873 each" -> "N≈31,873 each"  (factual typo; S14 L217, S30c
     L468 and L57 all give 31,873). Pure numeric token, word-count neutral.
  B. L174 Table-2 footnote: remove the revision-history sentence about seven
     weighted-median values "previously tabulated did not ... recomputed ...
     corrected", and the internal script id "(code/129)". This is audit-trail
     material in a journal-facing table footnote.
  C. L174 "no OR is given for the reverse pairs" -> "... for the OA→BMD
     reverse pairs" (L170 OA knee→fracture does carry an OR, so the original
     wording was ambiguous/contradictory).
  D. L560 AI declaration: drop "internal-consistency checking" and the audit
     phrasing "verified against the deposited code and outputs".
  E. L620-627 supplementary index: strip internal script ids "(code/NNN)"
     and the revision marker "this revision" from journal-facing entries.
  F. L631: stop pointing at the stale `Figures/` directory.

Hard constraints honoured:
  * No reported beta / SE / p / OR / CI is changed.
  * Edits B-F are all BELOW "## Figure legends" (L105), which code/148 excludes
    from the 4000-word body; edit A is word-count neutral.
  * docs/08 keeps its inline supplementary tables as the auditable source of
    truth (gates 182/310 depend on them); only the built docx strips them.

Every replacement is asserted to occur exactly once and is printed.
"""
from __future__ import annotations

import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "docs", "08_投稿稿_OA_BMD_fracture.md")

EDITS: list[tuple[str, str, str]] = [
    ("A. L49 sample-size typo",
     "ten regional BMD sites (Qian et al. 2025, N≈1.873 each)",
     "ten regional BMD sites (Qian et al. 2025, N≈31,873 each)"),

    ("B1. L174 revision history removed",
     "Weighted median is Bowden's estimator with weights (β_exposure/SE_outcome)², recomputed from the deposited harmonised instrument tables with 10,000-replicate bootstrap CIs (code/129); the IVW, MR-Egger, Cochran Q and I² columns reproduce from the same deposit to machine precision, and seven weighted-median values previously tabulated did not, so that column was recomputed from the deposited harmonised instrument tables (code/129) and corrected.",
     "Weighted median is Bowden's estimator with weights (β_exposure/SE_outcome)², computed from the harmonised instrument tables with 10,000-replicate bootstrap CIs; the IVW, MR-Egger, Cochran Q and I² columns were obtained from the same instrument tables."),

    ("C. L174 reverse-pair wording",
     "no OR is given for the reverse pairs, whose exposure is a binary liability rather than a 1-SD continuous exposure.",
     "no OR is given for the OA→BMD reverse pairs, whose exposure is a binary liability rather than a 1-SD continuous exposure."),

    ("D. L560 AI declaration",
     "An AI writing assistant was used during preparation for language editing, formatting and internal-consistency checking; all analyses, results, interpretations and conclusions are the authors' own and were verified against the deposited code and outputs.",
     "An AI writing assistant was used during preparation for language editing and formatting; all analyses, results, interpretations and conclusions are the authors' own."),

    ("E1. L620 S31",
     "conditioning on BMD (code/264):", "conditioning on BMD:"),
    ("E2. L621 S25",
     "competitive gene-set test (F3, this revision):",
     "competitive gene-set test (F3):"),
    ("E3. L621 S25 code",
     "total-body BMD enriches skeletal-system-development most strongly (code/171).",
     "total-body BMD enriches skeletal-system-development most strongly."),
    ("E4. L622 S26",
     "tissue-literal (cartilage and skeletal-muscle bone-proxy) arms of F4 (this revision):",
     "tissue-literal (cartilage and skeletal-muscle bone-proxy) arms of F4:"),
    ("E5. L622 S26 code",
     "consistent with a formation-side rather than cartilage-matrix or muscle-mediated interpretation (code/172a, 172b, 173).",
     "consistent with a formation-side rather than cartilage-matrix or muscle-mediated interpretation."),
    ("E6. L625 S29 code",
     "with F4-style Bayesian colocalisation at the 31 instrument loci (code/200, code/201).",
     "with F4-style Bayesian colocalisation at the 31 instrument loci."),
    ("E7. L626 S30 code",
     "(docs/Supplementary/Table_S29_ldsc_site_rg.md; code/207).",
     "(docs/Supplementary/Table_S30_ldsc_site_rg.md)."),
    ("E8. L627 S32 code",
     "(`docs/Supplementary/Table_S32_f10_ancestry_proxy.md`; code/314).",
     "(`docs/Supplementary/Table_S32_f10_ancestry_proxy.md`)."),

    ("F. L631 stale figure path",
     "provided as high-resolution PNG in the submission repository (`Figures/`);",
     "provided as high-resolution PNG in the submission repository;"),
]


def main() -> int:
    text = open(SRC, encoding="utf-8").read()
    orig = text
    leg = text.find("## Figure legends")
    body = text[:leg]
    before_body_words = len(re.findall(r"[A-Za-z][A-Za-z0-9'\-]*", body))

    ok = True
    for label, old, new in EDITS:
        n = text.count(old)
        if n != 1:
            print(f"FAIL  {label}: expected exactly 1 occurrence, found {n}")
            ok = False
            continue
        text = text.replace(old, new, 1)
        print(f"OK    {label}")

    if not ok:
        print("ABORTED: no file written (one or more patterns did not match).")
        return 1

    # Word-count guard: the only edit inside the counted body is a pure numeral.
    leg2 = text.find("## Figure legends")
    body2 = text[:leg2]
    after_body_words = len(re.findall(r"[A-Za-z][A-Za-z0-9'\-]*", body2))
    print(f"\nBody word count: {before_body_words} -> {after_body_words} "
          f"(delta {after_body_words - before_body_words:+d})")
    if before_body_words != after_body_words:
        print("ABORTED: body word count changed unexpectedly.")
        return 1

    # Reported-value guard: no numeric token may have been altered except the L49 typo.
    old_nums = set(re.findall(r"\d+\.\d+|\d+[eE]-\d+|\d+,\d+", orig))
    new_nums = set(re.findall(r"\d+\.\d+|\d+[eE]-\d+|\d+,\d+", text))
    removed = sorted(old_nums - new_nums)
    added = sorted(new_nums - old_nums)
    print(f"Numeric tokens removed: {removed}")
    print(f"Numeric tokens added  : {added}")
    # The corrected value "31,873" already occurs elsewhere in the manuscript
    # (Table S14 L217, Table S30c L468, L57), so it is NOT a newly added token.
    # The only legitimate numeric delta is the removal of the typo "1.873".
    if removed != ["1.873"] or added:
        print("ABORTED: unexpected numeric change.")
        return 1
    if text.count("N≈31,873 each") != 1 or "N≈1.873" in text:
        print("ABORTED: L49 correction not verified in place.")
        return 1

    open(SRC, "w", encoding="utf-8", newline="").write(text)
    print(f"\nWRITTEN {SRC} ({len(text)} chars)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
