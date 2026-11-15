#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
447_n3n7_docx_patch.py -- finish the two docx-only replacements that code/446
reported as MISSING. The docx builder strips markdown bold markers (the words
are real bold runs), so the anchors must be written without '**'.

The other two MISSING items in code/446 (S14 note, S17 index line) are expected:
the docx builder strips the supplementary block from "## Supplementary Table
S8." onward, and code/310 scopes those assertions to md/mirror only.

Idempotent. Usage: python3 code/447_n3n7_docx_patch.py
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOCXS = [
    os.path.join(ROOT, "submission_package", "Manuscript.docx"),
    os.path.join(ROOT, "submission_package", "Submission_MainText.docx"),
]

TRIPLES = [
 ("N3-L59-posterior-rule (docx)",
  "Colocalisation. Bayesian colocalisation (coloc.abf; Giambartolomei et al. 2014 [22]) tested the three a priori loci (WNT16, SOX6, GDF5) and all 31 BMD instrument loci (±500 kb; PP.H4≥0.8), with a BMD–fracture control.",
  "Colocalisation. Bayesian coloc.abf [22] tested WNT16, SOX6, GDF5 and 31 BMD instrument loci (±500 kb), with a BMD–fracture control; PP.H4≥0.8 is a posterior rule, not error-rate control."),
 ("N4/N5-S11-caption (docx)",
  "S11. Multi-estimator robustness grid for the ten MR pairs (IVW fixed/random, MR-Egger, weighted median, PRESSO outlier-removed, MRlap-corrected, MVMR).",
  "S11. Multi-estimator comparison grid for the ten MR pairs (IVW fixed/random, MR-Egger, weighted median, PRESSO outlier-removed, MRlap-corrected, MVMR); sign concordance is not magnitude concordance (Table S17-e/g)."),
]


def main():
    from docx import Document
    os.chdir(ROOT)
    for path in DOCXS:
        d = Document(path)
        changed = False
        print("[" + os.path.basename(path) + "]")
        for name, old, new in TRIPLES:
            hits = 0
            for p in d.paragraphs:
                full = "".join(r.text for r in p.runs)
                if old in full and new not in full:
                    merged = full.replace(old, new, 1)
                    if p.runs:
                        p.runs[0].text = merged
                        for r in p.runs[1:]:
                            r.text = ""
                    else:
                        p.add_run(merged)
                    hits += 1
                    changed = True
            if hits:
                print(f"  {name}: APPLIED ({hits} para)")
            elif new in "".join(p.text for p in d.paragraphs):
                print(f"  {name}: SKIP (already applied)")
            else:
                print(f"  {name}: MISSING anchor")
        if changed:
            d.save(path)
    # verify both carriers carry the new phrases
    print("[verify]")
    for path in DOCXS:
        txt = "".join(p.text for p in Document(path).paragraphs)
        print(f"  {os.path.basename(path)}: posterior-rule={('posterior rule, not error-rate control' in txt)}, "
              f"comparison-grid={('sign concordance is not magnitude concordance' in txt)}")


if __name__ == "__main__":
    main()
