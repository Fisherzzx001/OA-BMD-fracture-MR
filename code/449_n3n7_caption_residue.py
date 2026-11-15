#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
449_n3n7_caption_residue.py -- third-pass cleanup of N4/N5 wording residue.

code/448 (cross-carrier verification) found that the N4/N5 demotion from
"robustness grid" to "comparison grid" had reached the Figure-legend body but
NOT three derivative locations:

  1. docs/08 Supplementary Material index entries for Figures S11 and S20;
  2. the Figure S20 legend *heading* ("Multi-estimator robustness of ...");
  3. the STROBE-MR checklist item 13 (both copies), which names Figure S11.

All of these are outside the counted body (after `## Figure legends`, or in
supplementary files), so the edit is word-count neutral.  Idempotent.

Hard constraints: no reported beta/SE/p/OR/CI changed; no statistics re-run.

Usage: python3 code/449_n3n7_caption_residue.py
"""
import os
import shutil

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "docs", "08_投稿稿_OA_BMD_fracture.md")
MIRROR = os.path.join(ROOT, "Manuscript_Ready", "01_Manuscript.md")
DOCXS = [
    os.path.join(ROOT, "submission_package", "Manuscript.docx"),
    os.path.join(ROOT, "submission_package", "Submission_MainText.docx"),
]
STROBE = [
    os.path.join(ROOT, "docs", "Supplementary", "STROBE-MR_checklist.md"),
    os.path.join(ROOT, "Manuscript_Ready", "Supplementary", "STROBE-MR_checklist.md"),
]

# (name, old, new)
MD = [
    ("S20-legend-heading",
     "**Supplementary Figure S20.** Multi-estimator robustness of the ten MR pairs.",
     "**Supplementary Figure S20.** Multi-estimator comparison of the ten MR pairs (direction agreement only)."),
    ("S11-index-entry",
     "- **Figure S11.** Multi-estimator robustness grid for the ten MR pairs.",
     "- **Figure S11.** Multi-estimator comparison grid for the ten MR pairs (direction agreement, not robustness evidence)."),
    ("S20-index-entry",
     "- **Figure S20.** Multi-estimator robustness grid of the ten MR pairs (IVW, MR-Egger, weighted median, leave-one-out, MR-PRESSO).",
     "- **Figure S20.** Multi-estimator comparison grid of the ten MR pairs (IVW, MR-Egger, weighted median, leave-one-out, MR-PRESSO); direction agreement only."),
]

STROBE_TRIPLES = [
    ("STROBE-item-13",
     "multi-estimator robustness grid (Figure S11)",
     "multi-estimator comparison grid (Figure S11; direction agreement only)"),
]

# docx has no markdown bold markers and strips the supplementary index block
DOCX_TRIPLES = [
    ("S20-legend-heading (docx)",
     "Multi-estimator robustness of the ten MR pairs.",
     "Multi-estimator comparison of the ten MR pairs (direction agreement only)."),
]


def edit_text(path, triples, label):
    txt = open(path, encoding="utf-8").read()
    orig = txt
    out = []
    for name, old, new in triples:
        if new in txt and old not in txt:
            out.append(f"  {name}: SKIP (already applied)")
            continue
        n = txt.count(old)
        if n == 0:
            out.append(f"  {name}: MISSING anchor")
            continue
        txt = txt.replace(old, new, 1)
        out.append(f"  {name}: APPLIED ({n} match)")
    if txt != orig:
        open(path, "w", encoding="utf-8", newline="\n").write(txt)
    print(f"[{label}] {os.path.relpath(path, ROOT)}")
    print("\n".join(out))


def edit_docx(path, triples):
    from docx import Document
    d = Document(path)
    out = []
    changed = False
    for name, old, new in triples:
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
            out.append(f"  {name}: APPLIED ({hits} para)")
        elif new in "".join(p.text for p in d.paragraphs):
            out.append(f"  {name}: SKIP (already applied)")
        else:
            out.append(f"  {name}: MISSING anchor")
    if changed:
        d.save(path)
    print(f"[docx] {os.path.relpath(path, ROOT)}")
    print("\n".join(out))


def main():
    os.chdir(ROOT)
    edit_text(SRC, MD, "docs/08")
    shutil.copyfile(SRC, MIRROR)
    print(f"[mirror] byte copy -> {os.path.relpath(MIRROR, ROOT)}")
    for p in DOCXS:
        edit_docx(p, DOCX_TRIPLES)
    for p in STROBE:
        edit_text(p, STROBE_TRIPLES, "strobe")

    md = open(SRC, encoding="utf-8").read()
    st = [open(p, encoding="utf-8").read() for p in STROBE]
    print("[check] mirror == source:", open(SRC, "rb").read() == open(MIRROR, "rb").read())
    print("[check] STROBE dual copies identical:", st[0] == st[1])
    print("[check] no 'robustness grid' left in docs/08:", "robustness grid" not in md)
    print("[check] no 'Multi-estimator robustness of' in docs/08:",
          "Multi-estimator robustness of" not in md)
    print("[check] no 'robustness grid' in STROBE:", all("robustness grid" not in s for s in st))
    print("DONE")


if __name__ == "__main__":
    main()
