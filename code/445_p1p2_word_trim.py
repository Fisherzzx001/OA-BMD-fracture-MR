# -*- coding: utf-8 -*-
"""Word-budget trim after the P1/P2 remediation (code/444).

Gate 148 counts the markdown body (4000/4000) while gate 336 counts the docx
Word-like body (4005/4000) — the docx carries 5 extra numeric/citation tokens.
This script makes three meaning-preserving trims to reclaim ~10 words across
all three carriers, restoring a safety margin without dropping any scientific
content, reported value, or gate-locked string. Idempotent.
"""
from pathlib import Path
from docx import Document

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "docs/08_投稿稿_OA_BMD_fracture.md"
DOCX = [ROOT / "submission_package/Manuscript.docx",
        ROOT / "submission_package/Submission_MainText.docx"]

TRIMS = [
    # N11: "not independent loci" -> "not independent" (loci already implied by SNPs)
    ("trim-N11",
     "SNPs (LD-correlated, not independent loci), 94.4%",
     "SNPs (LD-correlated, not independent), 94.4%"),
    # N3-L99: compress exploratory-thresholding clause, keep every substantive point
    ("trim-N3",
     "Colocalisation (31 tests) applied exploratory Bayesian posterior thresholding (PP.H4≥0.8), not frequentist error control, and MVMR grids and the ten PRESSO distortion tests are sensitivity diagnostics; the colocalisation 9/31 and random-effects 2/5 findings therefore remain hypothesis-generating.",
     "Colocalisation (31 tests) is exploratory Bayesian thresholding (PP.H4≥0.8), not frequency control; MVMR grids and the ten PRESSO distortion tests are sensitivity diagnostics, so colocalisation 9/31 and random-effects 2/5 remain hypothesis-generating."),
    # N5: compress PRESSO clause
    ("trim-N5",
     "Outlier removal left signs unchanged (Table S17), not robustness evidence given the distortion below.",
     "Outlier removal left signs unchanged (Table S17), not robustness given the distortion below."),
    # N10: drop redundant tilde, shorten "absolute-scale"
    ("trim-N10",
     "95% CI 1.025–1.135 across sites, ±~20% absolute-scale)",
     "95% CI 1.025–1.135 across sites, ±20% scale)"),
]


def edit_text(path, trims):
    t = path.read_text(encoding="utf-8")
    out = []
    for name, old, new in trims:
        if new in t:
            out.append(f"  {name}: SKIP")
        elif old in t:
            t = t.replace(old, new, 1)
            out.append(f"  {name}: APPLIED")
        else:
            out.append(f"  {name}: MISSING")
    path.write_text(t, encoding="utf-8")
    return out


def edit_docx(path, trims):
    d = Document(path)
    out = []
    changed = False
    for name, old, new in trims:
        hits = 0
        for par in d.paragraphs:
            runs = par.runs
            full = "".join(r.text for r in runs)
            if old in full:
                merged = full.replace(old, new, 1)
                if runs:
                    runs[0].text = merged
                    for r in runs[1:]:
                        r.text = ""
                hits += 1
                changed = True
        out.append(f"  {name}: {hits} para")
    if changed:
        d.save(path)
    return out


def main():
    print("== docs/08 ==")
    for x in edit_text(SRC, TRIMS):
        print(x)
    for p in DOCX:
        print(f"== {p.name} ==")
        for x in edit_docx(p, TRIMS):
            print(x)


if __name__ == "__main__":
    main()
