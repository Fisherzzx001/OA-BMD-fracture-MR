# -*- coding: utf-8 -*-
"""
441 · R-a + R-b refinement — close the last two residual gaps from the Reviewer-1 re-review.

R-a (Loophole 2 strawman): in L93, after the rho-hat sentence, explicitly flag
    rho=-1 as a boundary strawman and state the two-module mixture analysis is
    uninformative about set-point versus modules.
R-b (recentring):
    - Abstract Objective (L23): lead with the descriptive / MR aim; the set-point
      becomes a secondary, explicitly post-hoc, unconfirmed hypothesis.
    - L93 opening: "the set-point we instrument" -> "the BMD-raising profile we
      instrument".

Edits the single source of truth docs/08, rebuilds the Manuscript_Ready mirror via
code/135, and syncs both delivered docx via run-level merge. Idempotent; preserves
every gate-locked string and every reported beta / SE / p / OR / CI. Does not touch
any reported value.
"""
import os
import subprocess

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOC = os.path.join(ROOT, "docs", "08_投稿稿_OA_BMD_fracture.md")
DOCX = [
    os.path.join(ROOT, "submission_package", "Manuscript.docx"),
    os.path.join(ROOT, "submission_package", "Submission_MainText.docx"),
]

# --- R-a: flag rho=-1 as a boundary strawman (append after the rho-hat sentence) ---
R_A_OLD = "so it indicates uncertainty, not refutation."
R_A_NEW = (
    "so it indicates uncertainty, not refutation. \u03c1=\u22121 is a boundary strawman; "
    "its rejection carries no discriminating weight, so the mixture analysis locates "
    "no coupling strength and remains uninformative about set-point versus modules."
)
R_A_DONE = "boundary strawman"

# --- R-b L23: common anchor (source 'then tested' / docx 'and tested' both contain it) ---
R_B_L23_OLD = (
    "tested whether one BMD-raising set-point (a hypothesis, not a demonstrated "
    "mechanism) generates opposing fracture and joint consequences."
)
R_B_L23_NEW = (
    "examined whether a single BMD-raising set-point could explain the opposing "
    "endpoints \u2014 a post-hoc hypothesis, not a demonstrated mechanism."
)
R_B_L23_DONE = "examined whether a single BMD-raising set-point"

# --- R-b L93 opening: drop the set-point lens ---
R_B_L93_OLD = "the set-point we instrument acts upstream of subchondral remodelling."
R_B_L93_NEW = "the BMD-raising profile we instrument acts upstream of subchondral remodelling."
R_B_L93_DONE = "the BMD-raising profile we instrument acts upstream"

REPLS = [
    (R_A_OLD, R_A_NEW, R_A_DONE, "R-a rho strawman"),
    (R_B_L23_OLD, R_B_L23_NEW, R_B_L23_DONE, "R-b L23 Abstract"),
    (R_B_L93_OLD, R_B_L93_NEW, R_B_L93_DONE, "R-b L93 opening"),
]


def edit_source():
    text = open(DOC, encoding="utf-8").read()
    applied = []
    for old, new, done, tag in REPLS:
        if done in text:
            print(f"[skip] {tag} (already applied)")
            continue
        if old in text:
            text = text.replace(old, new, 1)
            applied.append(tag)
        else:
            print(f"[skip] {tag} (anchor absent in source)")
    open(DOC, "w", encoding="utf-8").write(text)
    print("docs/08 edits:", applied)


def merge_replace(runs, old, new):
    full = "".join(r.text for r in runs)
    if old not in full:
        return 0
    newfull = full.replace(old, new, 1)
    if not runs:
        return 0
    runs[0].text = newfull
    for r in runs[1:]:
        r.text = ""
    return 1


def edit_docx(path):
    from docx import Document
    doc = Document(path)
    total = 0
    for old, new, done, tag in REPLS:
        n = 0
        for par in doc.paragraphs:
            pt = "".join(r.text for r in par.runs)
            if done in pt:
                break
            if old in pt:
                n += merge_replace(par.runs, old, new)
        if n:
            total += n
        else:
            print(f"  [skip docx {os.path.basename(path)}] {tag} (anchor absent)")
    doc.save(path)
    print(f"  {os.path.basename(path)}: {total} docx edits")


if __name__ == "__main__":
    edit_source()
    # Rebuild the Manuscript_Ready mirror from the edited source.
    r = subprocess.run(["python", "code/135_build_manuscript_ready.py"], cwd=ROOT)
    print("135 exit:", r.returncode)
    for p in DOCX:
        if os.path.isfile(p):
            edit_docx(p)
        else:
            print("MISSING", p)
    print("\nDONE — run code/310 / 331 / 334 / 336 to verify zero regression.")
