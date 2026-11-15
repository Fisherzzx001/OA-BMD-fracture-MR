# -*- coding: utf-8 -*-
"""
440 · Fix a pre-existing docx defect: the two delivered docx files carry a
"Data availability" heading but its body text (which contains the Zenodo
concept DOI and the OSF DOI) was never synced in, and a stray duplicate of the
post-hoc framework paragraph sits at the tail. This breaks 182's
"概念 DOI 存在" and "OSF 状态自洽" checks even though the markdown source
(docs/08 L597) is correct.

Fix (idempotent):
  - locate the "Data availability" heading paragraph;
  - delete any *later* duplicate "Post-hoc interpretive framework" paragraph;
  - insert the full Data availability statement (from docs/08 L597) directly
    after the heading.
"""
import os
import re

from docx import Document

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOC = os.path.join(ROOT, "docs", "08_投稿稿_OA_BMD_fracture.md")
DOCX = [
    os.path.join(ROOT, "submission_package", "Manuscript.docx"),
    os.path.join(ROOT, "submission_package", "Submission_MainText.docx"),
]

CONCEPT = "10.5281/zenodo.23400860"
OSF = "10.17605/OSF.IO/CSMP7"


def data_availability_text():
    for line in open(DOC, encoding="utf-8"):
        if line.strip().startswith("- **Data availability.**"):
            # strip the leading "- **Data availability.** " marker
            return line.strip()[len("- **Data availability.** "):].strip()
    raise RuntimeError("Data availability line not found in docs/08")


TEXT = data_availability_text()


def fix(path):
    doc = Document(path)
    paras = doc.paragraphs
    # find heading index
    h_idx = None
    for i, p in enumerate(paras):
        if p.text.strip() == "Data availability":
            h_idx = i
            break
    if h_idx is None:
        print(f"[skip] {os.path.basename(path)}: no 'Data availability' heading")
        return

    # delete later duplicate "Post-hoc interpretive framework" paragraphs
    removed = 0
    for p in list(paras):
        if "Post-hoc interpretive framework" in p.text and \
           paras.index(p) > h_idx:
            p._element.getparent().remove(p._element)
            removed += 1

    heading = doc.paragraphs[h_idx]
    # insert new content paragraph right after the heading
    new_p = doc.add_paragraph(TEXT)        # appended at end first
    new_p.style = doc.styles["Normal"]
    heading._p.addnext(new_p._p)            # then move after heading

    doc.save(path)
    # verify
    t = " ".join(pp.text for pp in Document(path).paragraphs)
    print(f"  {os.path.basename(path)}: removed {removed} dup, "
          f"conceptDOI={'Y' if CONCEPT in t else 'N'}, "
          f"osf={'Y' if re.search(OSF, t, re.I) else 'N'}")


if __name__ == "__main__":
    for p in DOCX:
        if os.path.isfile(p):
            fix(p)
        else:
            print(f"MISSING {p}")
    print("\nDONE — run code/182 to verify.")
