# -*- coding: utf-8 -*-
"""Regenerate submission_package/Cover_Letter.docx from docs/cover_letter_OAC.md.

The legacy build_cover_letter.py emits an obsolete title/date, so this converter
tracks the live markdown source of truth instead. Handles the md conventions used
by cover_letter_OAC.md: '# ' title, '- ' bullets, 'N. ' numbered, and inline
**bold** / *italic* (including inside the title).
"""
from docx import Document
from docx.shared import Pt
from docx.enum.text import WD_ALIGN_PARAGRAPH
import re, os, sys

SRC = "docs/cover_letter_OAC.md"
OUT = "submission_package/Cover_Letter.docx"

def add_runs(p, text):
    """Render inline **bold** and *italic* into a paragraph."""
    # Tokenise on **...** and *...*
    parts = re.split(r"(\*\*[^*]+\*\*|\*[^*]+\*)", text)
    for seg in parts:
        if not seg:
            continue
        if seg.startswith("**") and seg.endswith("**"):
            r = p.add_run(seg[2:-2]); r.bold = True
        elif seg.startswith("*") and seg.endswith("*"):
            r = p.add_run(seg[1:-1]); r.italic = True
        else:
            p.add_run(seg)

def main():
    with open(SRC, encoding="utf-8") as f:
        lines = f.read().splitlines()

    doc = Document()
    st = doc.styles["Normal"]
    st.font.name = "Times New Roman"
    st.font.size = Pt(11)

    for raw in lines:
        line = raw.rstrip()
        if not line.strip():
            continue
        if line.startswith("# "):
            p = doc.add_paragraph()
            p.style = doc.styles["Title"]
            add_runs(p, line[2:])
            continue
        if line.lstrip().startswith("- "):
            p = doc.add_paragraph(style="List Bullet")
            add_runs(p, line.lstrip()[2:])
            continue
        m = re.match(r"^\s*(\d+)\.\s+(.*)$", line)
        if m:
            p = doc.add_paragraph(style="List Number")
            add_runs(p, m.group(2))
            continue
        p = doc.add_paragraph()
        add_runs(p, line)

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    doc.save(OUT)
    print("Saved", OUT, "(%d paragraphs)" % len(doc.paragraphs))

if __name__ == "__main__":
    main()
