# -*- coding: utf-8 -*-
"""Build a journal-submission-ready Manuscript.docx from docs/08_投稿稿_OA_BMD_fracture.md.
Cleans draft scaffolding, resolves [REF] placeholders, restructures as a title-page + body docx.
"""
import re, os
from docx import Document
from docx.shared import Pt, RGBColor, Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

SRC = r"docs/08_投稿稿_OA_BMD_fracture.md"
OUT = (os.environ.get("OUTDIR", "submission_package") + r"/Manuscript.docx").replace("\\", "/")

# ---- author / affiliation block (from user input) ----
AUTHORS = "Zhenxin Zheng¹#, Song Chen¹#, Shenbo Huang¹, Chao Qian¹, Hang Yuan¹, Huricha Bao¹*, Qiuyan Jiang²*"
AFFIL1 = "¹ Department of Orthopedics, The Quzhou Affiliated Hospital of Wenzhou Medical University, Quzhou People's Hospital, Quzhou, Zhejiang Province, 324000, China"
AFFIL2 = "² Department of Surgery, The Quzhou Affiliated Hospital of Wenzhou Medical University, Quzhou People's Hospital, Quzhou, Zhejiang Province, 324000, China"
EQUAL = "# Zhenxin Zheng and Song Chen contributed equally to this work."
CORR = "* Corresponding authors: Huricha Bao, MD (baohuricha@126.com) and Qiuyan Jiang, MD (after9951023@gmail.com)."
FUND_COI = "Funding: None.   Competing interests: The authors declare no competing interests."

def read_text():
    with open(SRC, encoding="utf-8") as f:
        return f.read()

def clean(text):
    # drop draft-status lines
    text = re.sub(r"\*\*Manuscript draft — target journal.*?\*\*\n", "", text)
    text = re.sub(r"\*\*Status:.*?\*\*\n", "", text)
    # resolve [REF] placeholders
    repl = {
        "in 1972 [REF]": "in 1972 [11]",
        "inversely related [REF]": "inversely related [11]",
        "opposing conditions [REF]": "opposing conditions [11]",
        "[Bulik-Sullivan 2015, REF]": "[6]",
        "interrogate causality [REF]": "interrogate causality [7]",
        "[Hackinger 2017 HMG, REF]": "[4]",
        "[Zhou 2025 Front Genet, REF]": "[5]",
    }
    for k, v in repl.items():
        text = text.replace(k, v)
    if "[REF]" in text:
        import sys
        for m in re.finditer(r".{80}\[REF\].{30}", text):
            print("UNRESOLVED [REF]:", m.group(0), file=sys.stderr)
        raise SystemExit("Aborting: unresolved [REF] placeholder(s) remain in the manuscript.")
    # fix dangling cross-reference to the removed §5 checklist
    text = text.replace("(§5 and the Supplement, ", "(the Supplement, ")
    # append Foss & Byers 1972 as Ref 11 (real citation)
    text = text.replace(
        "PMID 30158200.",
        "PMID 30158200.\n11. Foss MV, Byers PD. Bone density, osteoarthrosis of the hip, and fracture. *Ann Rheum Dis*. 1972;31(3):259–264.",
    )
    # cut the §5 pre-submission checklist section (and everything after it)
    m = re.search(r"\n## §5 Pre-submission checklist.*$", text, re.S)
    if m:
        text = text[: m.start()]
    # extract data-availability sentence, then remove the §Authors block (now at end of text)
    da = ""
    mda = re.search(r"Summary statistics are from publicly published GWAS[^\n]*", text)
    if mda:
        da = mda.group(0)
    mauth = re.search(r"\n## Authors, affiliations, funding and competing interests\n", text)
    if mauth:
        text = text[:mauth.start()]
    return text.strip(), da

def parse_inline(p, text):
    """Add runs to paragraph p with **bold**, *italic*, `code`."""
    token_re = re.compile(r"(\*\*.+?\*\*|\*.+?\*|`.+?`)")
    last = 0
    for mt in token_re.finditer(text):
        if mt.start() > last:
            p.add_run(text[last:mt.start()])
        tok = mt.group(0)
        if tok.startswith("**"):
            r = p.add_run(tok[2:-2]); r.bold = True
        elif tok.startswith("`"):
            r = p.add_run(tok[1:-1]); r.font.name = "Consolas"; r.font.size = Pt(9)
        else:
            r = p.add_run(tok[1:-1]); r.italic = True
        last = mt.end()
    if last < len(text):
        p.add_run(text[last:])

def add_heading(doc, text, level):
    h = doc.add_heading(level=level)
    parse_inline(h, text)
    return h

def set_cell(cell, text, bold=False, size=8.5):
    cell.text = ""
    p = cell.paragraphs[0]
    parse_inline(p, text)
    for r in p.runs:
        r.font.size = Pt(size)
        if bold:
            r.bold = True
    # shade header row
    if bold:
        tcPr = cell._tc.get_or_add_tcPr()
        shd = OxmlElement("w:shd")
        shd.set(qn("w:val"), "clear"); shd.set(qn("w:color"), "auto"); shd.set(qn("w:fill"), "D9E2F3")
        tcPr.append(shd)

def add_table(doc, block):
    rows = [ln for ln in block if ln.strip().startswith("|")]
    # remove leading/trailing pipe and split
    data = []
    for ln in rows:
        cells = [c.strip() for c in ln.strip().strip("|").split("|")]
        data.append(cells)
    # data[0] header, data[1] separator, rest rows
    header = data[0]
    body = data[2:]
    ncol = len(header)
    t = doc.add_table(rows=1, cols=ncol)
    t.style = "Table Grid"
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    for j, h in enumerate(header):
        set_cell(t.rows[0].cells[j], h, bold=True)
    for r in body:
        cells = t.add_row().cells
        for j in range(ncol):
            set_cell(cells[j], r[j] if j < len(r) else "")
    # column widths
    try:
        for row in t.rows:
            for j, c in enumerate(row.cells):
                c.width = Inches(max(0.6, 6.5 / ncol))
    except Exception:
        pass
    doc.add_paragraph()

def build():
    text, da = clean(read_text())
    lines = text.split("\n")
    doc = Document()
    base = doc.styles["Normal"]
    base.font.name = "Times New Roman"
    base.font.size = Pt(11)
    base.paragraph_format.space_after = Pt(6)
    base.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY

    # ---- title page ----
    title = lines[0].lstrip("# ").strip()
    p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run(title); r.bold = True; r.font.size = Pt(15)
    p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.add_run(AUTHORS).font.size = Pt(11)
    p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.add_run(AFFIL1).font.size = Pt(9)
    p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.add_run(AFFIL2).font.size = Pt(9)
    p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.add_run(EQUAL).font.size = Pt(9); p.runs[0].italic = True
    p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.add_run(CORR).font.size = Pt(9)
    p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.add_run(FUND_COI).font.size = Pt(9)
    doc.add_paragraph()

    # ---- body ----
    i = 0
    n = len(lines)
    while i < n:
        s = lines[i].rstrip()
        if s.strip() == "" or s.strip() == "---" or s.strip() == "* * *":
            i += 1; continue
        if s.startswith("# "):
            i += 1; continue  # title already placed
        if s.startswith("## "):
            add_heading(doc, s[3:].strip(), 1); i += 1; continue
        if s.startswith("### "):
            add_heading(doc, s[4:].strip(), 2); i += 1; continue
        if s.strip().startswith("|"):
            block = []
            while i < n and lines[i].strip().startswith("|"):
                block.append(lines[i].rstrip()); i += 1
            add_table(doc, block); continue
        if s.startswith("- ") or s.startswith("  - ") or s.startswith("    - "):
            lvl = (len(s) - len(s.lstrip(" "))) // 2
            para = doc.add_paragraph(style="List Bullet" if lvl == 0 else "List Bullet 2")
            parse_inline(para, s.strip()[2:].strip())
            i += 1; continue
        # plain paragraph
        para = doc.add_paragraph()
        parse_inline(para, s.strip())
        i += 1

    # data availability at end
    if da:
        doc.add_paragraph()
        h = doc.add_heading(level=2); h.add_run("Data availability")
        p = doc.add_paragraph(); p.add_run(da)

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    doc.save(OUT)
    print("Saved", OUT)

if __name__ == "__main__":
    build()
