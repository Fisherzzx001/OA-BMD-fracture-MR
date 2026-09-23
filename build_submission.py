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

# ---- optional figure embedding (off by default) -------------------------------
# Set EMBED_FIG=1 and FIGDIR=<dir with Figure*.png> to also emit
# "Manuscript_with_Figures.docx", which inserts each figure after its legend.
# Journals normally want figures as separate files, so the default package keeps
# the text-only docx; this variant is for co-author / reviewer convenience only.
EMBED_FIG = os.environ.get("EMBED_FIG", "0") == "1"
FIGDIR = os.environ.get("FIGDIR", "Figures_v17_Nature")
FIG_MAP = {
    "Figure 1": "Figure1_rg_heatmap_13.png",
    "Figure 2": "Figure2_forest_key_results.png",
    "Figure 3": "Figure3_MR_forest.png",
    "Figure 4": "Figure4_h2_bar.png",
    "Figure 5": "Figure5_causal_network.png",
    "Figure 6": "Figure6_MR_scatter.png",
    "Figure 7": "Figure7_MR_funnel.png",
    "Figure 8": "Figure8_MR_power.png",
    "Figure 9": "Figure9_overlap_pleiotropy.png",
    "Figure 10": "Figure10_MRlap_overlap_correction.png",
    "Supplementary Figure S1": "SupplFig1_MR_scatter_remaining.png",
    "Supplementary Figure S2": "SupplFig2_leave_one_out.png",
    "Supplementary Figure S3": "SupplFig3_GEFOS_benchmark.png",
    "Supplementary Figure S4": "SupplFig4_STROBE_flow.png",
    "Supplementary Figure S5": "SupplFig5_OA_knee_N_sensitivity.png",
    "Supplementary Figure S6": "SupplFig6_fracture_N_decision.png",
    "Supplementary Figure S7": "SupplFig7_pleiotropy_diagnostics.png",
    "Graphical Abstract": "GraphicalAbstract.png",
}

def maybe_embed_figure(doc, raw_line):
    """If raw_line opens a figure legend, append the corresponding image."""
    if not EMBED_FIG:
        return False
    s = raw_line.strip()
    if not s.startswith("**"):
        return False
    key = None
    for k in FIG_MAP:
        if s.startswith("**" + k + ".**") or s.startswith("**" + k + ":**"):
            key = k
            break
    if key is None:
        return False
    path = os.path.join(FIGDIR, FIG_MAP[key]).replace("\\", "/")
    if not os.path.exists(path):
        print("  [warn] figure file missing, skipped:", path)
        return False
    from docx.shared import Inches as _Inches
    p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.add_run().add_picture(path, width=_Inches(6.2))
    return True

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
    # [REF] placeholders are resolved directly in the source markdown (docs/08).
    if "[REF]" in text:
        import sys
        for m in re.finditer(r".{80}\[REF\].{30}", text):
            print("UNRESOLVED [REF]:", m.group(0), file=sys.stderr)
        raise SystemExit("Aborting: unresolved [REF] placeholder(s) remain in the manuscript.")
    # fix dangling cross-reference to the removed §5 checklist
    text = text.replace("(§5 and the Supplement, ", "(the Supplement, ")
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
        maybe_embed_figure(doc, s)
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
    if EMBED_FIG:
        out_embed = OUT.replace(".docx", "_with_Figures.docx")
        os.replace(OUT, out_embed)
        print("Saved", out_embed)
        EMBED_FIG = False          # second pass: text-only canonical docx
        build()
