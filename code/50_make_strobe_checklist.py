# -*- coding: utf-8 -*-
"""Generate the STROBE-MR checklist as both .docx and .pdf (elsevier-submission ready)."""
import os
from docx import Document
from docx.shared import Pt, RGBColor, Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT

from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer
from reportlab.lib.enums import TA_LEFT

OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "docs", "Supplementary")
os.makedirs(OUT_DIR, exist_ok=True)

STUDY = ("Genetically higher bone mineral density raises knee osteoarthritis but "
         "protects against fracture: Mendelian randomisation of 13 skeletal phenotypes")
GUIDE = "Reporting guideline: STROBE-MR (Skrivankova et al., BMJ 2021;372:n223)."

# (#, item, status, location)
ITEMS = [
    (1, "Title — identifies the report as a Mendelian randomisation study",
     "Reported", "Title: \"…Mendelian randomisation of 13 skeletal phenotypes\""),
    (2, "Abstract — follows STROBE-MR abstract structure (background, objectives, methods, results, conclusions)",
     "Reported", "Abstract: background/objectives/methods/results/conclusions present (unstructured per OAC)"),
    (3, "Introduction – Background: rationale for the MR investigation",
     "Reported", "§1 Introduction (Foss & Byers 1972 → present debate; why MR)"),
    (4, "Introduction – Objectives: specific objectives/hypotheses",
     "Reported", "§1 final paragraph (three open questions; 13-phenotype design)"),
    (5, "Methods – Data sources: GWAS provenance, search, eligibility, date accessed",
     "Reported", "§2 Data sources (13 GWAS; Refs 1–3,10; accessions)"),
    (6, "Methods – GWAS description: phenotype definition, participants, QC, effect allele, units, population",
     "Reported", "§2 Data sources; Table 1 (h², N, intercepts); §3 heritability"),
    (7, "Methods – Exposure instrument: selection, number, pruning/clumping, population, proxies, palindromes",
     "Reported", "§2 MR (p<5×10⁻⁸; r²<0.001/10 Mb; 20–35 per BMD; 116 knee OA; HapMap3 EUR)"),
    (8, "Methods – Outcome instrument / outcome GWAS: same detail",
     "Reported", "§2 Data sources (knee-OA GO 2025 N≈1.21M; self-rep IEU ukb-b-14486; fracture GEFOS)"),
    (9, "Methods – Assumption 1 (relevance): instrument–exposure association strength (F, R²)",
     "Reported", "§2 MR (instrument counts; Table 3; S9 conditional F ≥10); §2 Effect-scale reconstruction"),
    (10, "Methods – Assumption 2 (independence): pleiotropy assessment plan",
     "Reported", "§2 MR (MR-Egger intercept, MR-PRESSO, leave-one-out); §3 robustness"),
    (11, "Methods – Assumption 3 (exclusion restriction): no confounders of instrument–outcome",
     "Reported", "§2 Sample overlap (Burgess bound); §2 MRlap; §4 limitations (residual overlap)"),
    (12, "Methods – Statistical analysis: MR estimators, heterogeneity, sensitivity, software",
     "Reported", "§2 MR/LDSC/Colocalisation/MVMR (IVW, MR-Egger, weighted median, MR-PRESSO, LOO; Python 3.13; MRlap R 4.6.1)"),
    (13, "Methods – Sample overlap: assessment and handling",
     "Reported", "§2 Sample overlap (Burgess approximation; bound 2.6% BMD→OA, 12% BMD→fracture; MRlap)"),
    (14, "Results – GWAS descriptives: N, cases, h², intercept, population",
     "Reported", "Table 1; §3 Genetic-correlation architecture (intercepts 1.01–1.07)"),
    (15, "Results – Exposure instrument: SNPs, F-statistic, R²",
     "Reported", "§2 MR (counts); Table 3 (n_snp); S9 (Fc≥10)"),
    (16, "Results – Outcome instrument: same",
     "Reported", "Table 3 (knee OA 69 SNPs); §2 Data sources"),
    (17, "Results – Assumption testing: heterogeneity (Cochran Q/I²), pleiotropy (MR-Egger, MR-PRESSO), sensitivity",
     "Reported", "§3 (I² 68–86%; MR-Egger intercept p>0.24; MR-PRESSO 9/10 significant, 0–3 outliers); Table 3"),
    (18, "Results – Main MR estimates: with CIs, p-values, direction",
     "Reported", "Table 3 (IVW β, SE, p, random-effects p, Egger, weighted median, I², LOO); §3"),
    (19, "Results – Supplementary/stratified analyses: MVMR, coloc, reverse, MRlap",
     "Reported", "§3 (positive control, reverse MR, coloc, MVMR S9, MRlap Fig 10); S8, S9, S10"),
    (20, "Results – Sample overlap results: empirical overlap diagnostics",
     "Reported", "§3 (Fig 9 ρ_GE; §2 MRlap median +20.5% away from null)"),
    (21, "Discussion – Key results: interpreted in context",
     "Reported", "§4 Principal findings; Resolving the OA–BMD paradox; Site and direction"),
    (22, "Discussion – Limitations, interpretation, generalisability, conclusions, funding, ethics",
     "Reported", "§4 Strengths and limitations; Conclusion; Declarations (Funding: None; Competing interests: none); Data availability"),
]

# ---------------- DOCX ----------------
doc = Document()
style = doc.styles["Normal"]
style.font.name = "Calibri"
style.font.size = Pt(10.5)

h = doc.add_heading("STROBE-MR Checklist", level=1)
p = doc.add_paragraph()
run = p.add_run("Study: " + STUDY)
run.bold = True
run.font.size = Pt(10.5)
doc.add_paragraph(GUIDE)

table = doc.add_table(rows=1, cols=4)
table.style = "Table Grid"
table.alignment = WD_TABLE_ALIGNMENT.CENTER
hdr = table.rows[0].cells
for i, t in enumerate(["#", "STROBE-MR item", "Status", "Location in manuscript"]):
    hdr[i].text = t
    for pa in hdr[i].paragraphs:
        for r in pa.runs:
            r.bold = True
for num, item, status, loc in ITEMS:
    row = table.add_row().cells
    row[0].text = str(num)
    row[1].text = item
    row[2].text = status
    row[3].text = loc
# column widths
widths = [Inches(0.4), Inches(2.9), Inches(0.9), Inches(2.6)]
for r in table.rows:
    for i, w in enumerate(widths):
        r.cells[i].width = w

doc.add_paragraph()
foot = doc.add_paragraph()
fr = foot.add_run("Overall: all 22 STROBE-MR items are addressed in the main manuscript and accompanying supplementary files. No item is Not Applicable.")
fr.italic = True

docx_path = os.path.join(OUT_DIR, "STROBE-MR_checklist.docx")
doc.save(docx_path)
print("WROTE", docx_path)

# ---------------- PDF ----------------
styles = getSampleStyleSheet()
cell_style = ParagraphStyle("cell", parent=styles["Normal"], fontSize=8.2, leading=10.5, alignment=TA_LEFT)
hdr_style = ParagraphStyle("hcell", parent=styles["Normal"], fontSize=8.5, leading=10.5, alignment=TA_LEFT, textColor=colors.white)
title_style = ParagraphStyle("title", parent=styles["Title"], fontSize=15)
sub_style = ParagraphStyle("sub", parent=styles["Normal"], fontSize=9.5, leading=12)

pdf_path = os.path.join(OUT_DIR, "STROBE-MR_checklist.pdf")
pdf = SimpleDocTemplate(pdf_path, pagesize=A4,
                        leftMargin=15*mm, rightMargin=15*mm, topMargin=15*mm, bottomMargin=15*mm)
elems = []
elems.append(Paragraph("STROBE-MR Checklist", title_style))
elems.append(Spacer(1, 4))
elems.append(Paragraph("<b>Study:</b> " + STUDY, sub_style))
elems.append(Paragraph(GUIDE, sub_style))
elems.append(Spacer(1, 8))

data = [[Paragraph("#", hdr_style), Paragraph("STROBE-MR item", hdr_style),
         Paragraph("Status", hdr_style), Paragraph("Location in manuscript", hdr_style)]]
for num, item, status, loc in ITEMS:
    data.append([Paragraph(str(num), cell_style), Paragraph(item, cell_style),
                 Paragraph(status, cell_style), Paragraph(loc, cell_style)])

pt = Table(data, colWidths=[8*mm, 62*mm, 20*mm, 85*mm], repeatRows=1)
pt.setStyle(TableStyle([
    ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1F4E79")),
    ("GRID", (0, 0), (-1, -1), 0.5, colors.grey),
    ("VALIGN", (0, 0), (-1, -1), "TOP"),
    ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#EEF3F8")]),
    ("LEFTPADDING", (0, 0), (-1, -1), 4),
    ("RIGHTPADDING", (0, 0), (-1, -1), 4),
    ("TOPPADDING", (0, 0), (-1, -1), 3),
    ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
]))
elems.append(pt)
elems.append(Spacer(1, 8))
elems.append(Paragraph("<i>Overall: all 22 STROBE-MR items are addressed in the main manuscript and "
                       "accompanying supplementary files. No item is Not Applicable.</i>", sub_style))
pdf.build(elems)
print("WROTE", pdf_path)
