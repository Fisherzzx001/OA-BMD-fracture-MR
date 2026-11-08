# -*- coding: utf-8 -*-
"""Generate the STROBE-MR checklist as both .docx and .pdf (elsevier-submission ready)."""
import io
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

STUDY = ("Genetically higher BMD raises knee osteoarthritis risk yet protects "
         "against fracture on pleiotropy-robust evidence")
GUIDE = ("Reporting guideline: STROBE-MR statement (Skrivankova VW et al. JAMA "
         "2021;326(16):1614-1621, doi:10.1001/jama.2021.18236) - 20 main items, 30 subitems.")

FOOTER = ("Overall: all 20 STROBE-MR items are addressed in the main manuscript and "
          "accompanying supplementary files. Two sub-items are not applicable to this "
          "design (11c, absolute-risk translation; 17b, timing of exposure) and are "
          "reported as such; no main item is Not Applicable.")

# (#, item, status, location)  -- official STROBE-MR numbering, verbatim item titles
ITEMS = [
    (1, "TITLE and ABSTRACT: Indicate Mendelian randomisation (MR) as the study's design in the title and/or the abstract if that is a main purpose of the study",
     "Reported", "Title ('...genetic correlation and bidirectional Mendelian randomisation across 13 skeletal phenotypes'); Abstract 'Design' paragraph"),
    (2, "INTRODUCTION - Background: scientific background and rationale; what is the exposure; is a potential causal relationship plausible; justify why MR is a helpful method",
     "Reported", "Section 1, paragraphs 1-2: OA-BMD debate since Foss and Byers 1972, confounding and reverse causation in observational designs, why MR is informative"),
    (3, "INTRODUCTION - Objectives: state specific objectives, including pre-specified causal hypotheses; state that MR intends to estimate causal effects under assumptions",
     "Reported", "Section 1, paragraphs 4-5 (three open questions; pre-specified predictions P1-P6 and falsifying tests F1-F10 in Table S11); MR as an assumption-bound estimator stated in Section 2 MR opening sentence and Section 4 limitations"),
    (4, "METHODS - Study design and data sources: (a) setting; (b) participants, eligibility, sample size; (c) measurement, quality control and selection of genetic variants; (d) methods of assessment and diagnostic criteria; (e) ethics approval and consent",
     "Reported", "Section 2 Data sources: 13 GWAS with accessions (Refs 14-17, 23, 26), N and case counts, Table 1 ancestry and heritability; variant QC, palindromic exclusion and clumping in Section 2 LDSC and MR; DXA, radiographic and self-report phenotype definitions in Section 2 Data sources and Section 3; original ethics approvals and no new approval required, and the summary-statistic access date (17 July 2026), in Data availability"),
    (5, "METHODS - Assumptions: explicitly state the three core IV assumptions (relevance, independence and exclusion restriction) as well as assumptions for additional or sensitivity analyses",
     "Reported", "Section 2 MR, opening sentence states all three assumptions verbatim. Assumptions of the added analyses are stated alongside them: MRlap (Section 2 Sample overlap) and MVMR including its no-measurement-error assumption and conditional F (Section 2 MVMR)"),
    (6, "METHODS - Statistical methods: main analysis: (a) scale, units, model; (b) how variants were handled and weighted; (c) MR estimator, covariates, and whether the same covariate set was used in both samples; (d) missing data; (e) multiple testing",
     "Reported", "Section 2 Effect-scale reconstruction (per-SD BMD exposure; log-odds reconstruction for binary outcomes via 1/sqrt[K(1-K)] = 3.389; no re-centring); Section 2 MR (inverse-variance weighting; random-effects multiplicative convention SE_random = SE_fixed x sqrt[max(1, Q/df)]); Section 2 MVMR (MVMR-IVW and MVMR-Egger; no individual-level covariates, so the covariate set is identical across samples); no missing data were imputed - complete-case harmonisation with attrition reported in Figure S9; multiple testing: Bonferroni alpha = 0.05/78 and Benjamini-Hochberg FDR across the 78-pair matrix, alpha = 0.005 in MR"),
    (7, "METHODS - Assessment of assumptions: describe any methods or prior knowledge used to assess the assumptions or justify their validity",
     "Reported", "Section 2 MR: Steiger directionality, per-SNP and aggregate F, R2_exposure, MR-Egger intercept, weighted median with bootstrap CIs, MR-PRESSO, leave-one-out, conservative 1-Mb clumping sensitivity. Section 2 Sample overlap: Burgess bound (<=2.6% BMD to OA, <=12% BMD to fracture) and MRlap. Section 2 Colocalisation: shared-causal-variant evidence at PP.H4 >= 0.8. Section 4 limitations: residual overlap and ancestry mismatch"),
    (8, "METHODS - Sensitivity analyses and additional analyses: describe any sensitivity analyses or additional analyses performed",
     "Reported", "Section 2 MR, Sample overlap and MVMR: internal positive control (BMD to fracture calibrated against GEFOS), reverse MR, MRlap sample-overlap correction, MVMR conditioning on height and separately on BMI, leave-one-out, MR-PRESSO outlier removal, colocalisation; results in Section 3"),
    (9, "METHODS - Software and pre-registration: (a) name statistical software and packages including version and settings; (b) state whether the protocol and details were pre-registered, and when and where",
     "Reported", "Section 2 opening paragraph: Python 3.13 with numpy, pandas, scipy and matplotlib; R 4.6.1 with MRlap 0.0.3.3. Data availability: the study pre-registration (P1-P6 as derived predictions, F1-F10 as registered tests with decision rules fixed before execution) is deposited with the public archive (Zenodo concept DOI 10.5281/zenodo.23400860)"),
    (10, "RESULTS - Descriptive data: (a) numbers of individuals at each stage and reasons for exclusion, consider a flow diagram; (b) summary statistics; (c) heterogeneity across contributing studies; (d) for two-sample MR, similarity of variant-exposure associations between samples and extent of sample overlap",
     "Reported", "Table 1 (SNP-based heritability, N, cases, LDSC intercept); Table S14 ancestry composition; Figure S9 instrument attrition flow; Section 3 paragraph 1; similarity of variant-exposure associations assessed by cross-trait rho_GE (Figure 9) and overlap bounded in Section 2 Sample overlap"),
    (11, "RESULTS - Main results: (a) variant-exposure and variant-outcome associations on an interpretable scale; (b) MR estimates with uncertainty on an interpretable scale, e.g. odds ratio per SD; (c) if relevant, translate relative risk into absolute risk; (d) consider plots to visualise results",
     "Reported (11c N/A)", "Table 2 (IVW beta, SE, p and random-effects p; MR-Egger slope and intercept; weighted median; OR per SD with 95% CI), Table 3, Figures 3-6 (forest, scatter, leave-one-out) and Section 3. Sub-item 11c was not translated to absolute risk: the source knee-OA GWAS defines no incidence window over which absolute risk could be computed"),
    (12, "RESULTS - Assessment of assumptions: (a) report the assessment of the validity of the assumptions; (b) report additional statistics such as heterogeneity across variants (I2, Q statistic or E-value)",
     "Reported", "Section 3 and Table 2: MR-Egger intercept p > 0.24 at all five sites; Cochran Q and I2 = 68-86% forward versus 78-83% in the positive control; Steiger z = 3.05-6.11 forward and non-significant in reverse; MR-PRESSO distortion significant for 9/10 pairs with 0-3 outliers each"),
    (13, "RESULTS - Sensitivity analyses and additional analyses: (a) robustness to violations of the assumptions; (b) other sensitivity analyses; (c) assessment of direction of causality, e.g. bidirectional MR; (d) comparison with non-MR estimates; (e) additional plots",
     "Reported", "Section 3: random-effects IVW; reverse MR; MVMR on height (Tables S9, S10) and BMI (Table S16); MRlap correction (Figure 10); comparison with the non-MR GEFOS estimate (OR 1.55) and with Hartley and Funck-Brentano (Figure S12); leave-one-out (Figure S2); multi-estimator robustness grid (Figure S11)"),
    (14, "DISCUSSION - Key results: summarise key results with reference to study objectives",
     "Reported", "Section 4, paragraph 1"),
    (15, "DISCUSSION - Limitations: taking into account validity of the IV assumptions, other sources of potential bias and imprecision; discuss direction and magnitude of any potential bias",
     "Reported", "Section 4 Strengths and limitations: ancestry mismatch versus European-only instruments; five prespecified sites deferred on compute; areal versus volumetric BMD; dependence of the significant count on the random-effects convention; reverse and OA-to-fracture power 0.14-0.26; residual overlap bounded at <=2.6% with its bias direction stated explicitly"),
    (16, "DISCUSSION - Interpretation: (a) meaning - cautious overall interpretation in context and against other studies; (b) mechanism - biological mechanisms that could drive the causal relationship and whether gene-environment equivalence is reasonable; (c) clinical relevance",
     "Reported", "Section 4: Systemic Osteogenic Set-point posed as a testable hypothesis, subchondral bone remodelling within the osteochondral unit, hip-versus-knee site contrast, explicit statement that these estimates do not support changing prescribing practice and that a joint cost requires trials powered for incident OA at the joint of interest; Conclusion"),
    (17, "DISCUSSION - Generalizability: (a) to other populations; (b) across other exposure periods/timings; (c) across other levels of exposure",
     "Reported (17b N/A)", "Section 4 limitations and Table S14: generalisation restricted to European-ancestry populations (87.31%-European knee-OA GWAS); F8 sex-stratified MR positive in both sexes. Sub-item 17b is not applicable because the germline instrument indexes a lifelong exposure with no definable timing; 17c is addressed by the per-SD scaling and by the MVMR attenuation ranges"),
    (18, "OTHER INFORMATION - Funding: sources of funding and role of funders, including for the original databases",
     "Reported", "Funding: None; Acknowledgements state no funding was received. Source GWAS consortia acknowledged"),
    (19, "OTHER INFORMATION - Data and data sharing: provide the data or report where and how they can be accessed, and provide the statistical code needed to reproduce the results",
     "Reported", "Data availability: all inputs are publicly published summary statistics (Refs 14-17, 23, 26) with accessions listed; every script and result table is archived at Zenodo concept DOI 10.5281/zenodo.23400860"),
    (20, "OTHER INFORMATION - Conflicts of interest: all authors declare all potential conflicts of interest",
     "Reported", "Declaration of competing interest (none declared)"),
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
fr = foot.add_run(FOOTER)
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
elems.append(Paragraph("<i>" + FOOTER + "</i>", sub_style))
pdf.build(elems)
print("WROTE", pdf_path)

# ---------------- MARKDOWN (single source of truth alongside the ITEMS list) ----------------
md_path = os.path.join(OUT_DIR, "STROBE-MR_checklist.md")
lines = [
    "# STROBE-MR Checklist (Skrivankova et al., *JAMA* 2021;326(16):1614-1621; Ref. 13)",
    "",
    "**Study:** " + STUDY + ".",
    "**Reporting guideline:** " + GUIDE.replace("Reporting guideline: ", ""),
    "Each item is marked **Reported**, with the location given as a section or "
    "display item of the main manuscript (`docs/08_\u6295\u7a3f\u7a3f_OA_BMD_fracture.md`); "
    "sub-items that cannot apply to this design are marked N/A.",
    "",
    "| # | STROBE-MR item | Status | Location in manuscript |",
    "|---|---|---|---|",
]
for num, item, status, loc in ITEMS:
    lines.append("| %d | %s | %s | %s |" % (num, item, status, loc))
lines += ["", "**" + FOOTER + "**", ""]
io.open(md_path, "w", encoding="utf-8").write("\n".join(lines))
print("WROTE", md_path)

