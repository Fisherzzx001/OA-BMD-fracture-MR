# -*- coding: utf-8 -*-
"""Generate Cover_Letter.docx for the OA/BMD/fracture manuscript."""
from docx import Document
from docx.shared import Pt
from docx.enum.text import WD_ALIGN_PARAGRAPH

import os
OUT = (os.environ.get("OUTDIR", "submission_package") + r"/Cover_Letter.docx").replace("\\", "/")
TITLE = ("Bone Mineral Density and Osteoarthritis Share Genetic Architecture, but "
         "Evidence for an Osteoarthritis-to-Fracture Causal Pathway Remains Limited: "
         "A Cross-Trait Analysis of 13 Skeletal Phenotypes")
JOURNAL = "Osteoarthritis and Cartilage"  # target journal per P6 ruling (docs/103_期刊选定_P6裁决); OAC cover letter also at docs/cover_letter_OAC.md

def line(doc, text, bold=False, size=11, space=6, align=None):
    p = doc.add_paragraph()
    if align is not None:
        p.alignment = align
    p.paragraph_format.space_after = Pt(space)
    r = p.add_run(text)
    r.bold = bold; r.font.size = Pt(size)
    return p

def main():
    doc = Document()
    st = doc.styles["Normal"]; st.font.name = "Times New Roman"; st.font.size = Pt(11)

    line(doc, "23 September 2026")
    line(doc, "")
    line(doc, f"Editor-in-Chief, {JOURNAL}")
    line(doc, f"{JOURNAL} Editorial Office")
    line(doc, "")
    line(doc, "Dear Editor,")
    line(doc, "")
    line(doc,
        f'Re: Submission of original research: "{TITLE}"',
        bold=True)
    line(doc, "")
    line(doc,
        "We submit the above original research article for consideration for "
        f"publication in {JOURNAL}. The work is not under consideration elsewhere, and all "
        "authors have approved the manuscript and its submission.")
    line(doc, "")
    line(doc,
        "People with osteoarthritis (OA) tend to have higher bone mineral density (BMD), yet "
        "whether this reflects shared genetic aetiology, causal direction, or confounding has "
        "remained unresolved. We addressed this in 13 skeletal phenotypes (two OA definitions, "
        "ten regional BMD sites, and any-fracture), using LDSC across all 78 phenotype pairs "
        "on a common SNP set, two-sample Mendelian randomisation for causal direction with an "
        "internal positive control (BMD\u2192fracture), and reverse MR to distinguish "
        "BMD\u2192OA from OA\u2192BMD. Instruments were selected by strict PLINK-equivalent "
        "clumping (r\u00b2<0.001 over 10 Mb, 1000G Phase 3 EUR), and the MR-PRESSO distortion "
        "test was applied to all pairs.")
    line(doc, "")
    line(doc,
        "Our principal findings are: (1) a three-block genetic architecture (BMD sites strongly "
        "inter-correlated; OA positively correlated with all BMD sites; fracture negatively "
        "correlated with all BMD sites, replicating the established relationship and validating "
        "the pipeline); (2) genetically higher BMD raises knee OA risk, whereas OA does not "
        "raise BMD and has no causal effect on fracture; and (3) the BMD\u2192fracture MR "
        "estimate (OR 1.40\u20131.42 per SD lower BMD) falls within about 10% of the independent "
        "GEFOS consortium benchmark (OR 1.55). These results resolve the long-standing "
        "OA\u2013BMD paradox, and they argue for site-specific assessment of OA outcomes when "
        "bone mass is increased.")
    line(doc, "")
    line(doc,
        "This manuscript fits " + JOURNAL + " because it addresses a core question for the "
        "osteochondral unit \u2014 whether higher BMD reflects shared aetiology with or causation "
        "of OA \u2014 and because it independently reproduces the established BMD\u2013fracture "
        "relationship and extends it to OA causality. All summary statistics come from "
        "publicly published, consortium-approved GWAS and were not redistributed; full "
        "provenance, accessions and reproducibility commands are provided in the Data "
        "availability statement and Supplementary material.")
    line(doc, "")
    line(doc,
        "Funding: None. Competing interests: The authors declare no competing interests. "
        "Ethics: analyses used only published summary-level data; no new human-subject "
        "research was conducted.")
    line(doc, "")
    line(doc, "We hope the manuscript can be assessed for publication in " + JOURNAL + ".")
    line(doc, "")
    line(doc, "Sincerely,")
    line(doc, "")
    line(doc, "Huricha Bao, MD (baohuricha@126.com)")
    line(doc, "Qiuyan Jiang, MD (after9951023@gmail.com)")
    line(doc, "Joint corresponding authors")
    line(doc,
        "On behalf of all authors: Zhenxin Zheng#, Song Chen#, Shenbo Huang, Chao Qian, "
        "Hang Yuan, Huricha Bao*, Qiuyan Jiang*", size=10)

    doc.save(OUT)
    print("Saved", OUT)

if __name__ == "__main__":
    main()
