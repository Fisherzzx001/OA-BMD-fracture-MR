# -*- coding: utf-8 -*-
"""Generate Cover_Letter.docx for the OA/BMD/fracture manuscript."""
from docx import Document
from docx.shared import Pt
from docx.enum.text import WD_ALIGN_PARAGRAPH

OUT = r"submission_package/Cover_Letter.docx"
TITLE = ("Bone Mineral Density and Osteoarthritis Share Genetic Architecture but Not "
         "Causal Pathways to Fracture: A Cross-Trait Analysis of 13 Skeletal Phenotypes")
JOURNAL = "Bone"  # adjust to Osteoporosis International if preferred

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

    line(doc, "17 September 2026")
    line(doc, "")
    line(doc, f"Editor-in-Chief, {JOURNAL}")
    line(doc, f"{JOURNAL} Editorial Office")
    line(doc, "")
    line(doc, "Dear Editor,")
    line(doc, "")
    line(doc,
        f"Re: Submission of original research — \u201c{TITLE}\u201d",
        bold=True)
    line(doc, "")
    line(doc,
        "We submit the above original research article for consideration for "
        f"publication in {JOURNAL}. The work is not under consideration elsewhere, and all "
        "authors have approved the manuscript and its submission.")
    line(doc, "")
    line(doc,
        "Osteoarthritis (OA) and bone mineral density (BMD) are clinically entangled, yet "
        "whether their association reflects shared genetic aetiology, causal direction, or "
        "confounding has remained unresolved. We addressed this by analysing 13 skeletal "
        "phenotypes — two OA definitions, ten regional BMD sites, and any-fracture — with "
        "(i) LDSC across all 78 phenotype pairs on a common SNP set, (ii) two-sample "
        "Mendelian randomisation for causal direction with a built-in positive control "
        "(BMD\u2192fracture), and (iii) reverse MR to discriminate BMD\u2192OA from OA\u2192BMD. "
        "Instruments used strict PLINK-equivalent clumping (r\u00b2<0.001 over 10 Mb, 1000G "
        "Phase 3 EUR), and the MR-PRESSO distortion test was applied to all pairs.")
    line(doc, "")
    line(doc,
        "Our principal findings are: (1) a coherent three-block genetic architecture "
        "(BMD sites strongly inter-correlated; OA uniformly positively correlated with all "
        "BMD sites; fracture negatively correlated with all BMD sites, replicating the "
        "canonical relationship as a pipeline validation); (2) genetically elevated BMD "
        "causally increases knee OA risk, whereas OA does not causally raise BMD and has no "
        "causal effect on fracture; and (3) the BMD\u2192fracture MR estimate (OR 1.40\u20131.42 "
        "per SD lower BMD) agrees within ~10% with the independent GEFOS consortium benchmark "
        "(OR 1.55), externally validating the pipeline. These results reconcile the long-standing "
        "OA\u2013BMD paradox and suggest that interventions increasing bone mass should be "
        "evaluated for site-specific OA consequences.")
    line(doc, "")
    line(doc,
        "We believe this manuscript fits " + JOURNAL + " because it centres on the genetics of "
        "bone density, fracture and their skeletal comorbidity with OA, and provides an "
        "independent external validation of the established BMD\u2013fracture relationship while "
        "extending it to OA causality. All summary statistics derive from publicly published, "
        "consortium-approved GWAS and were not redistributed; full provenance, accessions and "
        "reproducibility commands are provided in the Data availability statement and "
        "Supplementary material.")
    line(doc, "")
    line(doc,
        "Funding: None. Competing interests: The authors declare no competing interests. "
        "Ethics: analyses used only published summary-level data; no new human-subject "
        "research was conducted.")
    line(doc, "")
    line(doc, "We look forward to your editorial assessment.")
    line(doc, "")
    line(doc, "Sincerely,")
    line(doc, "")
    line(doc, "Huricha Bao, MD (baohuricha@126.com)")
    line(doc, "Qiuyan Jiang, MD (after9951023@gmail.com)")
    line(doc, "Joint corresponding authors")
    line(doc,
        "On behalf of all authors — Zhenxin Zheng#, Song Chen#, Shenbo Huang, Chao Qian, "
        "Hang Yuan, Huricha Bao*, Qiuyan Jiang*", size=10)

    doc.save(OUT)
    print("Saved", OUT)

if __name__ == "__main__":
    main()
