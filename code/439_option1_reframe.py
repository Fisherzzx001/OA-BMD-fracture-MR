# -*- coding: utf-8 -*-
"""
439 · Option-① reframing (descriptive architecture + triangulation; SOS transparently post-hoc)

Edits the single source of truth docs/08_投稿稿_OA_BMD_fracture.md and the two
delivered docx files in lock-step, preserving every gate-locked string and every
reported β/SE/p/OR/CI. Rerunnable / idempotent.

Source (docs/08) edits
  L93  Discussion   : reframe as descriptive multi-phenotype genetic architecture +
                      triangulation; SOS explicitly post-hoc / secondary; promote F6;
                      relocate `permutation-calibrated` + `not formation-side-exclusive`
                      (they only lived in the L95 register); KEEP the bone-forming triad
                      OR 3.00 [39] so OAC reference numbering and 149b stay green.
  L95  Discussion   : replace the in-body P1–P6/F1–F10 register with a concise
                      post-hoc pointer to Supplementary Material S11; retain the literal
                      "F2's decisive arm is unexecuted" (182 requires it).
  L99  Strengths    : promote F6 (cross-site homogeneity) as a genuine strength.
  S21  caption      : mark post-hoc interpretive framework, not central hypothesis.
  S22  caption      : make orthogonality explicitly neutral about a single set-point.
  S23  caption      : mark post-hoc; KEEP all theory-boundary locked strings intact.

Docx edits (submission_package/*.docx) mirror the same changes via run-level sync.
"""
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOC = os.path.join(ROOT, "docs", "08_投稿稿_OA_BMD_fracture.md")
DOCX = [
    os.path.join(ROOT, "submission_package", "Manuscript.docx"),
    os.path.join(ROOT, "submission_package", "Submission_MainText.docx"),
]

# ---------------------------------------------------------------------------
# New L93 — identical plain text for both source and docx (no markdown).
# Carries every relocated/locked phrase.
# ---------------------------------------------------------------------------
NEW_L93 = (
    "OA is an osteochondral-unit disorder [35–38]; the set-point we instrument acts "
    "upstream of subchondral remodelling. The 13-phenotype genetic architecture is, "
    "first and foremost, descriptive: BMD and the two OA definitions occupy two "
    "orthogonal genetic dimensions (PC1 bone-mass/fracture axis versus PC2 OA axis; "
    "Fig. S22, exploratory, OA h² 0.035), and triangulation across LDSC, MR and "
    "colocalisation shows one BMD-raising profile that is protective at the whole-bone "
    "endpoint yet pathogenic at the joint — the long-recognised bone–cartilage paradox "
    "[3, 4] (distinct in kind from the phenotype-level bone-forming triad — osteophytes, "
    "enthesophytes, high bone mass; OR 3.00, 95% CI 1.96–4.58 [39]), here quantified "
    "across ten sites rather than asserted. This descriptive separation is explicitly "
    "neutral about whether a single polygenic driver exists: orthogonality alone neither "
    "establishes nor refutes a single set-point (S22), and it is consistent also with two "
    "coupled modules. A post-hoc interpretive hypothesis, the Systemic Osteogenic "
    "Set-point (SOS; Figs. S21, S23), repackages this descriptive pattern plus one "
    "falsifiable single-set-point discriminator (a hypothesis, not a demonstrated "
    "mechanism): one latent set-point, read through a shared set of inherited variants, "
    "must fix opposite signs at whole-bone and joint endpoints. Its discriminating kernel "
    "is not supported by the data: the colocalisation split (9/31 fracture vs 0/31 "
    "knee-OA) is power-limited, so it reads neither as the predicted opposite-sign shared "
    "variant nor as distinct-variant refutation; a two-module mixture analysis (S23) gives "
    "ρ̂=−0.13 (95% CI −0.58 to +0.40, spanning zero), rejecting the strict single-set-point "
    "(ρ=−1) model yet unable to separate a weak set-point from two near-independent modules "
    "— the point estimate is negative (set-point direction) but the interval spans zero, so "
    "it indicates uncertainty, not refutation. The bone-versus-cartilage tissue enrichment "
    "(F3) is permutation-calibrated and mixed, not formation-side-exclusive (S25). Candidate "
    "cellular correlates span the osteochondral unit (RSPO3+ subchondral-bone precursors, "
    "NFATC1+ osteoclasts) but tissue enrichment is mixed, and bisphosphonate users had lower "
    "knee-replacement (RR 0.71, 0.66–0.77) [40], so areal BMD cannot by itself distinguish "
    "formation from resorption. Because areal BMD fuses mass, size and geometry, F2 tests "
    "whether geometry transmits the effect: a BMD-plus-hip-shape MVMR leaves the knee-OA "
    "conditional BMD effect positive in the better-powered modes (geometry does not absorb it), "
    "whereas hip-OA estimates are heterogeneous and sign-flipping and the volumetric arm is "
    "untested. A genuine strength of the architecture is its cross-site homogeneity: the five "
    "BMD→knee-OA estimates are statistically homogeneous (F6: Q=6.92, 4 df, p=0.14, τ²=0), "
    "favouring a systemic over a site-specific explanation and surviving covariance-corrected "
    "and conservative forms. The arms result (β=+0.075, p=1.7×10⁻⁶) supports the systemic "
    "clause, yet P2's ordering clause was negative (p=0.228). The directional pattern is "
    "supported; the osteogenic route and the set-point core are not."
)

# New L95 — markdown bold heading for the source; plain for the docx.
NEW_L95_MD = (
    "**Post-hoc interpretive framework (SOS predictions and falsifying tests).** The six "
    "predictions (P1–P6) and ten falsifying tests (F1–F10) cited above were derived post-hoc "
    "from the descriptive architecture and are disclosed as such in the preregistration (Data "
    "availability), which also discloses unexecuted arms (F2 volume, F10 cross-ancestry). They "
    "constitute an interpretive framework, not the central hypothesis of this study. The full "
    "register, with the current judgement on each item, is provided as Supplementary Material "
    "S11; in brief, the directional pattern is supported across tests, whereas the osteogenic "
    "route and the single-set-point core are not demonstrated. F2's decisive arm is unexecuted."
)
NEW_L95_DOCX = NEW_L95_MD.replace(
    "**Post-hoc interpretive framework (SOS predictions and falsifying tests).** ",
    "Post-hoc interpretive framework (SOS predictions and falsifying tests). ")

# L99 Strengths — promote F6 (append before "Limitations:").
L99_OLD = "directional concordance across IVW, MR-Egger, weighted median and MR-PRESSO. Limitations:"
L99_NEW = ("directional concordance across IVW, MR-Egger, weighted median and MR-PRESSO, "
           "and cross-site homogeneity of the BMD→knee-OA effect (F6: Q=6.92, 4 df, p=0.14, "
           "τ²=0) favouring a systemic over site-specific architecture. Limitations:")
L99_DONE = "favouring a systemic over site-specific architecture. Limitations:"

# S21 caption — mark post-hoc.
S21_OLD = "but not hip OA (Faber 2025 [28])."
S21_NEW = ("but not hip OA (Faber 2025 [28]). This figure presents a post-hoc interpretive "
           "framework (SOS), not the central hypothesis of this study; the descriptive "
           "architecture is reported in §3–§4.")
S21_DONE = "This figure presents a post-hoc interpretive framework"

# S22 caption — orthogonality explicitly neutral about a single set-point.
S22_OLD = "the two axes are orthogonal, so BMD and OA occupy distinct genetic dimensions."
S22_NEW = ("the two axes are orthogonal, so BMD and OA occupy distinct genetic dimensions; "
           "this orthogonality is a descriptive result and is explicitly neutral about whether "
           "a single set-point exists (S23) — it neither establishes nor refutes one.")
S22_DONE = "this orthogonality is a descriptive result"

# S23 caption — mark post-hoc (append at end; locked strings untouched).
S23_OLD = ("the ancestry-composition robustness proxy in S32 is a within-outcome consistency "
           "check and does not substitute for that replication.")
S23_NEW = (S23_OLD + " This framework was derived post-hoc from the descriptive architecture "
           "and is disclosed as such in the preregistration (Data availability); it is an "
           "interpretive hypothesis, not the central claim of this study.")
S23_DONE = "This framework was derived post-hoc from the descriptive architecture"


# ---------------------------------------------------------------------------
# Source edits
# ---------------------------------------------------------------------------
def edit_source():
    text = open(DOC, encoding="utf-8").read()
    lines = text.split("\n")
    applied = []

    # L93 (1-based line 93 -> index 92) full replace; guard on the opening phrase.
    if lines[92].startswith("OA is an osteochondral-unit disorder"):
        lines[92] = NEW_L93
        applied.append("L93 reframe")
    else:
        print("[skip] L93 (opening phrase absent)")

    # L95 (index 94) full replace; guard on the SOS-register phrase.
    if "SOS predictions and falsifying tests" in lines[94]:
        lines[94] = NEW_L95_MD
        applied.append("L95 -> post-hoc pointer")
    else:
        print("[skip] L95 (SOS-register phrase absent)")

    text = "\n".join(lines)

    # Substring edits (unique anchors) with idempotency guards.
    for old, new, done, tag in [
        (L99_OLD, L99_NEW, L99_DONE, "L99 promote F6"),
        (S21_OLD, S21_NEW, S21_DONE, "S21 post-hoc"),
        (S22_OLD, S22_NEW, S22_DONE, "S22 neutrality"),
        (S23_OLD, S23_NEW, S23_DONE, "S23 post-hoc"),
    ]:
        if done in text:
            print(f"[skip] {tag} (already applied)")
            continue
        if old in text:
            text = text.replace(old, new, 1)
            applied.append(tag)
        else:
            print(f"[skip] {tag} (anchor absent)")

    open(DOC, "w", encoding="utf-8").write(text)
    print("docs/08 edits applied:", applied)


# ---------------------------------------------------------------------------
# Docx edits
# ---------------------------------------------------------------------------
def merge_replace(runs_getter, old, new):
    runs = runs_getter()
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


def swap_paragraph(par, new):
    runs = par.runs
    if not runs:
        return 0
    runs[0].text = new
    for r in runs[1:]:
        r.text = ""
    return 1


def edit_docx(path):
    from docx import Document
    doc = Document(path)
    total = 0
    for par in doc.paragraphs:
        full = "".join(r.text for r in par.runs)
        if full.startswith("OA is an osteochondral-unit disorder"):
            total += swap_paragraph(par, NEW_L93)
        elif "SOS predictions and falsifying tests" in full:
            total += swap_paragraph(par, NEW_L95_DOCX)
    for old, new, done, tag in [
        (L99_OLD, L99_NEW, L99_DONE, "L99"),
        (S21_OLD, S21_NEW, S21_DONE, "S21"),
        (S22_OLD, S22_NEW, S22_DONE, "S22"),
        (S23_OLD, S23_NEW, S23_DONE, "S23"),
    ]:
        # locate the paragraph containing the marker to append only once
        n = 0
        for par in doc.paragraphs:
            ptext = "".join(r.text for r in par.runs)
            if done in ptext:
                break
            if old in ptext:
                n += merge_replace(lambda p=par: p.runs, old, new)
        if n:
            total += n
        else:
            print(f"  [skip docx {os.path.basename(path)}] {tag} (anchor absent)")
    doc.save(path)
    print(f"  {os.path.basename(path)}: {total} docx edits")


if __name__ == "__main__":
    edit_source()
    for p in DOCX:
        if os.path.isfile(p):
            edit_docx(p)
        else:
            print(f"MISSING {p}")
    print("\nDONE — run code/310 to verify.")
