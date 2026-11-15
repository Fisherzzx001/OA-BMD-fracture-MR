# -*- coding: utf-8 -*-
"""Reviewer-3 adversarial review (docs/151) — P0 fixes N1 + N2.

N1: resolve the mutually-exclusive multiplicity claims. L75 and L99 previously
    asserted "all five sites survive family-2 Bonferroni α=0.005" using FIXED-effect
    p-values, while the Abstract/Conclusion designate RANDOM-effects IVW as the
    primary inference (only 2/5 nominal, and 0/5 pass α=0.005: smallest random p=0.006).
    Fix: state both conventions explicitly — fixed 5/5, random 0/5 — no longer
    juxtaposing an inference with its own refutation. No p-value is altered.

N2: correct Table S27's instrument-strength description. The adjudication quoted
    "mean F=807, minimum F=168; SNP-level R2 29.4 vs 15.4", contradicted by the
    table's OWN 31 rows (per-SNP F=(beta/se)^2: min 31.2, median 41.8, mean 64.8;
    30/31 < 168) and by CSV instrument_strength_steiger.csv (BMD_total->OA_knee:
    F_min 31.1, F_agg 63.8, R2_exposure 0.0629). Fix: replace with recomputed,
    explicitly-defined values consistent with Tables 3/S13. (F/R2 are not among
    the reported beta/SE/p/OR/CI and are not gate-locked.)

Idempotent. Edits docs/08 (source), the S27 standalone supplement, and both docx
(N1 lives in the main body; N2's adjudication is inside the docx-stripped
supplement block, so docx needs only N1). Mirror is rebuilt separately via 135.
"""
from pathlib import Path
from docx import Document

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "docs/08_投稿稿_OA_BMD_fracture.md"
S27 = ROOT / "docs/Supplementary/Table_S27_hip_oa_replication.md"
DOCX = [ROOT / "submission_package/Manuscript.docx",
        ROOT / "submission_package/Submission_MainText.docx"]

# ---------------- N1 (main body; md + docx) ----------------
N1_L75_OLD = ("All significant findings survived the family-2 Bonferroni "
              "correction (α=0.005) across the ten primary MR pairs.")
N1_L75_NEW = ("All five primary BMD→knee-OA effects survived family-2 Bonferroni "
              "(α=0.005) under fixed effects; none survived under the random-effects "
              "primary convention (smallest p=0.006).")

N1_L99_OLD = ("family 2: the ten primary MR pairs (Bonferroni α=0.005, which all "
              "five BMD→knee-OA sites survive);")
N1_L99_NEW = ("family 2: the ten primary MR pairs (Bonferroni α=0.005; all five "
              "sites survive under fixed effects, none under random effects);")

# ---------------- N2 (supplement block; md + S27 standalone; NOT docx) ----------------
N2_L308_OLD = ("Instrument strength: mean F = 807, minimum F = 168; SNP-level R2 "
               "for the exposure 29.4 versus 15.4 for the outcome, Steiger-type "
               "directionality +0.312.")
N2_L308_NEW = ("Instrument strength (per-SNP F=(β_e/se_e)², recomputed from the 31 "
               "instruments tabulated below; exposure-side F is outcome-independent): "
               "minimum 31.2, median 41.8, mean 64.8, all >10; the matched "
               "total-body pair gives aggregate F 63.8 and R²_exposure 0.063 "
               "(Tables 3/S13; Qian 2025, N=31,873). Steiger-type directionality +0.312.")


def edit_text(path, triples):
    t = path.read_text(encoding="utf-8")
    out = []
    for name, old, new in triples:
        if new in t:
            out.append(f"  {name}: SKIP (already applied)")
        elif old in t:
            t = t.replace(old, new, 1)
            out.append(f"  {name}: APPLIED")
        else:
            out.append(f"  {name}: MISSING anchor")
    path.write_text(t, encoding="utf-8")
    return out


def edit_docx(path, triples):
    d = Document(path)
    out = []
    changed = False
    for name, old, new in triples:
        hits = 0
        for par in d.paragraphs:
            runs = par.runs
            full = "".join(r.text for r in runs)
            if old in full:
                merged = full.replace(old, new, 1)
                if runs:
                    runs[0].text = merged
                    for r in runs[1:]:
                        r.text = ""
                hits += 1
                changed = True
        if hits:
            out.append(f"  {name}: APPLIED ({hits} para)")
        elif new in "".join(p.text for p in d.paragraphs):
            out.append(f"  {name}: SKIP (already applied)")
        else:
            out.append(f"  {name}: MISSING anchor")
    if changed:
        d.save(path)
    return out


def main():
    print("== docs/08 (source) ==")
    for x in edit_text(SRC, [
        ("N1-L75", N1_L75_OLD, N1_L75_NEW),
        ("N1-L99", N1_L99_OLD, N1_L99_NEW),
        ("N2-L308", N2_L308_OLD, N2_L308_NEW),
    ]):
        print(x)
    print("== S27 standalone supplement ==")
    for x in edit_text(S27, [("N2-L308", N2_L308_OLD, N2_L308_NEW)]):
        print(x)
    for p in DOCX:
        print(f"== {p.name} ==")
        for x in edit_docx(p, [
            ("N1-L75", N1_L75_OLD, N1_L75_NEW),
            ("N1-L99", N1_L99_OLD, N1_L99_NEW),
        ]):
            print(x)


if __name__ == "__main__":
    main()
