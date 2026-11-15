#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
448_n3n7_verify_carriers.py -- independent cross-carrier verification of the
N3/N4/N5/N6/N7 (MAJOR) second-pass remediation (code/446 + code/447).

Why this exists: the seven gates check word counts, locked strings and
reported-value integrity, but they do NOT check that every remediation
sentence actually reached all four carriers.  The N1/N2 episode showed that
text-only gates can pass while a claim is internally wrong; this script
closes the narrower "did the edit land everywhere" gap.

Read-only: never writes any manuscript, mirror, supplementary or docx file.

Usage: python3 code/448_n3n7_verify_carriers.py
"""
import hashlib
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "docs", "08_投稿稿_OA_BMD_fracture.md")
MIRROR = os.path.join(ROOT, "Manuscript_Ready", "01_Manuscript.md")
DOCXS = [
    os.path.join(ROOT, "submission_package", "Manuscript.docx"),
    os.path.join(ROOT, "submission_package", "Submission_MainText.docx"),
]
S17 = [
    os.path.join(ROOT, "docs", "Supplementary", "Table_S17_sensitivity.md"),
    os.path.join(ROOT, "Manuscript_Ready", "Supplementary", "Table_S17_sensitivity.md"),
]
S14 = [
    os.path.join(ROOT, "docs", "Supplementary", "Table_S14_ancestry.md"),
    os.path.join(ROOT, "Manuscript_Ready", "Supplementary", "Table_S14_ancestry.md"),
]
STROBE = [
    os.path.join(ROOT, "docs", "Supplementary", "STROBE-MR_checklist.md"),
    os.path.join(ROOT, "Manuscript_Ready", "Supplementary", "STROBE-MR_checklist.md"),
]

# --- new text that MUST be present in md, mirror AND both docx -------------
NEW_ALL = [
    ("N4 abstract",              "pleiotropy-robust estimators bound the causal reading"),
    ("N3 posterior rule",        "PP.H4≥0.8 is a posterior rule, not error-rate control"),
    ("N3 unadjusted p",          "unadjusted exploratory paired Wilcoxon"),
    ("N4 legs CI",               "legs CI spans zero"),
    ("N4/N5 strengths",          "sign-only estimator agreement with attenuated magnitudes (S17e)"),
    ("N5 PRESSO",                "not robustness against systematic pleiotropy"),
    ("N6 exclusion bound L71",   "exclude only |β| above their 95% CI upper bounds (0.090–0.105"),
    ("N6 Steiger not MDE",       "not an effect-size MDE"),
    ("N6 conclusion bound",      "exclude only |β| above 0.090–0.105"),
    ("N7 negative lambda",       "eight of ten negative, |λ|≤0.037"),
    ("N6 Table 3 note",          "MDE, exclusion limits and post-hoc power are three different quantities"),
    ("N6 CI upper 0.1047",       "0.1047 (OA→BMD_femoral_neck)"),
    ("N4/N5 S11 caption",        "Multi-estimator comparison grid for the ten MR pairs"),
    ("N4/N5 S20 caption",        "Sign agreement is descriptive, not robustness evidence"),
    ("N4/N5 S20 heading",        "Multi-estimator comparison of the ten MR pairs (direction agreement only)"),
]

# --- new text that is md/mirror-only (docx strips supplementary/index) ----
NEW_MD_ONLY = [
    ("N7 S14 main-text note",    "Adjudication of the negative intercepts (N7)"),
    ("S17 index (e)(f)(g)",      "(g) posterior thresholds"),
    ("N4/N5 S11 index entry",    "- **Figure S11.** Multi-estimator comparison grid for the ten MR pairs"),
    ("N4/N5 S20 index entry",    "- **Figure S20.** Multi-estimator comparison grid of the ten MR pairs"),
]

# --- obsolete wording that MUST be gone from every carrier ----------------
OLD_FORBIDDEN_ALL = [
    ("N4 abstract old",          "the causal reading rests on pleiotropy-robust estimators"),
    ("N3 PP.H4 as test old",     "PP.H4≥0.8), with a BMD–fracture control"),
    ("N4/N5 strengths old",      "directional concordance across IVW"),
    ("N5 PRESSO old",            "removing outliers did not reverse any direction"),
    ("N6 exclusion old L71",     "excluded only |β|≥0.09 in SD of BMD per unit OA liability"),
    ("N6 Steiger old",           "minimum detectable effect exceeds the null"),
    ("N6 conclusion old",        "exclude only |β|≥0.09; the positive knee-OA"),
    ("N7 lambda old",            "with no overlap fraction supplied (configuration"),
    ("N4/N5 S11 caption old",    "Multi-estimator robustness grid"),
    ("N4/N5 S20 heading old",    "Multi-estimator robustness of the ten MR pairs"),
]

NEW_SUPP = {
    "S17": [
        ("S17e N4 consequence",  "Inferential consequence (N4)"),
        ("S17e legs zero",       "the legs interval (−0.003 to +0.132)"),
        ("S17f exclusion",       "0.0897 / 0.0937 / 0.1047"),
        ("S17f incommensurable", "Forward and reverse β are **not** on one scale"),
        ("S17g block",           "### S17-g. Three statistical distinctions behind the N3/N4/N5/N6 remediation"),
        ("S17g no-values-change","No reported β/SE/p/OR/CI was changed by this remediation."),
    ],
    "S14": [
        ("S14 N7 correction",    "Correction of an earlier over-reading (N7)"),
        ("S14 lambda r",         "Pearson r = −0.23"),
        ("S14 winner's curse",   "weak-instrument and winner's-curse correction"),
    ],
}
OLD_SUPP = {
    "S17": [("S17f old scale",   "Forward and reverse estimates are on the same as-deposited scale")],
    "S14": [("S14 old over-read","no positive overlap-driven inflation needs to be subtracted")],
}


def read_bytes(p):
    with open(p, "rb") as f:
        return f.read()


def read_text(p):
    return read_bytes(p).decode("utf-8")


def docx_text(p):
    from docx import Document
    d = Document(p)
    parts = [para.text for para in d.paragraphs]
    for t in d.tables:
        for row in t.rows:
            for c in row.cells:
                parts.append(c.text)
    return "\n".join(parts)


def sha256(p):
    return hashlib.sha256(read_bytes(p)).hexdigest()


def main():
    checks = []

    def ck(ok, label, detail=""):
        checks.append((bool(ok), label, detail))

    md = read_text(SRC)
    mir = read_text(MIRROR)
    docx = {p: docx_text(p) for p in DOCXS}

    # 1. byte-level identity of the two carriers that must be identical
    ck(read_bytes(SRC) == read_bytes(MIRROR), "docs/08 == Manuscript_Ready mirror (byte)")
    ck(sha256(S17[0]) == sha256(S17[1]), "Table S17 dual copies identical (sha256)")
    ck(sha256(S14[0]) == sha256(S14[1]), "Table S14 dual copies identical (sha256)")
    ck(sha256(STROBE[0]) == sha256(STROBE[1]), "STROBE checklist dual copies identical (sha256)")

    # 2. every new sentence reached all four carriers
    carriers = [("md", md), ("mirror", mir)] + [
        (os.path.basename(p), docx[p]) for p in DOCXS
    ]
    for name, needle in NEW_ALL:
        for cname, txt in carriers:
            ck(needle in txt, f"NEW {name} @ {cname}")

    # 3. md/mirror-only blocks
    for name, needle in NEW_MD_ONLY:
        ck(needle in md, f"NEW {name} @ md")
        ck(needle in mir, f"NEW {name} @ mirror")

    # 4. obsolete wording fully removed
    for name, needle in OLD_FORBIDDEN_ALL:
        for cname, txt in carriers:
            ck(needle not in txt, f"OLD-GONE {name} @ {cname}")

    # 5. supplementary edits in BOTH copies, old wording gone
    for key, paths in (("S17", S17), ("S14", S14)):
        for p in paths:
            t = read_text(p)
            tag = "docs" if os.sep + "docs" + os.sep in p else "ready"
            for name, needle in NEW_SUPP[key]:
                ck(needle in t, f"SUPP-NEW {name} @ {key}/{tag}")
            for name, needle in OLD_SUPP[key]:
                ck(needle not in t, f"SUPP-OLD-GONE {name} @ {key}/{tag}")

    # 5b. STROBE-MR item 13 must describe S11 as a comparison grid in BOTH copies
    for p in STROBE:
        t = read_text(p)
        tag = "docs" if os.sep + "docs" + os.sep in p else "ready"
        ck("multi-estimator comparison grid (Figure S11; direction agreement only)" in t,
           f"STROBE-NEW comparison grid @ {tag}")
        ck("robustness grid" not in t, f"STROBE-OLD-GONE robustness grid @ {tag}")

    # 6. hard constraints: no statistics re-run / reported values untouched
    ck("achieved power is computed post-hoc" in md, "310 lock string preserved")
    ck("All five primary BMD→knee-OA effects survived family-2 Bonferroni" in md,
       "N1 fixed-vs-random adjudication preserved")
    ck("minimum 31.2, median 41.8, mean 64.8" in md,
       "N2 recomputed instrument strength preserved")

    # ---- report -----------------------------------------------------------
    nfail = sum(1 for ok, _, _ in checks if not ok)
    print("=== 448 N3-N7 CROSS-CARRIER VERIFICATION ===")
    print(f"docs/08 SHA256 = {sha256(SRC)}")
    print(f"mirror  SHA256 = {sha256(MIRROR)}")
    for p in DOCXS:
        print(f"docx    {os.path.basename(p)}: {len(docx[p])} chars of text")
    print(f"checks: {len(checks) - nfail}/{len(checks)} PASS, {nfail} FAIL")
    if nfail:
        for ok, label, _ in checks:
            if not ok:
                print(f"  [FAIL] {label}")
        print("VERDICT: FAIL - remediation did not land on every carrier")
        return 1
    print("VERDICT: PASS - every N3-N7 edit landed on all required carriers; "
          "all obsolete wording removed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
