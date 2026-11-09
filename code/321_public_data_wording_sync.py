#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""321_public_data_wording_sync.py -- 2026-11-09 public-data recheck wording.

The 2026-11-09 GWAS Catalog recheck (`code/320`) found that non-European /
ancestry-stratified OA and BMD summary statistics DO exist but none is
fit for purpose, and that `GCST90566811` is the X-chromosome knee-OA GWAS,
not a European-restricted autosomal outcome.  Two blanket phrases in the
submission were therefore imprecise and one S11 sentence was factually wrong.

This script replaces those phrases only (no beta/SE/p/OR/CI is touched),
keeps the registered F10 status UNEXECUTED and the S23 "9/10" verdict, then
rebuilds every derived artefact from `docs/08` and verifies the sync.

    python3 code/321_public_data_wording_sync.py           # apply + sync
    python3 code/321_public_data_wording_sync.py --check   # verify only

Idempotent.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "docs" / "08_投稿稿_OA_BMD_fracture.md"
MIRROR = ROOT / "Manuscript_Ready" / "01_Manuscript.md"
S11 = ROOT / "docs" / "Supplementary" / "01_predictions_and_falsification.md"
S11_MIRROR = ROOT / "Manuscript_Ready" / "Supplementary" / "01_predictions_and_falsification.md"
S32 = ROOT / "docs" / "Supplementary" / "Table_S32_f10_ancestry_proxy.md"
S32_MIRROR = ROOT / "Manuscript_Ready" / "Supplementary" / "Table_S32_f10_ancestry_proxy.md"
GEN314 = ROOT / "code" / "314_f10_ancestry_proxy.py"
DOCX = [ROOT / "submission_package" / "Manuscript.docx",
        ROOT / "submission_package" / "Submission_MainText.docx"]

F10_OLD = ("F10's registered cross-ancestry replication not performed for want of "
           "public ancestry-stratified data, with an ancestry-composition robustness "
           "proxy reported in Table S32")
F10_NEW = ("F10's registered cross-ancestry replication not performed for want of "
           "ancestry-stratified summary statistics for this outcome study and "
           "ancestry-matched BMD instruments, with an ancestry-composition robustness "
           "proxy reported in Table S32")

S23_OLD = ("while direct non-European replication stays blocked by the absence of "
           "non-European BMD exposure summary statistics")
S23_NEW = ("while direct non-European replication stays blocked by the lack of "
           "fit-for-purpose non-European BMD instruments and site-matched "
           "ancestry-stratified knee-OA statistics")

F2_OLD = ("because no public summary-statistics file for a volumetric BMD trait "
          "was found in the GWAS Catalog")
F2_NEW = ("because no fit-for-purpose public common-SNP summary-statistics file "
          "for a volumetric BMD trait was found in the GWAS Catalog")

CLOSING_OLD = ("Closing F10 still requires non-European BMD exposure summary "
               "statistics or ancestry-stratified OA GWAS; F2's volumetric arm "
               "likewise awaits volumetric BMD/hip-osteophyte data from the "
               "Faber/Tobias releases")
CLOSING_NEW = ("Closing F10 still requires fit-for-purpose non-European BMD "
               "instruments and site-matched ancestry-stratified knee-OA summary "
               "statistics; F2's volumetric arm awaits volumetric BMD data and F5 "
               "awaits a hip-osteophyte GWAS from the Faber/Tobias releases")

S11_F10_OLD = ("The instruments are European-only, and direct non-European MR "
               "replication remains blocked by the absence of non-European BMD "
               "exposure summary statistics. Two routes to a true cross-ancestry "
               "replication are blocked by data availability, not design: the GWAS "
               "Catalog exposes no ancestry-stratified summary statistics for this "
               "study, and the available non-European BMD exposure (GCST90321120, "
               "5,862 South Asian participants) has no deposited summary statistics. "
               "A European-restricted outcome file (GCST90566811, 89,555 cases) has "
               "been staged so that the sensitivity of the estimate to the "
               "non-European outcome share can be quantified directly")
S11_F10_NEW = ("The instruments are European-only, and direct non-European MR "
               "replication remains blocked by the absence of fit-for-purpose data "
               "rather than by an absence of any non-European GWAS. A 2026-11-09 "
               "recheck of the GWAS Catalog (code/320) confirms that no "
               "ancestry-stratified summary statistics are deposited for this "
               "knee-OA outcome study [14], while the registered test requires de "
               "novo ancestry-matched BMD instruments plus a site-matched knee-OA "
               "endpoint. Public alternatives now exist but none is fit for purpose: "
               "Pan-UKB Central/South Asian heel BMD (GCST90692328, N=2,633; heel "
               "quantitative ultrasound rather than total-body DXA areal BMD), "
               "pooled multi-ancestry total-body BMD files without "
               "ancestry-specific instruments, MVP+UKB ancestry-stratified all-site "
               "OA (GCST90134280/81/82/87) and Genes & Health British South Asian OA "
               "excluding spine (GCST90727231; exome-wide rare-variant sequencing, "
               "GRCh38) — none is a knee-OA outcome with common-SNP statistics "
               "harmonisable to ancestry-matched BMD instruments. The previously "
               "staged GCST90566811 is the X-chromosome knee-OA GWAS (0/31 autosomal "
               "BMD instruments matched), not a European-restricted autosomal "
               "outcome, and is not used. The executed ancestry-composition robustness "
               "proxy in Table S32 quantifies dependence on the non-European outcome "
               "share within the mixed-ancestry GWAS and is explicitly not a "
               "cross-ancestry replication")

S11_F2_A = ("Volume arm: a volumetric or area-conditioned BMD instrument, for which "
            "no public summary-statistics file could be located")
S11_F2_A_NEW = ("Volume arm: a volumetric or area-conditioned BMD instrument, for "
                "which no fit-for-purpose public common-SNP summary-statistics file "
                "could be located")
S11_F2_B = ("No volumetric-BMD, bone-area or bone-width summary statistics exist in "
            "the GWAS Catalog (and Tobias 2024 femoral neck width was never "
            "deposited), so the volume arm cannot be run.")
S11_F2_B_NEW = ("No fit-for-purpose volumetric-BMD, bone-area or bone-width "
                "common-SNP summary statistics exist in the GWAS Catalog (the 2026 "
                "bone-area PheWAS GCST90826226–GCST90826232 scores copy-number "
                "burden from read depth rather than SNP genotypes; the trabecular "
                "vBMD GWAS GCST007014 has no deposited file; and Tobias 2024 femoral "
                "neck width was never deposited), so the volume arm cannot be run.")
S11_F5_OLD = ("No subtype-stratified OA GWAS exists in the GWAS Catalog, so the "
              "registered hip-osteophyte test is substituted by a publicly available "
              "site-stratified proxy")
S11_F5_NEW = ("No subtype-stratified common-SNP OA GWAS and no hip-osteophyte "
              "summary statistics exist in the GWAS Catalog (the Backman 2021 "
              "osteophyte files GCST90080462/64 are exome-wide rare-variant burden "
              "at knee or any site, with 0/31 matches to the locked BMD "
              "instruments), so the registered hip-osteophyte test is substituted "
              "by a publicly available site-stratified proxy")
S11_HB_1 = ("the volume arm of F2 is untested because no volumetric BMD, bone-area "
            "or bone-width summary statistics are public")
S11_HB_1_NEW = ("the volume arm of F2 is untested because no fit-for-purpose "
                "common-SNP volumetric BMD, bone-area or bone-width summary "
                "statistics are public")
S11_HB_2 = ("F10 is partial because no ancestry-stratified outcome summary statistics "
            "are deposited for re-running MR and no non-European BMD exposure "
            "summary statistics are available")
S11_HB_2_NEW = ("F10 is partial because no ancestry-stratified summary statistics "
                "are deposited for this knee-OA outcome and no fit-for-purpose "
                "non-European BMD instrument set is available (public non-European "
                "and ancestry-stratified GWAS exist but differ in phenotype, "
                "genotyping and anatomical site)")

# code/314 regenerates the standalone S32 file, so its template must carry the
# same correction.  Here the phrase is split over adjacent Python string literals.
GEN314_OLD = (
    '        "execution. (iv) Closing F10 still requires non-European BMD exposure summary statistics or "\n'
    '        "ancestry-stratified OA GWAS; F2\'s volumetric arm likewise awaits volumetric BMD/hip-osteophyte "\n'
    '        "data from the Faber/Tobias releases.",'
)
GEN314_NEW = (
    '        "execution. (iv) Closing F10 still requires fit-for-purpose non-European BMD "\n'
    '        "instruments and site-matched ancestry-stratified knee-OA summary statistics; "\n'
    '        "F2\'s volumetric arm awaits volumetric BMD data and F5 awaits a hip-osteophyte "\n'
    '        "GWAS from the Faber/Tobias releases.",'
)

# (path, old, new, expected_occurrences)
EDITS: list[tuple[Path, str, str, int]] = [
    (SRC, F10_OLD, F10_NEW, 1),
    (SRC, S23_OLD, S23_NEW, 1),
    (SRC, F2_OLD, F2_NEW, 2),                 # S27 Part A interpretation + F2 conclusion
    (SRC, CLOSING_OLD, CLOSING_NEW, 1),       # inline Table S32
    (S32, CLOSING_OLD, CLOSING_NEW, 1),       # standalone Table S32
    (GEN314, GEN314_OLD, GEN314_NEW, 1),      # generator template (split literals)
    (S11, S11_F10_OLD, S11_F10_NEW, 1),
    (S11, S11_F2_A, S11_F2_A_NEW, 1),
    (S11, S11_F2_B, S11_F2_B_NEW, 1),
    (S11, S11_F5_OLD, S11_F5_NEW, 1),
    (S11, S11_HB_1, S11_HB_1_NEW, 1),
    (S11, S11_HB_2, S11_HB_2_NEW, 1),
]

# strings that must no longer occur in any submission artefact
FORBIDDEN = [
    "public ancestry-stratified data",
    "no non-European BMD exposure summary statistics are available",
    "A European-restricted outcome file (GCST90566811",
    "the GWAS Catalog exposes no ancestry-stratified summary statistics for this study",
    "absence of non-European BMD exposure summary statistics",
]
MARKERS = [
    "ancestry-matched BMD instruments",
    "fit-for-purpose non-European BMD instruments",
    "X-chromosome knee-OA GWAS (0/31 autosomal",
    "9/10 executed tests",
]


def replace(path: Path, old: str, new: str, expected: int) -> str:
    txt = path.read_text(encoding="utf-8")
    n = txt.count(old)
    if n == expected:
        path.write_text(txt.replace(old, new), encoding="utf-8", newline="\n")
        return f"applied ({n})"
    if n == 0 and txt.count(new) == expected:
        return "already applied"
    raise RuntimeError(f"{path.name}: expected {expected} occurrence(s) of old "
                       f"phrase, found {n}; new phrase count={txt.count(new)}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()

    if not args.check:
        print("[1/4] wording edits")
        for path, old, new, n in EDITS:
            print(f"  {path.relative_to(ROOT)}: {replace(path, old, new, n)}")

        print("[2/4] mirrors (byte-identical copies)")
        shutil.copyfile(SRC, MIRROR)
        shutil.copyfile(S11, S11_MIRROR)
        shutil.copyfile(S32, S32_MIRROR)

        print("[3/4] rebuild two docx from docs/08")
        for script in ("build_submission.py", "build_maintext.py"):
            r = subprocess.run([sys.executable, script], cwd=ROOT,
                               capture_output=True, text=True)
            if r.returncode != 0:
                raise RuntimeError(f"{script} failed:\n{r.stdout}\n{r.stderr}")
            print(f"  rebuilt via {script}")
    else:
        print("[check] verifying only")

    print("[4/4] verification")
    artefacts = [SRC, MIRROR, S11, S11_MIRROR, S32, S32_MIRROR, GEN314]
    try:
        from docx import Document
        docx_text = ["\n".join(p.text for p in Document(p).paragraphs) for p in DOCX]
    except Exception as exc:                                   # noqa: BLE001
        raise RuntimeError(f"cannot read docx: {exc}")
    all_text = [p.read_text(encoding="utf-8") for p in artefacts] + docx_text

    src = SRC.read_text(encoding="utf-8")
    for m in MARKERS[:2]:
        assert m in src, f"docs/08 missing marker: {m}"
    assert "X-chromosome knee-OA GWAS (0/31 autosomal" in S11.read_text(encoding="utf-8"), \
        "S11 missing corrected GCST90566811 statement"
    assert "9/10" in src and "registered cross-ancestry replication still unexecuted" in src, \
        "F10 registered status / 9-10 verdict must be unchanged"

    for name, txt in zip([str(p.relative_to(ROOT)) for p in artefacts]
                         + [p.name for p in DOCX], all_text):
        bad = [f for f in FORBIDDEN if f in txt]
        assert not bad, f"{name} still contains stale wording: {bad}"
    assert SRC.read_bytes() == MIRROR.read_bytes(), "mirror not byte-identical to docs/08"
    assert S11.read_bytes() == S11_MIRROR.read_bytes(), "S11 mirror differs"
    assert S32.read_bytes() == S32_MIRROR.read_bytes(), "S32 mirror differs"
    for p, txt in zip(DOCX, docx_text):
        assert "ancestry-matched BMD instruments" in txt, f"{p.name} missing Data-availability edit"
        print(f"  {p.name}: wording OK")
    print("  docs/08 / S11 / S32 mirrors byte-identical: OK")
    print("\nSYNC COMPLETE -- wording precise, F10 status UNEXECUTED unchanged, "
          "no reported value touched")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
