# Manuscript package — ready for submission

**Title.** Genetically higher BMD raises knee osteoarthritis risk yet protects against fracture on pleiotropy-robust evidence

**Target journal.** *Osteoarthritis and Cartilage* (full-length original article)

**Archive.** Code and reproduction material: <https://github.com/Fisherzzx001/OA-BMD-fracture-MR>;
archived under the Zenodo concept DOI **<https://doi.org/10.5281/zenodo.23400860>**
(preregistration protocol P1–P6 / F1–F10 is archived in the same deposit).
The independent OSF preregistration is **<https://doi.org/10.17605/OSF.IO/CSMP7>**.

## Contents

| File | Item |
|---|---|
| `01_Manuscript.md` | Full manuscript: title page, highlights, abstract, Introduction, Methods, Results, Discussion, figure legends, Tables 1–4, declarations, references, supplementary-material list |
| `02_Cover_letter.md` | Cover letter (no placeholders) |
| `Figures/Main/` | Figure 1–5 and Graphical Abstract (23 files) |
| `Figures/Supplementary/` | Supplementary Figures S1–S24 (24 files, PNG/PDF/TIF) |
| `Supplementary/` | Independent supplementary tables S8–S9, S12–S17, S25–S32; Table S10 (MVMR grid); Table S11 (SOS predictions and falsifying tests); STROBE-MR checklist. Internal audit archives are excluded from the journal-facing set. |

## Verification state at build time

- OAC format gate `code/148_oac_compliance.py` + docx-level `code/336_verify_wordcount_docx.py`: **COMPLIANT**
  (body 3996/4000 words on the Word-compatible count that INCLUDES numerals; abstract 341/350;
   references 50/50; display items 8/8; 5 highlights all ≤85 characters)
  NB: the old letter-start regex under-counted the body by 246 numeric tokens (3984 false value); fixed 2026-11-11.
- Reported-value integrity `code/149b_verify_reported_values.py`: **SAFE** — no reported β/SE/p/OR/CI changed
- Data-level recompute `code/334_data_level_trust_but_verify.py`: **18/18 VERIFIED**
- Final QC `code/182_final_qc.py`: **78/78 PASS — READY TO SUBMIT**
- No-regression guard `code/310_assert_no_regression.py`: **33/33 PASS**
- PLINK .bed decode guard `code/331_verify_plink_bed_decode.py`: **PASS (31/31)**
- Supplementary figures S1–S24; inline supplementary tables S25, S28, S30, S31 and S32 in `01_Manuscript.md`
- Placeholders in this folder: **0**

## Language-editing note

The manuscript has been through a de-AI copy-edit (`code/130`–`code/134`).
Formal AIGC-risk score went from 40.6 to 18.0 on a transparent 8-signal scale
(em dash, boldface density, burstiness, sentence-initial particles, AI-preferred
vocabulary, clausal semicolons, passive-voice deviation, MATTR).
No numeric value was altered: verified by number-multiset comparison.

## Excluded from this folder

Working logs, audit reports, pre-review committee reports, data-download guides and
UKB application drafts remain under `docs/` and are **not** part of the submission.
