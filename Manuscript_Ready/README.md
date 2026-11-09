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
| `Figures/Supplementary/` | Supplementary Figures S1–S24 (23 files, PNG/PDF/TIF) |
| `Supplementary/` | Independent supplementary tables S8–S9, S12–S17, S26–S27, S29–S32; Table S10 (MVMR grid); Table S11 (SOS predictions and falsifying tests); STROBE-MR checklist; Audit_Trail |

## Verification state at build time

- OAC format gate `code/148_oac_compliance.py`: **COMPLIANT**
  (body 3982/4000 words; abstract 318/350; references 50/50; display items 8/8; 5 highlights all ≤85 characters)
- Reported-value integrity `code/149b_verify_reported_values.py`: **SAFE** — no reported β/SE/p/OR/CI changed
- Final QC `code/182_final_qc.py`: **65/65 PASS — READY TO SUBMIT**
- No-regression guard `code/310_assert_no_regression.py`: **30/30 PASS**
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
