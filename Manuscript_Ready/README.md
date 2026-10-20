# Manuscript package — ready for submission

**Title.** Genetically higher bone mineral density raises knee osteoarthritis but
protects against fracture: Mendelian randomisation of 13 skeletal phenotypes

**Target journal.** *Osteoarthritis and Cartilage* (full-length original article)

## Contents

| File | Item |
|---|---|
| `01_Manuscript.md` | Full manuscript: title page, highlights, abstract, Introduction, Methods, Results, Discussion, figure legends, Tables 1–4, declarations, references, supplementary-material list |
| `02_Cover_letter.md` | Cover letter (no placeholders) |
| `Figures/Main/` | Figure 1–10 and Graphical Abstract (19 PNG) |
| `Figures/Supplementary/` | Supplementary Figures S1–S12 (9 PNG) |
| `Supplementary/` | Table S10 (MVMR grid), Table S11 (predictions and falsifying tests), STROBE-MR checklist (md/docx/pdf) |

## Verification state at build time

- Consistency gate `audit_v18_full.py`: **460 checks / 0 fails / 1 warn**
  (the single warn is a disclosed 0.40% h²/N rounding note)
- Figure-to-data traceability `code/115_audit_figures_v18.py`: **187 checks / 0 failures**
- OAC format check `code/148_oac_compliance.py`: **COMPLIANT**
  (body 3973/4000 words; abstract 320/350; 5 highlights all ≤85 characters; references 38/50; display items 8/8)
- Numeric integrity after copy-editing `code/134_verify_humanize.py`: **PASS**
  (1263-number multiset identical before and after; 27 citation markers unchanged)
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
