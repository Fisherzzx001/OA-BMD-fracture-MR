# Submission package — *Osteoarthritis and Cartilage*

This file tracks what is required to submit the manuscript and what has been prepared.

## Manuscript compliance (verified by `code/49_oac_wordcount.py`)
- IMRaD body: **2,375 / 4,000 words** ✅
- Abstract: **210 / 250 words** ✅
- Highlights: **5 bullets**, each ≤ 85 characters ✅
- References: Vancouver, continuous 1–36 ✅

## Files ready for upload
| Item | File | Notes |
| --- | --- | --- |
| Manuscript (full text) | `docs/08_投稿稿_OA_BMD_fracture.md` | convert to .docx/.pdf at upload |
| Cover letter | `docs/cover_letter_OAC.md` | fill reviewer names + corresponding-author email |
| STROBE-MR checklist | `docs/Supplementary/STROBE-MR_checklist.docx` (or `.pdf`) | 22 items, all Reported |
| Main figures | `Figures/Figure1.png` … `Figure10.png` (+ `.pdf` vectors) | 10 main figures |
| Graphical Abstract | `Figures/GraphicalAbstract.png` (+ `.pdf`) | |
| Supplementary figures | `Figures/SupplFig1.png` … `SupplFig8.png` (+ `.pdf`) | S1–S8 |
| Supplementary tables | `docs/Supplementary/01_predictions_and_falsification.md`, `02_mvmr_height_grid.md` | S10 (MVMR grid), S11 (SOS P/F) |
| Colocalisation table | `docs/08` Table S8 (in-manuscript) | S8 |
| CRediT author contributions | `docs/08` Author contributions table | |

## Pre-submission actions still requiring the author's account
1. **Public GitHub repository** — create `OA-BMD-fracture-MR` and push the repository (scripts + small derived outputs + figures). The bulky regenerable intermediates (`code/output/genetics/*_clump_input.txt`, raw extracts) are excluded from git and archived in Zenodo instead (see `.gitignore`).
2. **Zenodo DOI** — after the GitHub repo exists, create a Zenodo upload linked to the repo (or upload the full `code/output/` archive directly) to obtain a DOI.
3. **Backfill Data availability** — replace the two `to be inserted here` placeholders in `docs/08` with the repository URL and the Zenodo DOI.

## Commands (run by the author, or by the assistant once a token path is provided)
```bash
cd OA_B_PROJECT
# 1) create public repo via GitHub API (token used once, then stripped)
curl -H "Authorization: Bearer <TOKEN>" -d "{\"name\":\"OA-BMD-fracture-MR\",\"private\":false}" https://api.github.com/user/repos
# 2) push
git remote add origin https://<USER>:<TOKEN>@github.com/<USER>/OA-BMD-fracture-MR.git
git push -u origin master
git remote set-url origin https://github.com/<USER>/OA-BMD-fracture-MR.git   # strip token
# 3) Zenodo -> New upload -> link GitHub repo -> get DOI -> backfill docs/08
```
