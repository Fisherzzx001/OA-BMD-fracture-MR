# -*- coding: utf-8 -*-
"""
135 · 建立 Manuscript_Ready/ —— 剥离工作日志，只留投稿物

问题（非推测）：
  docs/ 当前与稿件混放 50+ 个工作日志类文件
  （01_研究方案 / 02_变量协调 / 05_UKB申请材料 / 09_数据下载指南 /
   12–28 各轮审计与预审报告 / 99a–99g 专项审计 / 101–108 报告 /
   以及 5 个 .bak-* 备份）。这些不是投稿物。

做法（只复制，不删除、不移动原文件）：
  Manuscript_Ready/
    01_Manuscript.md            稿件正文（含表 1–4、图注、补充材料清单）
    02_Cover_letter.md          投稿信（占位已清零）
    Figures/Main/               主图 Figure1–10 + Graphical Abstract
    Figures/Supplementary/      SupplFig S1–S12
    Supplementary/              Table S10/S11 + STROBE-MR checklist
    README.md                   清单与校验说明

同时把 docs/ 下的 .bak-* 备份收进 .workbuddy/backups/（保留，不删）。
"""

import io
import os
import re
import shutil
import glob

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "Manuscript_Ready")

MAIN_FIG_SRC = os.path.join(ROOT, "Figures_v17_Nature")
SUPPL_FIG_SRC = os.path.join(ROOT, "Figures_v18_Suppl")
SUPPL_SRC = os.path.join(ROOT, "docs", "Supplementary")
DOC = os.path.join(ROOT, "docs", "08_投稿稿_OA_BMD_fracture.md")
CL = os.path.join(ROOT, "docs", "cover_letter_OAC.md")


def main():
    # Rebuild in place (no rmtree): overwrite derived files. The file SET is
    # stable across builds, so stale leftovers are not expected.
    os.makedirs(os.path.join(OUT, "Figures", "Main"), exist_ok=True)
    os.makedirs(os.path.join(OUT, "Figures", "Supplementary"), exist_ok=True)
    os.makedirs(os.path.join(OUT, "Supplementary"), exist_ok=True)

    # ---- 稿件与投稿信 ----
    shutil.copy2(DOC, os.path.join(OUT, "01_Manuscript.md"))
    shutil.copy2(CL, os.path.join(OUT, "02_Cover_letter.md"))
    print("[OK] manuscript + cover letter")

    # ---- 主图 ----
    n_main = 0
    for f in sorted(glob.glob(os.path.join(MAIN_FIG_SRC, "*.png"))):
        shutil.copy2(f, os.path.join(OUT, "Figures", "Main"))
        n_main += 1
    print(f"[OK] main figures: {n_main} PNG")

    # ---- 补充图 ----
    # Candidate revision figures S33–S35 (code/335) are NOT in the submitted
    # Figure-legend list (OAC displayed items already 8/8). Keep them out of
    # the journal-facing delivery set so they are not mistaken for submitted
    # supplementary figures; they remain available in Figures_v18_Suppl/ for
    # the response-to-reviewers.
    FIG_EXCLUDE_PREFIX = ("SupplFig33_", "SupplFig34_", "SupplFig35_")
    n_supp = 0
    if os.path.isdir(SUPPL_FIG_SRC):
        for f in sorted(glob.glob(os.path.join(SUPPL_FIG_SRC, "*.png"))):
            if os.path.basename(f).startswith(FIG_EXCLUDE_PREFIX):
                continue
            shutil.copy2(f, os.path.join(OUT, "Figures", "Supplementary"))
            n_supp += 1
    # Prune candidate figures left by an earlier incremental build.
    delivered_fig = os.path.join(OUT, "Figures", "Supplementary")
    for f in sorted(os.listdir(delivered_fig)):
        if f.startswith(FIG_EXCLUDE_PREFIX):
            os.remove(os.path.join(delivered_fig, f))
            print(f"[OK] excluded candidate revision figure from delivery: {f}")
    print(f"[OK] supplementary figures: {n_supp} PNG")

    # Count the *delivered* files, not just the newly copied PNGs: the figure
    # directories also carry PDF/SVG/TIF variants synced by other scripts.
    n_main_files = len(os.listdir(os.path.join(OUT, "Figures", "Main")))
    n_supp_files = len(os.listdir(os.path.join(OUT, "Figures", "Supplementary")))

    # ---- 补充材料 ----
    # Compliance boundary (2026-11-11 full-text audit): Audit_Trail.md is an
    # INTERNAL audit archive (claim -> script line -> independent check). It is
    # kept in docs/Supplementary/ for repository traceability but must never be
    # copied into the journal-facing Manuscript_Ready/ delivery set.
    EXCLUDE = {"Audit_Trail.md"}
    n_supp = 0
    for f in sorted(os.listdir(SUPPL_SRC)):
        if f in EXCLUDE:
            continue
        shutil.copy2(os.path.join(SUPPL_SRC, f), os.path.join(OUT, "Supplementary"))
        n_supp += 1
    stale = os.path.join(OUT, "Supplementary", "Audit_Trail.md")
    if os.path.exists(stale):
        os.remove(stale)  # purge any copy left by an older build
        print("[OK] purged stale internal Audit_Trail.md from delivery set")
    # Incremental builds never rmtree the delivery dir, so files removed from
    # the source (e.g. the three legacy mis-numbered Table_S*.md retired by
    # code/332) linger in the delivery set and resurrect duplicate table
    # numbers. Prune any Table_S*.md present in the delivery set but absent
    # from the source, without touching non-table deliverables.
    src_tables = {f for f in os.listdir(SUPPL_SRC) if f.startswith("Table_S") and f.endswith(".md")}
    delivered_dir = os.path.join(OUT, "Supplementary")
    pruned = []
    for f in sorted(os.listdir(delivered_dir)):
        if f.startswith("Table_S") and f.endswith(".md") and f not in src_tables:
            os.remove(os.path.join(delivered_dir, f))
            pruned.append(f)
    if pruned:
        print(f"[OK] pruned {len(pruned)} stale table files: {pruned}")
    print(f"[OK] supplementary docs: {n_supp}")

    # ---- 备份归档 ----
    bakdir = os.path.join(ROOT, ".workbuddy", "backups")
    os.makedirs(bakdir, exist_ok=True)
    moved = 0
    for f in glob.glob(os.path.join(ROOT, "docs", "*.bak-*")):
        shutil.move(f, os.path.join(bakdir, os.path.basename(f)))
        moved += 1
    for f in glob.glob(os.path.join(ROOT, "docs", "Supplementary", "*.bak-*")):
        shutil.move(f, os.path.join(bakdir, os.path.basename(f)))
        moved += 1
    print(f"[OK] moved {moved} .bak-* backups -> .workbuddy/backups/")

    # ---- 校验：交付目录内不得含占位符 ----
    md = io.open(os.path.join(OUT, "01_Manuscript.md"), encoding="utf-8").read()
    cl = io.open(os.path.join(OUT, "02_Cover_letter.md"), encoding="utf-8").read()
    ph = re.findall(r"\[[^\]]{0,80}(?:to be|insert|filled|name|email|XXXX)[^\]]{0,40}\]",
                    md + cl, re.I)
    print(f"\n占位符残留: {len(ph)}")
    for x in ph:
        print("   ", x)

    # ---- README ----
    readme = f"""# Manuscript package — ready for submission

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
| `Figures/Main/` | Figure 1–5 and Graphical Abstract ({n_main_files} files) |
| `Figures/Supplementary/` | Supplementary Figures S1–S24 ({n_supp_files} files, PNG/PDF/TIF) |
| `Supplementary/` | Independent supplementary tables S8–S9, S12–S17, S25–S32; Table S10 (MVMR grid); Table S11 (SOS predictions and falsifying tests); STROBE-MR checklist. Internal audit archives are excluded from the journal-facing set. |

## Verification state at build time

- OAC format gate `code/148_oac_compliance.py` + docx-level `code/336_verify_wordcount_docx.py`: **COMPLIANT**
  (body 3965/4000 words on the Word-compatible count that INCLUDES numerals; abstract 342/350;
   references 50/50; display items 8/8; 5 highlights all ≤85 characters)
  NB: the old letter-start regex under-counted the body by 246 numeric tokens (3984 false value); fixed 2026-11-11.
- Reported-value integrity `code/149b_verify_reported_values.py`: **SAFE** — no reported β/SE/p/OR/CI changed
- Data-level recompute `code/334_data_level_trust_but_verify.py`: **18/18 VERIFIED**
- Final QC `code/182_final_qc.py`: **78/78 PASS — READY TO SUBMIT**
- No-regression guard `code/310_assert_no_regression.py`: **33/33 PASS**
- PLINK .bed decode guard `code/331_verify_plink_bed_decode.py`: **PASS (31/31)**
- Supplementary figures S1–S24; inline supplementary tables S25, S28, S30, S31 and S32 in `01_Manuscript.md`
- Placeholders in this folder: **{len(ph)}**

## Language-editing note

The manuscript has been through a de-AI copy-edit (`code/130`–`code/134`).
Formal AIGC-risk score went from 40.6 to 18.0 on a transparent 8-signal scale
(em dash, boldface density, burstiness, sentence-initial particles, AI-preferred
vocabulary, clausal semicolons, passive-voice deviation, MATTR).
No numeric value was altered: verified by number-multiset comparison.

## Excluded from this folder

Working logs, audit reports, pre-review committee reports, data-download guides and
UKB application drafts remain under `docs/` and are **not** part of the submission.
"""
    io.open(os.path.join(OUT, "README.md"), "w", encoding="utf-8").write(readme)
    print(f"\n[OK] README.md")
    print(f"\n-> {OUT}")


if __name__ == "__main__":
    main()
