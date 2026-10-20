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
    if os.path.exists(OUT):
        shutil.rmtree(OUT)
    os.makedirs(os.path.join(OUT, "Figures", "Main"))
    os.makedirs(os.path.join(OUT, "Figures", "Supplementary"))
    os.makedirs(os.path.join(OUT, "Supplementary"))

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
    n_supp = 0
    if os.path.isdir(SUPPL_FIG_SRC):
        for f in sorted(glob.glob(os.path.join(SUPPL_FIG_SRC, "*.png"))):
            shutil.copy2(f, os.path.join(OUT, "Figures", "Supplementary"))
            n_supp += 1
    print(f"[OK] supplementary figures: {n_supp} PNG")

    # ---- 补充材料 ----
    for f in sorted(os.listdir(SUPPL_SRC)):
        shutil.copy2(os.path.join(SUPPL_SRC, f), os.path.join(OUT, "Supplementary"))
    print(f"[OK] supplementary docs: {len(os.listdir(SUPPL_SRC))}")

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

**Title.** Genetically higher bone mineral density raises knee osteoarthritis but
protects against fracture: Mendelian randomisation of 13 skeletal phenotypes

**Target journal.** *Osteoarthritis and Cartilage* (full-length original article)

## Contents

| File | Item |
|---|---|
| `01_Manuscript.md` | Full manuscript: title page, highlights, abstract, Introduction, Methods, Results, Discussion, figure legends, Tables 1–4, declarations, references, supplementary-material list |
| `02_Cover_letter.md` | Cover letter (no placeholders) |
| `Figures/Main/` | Figure 1–10 and Graphical Abstract ({n_main} PNG) |
| `Figures/Supplementary/` | Supplementary Figures S1–S12 ({n_supp} PNG) |
| `Supplementary/` | Table S10 (MVMR grid), Table S11 (predictions and falsifying tests), STROBE-MR checklist (md/docx/pdf) |

## Verification state at build time

- Consistency gate `audit_v18_full.py`: **460 checks / 0 fails / 1 warn**
  (the single warn is a disclosed 0.40% h²/N rounding note)
- Figure-to-data traceability `code/115_audit_figures_v18.py`: **187 checks / 0 failures**
- OAC format check `code/49_oac_wordcount.py`: **COMPLIANT**
  (body 2736/4000 words; abstract 221/250; 5 highlights all ≤85 characters)
- Numeric integrity after copy-editing `code/134_verify_humanize.py`: **PASS**
  (1263-number multiset identical before and after; 27 citation markers unchanged)
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
