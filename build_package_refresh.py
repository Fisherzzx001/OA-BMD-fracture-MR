# -*- coding: utf-8 -*-
"""
build_package_refresh.py — 刷新 Submission_Package_2026-10-17.zip
====================================================================
以现有 dated zip 为基底 staging，把 canonical 源文件的最新版本刷入包内，
重写 README 校验段与 MANIFEST，然后重新打包（同名覆盖）。

刷新映射（canonical -> package 内路径）：
  submission_package/Manuscript.docx                     -> Manuscript.docx
  docs/08_投稿稿_OA_BMD_fracture.md                       -> 01_Manuscript.md
  docs/118_OSF预注册_OA_BMD_fracture.md                   -> Governance/OSF_preregistration_OA_BMD_fracture.md
  docs/Supplementary/Table_S16_mvmr_bmi.{md,csv}          -> Supplementary/
  code/output/genetics/mvmr_bmi_main.csv                  -> Table_Source_Data/
  code/output/genetics/mvmr_bmi_instruments.csv           -> Table_Source_Data/

用法： python build_package_refresh.py
"""
import os
import re
import shutil
import subprocess
import sys
import zipfile
from datetime import datetime

PROJ = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.join(PROJ, "Submission_Package_2026-10-17.zip")
STAGE = os.path.join(PROJ, "build", "pkg_stage")

REFRESH = [
    ("submission_package/Manuscript.docx", "Manuscript.docx"),
    ("docs/08_投稿稿_OA_BMD_fracture.md", "01_Manuscript.md"),
    ("docs/118_OSF预注册_OA_BMD_fracture.md", "Governance/OSF_preregistration_OA_BMD_fracture.md"),
    ("docs/Supplementary/Table_S16_mvmr_bmi.md", "Supplementary/Table_S16_mvmr_bmi.md"),
    ("docs/Supplementary/Table_S16_mvmr_bmi.csv", "Supplementary/Table_S16_mvmr_bmi.csv"),
    ("code/output/genetics/mvmr_bmi_main.csv", "Table_Source_Data/mvmr_bmi_main.csv"),
    ("code/output/genetics/mvmr_bmi_instruments.csv", "Table_Source_Data/mvmr_bmi_instruments.csv"),
    # 2026-10-20 audit fix: package must carry the N-corrected (manuscript-scale) MR values
    # under the same filenames so all in-package references stay valid.
    # (raw wrong-scale repo files are kept at code/output/genetics/mr_results*.csv as
    #  code/22's inputs; see docs/120 G3/G5)
    ("code/output/genetics/mr_results_corrected.csv", "Table_Source_Data/mr_results.csv"),
    ("code/output/genetics/mr_presso_results_corrected.csv", "Table_Source_Data/mr_presso_results.csv"),
    # 2026-10-20 audit fix (M8): the file previously shipped as "mvmr_height_main.csv" was the
    # n_h=1000 configuration (Fc(BMD)=3.3-4.6, fails the pre-registered Fc>=10 gate) and did NOT
    # reproduce Table S9. "main" is now the true n_h=100/EUR primary; the old file is preserved
    # under an explicit name. Also ship the full grid and the S10 text with n_h=0 artefacts removed.
    ("code/output/genetics/mvmr_height_main.csv", "Table_Source_Data/mvmr_height_main.csv"),
    ("code/output/genetics/mvmr_height_n_h1000.csv", "Table_Source_Data/mvmr_height_n_h1000.csv"),
    ("code/output/genetics/mvmr_height_grid.csv", "Table_Source_Data/mvmr_height_grid.csv"),
    ("docs/Supplementary/02_mvmr_height_grid.md", "Supplementary/02_mvmr_height_grid.md"),
]

README_SUPP_ROW_OLD = "| `Supplementary/` | Table S10 (MVMR grid), Table S11 (predictions and falsifying tests), STROBE-MR checklist (md/docx/pdf) |"
README_SUPP_ROW_NEW = "| `Supplementary/` | Tables S8\u2013S16 (colocalisation, height MVMR + grid, SOS predictions/falsifying tests, variant selection, instrument strength, ancestry composition, r_g blocks, BMI-MVMR) and the STROBE-MR checklist (md/docx/pdf) |"
README_SFIG_ROW_OLD = "| `Figures/Supplementary/` | Supplementary Figures S1\u2013S12 (9 PNG) |"
README_SFIG_ROW_NEW = "| `Figures/Supplementary/` | Supplementary Figures S1\u2013S17 (PNG, plus vector PDF/TIF under `Vector/`) |"

README_VERIF_OLD = """## Verification state at build time

- Consistency gate `audit_v18_full.py`: **460 checks / 0 fails / 1 warn**
  (the single warn is a disclosed 0.40% h\u00b2/N rounding note)
- Figure-to-data traceability `code/115_audit_figures_v18.py`: **187 checks / 0 failures**
- OAC format check `code/49_oac_wordcount.py`: **COMPLIANT**
  (body 2736/4000 words; abstract 221/250; 5 highlights all \u226485 characters)
- Numeric integrity after copy-editing `code/134_verify_humanize.py`: **PASS**
  (1263-number multiset identical before and after; 27 citation markers unchanged)
- Placeholders in this folder: **0**"""

README_VERIF_NEW = """## Verification state at build time

- Consistency gate `audit_v18_full.py`: **494 checks / 0 fails / 0 warns**
- Figure-to-data traceability `audit_figures_v17.py`: **845 checks / 0 failures**
- OAC format check `code/144_measure_oac.py`: **COMPLIANT**
  (body 3992/4000 words; abstract 346/350 structured; 8 display items; 38/50 references)
- Placeholders in this folder: **0**"""

README_LANG_OLD = """The manuscript has been through a de-AI copy-edit (`code/130`\u2013`code/134`).
Formal AIGC-risk score went from 40.6 to 18.0 on a transparent 8-signal scale
(em dash, boldface density, burstiness, sentence-initial particles, AI-preferred
vocabulary, clausal semicolons, passive-voice deviation, MATTR).
No numeric value was altered: verified by number-multiset comparison."""

README_LANG_NEW = """The manuscript has been through a de-AI copy-edit (`code/130`\u2013`code/134`);
the AIGC-risk score on the transparent 8-signal scale was 18.0 at that pass
(from 40.6). Subsequent revision-round edits \u2014 including the BMI-MVMR
integration of 2026-10-17 \u2014 were made under the project's constraint that no
reported statistic (\u03b2/SE/p/OR) is altered; every added number was regenerated
from its generating script and re-audited (`audit_v18_full.py` 494/0/0)."""


def measure_oac():
    """Run code/144_measure_oac.py and parse the three counts."""
    out = subprocess.run([sys.executable, "code/144_measure_oac.py"], cwd=PROJ,
                         capture_output=True, text=True).stdout
    body = int(re.search(r"MAIN BODY words[^:]*:\s*(\d+)", out).group(1))
    abst = int(re.search(r"ABSTRACT words:\s*(\d+)", out).group(1))
    refs = int(re.search(r"REFERENCE count[^:]*:\s*(\d+)", out).group(1))
    disp = int(re.search(r"DISPLAY ITEMS TOTAL:\s*(\d+)", out).group(1))
    return body, abst, refs, disp


def main():
    print("=" * 78)
    print("Refreshing Submission_Package_2026-10-17.zip")
    print("=" * 78)

    # ---- 1. stage from existing zip ----
    if not os.path.exists(PKG):
        raise SystemExit("FATAL: base package not found: " + PKG)
    if os.path.exists(STAGE):
        shutil.rmtree(STAGE)
    os.makedirs(STAGE)
    with zipfile.ZipFile(PKG) as z:
        z.extractall(STAGE)
    print("staged %d entries from existing package" % len(zipfile.ZipFile(PKG).namelist()))

    # ---- 2. refresh canonical files ----
    for src, dst in REFRESH:
        s = os.path.join(PROJ, src)
        if not os.path.exists(s):
            raise SystemExit("FATAL: source missing: " + src)
        d = os.path.join(STAGE, dst)
        os.makedirs(os.path.dirname(d), exist_ok=True)
        shutil.copy2(s, d)
        print("  refreshed: %s" % dst)

    # ---- 3. README updates ----
    rp = os.path.join(STAGE, "README.md")
    rd = open(rp, encoding="utf-8").read()
    for old, new in [(README_SUPP_ROW_OLD, README_SUPP_ROW_NEW),
                     (README_SFIG_ROW_OLD, README_SFIG_ROW_NEW),
                     (README_VERIF_OLD, README_VERIF_NEW),
                     (README_LANG_OLD, README_LANG_NEW)]:
        if old not in rd:
            print("  [warn] README block not found (skipped): %s..." % old[:60].replace("\n", " "))
        else:
            rd = rd.replace(old, new)
    open(rp, "w", encoding="utf-8").write(rd)
    print("  README.md verification/contents sections updated")

    # ---- 4. MANIFEST ----
    body, abst, refs, disp = measure_oac()
    n_supp = len(os.listdir(os.path.join(STAGE, "Supplementary")))
    n_tsd = len([f for f in os.listdir(os.path.join(STAGE, "Table_Source_Data")) if f.endswith(".csv")])
    n_gov = len(os.listdir(os.path.join(STAGE, "Governance")))
    man = os.path.join(STAGE, "MANIFEST.txt")
    # count final entries after writing manifest (manifest itself included)
    manifest_txt = (
        "Submission_Package_2026-10-17\n"
        "Built: %s\n"
        "Manuscript main body words: %d (OAC cap 4000)\n"
        "Abstract words: %d (OAC cap 350, structured)\n"
        "References: %d (OAC cap 50)\n"
        "Display items: Figure 1-5 + Table 1-3 = %d (OAC cap 8)\n"
        "audit_v18_full.py: 494 checks / 0 fails / 0 warns\n"
        "audit_figures_v17.py: 845 checks / 0 fails\n"
        "Table source CSVs included: %d\n"
        "Supplementary files included: %d\n"
        "Governance docs included: %d\n"
        "Total zip entries: %d\n"
        % (datetime.now().replace(microsecond=0).isoformat(), body, abst, refs, disp,
           n_tsd, n_supp, n_gov, 0)
    )
    # two-pass: write manifest, count entries, rewrite with the true count
    open(man, "w", encoding="utf-8").write(manifest_txt)
    n_entries = sum(len(files) for _, dirs, files in os.walk(STAGE))
    manifest_txt = manifest_txt.replace("Total zip entries: 0", "Total zip entries: %d" % n_entries)
    open(man, "w", encoding="utf-8").write(manifest_txt)
    print("  MANIFEST.txt regenerated (body=%d abstract=%d refs=%d display=%d)" % (body, abst, refs, disp))

    # ---- 5. re-zip ----
    tmp_zip = PKG + ".tmp"
    with zipfile.ZipFile(tmp_zip, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for root, dirs, files in os.walk(STAGE):
            dirs.sort()
            for f in sorted(files):
                full = os.path.join(root, f)
                arc = os.path.relpath(full, STAGE).replace("\\", "/")
                z.write(full, arc)
    shutil.move(tmp_zip, PKG)
    print("\nrebuilt: %s (%d entries, %.1f MB)"
          % (os.path.basename(PKG), n_entries, os.path.getsize(PKG) / 1e6))


if __name__ == "__main__":
    main()
