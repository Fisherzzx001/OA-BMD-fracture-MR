#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
315_sync_f10_s32.py -- sync all derived submission artefacts after the F10 S32 addition.

docs/08_投稿稿_OA_BMD_fracture.md is the single source of truth.  This script:
  1. copies docs/08 -> Manuscript_Ready/01_Manuscript.md (byte-identical mirror)
  2. rebuilds submission_package/Manuscript.docx and Submission_MainText.docx
     from docs/08 via the repository's own builders (build_submission.py /
     build_maintext.py); a full rebuild is used instead of run-level patching
     because S32 adds a whole new supplementary section.
  3. copies the standalone Table_S32_f10_ancestry_proxy.md into
     Manuscript_Ready/Supplementary/
  4. verifies that every derived artefact contains the F10/S32 markers and that
     the mirror is byte-identical to docs/08.

Idempotent.  No reported beta/SE/p/OR/CI is touched.
"""
import os
import shutil
import subprocess
import sys
from docx import Document

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "docs", "08_投稿稿_OA_BMD_fracture.md")
MIRROR = os.path.join(ROOT, "Manuscript_Ready", "01_Manuscript.md")
SUPP_SRC = os.path.join(ROOT, "docs", "Supplementary", "Table_S32_f10_ancestry_proxy.md")
SUPP_DST = os.path.join(ROOT, "Manuscript_Ready", "Supplementary",
                        "Table_S32_f10_ancestry_proxy.md")
DOCX = [os.path.join(ROOT, "submission_package", "Manuscript.docx"),
        os.path.join(ROOT, "submission_package", "Submission_MainText.docx")]

# markers that must be present in every derived artefact after the F10 update
MARKERS = [
    "Supplementary Table S32",
    "ancestry-composition proxy",
    "Q_between=0.29",
    "registered cross-ancestry replication unexecuted",
]


def main():
    os.chdir(ROOT)

    # 1. mirror (byte-identical copy)
    shutil.copyfile(SRC, MIRROR)
    print("[1/4] mirror <- docs/08")

    # 2. rebuild both docx from docs/08 using the repository builders
    for script in ("build_submission.py", "build_maintext.py"):
        r = subprocess.run([sys.executable, script], capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"{script} failed:\n{r.stdout}\n{r.stderr}")
        print(f"[2/4] rebuilt via {script}")

    # 3. standalone supplementary table
    os.makedirs(os.path.dirname(SUPP_DST), exist_ok=True)
    shutil.copyfile(SUPP_SRC, SUPP_DST)
    print("[3/4] standalone Table S32 copied into Manuscript_Ready/Supplementary/")

    # 4. verification
    a = open(SRC, "rb").read()
    b = open(MIRROR, "rb").read()
    assert a == b, "mirror is not byte-identical to docs/08"
    print("[4/4] mirror byte-identical to docs/08: OK")

    for path in DOCX:
        text = "\n".join(p.text for p in Document(path).paragraphs)
        missing = [m for m in MARKERS if m not in text]
        assert not missing, f"{os.path.basename(path)} missing markers: {missing}"
        print(f"      {os.path.basename(path)}: {len(MARKERS)}/{len(MARKERS)} F10/S32 markers OK")

    supp = open(SUPP_DST, encoding="utf-8").read()
    assert "Table S32" in supp and "1.0812" in supp
    print("\nSYNC COMPLETE -- all derived artefacts in sync with docs/08")


if __name__ == "__main__":
    main()
