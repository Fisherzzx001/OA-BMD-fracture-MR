# -*- coding: utf-8 -*-
"""
311 · 将 docs/08 本轮 S23 图注末尾的 F10 跨祖先可移植性补强（以及 S11 F10 行 /
Honest boundaries 的同源修订）同步到派生稿与镜像。

派生稿：
  - submission_package/Manuscript.docx
  - submission_package/Submission_MainText.docx
  - Manuscript_Ready/01_Manuscript.md（镜像，直接 copy docs/08）
  - Manuscript_Ready/Supplementary/01_predictions_and_falsification.md（镜像）

做法：S23 末尾为纯文本图注段落，无 intra-paragraph 格式，故对 docx 采用
"整段文本替换 + 折叠为单 run" 的定点替换；镜像直接 shutil.copy2。
幂等：若 NEW 短语已存在则 skip，不报错。
"""
import os
import shutil

from docx import Document

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOC = os.path.join(ROOT, "docs", "08_投稿稿_OA_BMD_fracture.md")
S11 = os.path.join(ROOT, "docs", "Supplementary", "01_predictions_and_falsification.md")
MIRROR_MD = os.path.join(ROOT, "Manuscript_Ready", "01_Manuscript.md")
MIRROR_S11 = os.path.join(ROOT, "Manuscript_Ready", "Supplementary", "01_predictions_and_falsification.md")
DOCX = [
    os.path.join(ROOT, "submission_package", "Manuscript.docx"),
    os.path.join(ROOT, "submission_package", "Submission_MainText.docx"),
]

OLD_END = "could close the two remaining gaps."
NEW_TAIL = ("could close the two remaining gaps. The F10 gap is attenuated, not removed: "
            "the outcome GWAS [14] performed GWAS meta-analyses within four ancestry groups "
            "(East Asian, African, South Asian, Hispanic) and detected no ancestry-specific study-wide "
            "significant associations, indicating the osteoarthritis genetic architecture is not "
            "ancestry-specific at the GWAS level and supporting the cross-ancestry transportability of "
            "the European-instrument estimates, while direct non-European replication stays blocked by "
            "the absence of non-European BMD exposure summary statistics.")


def sync_mirror():
    shutil.copy2(DOC, MIRROR_MD)
    os.makedirs(os.path.dirname(MIRROR_S11), exist_ok=True)
    shutil.copy2(S11, MIRROR_S11)
    print("[mirror] copied docs/08 -> Manuscript_Ready/01_Manuscript.md")
    print("[mirror] copied S11 -> Manuscript_Ready/Supplementary/01_predictions_and_falsification.md")


def sync_docx(path):
    d = Document(path)
    hit = 0
    for para in d.paragraphs:
        full = "".join(r.text for r in para.runs)
        if "attenuated, not removed" in full:
            print(f"[skip] {path}: F10 note already present")
            d.save(path)
            return
        if OLD_END in full:
            new_full = full.replace(OLD_END, NEW_TAIL, 1)
            # 折叠为单 run，保留段落样式
            if para.runs:
                para.runs[0].text = new_full
                for r in para.runs[1:]:
                    r.text = ""
            else:
                para.add_run(new_full)
            hit += 1
    if hit == 0:
        print(f"[warn] {path}: OLD_END not found (no change)")
        return
    d.save(path)
    print(f"[ok] {path}: S23 F10 note inserted (replacements={hit})")


def main():
    # docx 先同步，再复制镜像（镜像为 md，不受 docx 影响）
    for f in DOCX:
        sync_docx(f)
    sync_mirror()
    print("DONE")


if __name__ == "__main__":
    main()
