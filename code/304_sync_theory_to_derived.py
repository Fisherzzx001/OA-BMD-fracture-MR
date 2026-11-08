# -*- coding: utf-8 -*-
"""
304 · 将 docs/08 本轮 L135（S23 图注）理论边界深化同步到派生稿

派生稿：
  - submission_package/Manuscript.docx
  - submission_package/Submission_MainText.docx
  - Manuscript_Ready/01_Manuscript.md（镜像）

做法：沿用 code/303 的 run 级定点替换（保留版式与元数据），并带
"恰好替换 1 处" 的防错守卫；镜像直接 shutil.copy2 docs/08。

幂等：重跑时 OLD 已不存在 → 命中 0 → 视为 [skip]，不报错、不改文件。
"""
import io
import os
import shutil
import sys

from docx import Document

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOC = os.path.join(ROOT, "docs", "08_投稿稿_OA_BMD_fracture.md")
MIRROR = os.path.join(ROOT, "Manuscript_Ready", "01_Manuscript.md")
DOCX = [
    os.path.join(ROOT, "submission_package", "Manuscript.docx"),
    os.path.join(ROOT, "submission_package", "Submission_MainText.docx"),
]

OLD = ("Theory boundary. The discriminating kernel is not a label for the observed pattern: "
       "one latent set-point, read through a shared set of inherited variants, must fix opposite signs "
       "at whole-bone and joint endpoints, with sign set by the set-point rather than by endpoint-specific "
       "pathways. The descriptive bone–cartilage paradox [3, 4] and the phenotype-level bone-forming "
       "triad [39] are rivals because neither requires one variant set nor sign opposition. The kernel is falsifiable:")

NEW = ("Theory boundary. The discriminating kernel is not a label for the observed pattern; it is "
       "non-trivial because BMD, OA and fracture occupy distinct genetic dimensions (S22: the PC1 "
       "bone-mass/fracture axis is orthogonal to the PC2 OA axis), so a single set-point cannot be "
       "reduced to a generic bone–OA correlation. One latent set-point, read through a shared set of "
       "inherited variants, must fix opposite signs at whole-bone and joint endpoints, with sign set by "
       "the set-point rather than by endpoint-specific pathways. The descriptive bone–cartilage paradox "
       "[3, 4] and the phenotype-level bone-forming triad [39] are rivals because neither requires one "
       "variant set nor sign opposition; what SOS adds is the single-set-point discriminator — one "
       "variant set fixing sign-opposing causal effects — rather than a restatement of the descriptive "
       "pattern. The kernel is falsifiable:")


def docx_text(doc):
    parts = [p.text for p in doc.paragraphs]
    for tbl in doc.tables:
        for row in tbl.rows:
            for cell in row.cells:
                parts += [p.text for p in cell.paragraphs]
    return "\n".join(parts)


def replace_in_docx(path):
    doc = Document(path)
    pre = docx_text(doc)
    if OLD not in pre:
        print(f"  [skip] {os.path.basename(path)}: OLD 已不存在（已同步）")
        return 0
    n = 0
    for par in doc.paragraphs:
        for run in par.runs:
            if OLD in run.text:
                run.text = run.text.replace(OLD, NEW, 1)
                n += 1
    for tbl in doc.tables:
        for row in tbl.rows:
            for cell in row.cells:
                for par in cell.paragraphs:
                    for run in par.runs:
                        if OLD in run.text:
                            run.text = run.text.replace(OLD, NEW, 1)
                            n += 1
    if n != 1:
        raise SystemExit(f"ABORT: {os.path.basename(path)} 替换计数={n}（应为 1），文件未改动")
    doc.save(path)
    post = docx_text(Document(path))
    assert NEW in post and OLD not in post, "post-condition failed"
    print(f"  [ok] {os.path.basename(path)}: 1 处替换已写入")
    return n


def main():
    print("304 · 同步 L135 理论边界深化到派生稿")
    print("-" * 60)
    for p in DOCX:
        replace_in_docx(p)
    # 镜像
    if not os.path.isfile(DOC):
        raise SystemExit("docs/08 缺失")
    shutil.copy2(DOC, MIRROR)
    print(f"  [ok] mirror: {os.path.relpath(MIRROR, ROOT)} <- docs/08")
    print("DONE.")


if __name__ == "__main__":
    main()
