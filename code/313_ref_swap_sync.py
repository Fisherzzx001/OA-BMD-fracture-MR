# -*- coding: utf-8 -*-
"""
313 · 将 docs/08 本轮参考文献等价替换（[11] Zhu 2025 -> Yu 2025）同步到派生稿与镜像。

做法（定点 run 级替换，保留段落样式，幂等；用 python-docx 拼接 runs 做匹配/校验，
避免原始 XML 因 run 拆分导致的子串误判）：
  - submission_package/Manuscript.docx
  - submission_package/Submission_MainText.docx
      * 正文 Introduction 第 37 段：[11] 结论句由 Zhu 的位点关节依赖性改为
        Yu 2025 的 "total BMD 显著正向中介 knee-OA 风险 (mediation 3.59%)"。
      * 文献表第 11 条：Zhu 2025 (Clin Rheumatol) -> Yu 2025 (Aging Clin Exp Res)。
  - Manuscript_Ready/01_Manuscript.md：直接 shutil.copy2(docs/08)（镜像与真源一致）。

等价替换：rcount 仍为 50（未新增编号）；[7] Alhassan 保持不动，R3-3 门禁
（"[7, 8]" 共存）不受影响；已报数值零改动。
"""
import os
import shutil

from docx import Document

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOC = os.path.join(ROOT, "docs", "08_投稿稿_OA_BMD_fracture.md")
MIRROR_MD = os.path.join(ROOT, "Manuscript_Ready", "01_Manuscript.md")
DOCX = [
    os.path.join(ROOT, "submission_package", "Manuscript.docx"),
    os.path.join(ROOT, "submission_package", "Submission_MainText.docx"),
]

OLD_CLAUSE = ("and a 2025 two-sample MR of BMD against OA progression reports the same "
              "site-by-joint dependence, femoral-neck BMD raising hip-replacement risk "
              "(OR 1.27, 95% CI 1.13–1.43) but not knee-replacement risk [11].")
NEW_CLAUSE = ("and a 2025 two-sample MR reported total BMD as a significant positive "
              "mediator of knee-OA risk (mediation proportion 3.59%) [11].")

OLD_BIB = ('11. Zhu R, Xing X, Bian J, Zhang X, Ge L, Cai G. Causal association between '
           'bone mineral density and the risk of joint replacement in patients with '
           'osteoarthritis: a Mendelian randomization study. Clin Rheumatol. '
           '2025;44(2):823–830. doi:10.1007/s10067-024-07289-5. PMID 39738847.')
NEW_BIB = ('11. Yu H, Wang J, Xu X, Li H, Guo J. Revealing the mediating mechanisms '
           'between BMI and osteoarthritis: a Mendelian randomization and mediation '
           'analysis. Aging Clin Exp Res. 2025;37(1):119. '
           'doi:10.1007/s40520-025-03035-2.')


def _joined(para):
    return "".join(r.text for r in para.runs)


def _collapse(para, new_text):
    """Replace paragraph content with a single run carrying new_text; keep <w:pPr>."""
    p = para._p
    # drop all run-level elements (w:r, w:hyperlink, w:bookmark* handled by parent)
    for tag in ("w:r", "w:hyperlink", "w:smartTag"):
        for el in p.findall(_qn(tag)):
            p.remove(el)
    from docx.oxml.ns import qn
    new_r = p.makeelement(qn("w:r"), {})
    new_t = p.makeelement(qn("w:t"), {qn("xml:space"): "preserve"})
    new_t.text = new_text
    new_r.append(new_t)
    p.append(new_r)


def _qn(tag):
    from docx.oxml.ns import qn
    return qn(tag)


def sync_docx(path):
    d = Document(path)
    n_clause, n_bib = 0, 0
    for para in d.paragraphs:
        full = _joined(para)
        # clause
        if NEW_CLAUSE in full:
            pass  # 已同步
        elif OLD_CLAUSE in full:
            _collapse(para, full.replace(OLD_CLAUSE, NEW_CLAUSE, 1))
            n_clause += 1
        # bib
        if NEW_BIB in full:
            pass
        elif OLD_BIB in full:
            _collapse(para, full.replace(OLD_BIB, NEW_BIB, 1))
            n_bib += 1
    if n_clause == 0 and n_bib == 0:
        print(f"[skip] {os.path.basename(path)}: 已是 Yu 2025 状态")
    else:
        d.save(path)
        print(f"[ok] {os.path.basename(path)}: 正文[11]句 x{n_clause}, 文献[11] x{n_bib}")
    # 用 python-docx 重新打开校验（拼接 runs，避免 XML run 拆分误判）
    d2 = Document(path)
    joined_all = "\n".join(_joined(p) for p in d2.paragraphs)
    assert OLD_CLAUSE not in joined_all, f"FAIL {path}: OLD_CLAUSE 残留"
    assert OLD_BIB not in joined_all, f"FAIL {path}: OLD_BIB 残留(Zhu)"
    assert NEW_CLAUSE in joined_all, f"FAIL {path}: NEW_CLAUSE 缺失"
    assert NEW_BIB in joined_all, f"FAIL {path}: NEW_BIB 缺失(Yu)"
    print(f"[verify] {os.path.basename(path)}: 校验通过")


def sync_mirror():
    os.makedirs(os.path.dirname(MIRROR_MD), exist_ok=True)
    shutil.copy2(DOC, MIRROR_MD)
    print(f"[mirror] copied docs/08 -> {os.path.relpath(MIRROR_MD, ROOT)}")


def main():
    for f in DOCX:
        sync_docx(f)
    sync_mirror()
    print("DONE")


if __name__ == "__main__":
    main()
