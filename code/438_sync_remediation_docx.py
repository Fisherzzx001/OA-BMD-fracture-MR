# -*- coding: utf-8 -*-
"""
438 · 将 docs/08 本轮 M1–M6 整改编辑同步到两份投递 docx

docs/08 是唯一真源（已通过 148/149b/334）。两份 docx 是增量派生稿，
由 python-docx 定点替换维护。本脚本把本轮对 docs/08 的 8 处编辑（L11/L23/L27/
L55/L65/L69/L99/ref[50]）原样移植到：
  - submission_package/Manuscript.docx
  - submission_package/Submission_MainText.docx

替换策略：对每个 (OLD, NEW) 在段落与表格单元格内做 run 级合并替换
（拼接该段落全部 run 文本 → 替换 → 写回首个 run、清空其余，保留 run[0] 格式），
对 L69 做整段替换（锚定句首，整段覆盖为 docs/08 原文）。

幂等：重跑时 OLD 已不存在 → 命中 0 → 视为 [skip]，不报错、不改文件。
前置：需 python-docx（docx）。
"""
import os
import sys

from docx import Document

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOCX = [
    os.path.join(ROOT, "submission_package", "Manuscript.docx"),
    os.path.join(ROOT, "submission_package", "Submission_MainText.docx"),
]

# (OLD, NEW) —— OLD 取 docx 当前（v2.4）原文，NEW 取 docs/08 整改后原文
REPL = [
    # L11 Highlights
    ("Higher genetically predicted BMD associates with higher knee OA, lower fracture.",
     "Genetically predicted higher BMD raises knee-OA risk (2/5 sites) and lowers fracture."),
    # L23 Abstract Objective
    ("tested the hypothesis that one BMD-raising set-point generates",
     "tested whether one BMD-raising set-point (a hypothesis, not a demonstrated mechanism) generates"),
    # L27 Abstract Results
    ("raised knee-OA risk at all five analysed sites",
     "raised knee-OA risk directionally at all five analysed sites"),
    # L55 Sample overlap — MRlap 透明度
    ("applied MRlap [21] (corrects weak-instrument bias and winner's curse while accounting for overlap) to all ten",
     "applied MRlap [21] (corrects weak-instrument bias and winner's curse; overlap is accommodated through its "
     "internally estimated cross-trait LDSC intercept, with no overlap fraction supplied — see §S14) to all ten"),
    # L65 EUR-only r_g 措辞
    ("preserves all five positive r_g (Table S33)",
     "preserves the positive r_g at all five sites (Table S33)"),
    # L99 Limitations — 多重性控制
    ("so the colocalisation 9/31 and random-effects 2/5 findings are descriptive rather than multiplicity-adjusted",
     "so the colocalisation 9/31 and random-effects 2/5 findings are exploratory (hypothesis-generating), "
     "not multiplicity-adjusted"),
    # ref [50] — 机制生物学对话（R3）
    ("50. Dequeker J, Mokassa L, Aerssens J. Bone density and osteoarthritis. "
     "J Rheumatol Suppl. 1995;43:98–100. PMID 7752151",
     "50. Karsdal MA, Bay-Jensen AC, Lories RJ, Spector T, Pastoureau P, Christiansen C, Henriksen K, et al. "
     "The coupling of bone and cartilage turnover in osteoarthritis: opportunities for bone antiresorptives "
     "and anabolics as potential treatments? Ann Rheum Dis . 2014;73(2):336–348. "
     "doi:10.1136/annrheumdis-2013-204111. PMID 24285494"),
]

# L69 整段替换（锚定句首，整段覆盖为 docs/08 原文）
L69_ANCHOR = "MR supported a positive effect of BMD on knee OA across all five exposures"
L69_NEW = (
    "MR supported a positive effect of BMD on knee OA across all five exposures (Table 2): "
    "total-body β=+0.078 (p=3.9×10⁻¹⁰), lumbar-spine +0.095 (p=7.7×10⁻⁹), arms +0.075 (p=1.7×10⁻⁶), "
    "legs +0.060 (p=2.3×10⁻⁵), femoral-neck +0.058 (p=6.4×10⁻⁴) — OR 1.060–1.099 per SD higher BMD, "
    "against OR 1.40 (1.20–1.63) per SD in the largest one-sample estimate [12], a 3.6–5.8-fold log-odds "
    "difference (adjudicated §4). MR-Egger intercepts [19] were all p>0.24 and weighted medians matched the "
    "IVW sign at all five sites, with bootstrap 95% CIs excluding the null at three of five. "
    "MR-Egger non-significance (all p>0.24) is underpowered for pleiotropy, so it is read as not detected "
    "rather than absent; the positive control is itself only partly calibrated (femoral-neck MR-Egger p=0.018), "
    "so pleiotropy-robust readings are contingent on a control that carries directional pleiotropy. "
    "Leave-one-out produced no sign reversals. Heterogeneity was high (I²=68–86%); under random-effects IVW "
    "only two of five exposures remained significant (total-body p=0.006, lumbar-spine p=0.017), "
    "so the count depends on the random-effects convention (Table S17). "
    "A EUR-only knee-OA MR (FinnGen) replicated the positive direction at all five sites (OR/SD 1.06–1.12; S35). "
    "All five effects remained positively causal after MR-PRESSO outlier removal (Table S17). "
    "The positive control was comparably heterogeneous (I² 78–83%), so high I² here reflects the BMD instruments."
)


def merge_replace(container_runs_getter, old, new):
    """替换一个容器内（段落或单元格）的 run 文本。命中返回 1，未命中返回 0。"""
    runs = container_runs_getter()
    full = "".join(r.text for r in runs)
    if old not in full:
        return 0
    newfull = full.replace(old, new)
    if not runs:
        return 0
    runs[0].text = newfull
    for r in runs[1:]:
        r.text = ""
    return 1


def process(path):
    doc = Document(path)
    print(f"\n== {os.path.relpath(path, ROOT)} ==")
    total = 0
    for old, new in REPL:
        n = 0
        for par in doc.paragraphs:
            n += merge_replace(lambda p=par: p.runs, old, new)
        for tbl in doc.tables:
            for row in tbl.rows:
                for cell in row.cells:
                    for par in cell.paragraphs:
                        n += merge_replace(lambda p=par: p.runs, old, new)
        state = "REPLACED" if n else "skip (not found)"
        print(f"  [{state}] {n}  {old[:54]!r}")
        total += n

    # L69 整段替换
    n69 = 0
    for par in doc.paragraphs:
        if L69_ANCHOR in "".join(r.text for r in par.runs):
            if par.runs:
                par.runs[0].text = L69_NEW
                for r in par.runs[1:]:
                    r.text = ""
                n69 += 1
    print(f"  [{'REPLACED' if n69 else 'skip'}] {n69}  L69 full-paragraph swap")
    total += n69

    doc.save(path)
    print(f"  saved ({total} edits total)")
    return total


def main():
    for p in DOCX:
        if not os.path.isfile(p):
            print(f"MISSING {p}", file=sys.stderr)
            sys.exit(1)
        process(p)
    print("\nDONE — verify with code/310 and a docx text diff.")


if __name__ == "__main__":
    main()
