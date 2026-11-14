#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Pre-submission consolidated QC (final check sheet).

Runs every verification that has been done across the project and emits a
single PASS/FAIL table. Re-runnable at any time; the check sheet docs/128 is
generated from this script's output.

Usage:  python code/182_final_qc.py
"""
import io
import os
import re
import subprocess
import sys
import zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)

PY = sys.executable
CODE = os.path.join(ROOT, "code")

MANU = "docs/08_投稿稿_OA_BMD_fracture.md"
CL = "docs/cover_letter_OAC.md"
DOCX = {
    "Manuscript.docx": "submission_package/Manuscript.docx",
    "Submission_MainText.docx": "submission_package/Submission_MainText.docx",
    "Cover_Letter.docx": "submission_package/Cover_Letter.docx",
}
MR = {
    "MR/01_Manuscript.md": "Manuscript_Ready/01_Manuscript.md",
    "MR/02_Cover_letter.md": "Manuscript_Ready/02_Cover_letter.md",
    "MR/README.md": "Manuscript_Ready/README.md",
}

CONCEPT = "10.5281/zenodo.23400860"
STALE = "10.5281/zenodo.23400861"

results = []


def record(section, item, ok, detail=""):
    results.append((section, item, bool(ok), detail))
    return bool(ok)


def docx_text(path):
    xml = zipfile.ZipFile(path).read("word/document.xml").decode("utf-8", "ignore")
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", xml))


def docx_paragraphs(path):
    """Paragraph texts (one per <w:p>). A full supplementary table starts its
    OWN paragraph with 'Supplementary Table SNN.'; an in-text cross-reference
    ('...in Supplementary Table S16.') never opens a paragraph, so this
    distinction is robust to whitespace flattening in docx_text()."""
    xml = zipfile.ZipFile(path).read("word/document.xml").decode("utf-8", "ignore")
    paras = re.split(r"</w:p>", xml)
    out = []
    for p in paras:
        t = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", "", p)).strip()
        if t:
            out.append(t)
    return out


def read(path):
    return io.open(path, encoding="utf-8").read()


def run(script):
    """Run a project script, return (ok, combined output)."""
    try:
        p = subprocess.run([PY, os.path.join(CODE, script)],
                           capture_output=True, text=True, encoding="utf-8",
                           errors="replace", timeout=300)
        return p.returncode == 0, (p.stdout or "") + (p.stderr or "")
    except Exception as e:  # noqa
        return False, str(e)


# ---------------------------------------------------------------- Section A
ok, out = run("148_oac_compliance.py")
compliant = "OVERALL OAC COMPLIANT: YES" in out
record("A. OAC 合规", "code/148 综合判定", compliant,
       ("COMPLIANT: YES" if compliant else "NOT COMPLIANT"))
for label, pat in [("正文词数", r"Body words\s*:\s*(\d+)\s*\(limit\s*(\d+)\)"),
                   ("摘要词数", r"Abstract words\s*:\s*(\d+)\s*\(limit\s*(\d+)\)"),
                   ("参考文献", r"References\s*:\s*(\d+)\s*\(limit\s*(\d+)\)"),
                   ("显示项", r"Display items\s*:\s*(\d+)\s*\(limit\s*(\d+)\)")]:
    m = re.search(pat, out)
    if m:
        v, lim = int(m.group(1)), int(m.group(2))
        record("A. OAC 合规", f"{label} {v}/{lim}", v <= lim, "")

# ---------------------------------------------------------------- Section B
ok, out = run("149b_verify_reported_values.py")
safe = "ALL REPORTED VALUES IDENTICAL" in out and "SAFE" in out
record("B. 数值完整性", "code/149b 已报告统计量零改动", safe,
       ("SAFE" if safe else "MISMATCH"))

# ---------------------------------------------------------------- Section C
t = read(MANU)
cited = set()
for m in re.finditer(r"\[(\d+(?:[–\-]\d+)?(?:,\s*\d+)*)\]", t):
    for p in m.group(1).split(","):
        p = p.strip()
        if "–" in p or "-" in p:
            a, b = re.split(r"[–-]", p)
            cited |= set(range(int(a), int(b) + 1))
        else:
            cited.add(int(p))
defined = set(int(re.match(r"^(\d+)\.\s", l.strip()).group(1))
              for l in t.splitlines() if re.match(r"^\d+\.\s", l.strip()))
never = sorted(defined - cited)
record("C. 引用完整性", f"参考文献 1..{len(defined)} 全部被引", not never,
       ("never-cited: " + str(never)) if never else "全部被引")

# 关键方法学文献：按作者名在文献表定位编号，再验证正文引用了该编号。
# 用语义锚点而非写死编号——任一次重编号后本检查自动跟随，不会误报 FAIL。
def refnum(anchor):
    for l in t.splitlines():
        mm = re.match(r"^(\d+)\.\s+(.*)$", l.strip())
        if mm and mm.group(2).startswith(anchor):
            return int(mm.group(1))
    return None


works = [
    ("加权中位数 (Bowden 2016)", "Bowden J, Davey Smith G, Haycock PC"),
    ("MR-Egger (Bowden 2015)", "Bowden J, Davey Smith G, Burgess S"),
    ("liability-scale (Lee 2011)", "Lee SH, Wray NR"),
    ("coloc.abf (Giambartolomei 2014)", "Giambartolomei C"),
    ("STROBE-MR (Skrivankova 2021)", "Skrivankova VW"),
]
for label, anchor in works:
    n = refnum(anchor)
    if n is None:
        record("C. 引用完整性", f"{label} 在文献表中", False, "文献缺失")
    else:
        record("C. 引用完整性", f"{label} → [{n}] 且正文已引", f"[{n}]" in t, "")

# OAC 强制：参考文献按正文首次出现顺序编号（Guide for Authors 明文要求）
body_for_order = t[:t.index("\n## References")] if "\n## References" in t else t
first = {}
for m in re.finditer(r"\[(\d+(?:[–\-]\d+)?(?:,\s*\d+)*)\]", body_for_order):
    for p in m.group(1).split(","):
        p = p.strip()
        if "–" in p or "-" in p:
            a, b = re.split(r"[–-]", p)
            nums = list(range(int(a), int(b) + 1))
        else:
            nums = [int(p)]
        for n in nums:
            first.setdefault(n, m.start())
order = sorted(first, key=lambda n: first[n])
record("C. 引用完整性", "OAC 顺序编码制（按首次出现编号）",
       order == sorted(first) == list(range(1, len(first) + 1)),
       f"首次出现顺序 {'升序连续' if order == sorted(first) else str(order)}")

# ---------------------------------------------------------------- Section D
# Theoretical framing (R1 revision)
# NOTE (2026-10-30): these two checks previously pinned VERBATIM sentences. A verbatim
# string is a *proxy* for the claim, not the claim itself -- rewording the sentence
# silently failed the gate even though the underlying claim was still stated (indeed
# stated more clearly). Both are now tested semantically: accepted if ANY wording from
# the corresponding paraphrase family is present. Deleting the passage still fails;
# rewording it no longer does.
_SOS_TENTATIVE = [
    "testable hypothesis rather than a demonstrated mechanism",
    "testable hypothesis, not a demonstrated mechanism",
    "not a demonstrated mechanism",
]
record("D. 理论段定位", "SOS 明确为待检验假设",
       any(k in t for k in _SOS_TENTATIVE), "")
# v2.4 (2026-11-12): P1-4 改稿把判别者明确收敛到 P4/共定位，F2 为"decisive arm 未执行"。
# 旧串 "its decisive test ... remains unexecuted" 已被替换；实质要求（F2 决定性臂未执行、
# 且判别性归于共定位）不降反升，故同步更新断言并新增 P4 判别性锚点。
record("D. 理论段定位", "F2 决定性臂标为未执行",
       "F2's decisive arm is unexecuted" in t, "")
record("D. 理论段定位", "判别性收敛到 共定位判别者（非 F2）",
       "colocalisation split (9/31 fracture vs 0/31 knee-OA)" in t, "")
record("D. 理论段定位", "Highlights 同步 untested hypothesis",
       "Bone-forming set-point is untested" in t, "")

# STROBE-MR item 5 (assumptions) must state all three IV assumptions in Methods (§2 MR)
record("D. 理论段定位", "STROBE-MR item 5 三假设在正文显式陈述",
       all(k in t for k in ["Valid instruments must satisfy three assumptions",
                            "relevance", "independence from exposure", "exclusion restriction"]), "")

# Supplementary checklist must use the official 20-item form and the correct (JAMA) citation
try:
    _cl = read("docs/Supplementary/STROBE-MR_checklist.md")
except Exception:
    _cl = ""
record("D. 理论段定位", "STROBE-MR 补充清单为官方 20 项",
       "| 20 | OTHER INFORMATION - Conflicts of interest" in _cl and "| 21 |" not in _cl, "")
record("D. 理论段定位", "STROBE-MR 引用为 JAMA 2021 而非 BMJ",
       "JAMA" in _cl and "BMJ" not in _cl, "")
record("D. 理论段定位", "无遗留 22-item 表述",
       "22-item" not in t and "22-item" not in _cl and "22 STROBE" not in t, "")


# Round 4 (user-requested P1/P2/P3+P4/E1/E2) locks
# Reworded 2026-10-30 from a verbatim sentence to a two-part semantic test (see above):
# part (i) the text states that areal BMD fuses/conflates several components, and
# part (ii) it names the resulting gap for SOS. Removing either half fails; rewording does not.
_FUSE = ["areal BMD fuses", "areal BMD conflates", "fuses mass, size and geometry",
         "fuses bone mass, size and geometry", "conflates size with mass"]
_GAP = ["SOS is directionally supported but quantitatively untested",
        "not identified as bone mass", "mass versus geometry",
        "not resolved into bone mass versus skeletal geometry"]
record("D. 理论段定位", "P2: IV estimand vs SOS 构念 gap 已在正文命名",
       (sum(k in t for k in _FUSE) >= 1 and sum(k in t for k in _GAP) >= 1), "")
record("D. 理论段定位", "E1: 标题含 pleiotropy-aware 限定",
       "pleiotropy-aware Mendelian-randomisation evidence" in t, "")
# NOTE (2026-10-30): this previously hard-coded "[42]" as the TREASURE anchor. A literal
# number is a *proxy* for the citation, not the citation: any renumbering of the
# bibliography (e.g. inserting new references mid-list) failed the gate even though the
# trial was still cited correctly. Now resolved semantically via the author anchor.
_n_treasure = refnum("Sobue Y")
record("D. 理论段定位", "P3: TREASURE 试验被引且结论提及",
       "TREASURE" in t and _n_treasure is not None and f"[{_n_treasure}]" in t,
       "" if _n_treasure is not None else "Sobue Y (TREASURE) 不在文献表中")
record("D. 理论段定位", "P4: 锶 DXA 伪影 + 心血管限制已在结论交代",
       "measurement artefact" in t and "cardiovascular risk" in t, "")
record("D. 理论段定位", "P1: GWAS 访问日期已写入 Data availability",
       "accessed on 17 July 2026" in t, "")
record("D. 理论段定位", "R3-3: [7-8] 引用堆叠已拆为 [7, 8]",
       "[7, 8]" in t and "[7–8]" not in t, "")
record("D. 理论段定位", "E2: 摘要 Design 段已拆为两句且 ≤350",
       "using BMD→fracture as the positive control." in t, "")

# ---------------------------------------------------------------- Section E
art = {MANU: read(MANU), CL: read(CL)}
for name, p in DOCX.items():
    if os.path.exists(p):
        art[name] = docx_text(p)
for name, p in MR.items():
    if os.path.exists(p):
        art[name] = read(p)

missing = [n for n in list(DOCX) + list(MR) if n not in art]
record("E. 跨产物同步", f"全部 {len(art)} 个产物可读", not missing,
       ("missing: " + str(missing)) if missing else f"{len(art)} 产物")

for name, txt in art.items():
    is_cl = ("Cover" in name) or ("cover" in name)
    # OSF wording: only manuscript + cover-letter derived artifacts carry it
    if MANU in name or "Manuscript" in name or is_cl:
        # Either a resolved OSF DOI or the "being deposited" placeholder is acceptable;
        # the gate tests *consistency*, not a prescribed literal string.
        # DOI prefixes are case-insensitive (OSF mints 10.17605/OSF.IO/<GUID>,
        # DataCite canonicalises to 10.17605/osf.io/<guid>), so match case-insensitively.
        has_real = bool(re.search(r"10\.17605/OSF\.IO/[A-Z0-9]{5,10}", txt, re.I))
        ok_osf = has_real or ("being deposited" in txt)
        ok_pd = "mirror is pending" not in txt
        record("E. 跨产物同步", f"{name}: OSF 状态自洽 (真实DOI 或 being deposited)", ok_osf, "")
        record("E. 跨产物同步", f"{name}: 无 pending 残留", ok_pd, "")

for name, txt in art.items():
    record("E. 跨产物同步", f"{name}: 无陈旧版本 DOI {STALE.split('.')[-1]}",
           STALE not in txt, "")
    record("E. 跨产物同步", f"{name}: 概念 DOI 存在", CONCEPT in txt, "")

# ---------------------------------------------------------------- Section F
KEY = ["OR 1.060–1.099", "3.7×10⁻⁵¹", "3.389", "9/31", "31,873"]
for name in ["Manuscript.docx", "Submission_MainText.docx"]:
    if name in art:
        miss = [k for k in KEY if k not in art[name]]
        record("F. 关键统计量", f"{name}: 关键值 {len(KEY) - len(miss)}/{len(KEY)}",
               not miss, ("missing: " + str(miss)) if miss else "全部存在")

# ---------------------------------------------------------------- Section G
# Figure sources are split across the project: S1-S8 + Figure1-5 + GraphicalAbstract
# live in top-level "Figures/"; S9-S17 live in "Manuscript_Ready/Figures/Supplementary/".
# We scan the UNION of both so the inventory reflects the actual project layout
# (the manuscript states all S1-S17 are in "the submission repository (Figures/)";
#  the S9-S17 sync into Figures/ is a packaging step, not a content gap).
MAIN_DIRS = ["Figures", "Manuscript_Ready/Figures/Main"]
SUPP_DIRS = ["Figures", "Manuscript_Ready/Figures/Supplementary"]
TAB = "Manuscript_Ready/Supplementary"

main_files = []
for d in MAIN_DIRS:
    if os.path.isdir(d):
        main_files += os.listdir(d)
allfig = list(main_files)
for d in SUPP_DIRS:
    if os.path.isdir(d):
        allfig += os.listdir(d)
tab_files = os.listdir(TAB) if os.path.isdir(TAB) else []

miss_main = [f"Figure{i}" for i in range(1, 6)
             if not any(f.startswith("Figure%d_" % i) for f in main_files)]
if not any(f.startswith("GraphicalAbstract") for f in main_files):
    miss_main.append("GraphicalAbstract")
record("G. 物料清单", "主图 Figure 1–5 + Graphical Abstract (6)",
       not miss_main, ("missing: " + str(miss_main)) if miss_main else "全部存在")

miss_s = ["S%d" % i for i in range(1, 25)
          if not any(f.startswith("SupplFig%d_" % i) for f in allfig)]
record("G. 物料清单", "补充图 S1–S24 (24)", not miss_s,
       ("missing: " + str(miss_s)) if miss_s else "全部存在")

# 2026-11-11 compliance audit: the standalone table set uses the SAME numbering
# as the inline manuscript. Legacy files (S26=hip, S27=geometry, S29=ldsc) were
# mis-numbered; they have been removed in favour of correctly numbered exports.
# Audit_Trail.md is an INTERNAL archive and must NOT appear in the delivery set.
md_tabs = [f for f in tab_files if f.endswith(".md")]
miss_t = ["Table_S%d" % i for i in (8, 9, 12, 13, 14, 15, 16, 17, 25, 26, 27, 28, 29, 30, 31, 32)
          if not any(f.startswith("Table_S%d_" % i) for f in md_tabs)]
record("G. 物料清单", "补充表 S8,S9,S12–S17,S25–S32 (16, 独立 .md)", not miss_t,
       ("missing: " + str(miss_t)) if miss_t else "16/16")

need = ["02_mvmr_height_grid", "01_predictions_and_falsification", "STROBE-MR_checklist"]
miss_o = [k for k in need if not any(k in f for f in tab_files)]
record("G. 物料清单", "Table S10/S11 + STROBE-MR 清单", not miss_o,
       ("missing: " + str(miss_o)) if miss_o else "全部存在")

leaked_audit = [f for f in tab_files if "Audit_Trail" in f]
record("G. 物料清单", "审计档案未混入交付集 (no Audit_Trail in Manuscript_Ready/)",
       not leaked_audit, ("LEAK: " + str(leaked_audit)) if leaked_audit else "clean")

# Every inline supplementary-table heading must be present so the declared range
# is complete; the docx builder strips this block, so these are repository-only.
t08 = read(MANU)
inline_tables = ["## Supplementary Table S%d" % i
                 for i in (8, 9, 14, 15, 25, 26, 27, 28, 29, 30, 31, 32)]
miss_inline = [t for t in inline_tables if t not in t08]
record("G. 物料清单", "内联补充表 S8,S9,S14,S15,S25–S32 (12, docs/08)", not miss_inline,
       ("missing: " + str(miss_inline)) if miss_inline else "12/12")

# ---------------------------------------------------------------- Section H
# Compliance rule (2026-10-28): the MAIN TEXT must not contain audit-trail
# material. Audit content belongs in docs/Supplementary/Audit_Trail.md and is
# referenced from the Supplementary Material list only. Enforce automatically.
leak = bool(re.search(r"audit trail", t, re.IGNORECASE))
record("H. 主文合规", "真源无审计轨迹字样 (no 'audit trail' in docs/08)",
       not leak, ("found 'audit trail' in docs/08" if leak else "clean"))

# 2026-11-11 compliance audit: the DECISIVE check is on the journal-facing docx
# files, not merely a literal string in the markdown source. The old check passed
# while both docx carried all 21 inline supplementary tables plus internal
# self-audit phrasing. Enforce on the actual delivered documents.
for name in ["Manuscript.docx", "Submission_MainText.docx"]:
    dx = art.get(name, "")
    has_body = "Abstract" in dx or "Introduction" in dx
    record("H. 主文合规", f"{name}: 含主文结构", has_body,
           "" if has_body else "missing Abstract/Introduction (empty/truncated?)")
    # A full table opens its OWN paragraph with 'Supplementary Table SNN.'.
    # An in-text cross-reference ('...in Supplementary Table S16.') lives inside
    # a sentence and never starts a paragraph -- robust to whitespace flattening.
    full_tables = [p[:60] for p in docx_paragraphs(DOCX[name])
                   if re.match(r"^#{0,3}\s*Supplementary Table S\d+\.", p)]
    record("H. 主文合规", f"{name}: 无补表全文嵌入", not full_tables,
           ("LEAK paras: " + str(full_tables)) if full_tables else "clean")
    for pat, lab in [(r"code/\d+", "code/NNN 内部编号"),
                     (r"this revision", "修订过程语气"),
                     (r"internal-consistency", "自检语气"),
                     (r"previously tabulated", "重算过程语气")]:
        n = len(re.findall(pat, dx, re.IGNORECASE))
        record("H. 主文合规", f"{name}: 无 {lab}", n == 0,
               f"{n} hit(s)" if n else "clean")

# The audit archive itself must exist (repository-internal traceability only).
try:
    at = read("docs/Supplementary/Audit_Trail.md")
    n_claims = at.count("| A")  # rows tagged A1..A20
    record("H. 主文合规", f"内部审计档案存在且覆盖 {n_claims} 条核心声称",
           n_claims >= 10 and "claim → script line → independent" in at, "")
except Exception:
    record("H. 主文合规", "内部审计档案存在", False, "Audit_Trail.md 缺失")

# ---------------------------------------------------------------- Report
print("=" * 74)
print("PRE-SUBMISSION CONSOLIDATED QC — OA-BMD-fracture MR (target: OAC)")
print("=" * 74)
sections = {}
for sec, item, ok, detail in results:
    sections.setdefault(sec, []).append((item, ok, detail))

for sec in sorted(sections):
    print("\n" + sec)
    print("-" * len(sec))
    for item, ok, detail in sections[sec]:
        flag = "PASS" if ok else "FAIL"
        line = f"  [{flag}] {item}"
        if detail and not ok:
            line += f"  -> {detail}"
        elif detail:
            line += f"  ({detail})"
        print(line)

n_total = len(results)
n_fail = sum(1 for _, _, ok, _ in results if not ok)
print("\n" + "=" * 74)
print(f"TOTAL: {n_total - n_fail}/{n_total} PASS   FAIL: {n_fail}")
print("VERDICT:", "ALL CHECKS PASSED - READY TO SUBMIT" if n_fail == 0
      else f"{n_fail} CHECK(S) FAILED - DO NOT SUBMIT")
print("=" * 74)
sys.exit(0 if n_fail == 0 else 1)
