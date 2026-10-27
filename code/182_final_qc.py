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
record("D. 理论段定位", "SOS 明确为待检验假设",
       "which we pose as a testable hypothesis rather than a demonstrated mechanism" in t, "")
record("D. 理论段定位", "F2 标为决定性检验且未执行",
       "its decisive test" in t and "remains unexecuted" in t, "")
record("D. 理论段定位", "Highlights 同步 testable hypothesis",
       "is a testable hypothesis consistent with data." in t, "")

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
        ok_bd = "being deposited" in txt
        ok_pd = "mirror is pending" not in txt
        record("E. 跨产物同步", f"{name}: OSF 措辞 being deposited", ok_bd, "")
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
MAIN = "Manuscript_Ready/Figures/Main"
SUPP = "Manuscript_Ready/Figures/Supplementary"
TAB = "Manuscript_Ready/Supplementary"

main_files = os.listdir(MAIN) if os.path.isdir(MAIN) else []
supp_files = os.listdir(SUPP) if os.path.isdir(SUPP) else []
tab_files = os.listdir(TAB) if os.path.isdir(TAB) else []

miss_main = [f"Figure{i}" for i in range(1, 11)
             if not any(f.startswith("Figure%d_" % i) for f in main_files)]
if not any(f.startswith("GraphicalAbstract") for f in main_files):
    miss_main.append("GraphicalAbstract")
record("G. 物料清单", "主图 Figure 1–10 + Graphical Abstract (11)",
       not miss_main, ("missing: " + str(miss_main)) if miss_main else "全部存在")

allfig = main_files + supp_files
miss_s = ["S%d" % i for i in range(1, 18)
          if not any(f.startswith("SupplFig%d_" % i) for f in allfig)]
record("G. 物料清单", "补充图 S1–S17 (17)", not miss_s,
       ("missing: " + str(miss_s)) if miss_s else "全部存在")

miss_t = ["Table_S%d" % i for i in (8, 9, 12, 13, 14, 15, 16)
          if not any(f.startswith("Table_S%d_" % i) for f in tab_files)]
record("G. 物料清单", "补充表 S8,S9,S12–S16 (7)", not miss_t,
       ("missing: " + str(miss_t)) if miss_t else "全部存在")

need = ["02_mvmr_height_grid", "01_predictions_and_falsification", "STROBE-MR_checklist"]
miss_o = [k for k in need if not any(k in f for f in tab_files)]
record("G. 物料清单", "Table S10/S11 + STROBE-MR 清单", not miss_o,
       ("missing: " + str(miss_o)) if miss_o else "全部存在")

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
