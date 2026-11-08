# -*- coding: utf-8 -*-
"""
310 · 防回归综合断言（no-regression guard）

在每一次整改后运行，确保：
  (A) 三连门禁全绿：code/148 OAC COMPLIANT / code/149b SAFE / code/182 全 PASS
  (B) docs/08 本轮理论边界深化已落地（S22 正交维度前提 + 单一设定点判别式 novelty）
  (C) 派生稿同步：Manuscript_Ready 镜像 + 两份 docx 均含新理论串
  (D) 已报数值零改动：5 个 KEY 串在 docs/08 与两份 docx 中均存在
  (E) 物料清单完整：S1–S24 图文件、独立补充表文件、内联 S25/S28/S30/S31

退出码 0 = 全部通过；1 = 任一断言失败。
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
MIRROR = "Manuscript_Ready/01_Manuscript.md"
DOCX = {
    "Manuscript.docx": "submission_package/Manuscript.docx",
    "Submission_MainText.docx": "submission_package/Submission_MainText.docx",
}
KEY = ["OR 1.060–1.099", "3.7×10⁻⁵¹", "3.389", "9/31", "31,873"]

results = []


def record(item, ok, detail=""):
    results.append((item, bool(ok), detail))
    print(f"  [{'PASS' if ok else 'FAIL'}] {item}" + (f"  -> {detail}" if (detail and not ok) else ""))


def read(path):
    return io.open(path, encoding="utf-8").read()


def docx_text(path):
    xml = zipfile.ZipFile(path).read("word/document.xml").decode("utf-8", "ignore")
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", xml))


def run(script):
    try:
        p = subprocess.run([PY, os.path.join(CODE, script)], capture_output=True,
                           text=True, encoding="utf-8", errors="replace", timeout=400)
        return p.returncode == 0, (p.stdout or "") + (p.stderr or "")
    except Exception as e:
        return False, str(e)


print("=" * 70)
print("310 · NO-REGRESSION ASSERTIONS")
print("=" * 70)

# ---------- (A) 三连门禁 ----------
print("\n(A) 三连门禁")
ok, out = run("148_oac_compliance.py")
record("code/148 OAC COMPLIANT", "OVERALL OAC COMPLIANT: YES" in out,
       "" if "OVERALL OAC COMPLIANT: YES" in out else "148 未 COMPLIANT")
ok, out = run("149b_verify_reported_values.py")
record("code/149b 已报数值 SAFE", "ALL REPORTED VALUES IDENTICAL" in out and "SAFE" in out,
       "" if ("ALL REPORTED VALUES IDENTICAL" in out and "SAFE" in out) else "149b 非 SAFE")
ok, out = run("182_final_qc.py")
record("code/182 综合 QC 全 PASS", "ALL CHECKS PASSED" in out,
       "" if "ALL CHECKS PASSED" in out else "182 有 FAIL")

# ---------- (B) docs/08 理论边界深化 ----------
print("\n(B) docs/08 理论边界深化（L135 S23 图注）")
t = read(MANU)
record("docs/08: S22 正交维度前提", "orthogonal to the PC2 OA axis" in t)
record("docs/08: 单一设定点判别式 novelty", "single-set-point discriminator" in t)
record("docs/08: 内核判别式语料保留（shared set of inherited variants, opposite signs）",
       "read through a shared set of inherited variants, must fix opposite signs at whole-bone and joint endpoints" in t)

# ---------- (C) 派生稿同步 ----------
print("\n(C) 派生稿同步（镜像 + 两份 docx）")
m = read(MIRROR)
record("Manuscript_Ready 镜像: 含新理论串", "orthogonal to the PC2 OA axis" in m)
record("Manuscript_Ready 镜像: 含 novelty 串", "single-set-point discriminator" in m)
for name, path in DOCX.items():
    dt = docx_text(path)
    record(f"{name}: 含新理论串", "orthogonal to the PC2 OA axis" in dt)
    record(f"{name}: 含 novelty 串", "single-set-point discriminator" in dt)

# ---------- (D) 已报数值零改动 ----------
print("\n(D) 已报数值零改动（5 KEY 串）")
for k in KEY:
    record(f"docs/08 含 KEY: {k}", k in t)
for name, path in DOCX.items():
    dt = docx_text(path)
    record(f"{name} 含 KEY: {k}", k in dt) if False else None
for k in KEY:
    for name, path in DOCX.items():
        dt = docx_text(path)
        record(f"{name} 含 KEY: {k}", k in dt)

# ---------- (E) 物料清单完整 ----------
print("\n(E) 物料清单完整")
main_files, supp_files = [], []
for d in ["Figures", "Manuscript_Ready/Figures/Main"]:
    if os.path.isdir(d):
        main_files += os.listdir(d)
for d in ["Figures", "Manuscript_Ready/Figures/Supplementary"]:
    if os.path.isdir(d):
        supp_files += os.listdir(d)
allfig = main_files + supp_files
miss_s = [f"S{i}" for i in range(1, 25)
          if not any(f.startswith(f"SupplFig{i}_") for f in allfig)]
record("S1–S24 图文件齐全", not miss_s, str(miss_s) if miss_s else "24/24")
tab_files = os.listdir("Manuscript_Ready/Supplementary") if os.path.isdir("Manuscript_Ready/Supplementary") else []
miss_t = [f"Table_S{i}" for i in (8, 9, 12, 13, 14, 15, 16, 17, 26, 27, 29)
          if not any(f.startswith(f"Table_S{i}_") for f in tab_files)]
record("独立补充表文件齐全", not miss_t, str(miss_t) if miss_t else "11/11")
inline = ["## Supplementary Table S25", "## Supplementary Table S28",
          "## Supplementary Table S30", "## Supplementary Table S31"]
miss_inline = [x for x in inline if x not in t]
record("内联补充表 S25/S28/S30/S31 在 docs/08", not miss_inline,
       str(miss_inline) if miss_inline else "4/4")

# ---------- 汇总 ----------
n_fail = sum(1 for _, ok, _ in results if not ok)
print("\n" + "=" * 70)
print(f"TOTAL: {len(results) - n_fail}/{len(results)} PASS   FAIL: {n_fail}")
print("VERDICT:", "NO REGRESSION - ALL ASSERTIONS PASSED" if n_fail == 0
      else f"{n_fail} ASSERTION(S) FAILED")
print("=" * 70)
sys.exit(0 if n_fail == 0 else 1)
