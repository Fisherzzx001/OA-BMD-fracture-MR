# -*- coding: utf-8 -*-
"""
310 · 防回归综合断言（no-regression guard）

在每一次整改后运行，确保：
  (A) 三连门禁全绿：code/148 OAC COMPLIANT / code/149b SAFE / code/182 全 PASS
  (B) docs/08 本轮理论边界深化已落地（S22 正交维度前提 + 单一设定点判别式 novelty）
  (C) 派生稿同步：Manuscript_Ready 镜像 + 两份 docx 均含新理论串
  (D) 已报数值零改动：5 个 KEY 串在 docs/08 与两份 docx 中均存在
  (E) 物料清单完整：S1–S24 图文件、16 份独立补表 .md、12 个内联补表标题；
      审计档案未混入交付集；两份 docx 无补表全文/内部编号

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


def docx_paragraphs(path):
    """Paragraph texts, one per <w:p>. A full supp. table opens its own
    paragraph 'Supplementary Table SNN.'; an in-text cross-reference never
    does, so this is robust to whitespace flattening."""
    xml = zipfile.ZipFile(path).read("word/document.xml").decode("utf-8", "ignore")
    out = []
    for p in re.split(r"</w:p>", xml):
        t = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", "", p)).strip()
        if t:
            out.append(t)
    return out


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
# v2.4 (2026-11-12) P1-5 改稿：S22 正交串去掉定冠词 "the"，并明确正交性对
# "单一设定点" 中性（R1-2/C 证据两用修复）。旧断言 "orthogonal to the PC2 OA axis"
# 已不存在；同步到新串，并把 R1-2 修复核心固化为实质检查（比旧门禁更严格）。
t = read(MANU)
S22_NEW = "PC1 bone-mass/fracture axis orthogonal to PC2 OA axis"
record("docs/08: S22 正交维度前提（新措辞）", S22_NEW in t)
record("docs/08: 正交性对单设定点中性（R1-2 修复）",
       "orthogonality alone is neutral about a single set-point" in t)
record("docs/08: 兼容双耦合模块（R1-2 修复）",
       "consistent also with two coupled modules" in t)
record("docs/08: 单一设定点判别式 novelty", "single-set-point discriminator" in t)
record("docs/08: 内核判别式语料保留（shared set of inherited variants, opposite signs）",
       "read through a shared set of inherited variants, must fix opposite signs at whole-bone and joint endpoints" in t)

# ---------- (C) 派生稿同步 ----------
print("\n(C) 派生稿同步（镜像 + 两份 docx）")
m = read(MIRROR)
record("Manuscript_Ready 镜像: 含新理论串", S22_NEW in m)
record("Manuscript_Ready 镜像: 含 R1-2 中性串",
       "orthogonality alone is neutral about a single set-point" in m)
record("Manuscript_Ready 镜像: 含 novelty 串", "single-set-point discriminator" in m)
for name, path in DOCX.items():
    dt = docx_text(path)
    record(f"{name}: 含新理论串", S22_NEW in dt)
    record(f"{name}: 含 R1-2 中性串",
           "orthogonality alone is neutral about a single set-point" in dt)
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
md_tabs = [f for f in tab_files if f.endswith(".md")]
miss_t = [f"Table_S{i}" for i in (8, 9, 12, 13, 14, 15, 16, 17,
                                  25, 26, 27, 28, 29, 30, 31, 32)
          if not any(f.startswith(f"Table_S{i}_") for f in md_tabs)]
record("独立补充表 .md 齐全 (16)", not miss_t, str(miss_t) if miss_t else "16/16")
leaked = [f for f in tab_files if "Audit_Trail" in f]
record("审计档案未混入交付集", not leaked, str(leaked) if leaked else "clean")
inline = ["## Supplementary Table S%d" % i
          for i in (8, 9, 14, 15, 25, 26, 27, 28, 29, 30, 31, 32)]
miss_inline = [x for x in inline if x not in t]
record("内联补充表 12 个标题在 docs/08", not miss_inline,
       str(miss_inline) if miss_inline else "12/12")
# Decisive compliance check on the delivered docx (not just the markdown).
for name, path in DOCX.items():
    dx = docx_text(path)
    # A full table opens its own paragraph 'Supplementary Table SNN.'; an
    # in-text cross-reference lives mid-sentence and never starts a paragraph.
    full = [p[:60] for p in docx_paragraphs(path)
            if re.match(r"^#{0,3}\s*Supplementary Table S\d+\.", p)]
    internal = re.findall(r"code/\d+|this revision|internal-consistency", dx, re.IGNORECASE)
    record(f"{name} 无补表全文/内部编号", (not full) and (not internal),
           f"full_tables={len(full)} internal={len(internal)}")

# ---------- (F) 第七轮极限审稿 P0-1/P0-3/P0-4/P1-1 实质整改 ----------
print("\n(F) docs/148 整改落地（P0-1 尺度披露 / P0-2 唯一权威表 / P0-3 逐SNP表 / P1-1 组织富集 / P0-4 协方差重建）")
m = read(MIRROR)
# P0-1: 尺度换算必须写明，且 "no rescaling" 必须彻底消失（正文真源+镜像+docx）
for label, text in [("docs/08", t), ("镜像", m)] + [(n, docx_text(p)) for n, p in DOCX.items()]:
    record(f"P0-1 {label}: 写明 √(400,000/31,873)=3.543", "3.543" in text)
    record(f"P0-1 {label}: 无 'no rescaling' 旧披露", "no rescaling" not in text)
# P1-1: 置换校准 + 软骨≥骨 + 撤回 formation-side 独占（正文真源+镜像+docx）
for label, text in [("docs/08", t), ("镜像", m)] + [(n, docx_text(p)) for n, p in DOCX.items()]:
    record(f"P1-1 {label}: 置换校准措辞", "permutation-calibrated" in text)
    record(f"P1-1 {label}: 撤回 formation-side 独占", "not formation-side-exclusive" in text)
    record(f"P1-1 {label}: 无旧串 'formation-side rather than'", "formation-side rather than" not in text)
# P1-1: S25 表新增 perm p 列 + contrast + 67% 污染披露（markdown 层，图注区）
record("P1-1 docs/08: S25 perm p 列", "perm p (length/chr-matched)" in t)
record("P1-1 docs/08: S25 cartilage-minus-bone contrast", "Cartilage-minus-bone contrast" in t)
record("P1-1 docs/08: 67% 软骨污染披露", "67%" in t and "shared with the cartilage sets" in t)
# P0-2: 唯一权威结果表声明 + 打包排除 raw/backup
auth = "docs/RESULTS_AUTHORITATIVE.md"
record("P0-2 权威数据声明存在", os.path.isfile(auth))
if os.path.isfile(auth):
    at = read(auth)
    record("P0-2 声明唯一权威表", "mr_results_corrected.csv" in at and "AUTHORITATIVE" in at)
    record("P0-2 记录 raw/backup 排除", "EXCLUDED" in at and "_backup_pre_G6" in at)
z = read("code/150_build_submission_zip.py")
# TABLE_CSVS must list the authoritative corrected table; raw mr_results.csv must
# NOT appear as an active entry (it may appear only in comments/docstrings).
record("P0-2 打包脚本用 corrected 权威表",
       '"mr_results_corrected.csv"' in z and "AUTHORITATIVE MR table" in z)
record("P0-2 打包脚本不交付 raw mr_results.csv",
       '    "mr_results.csv",' not in z and '"mr_results.csv", "mr_power.csv"' not in z)
record("P0-2 打包脚本交付 harmonised", "HARMONISED_CSVS" in z and "harmonised" in z)
# P0-3: 12 个 harmonised 逐 SNP 表存在且可复现 Table 2
harm_dir = os.path.join(CODE, "output", "genetics", "harmonised")
need_harm = ["BMD_total__fracture", "BMD_total__OA_knee", "BMD_femur_neck__fracture",
             "OA_knee__BMD_total", "OA_knee__fracture"]
record("P0-3 harmonised 逐 SNP 表齐全 (12)",
       os.path.isdir(harm_dir) and len([f for f in os.listdir(harm_dir) if f.endswith(".csv")]) == 12
       and all(os.path.isfile(os.path.join(harm_dir, p + ".csv")) for p in need_harm))
# P0-4: R1 summary 含 bootstrap p，且不再把 MRlap 截距当抽样协方差
import json
r1p = os.path.join(CODE, "output", "genetics", "_r1_model_comparison_summary.json")
if os.path.isfile(r1p):
    r1 = json.load(open(r1p, encoding="utf-8"))
    record("P0-4 R1 用 bootstrap 校准 LRT", "p_LRT_bootstrap" in r1 and r1["p_LRT_bootstrap"] is not None)
    record("P0-4 MRlap λ 降格为参照（非抽样相关）",
           "lambda_lcsc_intercept" in r1 and "NOT used as sampling correlation" in r1.get("covariance_note", ""))
    record("P0-4 结论限定: ρ=−1 被拒但耦合强度未定位",
           "rho_hat" in r1 and abs(r1["rho_hat"]) < 0.5 and r1["rho_compat_CI"]["rho_eq_0_within_CI"])
else:
    record("P0-4 R1 summary 存在", False)
r1src = read("code/340_r1_single_vs_dual_module.py")
record("P0-4 code/340 协方差=共享分母+ρ_s (非 MRlap λ)",
       "cov_sampling = (bo_fx * bo_oa / be ** 4) * se_be ** 2" in r1src
       and "rho_s * (so_fx * so_oa / np.abs(be) ** 2)" in r1src)

# ---------- 汇总 ----------
n_fail = sum(1 for _, ok, _ in results if not ok)
print("\n" + "=" * 70)
print(f"TOTAL: {len(results) - n_fail}/{len(results)} PASS   FAIL: {n_fail}")
print("VERDICT:", "NO REGRESSION - ALL ASSERTIONS PASSED" if n_fail == 0
      else f"{n_fail} ASSERTION(S) FAILED")
print("=" * 70)
sys.exit(0 if n_fail == 0 else 1)
