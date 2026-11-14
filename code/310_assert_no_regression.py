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

# ---------- (G) 第八轮 §2 整改 P1-2/3/4 + P2 ----------
print("\n(G) docs/148 §2 整改落地（P1-2 双分母 / P1-3 EUR-only / P1-4 K双尺度 / P2-1..5）")
# 正文新串须同时出现在 docs/08、镜像与两份 docx（三载体；防止只改真源漏派生）。
msrc = read("docs/08_投稿稿_OA_BMD_fracture.md")
mirror = read("Manuscript_Ready/01_Manuscript.md")
carriers = {"md": msrc, "mirror": mirror,
            "Manuscript.docx": docx_text(md_path:=("submission_package/Manuscript.docx")),
            "MainText.docx": docx_text("submission_package/Submission_MainText.docx")}
BODY_NEW = {
    "P1-2 并列 34% BMD-referenced 分母": "≤34% on the BMD-referenced basis",
    "P1-2 保留 8.4%-of-fracture 参照": "≤4.5% of the positive control on the fracture-referenced overlap",
    "P1-3 EUR-only 方向复制指针": "an EUR-only knee-OA GWAS preserves all five positive r_g",
    "P2-6 GDF5 腰椎反例进正文": "one lumbar GDF5 locus reverses this endpoint split",
    "P2-5 hip 改 cross-site 同研究": "restoring P6's pre-registered direction",
    "P2-3 阳性对照 only partly calibrates": "the control only partly calibrates the pipeline",
    "P2-3 directional pleiotropy flagged": "carries directional pleiotropy",
    "P2-1 Highlights 因果降调": "associates with higher knee OA, lower fracture",
    "P2-4 achieved power 标 post-hoc": "achieved power is computed post-hoc",
}
def flat(s):
    return re.sub(r"\s+", " ", s)
for item, needle in BODY_NEW.items():
    hit = {k: (flat(needle) in flat(v)) for k, v in carriers.items()}
    record(item, all(hit.values()), detail=str(hit))
# S33 index lives in the Supplementary-Material index block, which
# build_submission.py strips from BOTH docx (cut from "## Supplementary Table
# S8." onward). It is therefore correctly absent from the docx carriers. Scope
# this assertion to the two markdown carriers only (true source + byte mirror).
record("Table S33 索引行存在",
       all(flat("**Table S33.** Section-2 review remediation audit") in flat(v)
           for v in (msrc, mirror)),
       detail="md/mirror only; docx strips the supplementary index")
# 旧强断言串必须消失（四载体）
for item, old in [
    ("P2-1 Highlights 旧 raises 已去", "- Higher genetically predicted BMD raises knee OA risk"),
    ("P2-5 prospective replication 旧串已去", "Our prospective hip estimate resolves F1"),
    ("P2-3 calibrates 旧绝对串已去", "so the control calibrates the pipeline despite its non-clean diagnostic"),
]:
    record(item, all(old not in v for v in carriers.values()))

# 补表实质断言（独立补表文件，不计正文词数）
s14 = read("docs/Supplementary/Table_S14_ancestry.md")
record("P1-2 S14 修正 12% 误标并说明几何平均",
       "mislabeled the geometric-mean figure" in s14 and "geometric-mean overlap" in s14
       and "~34–36% on the literal 70%-of-BMD basis" in s14)
record("P1-2 S14 MRlap 校正对照 vs GEFOS (1.475/1.521)",
       "1.475" in s14 and "1.521" in s14 and "not a strictly two-sample comparison" in s14)
s17 = read("docs/Supplementary/Table_S17_sensitivity.md")
record("P1-4 S17-a 39SNP bootstrap CI 落表",
       "bootstrap 95% CI of 0.085457–0.087291" in s17 and "3.385–3.421" in s17)
record("P1-4 S17-a 双尺度披露（OR<1% vs log-odds~9%）",
       "Dual reporting scale" in s17 and "±9% K-sensitivity" in s17)
record("P1-4 S17-a 三 K 口径调和",
       "K=0.0964" in s17 and "K=0.096 is" in s17 and "K=0.12 is an external population-prevalence" in s17)
s33 = flat(read("docs/Supplementary/Table_S33_review_remediation.md"))
record("P1-2 S33 Part A 三分母齐全 (4.10/11.83/34.10)",
       "4.10%" in s33 and "11.83%" in s33 and "34.10%" in s33 and "Label correction" in s33)
record("P1-2 S33 Part B MRlap relative shift",
       "−0.1222→−0.1409" in s33 and "1.475" in s33 and "1.521" in s33)
record("P1-3 S33 Part C EUR-only 5/5 方向一致 + GO×FinnGen r_g",
       "Sign agrees" in s33 and "0.137" in s33 and "not publicly reconstructible" in s33)
record("P2-2 S33 Part D 反向 79→69 + Egger 翻转",
       "79 harmonised instruments" in s33 and "p=0.0392" in s33 and "p=0.376" in s33
       and "excluded from the deliverable" in s33)
s27s = read("docs/Supplementary/Table_S27_hip_oa_replication.md")
record("P2-5 S27 标题/导语 cross-site 且非独立队列",
       "Within-study cross-site" in s27s and "not an independent cohort" in s27s)

# 三份审计 JSON 存在且关键字段齐全（可复跑数值证据）
def jget(d, *path):
    for k in path:
        d = d[k]
    return d
for item, path, need in [
    ("P1-2 JSON 三分母+MRlap转移", "code/output/genetics/P1_2_overlap_dual_denominator.json",
     lambda j: ("bias_on_8.4pct_pct" in jget(j, "burgess_bias", "BMD_total->fracture")
                and "bias_on_geometric_mean_pct" in jget(j, "burgess_bias", "BMD_total->fracture")
                and "bias_on_70pct_pct" in jget(j, "burgess_bias", "BMD_total->fracture")
                and "BMD_total" in j["mrlap_transferred_positive_control"])),
    ("P1-3 JSON EUR-only 方向+祖先错配", "code/output/genetics/P1_3_ancestry_ldsc_sensitivity.json",
     lambda j: all(k in j for k in ["eur_only_independent_replication", "GO_vs_FinnGen_kneeOA", "h2_liab_sensitivity"])),
    ("P1-4 JSON K bootstrap+双尺度", "code/output/genetics/P1_4_K_bootstrap.json",
     lambda j: all(k in j["K_estimation"] for k in ["K_boot_CI", "scale_factor_boot_CI"])
                and "forward_beta_propagation" in j
                and "S17a_plusminus10pct_K_on_both_scales" in j),
]:
    if os.path.isfile(path):
        j = json.load(open(path, encoding="utf-8"))
        record(item, need(j) and all(c["pass"] for c in j["checks"]))
    else:
        record(item, False)

# ---------- 汇总 ----------
n_fail = sum(1 for _, ok, _ in results if not ok)
print("\n" + "=" * 70)
print(f"TOTAL: {len(results) - n_fail}/{len(results)} PASS   FAIL: {n_fail}")
print("VERDICT:", "NO REGRESSION - ALL ASSERTIONS PASSED" if n_fail == 0
      else f"{n_fail} ASSERTION(S) FAILED")
print("=" * 70)
sys.exit(0 if n_fail == 0 else 1)
