# -*- coding: utf-8 -*-
"""302 · R2 预审「建议改 3-6」四项文字整改（2026-11-08）

落地范围（只改表述，不改任何已报告 beta/SE/p/OR/CI）：
  #3 结论段：把 F2 条件 MVMR 的弱工具限制从 Table S28 提到正文结论段
  #4 摘要：明示双向 MR 的反向检验功效不足、未提供方向性证据
  #5 §2 Methods：给出 Burgess 偏倚上界的定量值，并说明最悲观重叠校正后仍为正
  #6 Table 3 表注：补随机效应 IVW 下实际显著比例（2/5），与 power 列并置消歧

同时：
  * docs/Supplementary/Table_S14_ancestry.md 补 Burgess 上界段（使正文引用的
    "Table S14" 在独立补充表载体上也可追溯）；
  * 两份 submission_package/*.docx 做 run 级定点替换（保留既有版式与元数据），
    避免全量重建带来的未知差异；
  * Manuscript_Ready/ 的稿件与 S14 表镜像同步。

幂等：new 已存在则跳过；old 与 new 均不存在则报错退出。
用法： python3 code/302_r2_minor_rectifications.py
"""
import hashlib
import os
import shutil
import sys

from docx import Document

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)

MANU = "docs/08_投稿稿_OA_BMD_fracture.md"
DOCX = ["submission_package/Manuscript.docx",
        "submission_package/Submission_MainText.docx"]
SUPPL_S14 = "docs/Supplementary/Table_S14_ancestry.md"
MIRROR_MANU = "Manuscript_Ready/01_Manuscript.md"
MIRROR_S14 = "Manuscript_Ready/Supplementary/Table_S14_ancestry.md"
BACKUP_DIR = ".workbuddy/backups/r2_minor_2026-11-08"

REPLS = [
    # (id, old, new)
    ("R2-4 abstract",
     "Reverse estimates were non-significant and underpowered (power 0.14–0.26).",
     "Reverse estimates were non-significant and underpowered (power 0.14–0.26), "
     "so the bidirectional MR provided no directional evidence."),
    ("R2-5 methods",
     "the internally estimated per-pair cross-trait intercepts and the default-setting "
     "sensitivity are in Table S14.",
     "the internally estimated per-pair cross-trait intercepts and the default-setting "
     "sensitivity are in Table S14. Worst-case Burgess biases from the tabulated overlaps "
     "are ≤0.7% of the BMD→knee-OA estimate and ≤4.5% of the positive-control estimate "
     "(Table S14); MRlap retained all ten directions (median +20.5% away from the null), "
     "so the positive joint-endpoint effect survives overlap correction."),
    ("R2-3 conclusion",
     "and not transmitted by bone geometry (F2);",
     "and not transmitted by bone geometry (F2), albeit on weak conditional instruments "
     "(conditional F below 4 in all but one mode), so this is direction-indicative rather "
     "than magnitude-precise (Table S28, Part B);"),
    ("R2-6 Table 3 note",
     "K = strict instrument count.",
     "K = strict instrument count. Achieved power refers to fixed-effect IVW; under "
     "the primary random-effects IVW only two of five forward BMD→knee-OA estimates "
     "remain significant (Table S17-b; §3)."),
]

BURGESS_NOTE = (
    "\n**Burgess overlap-bias bounds (bias = overlap fraction × r_g; §2 of the manuscript).** "
    "The tabulated counts are BMD 31,873 participants (100% UKB), OA_selfrep 462,933 "
    "(100% UKB), fracture 22,330 UKB participants (8.4% of 264,973; 70% of BMD), and "
    "OA_knee UKB participation not separately reported (overlap upper bound "
    "31,873 / 1,208,397 = 2.6%). Hence BMD→OA_knee bias ≤ ~0.3–0.7% "
    "(2.6% × r_g 0.12–0.25), BMD→OA_selfrep ~0.8–1.7% (6.9% × 0.12–0.25) and "
    "BMD→fracture ~3.4–4.5% on the conservative 8.4%-of-fracture basis (up to ~12% on "
    "the 70%-of-BMD basis); all act in the direction of r_g. MRlap accounts for overlap "
    "directly and retained all ten directions (median +20.5% away from the null), so the "
    "positive BMD→knee-OA effect survives the worst-case overlap correction.\n"
)
BURGESS_ANCHOR = "**Burgess overlap-bias bounds"


def md5(path):
    return hashlib.md5(open(path, "rb").read()).hexdigest()


def backup(paths):
    # Never overwrite an existing snapshot: an idempotent rerun must not
    # replace the pre-rectification backup with the already-patched file.
    os.makedirs(BACKUP_DIR, exist_ok=True)
    n = 0
    for p in paths:
        dst = os.path.join(BACKUP_DIR, os.path.basename(p))
        if not os.path.exists(dst):
            shutil.copy2(p, dst)
            n += 1
    print("[backup] %d new snapshot(s) -> %s" % (n, BACKUP_DIR))


def apply_text(path, repls):
    txt = open(path, encoding="utf-8").read()
    changed = False
    for label, old, new in repls:
        if new in txt:
            print("  [skip] %-20s already applied" % label)
            continue
        n = txt.count(old)
        if n != 1:
            raise SystemExit("ABORT: %s -> %s anchor count %d (expected 1) in %s"
                             % (label, old[:40], n, path))
        txt = txt.replace(old, new, 1)
        changed = True
        print("  [ok]   %-20s applied" % label)
    if changed:
        open(path, "w", encoding="utf-8").write(txt)
    return txt


def docx_text(doc):
    parts = [p.text for p in doc.paragraphs]
    for tb in doc.tables:
        for row in tb.rows:
            for cell in row.cells:
                parts += [p.text for p in cell.paragraphs]
    return "\n".join(parts)


def apply_docx(path, repls):
    doc = Document(path)
    full = docx_text(doc)
    todo = []
    for label, old, new in repls:
        if new in full:
            print("  [skip] %-20s already applied" % label)
        elif old not in full:
            raise SystemExit("ABORT: neither old nor new found for %s in %s" % (label, path))
        else:
            todo.append((label, old, new))
    if not todo:
        return
    counts = {old: 0 for _, old, _ in todo}
    def scan(par):
        for run in par.runs:
            for label, old, new in todo:
                if old in run.text:
                    run.text = run.text.replace(old, new, 1)
                    counts[old] += 1
    for par in doc.paragraphs:
        scan(par)
    for tb in doc.tables:
        for row in tb.rows:
            for cell in row.cells:
                for par in cell.paragraphs:
                    scan(par)
    for label, old, new in todo:
        if counts[old] != 1:
            raise SystemExit("ABORT: run-level replacement count %d for %s in %s"
                             % (counts[old], label, path))
        print("  [ok]   %-20s applied" % label)
    doc.save(path)
    after = docx_text(Document(path))
    for label, old, new in todo:
        assert new in after, (label, path)
        # #5/#6 append to the anchor, so old survives as a prefix of new.
        if old not in new:
            assert old not in after, (label, path)
    print("  [verify] %s: all new sentences present; replaced anchors absent" % path)


def main():
    backup([MANU] + DOCX + [SUPPL_S14])

    print("\n[1] main manuscript:", MANU)
    apply_text(MANU, REPLS)

    print("\n[2] independent supplementary S14:", SUPPL_S14)
    s14 = open(SUPPL_S14, encoding="utf-8").read()
    if BURGESS_ANCHOR in s14:
        print("  [skip] Burgess bounds note already present")
    else:
        if not s14.endswith("\n"):
            s14 += "\n"
        open(SUPPL_S14, "w", encoding="utf-8").write(s14 + BURGESS_NOTE)
        print("  [ok]   Burgess bounds note appended")

    print("\n[3] derived docx (run-level, layout-preserving):")
    for p in DOCX:
        print(" ", p)
        apply_docx(p, REPLS)

    print("\n[4] Manuscript_Ready mirrors:")
    shutil.copy2(MANU, MIRROR_MANU)
    shutil.copy2(SUPPL_S14, MIRROR_S14)
    for a, b in [(MANU, MIRROR_MANU), (SUPPL_S14, MIRROR_S14)]:
        same = md5(a) == md5(b)
        print("  [%s] %s <-> %s" % ("ok" if same else "FAIL", a, b))
        if not same:
            raise SystemExit("ABORT: mirror md5 mismatch")

    print("\nDONE. Next: python3 code/148_oac_compliance.py ; "
          "python3 code/149b_verify_reported_values.py ; python3 code/182_final_qc.py")


if __name__ == "__main__":
    main()
