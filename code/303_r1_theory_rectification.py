# -*- coding: utf-8 -*-
"""303 · R1「理论贡献质询师」三漏洞审视后的理论建构整改（2026-11-08）

依据：最近一轮全面预审报告（综合 7.2/10）中 `## 【Reviewer 1】理论贡献质询师`
的「三个最致命的理论漏洞」与「深化理论建构的方向」。

审视结论（重读 docs/08 真源逐条核验，不复用旧结论）：
  R1-1 机制内核未获数据支持  -> ②③ 已加「方向模式被支持、骨形成机制未被证实」，
      但 L93 仍写 SOS 是骨形成三联征的 "genetic counterpart"（等同口吻），
      未显式辨析 SOS 与 rivals 的差别，也未给出可区分 rivals 的判别式内核 -> 残留。
  R1-2 P2 排序子句阴性       -> ③ 已改为 "partially supported"（L95 与 L93 双处）-> 闭环。
  R1-3 "survives 9/10" 完备性 -> ② 已在 S23 图注加限定并删正文 "survived F6–F8"，
      但正文结论未同步显眼边界 -> 残留。
  深化方向（把 SOS 从描述性模式升级为带边界的可机制化断言：单一设定点 =
      同一组遗传变异、由其符号决定整骨与关节两端点反向；F2 体积臂 + F10
      跨祖先为关闭漏洞的优先后续）-> 全稿未落地 -> 本轮主载荷。

落地（只增补/改写理论表述，不改任何已报告 beta/SE/p/OR/CI，不新增参考文献编号）：
  V1 正文 §4 假说段：把 SOS 写成带判别式内核的可证伪断言，并显式区别于
      描述性骨-软骨悖论 [3,4] 与表型层面的骨形成三联征 [39]。
  V2 正文 §4 假说段末：补「存活计数只覆盖已执行检验，不等于登记册完备性」。
  V3 S23 图注（Figure legends 段，不计 4000 词正文）：完整理论边界 ——
      内核定义 / rivals 分歧 / 三条可证伪条件 / 现有 F4,F6,F7 只是间接探测 /
      F2 体积臂与 F10 为优先后续检验。
  V4 S30 note vii（不计正文词数）：明确 SOS 与骨形成三联征是「种类之别而非程度之别」。
  V6/V7 收尾核验发现的**内部矛盾**：L93 前半原称 SOS 是骨形成三联征的
      "genetic counterpart"（等同口吻），与本轮新增的 rivals / "in kind, not degree"
      判别式自相矛盾 -> 改为 "a testable hypothesis distinct in kind from … and from …"。

派生稿：两份 submission_package/*.docx 采用 run 级定点替换（保留版式与元数据）；
Manuscript_Ready/01_Manuscript.md 复制真源并校验 md5。
配套：code/300_fig_S23_sos_framework.py 的 verdict 文案同步改为判别式内核，
      重跑生成 S23 三目录图片（本脚本不负责图片，见整改记录 §7）。

幂等：new 已存在则跳过；old 与 new 均不存在则报错退出。
用法： python3 code/303_r1_theory_rectification.py
"""
import hashlib
import os
import shutil

from docx import Document

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)

MANU = "docs/08_投稿稿_OA_BMD_fracture.md"
DOCX = ["submission_package/Manuscript.docx",
        "submission_package/Submission_MainText.docx"]
MIRROR_MANU = "Manuscript_Ready/01_Manuscript.md"
BACKUP_DIR = ".workbuddy/backups/r1_theory_2026-11-08"

# (id, old, new) -- anchors verified unique in MANU and in a single run of each DOCX.
REPLS = [
    # First remove the internal contradiction: L93 opened by calling SOS the
    # "genetic counterpart" of the bone-forming triad, which would contradict
    # the rivals / "in kind, not degree" framing added below.
    ("R1-V6 premise",
     "and read as the genetic counterpart of the bone-forming triad",
     "distinct in kind from the bone-forming triad"),
    ("R1-V7 parallel",
     "and of the broader osteoporosis–osteoarthritis",
     "and from the broader osteoporosis–osteoarthritis"),

    ("R1-V1 kernel",
     "Unlike the descriptive bone–cartilage paradox [3, 4], SOS is a falsifiable causal "
     "assertion, not a restatement of the pattern: the arms result",
     "Unlike the descriptive bone–cartilage paradox [3, 4] or the phenotype-level "
     "bone-forming triad [39], SOS is a falsifiable causal assertion, not a restatement "
     "of the pattern. Its discriminating kernel is one latent set-point: the same inherited "
     "variants must fix opposite signs at whole-bone and joint endpoints, and distinct "
     "variant sets or concordant signs would refute it (Fig. S23). The arms result"),

    ("R1-V3 completeness",
     "The reading generates P1–P6 and F1–F10 (Table S11).",
     "The reading generates P1–P6 and F1–F10 (Table S11); survival counts executed "
     "tests only, not register completeness (Fig. S23)."),

    ("R1-V4 figure legend",
     "the osteogenic set-point mechanism is supported indirectly, not demonstrated (§4).",
     "the osteogenic set-point mechanism is supported indirectly, not demonstrated (§4). "
     "Theory boundary. The discriminating kernel is not a label for the observed pattern: "
     "one latent set-point, read through a shared set of inherited variants, must fix "
     "opposite signs at whole-bone and joint endpoints, with sign set by the set-point "
     "rather than by endpoint-specific pathways. The descriptive bone–cartilage paradox "
     "[3, 4] and the phenotype-level bone-forming triad [39] are rivals because neither "
     "requires one variant set nor sign opposition. The kernel is falsifiable: "
     "endpoint-specific instrument sets, loss of sign opposition under volumetric BMD "
     "instruments, or ancestry-specific reversal would refute it. F4 (cross-trait "
     "colocalisation), F6 (cross-site homogeneity) and F7 (joint-locus exclusion) probe "
     "it only indirectly; the F2 volumetric arm (true vBMD instruments) and F10 "
     "cross-ancestry test are the priority subsequent tests that could close the two "
     "remaining gaps."),

    ("R1-V5 note vii",
     "the 2025 MR evidence [11, 12] and the osteophyte-rich endotypes [46] are consistent "
     "with, but do not demonstrate, its formation-side core (§4).",
     "the 2025 MR evidence [11, 12] and the osteophyte-rich endotypes [46] are consistent "
     "with, but do not demonstrate, its formation-side core (§4). SOS therefore differs "
     "from the bone-forming triad [39] in kind, not degree: the triad describes "
     "co-occurring phenotypes, whereas SOS requires one variant set with sign-opposing "
     "causal effects at two endpoints — a constraint neither the triad nor the "
     "descriptive paradox [3, 4] entails (Fig. S23)."),
]


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
            print("  [skip] %-22s already applied" % label)
            continue
        n = txt.count(old)
        if n != 1:
            raise SystemExit("ABORT: %s -> %s anchor count %d (expected 1) in %s"
                             % (label, old[:40], n, path))
        txt = txt.replace(old, new, 1)
        changed = True
        print("  [ok]   %-22s applied" % label)
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
            print("  [skip] %-22s already applied" % label)
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
        print("  [ok]   %-22s applied" % label)
    doc.save(path)
    after = docx_text(Document(path))
    for label, old, new in todo:
        assert new in after, (label, path)
        # V3/V4/V5 append to the anchor, so old survives as a prefix of new.
        if old not in new:
            assert old not in after, (label, path)
    print("  [verify] %s: all new passages present; replaced anchors absent" % path)


def main():
    backup([MANU] + DOCX)

    print("\n[1] main manuscript:", MANU)
    apply_text(MANU, REPLS)

    print("\n[2] derived docx (run-level, layout-preserving):")
    for p in DOCX:
        print(" ", p)
        apply_docx(p, REPLS)

    print("\n[3] Manuscript_Ready mirror:")
    shutil.copy2(MANU, MIRROR_MANU)
    same = md5(MANU) == md5(MIRROR_MANU)
    print("  [%s] %s <-> %s" % ("ok" if same else "FAIL", MANU, MIRROR_MANU))
    if not same:
        raise SystemExit("ABORT: mirror md5 mismatch")

    print("\nDONE. Next: re-run code/300 (figure verdict), then "
          "code/148_oac_compliance.py ; code/149b_verify_reported_values.py ; "
          "code/182_final_qc.py")


if __name__ == "__main__":
    main()
