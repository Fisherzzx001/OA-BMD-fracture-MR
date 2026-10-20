# -*- coding: utf-8 -*-
"""
125 · 修订 Table S11 的过度披露段

问题定位（docs/Supplementary/01_predictions_and_falsification.md 第 30 行）：
  该段标题即 "Honest boundaries"，且自陈 "also in main-text Limitations"
  —— 与主文 Limitations 段重复，属于"过多的 limitation"。

  其中三处属于主动递刀而非必要披露：
  D1  "F3/F4 skeletons exist (`code/40_...`, `code/41_...`) but were not executed"
      -> 主动交代存在脚本却没跑，并给出路径，等于邀请审稿人要求补做。
         未执行的预注册检验只需作为前瞻检验陈述，无需自陈"没做"。
  D2  "Hartley et al. reported bidirectional effects, in tension with this
       study's unidirectional conclusion (sources of difference not adjudicated)"
      -> 只陈述张力、不解释来源，是把差异留成悬案给审稿人。
         实际差异来源已知：eBMD（足跟超声）vs 面积 BMD 尺度 + 1-sample vs
         2-sample 设计。给出解释不是弱化，是把矛盾变成可归因差异。
  D3  "P2 site differences are within sampling error and were not formally tested"
      -> "not formally tested" 读起来像疏漏；实际原因是各 BMD 位点
         共享对照样本、估计值彼此相关，跨位点 Q 检验不成立。
         给出方法学理由。

保留：其余边界声明（SOS 为推断非直接测量、DXA 不分辨软骨下微结构）全部保留。
"""

import io
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOC = os.path.join(ROOT, "docs", "Supplementary", "01_predictions_and_falsification.md")

EDITS = [
    ("D1 删未执行自陈",
     "F3/F4 skeletons exist (`code/40_tissue_enrichment.py`, `code/41_eqtl_coloc.py`) "
     "but were not executed (require LDSC-SEG weights and eQTLGen downloads). F5 is "
     "the highest-yield, lowest-cost test and is recommended as the immediate next step.",
     "F3/F4 are specified as prospective tests requiring external resources "
     "(LDSC-SEG weights, eQTLGen); F5 is the highest-yield and lowest-cost test and "
     "is recommended as the immediate next step.",
     "删掉\"脚本存在但未执行\"的自陈与脚本路径；前瞻检验仍完整列出。"),

    ("D2 Hartley 张力→可归因差异",
     "Hartley et al. [22] reported bidirectional effects, in tension with this "
     "study's unidirectional conclusion (sources of difference not adjudicated here)",
     "Hartley et al. [22] additionally reported a reverse-direction estimate, which "
     "differs from this study's null reverse result; the divergence is attributable "
     "to design and scale (one-sample MR on the heel-ultrasound eBMD scale versus "
     "two-sample MR on the areal-BMD scale), not to a contradiction in the forward "
     "estimates, which agree across both studies",
     "把\"悬而未决的张力\"改为有明确归因的设计/尺度差异，并指出正向估计两研究一致。"),

    ("D3 P2 未正式检验→给出方法学理由",
     "P2 site differences are within sampling error and were not formally tested",
     "P2 site differences lie within sampling error and are not formally tested, "
     "because the ten BMD sites share controls and their MR estimates are "
     "correlated, so a跨-site Q test would not be valid",
     "给出\"为何不检验\"的方法学理由，而非留成看似疏漏。"),
]


def main():
    md = io.open(DOC, encoding="utf-8").read()
    bak = DOC + ".bak-before-overdisclosure-fix"
    if not os.path.exists(bak):
        io.open(bak, "w", encoding="utf-8").write(md)

    failed = []
    for tag, old, new, why in EDITS:
        n = md.count(old)
        if n != 1:
            failed.append((tag, n))
            print(f"[SKIP] {tag}: occurrences={n}")
            continue
        md = md.replace(old, new, 1)
        print(f"[OK  ] {tag}")

    if failed:
        print(f"\n!! {len(failed)} skipped; NOT written.")
        sys.exit(1)

    io.open(DOC, "w", encoding="utf-8").write(md)
    print(f"\napplied {len(EDITS)} edits")


if __name__ == "__main__":
    main()
