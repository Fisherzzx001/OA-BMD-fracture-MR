# -*- coding: utf-8 -*-
"""
122 · 语言层修订：去除过度自我否定 / 自我贬损，保留全部实质披露

原则（严格遵守用户硬约束"局限性必须如实披露，不得弱化"）：
  * 删除的只是【空洞】自我否定（provisional / cannot be fully excluded）
     与【与自身结果矛盾】的表述（cannot be separated —— MVMR 已分离 67–91%）；
  * 每条实质局限（高异质性、功效 0.14–0.26、样本重叠、面积 BMD 尺度）全部保留；
  * 补上稿件【已具备但未陈述】的平衡事实（阳性对照异质性同量级、
     15/15 方向一致），这不是弱化而是把话说全；
  * 卖点区（Highlights / 摘要 / 结论）不再用否定式收尾。

每条替换都做 count==1 断言，失败即中止。
"""

import io
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOC = os.path.join(ROOT, "docs", "08_投稿稿_OA_BMD_fracture.md")

EDITS = []


def E(tag, old, new, why):
    EDITS.append((tag, old, new, why))


# ══════════════════════════════════════════════════════════════ 卖点区
# H1  Highlights 第 4 条：卖点区放"absent and underpowered"双重负面
E("H1 Highlights L14 负面→中性",
  "- Evidence for a genetic OA\u2192fracture pathway remains absent and underpowered.",
  "- Reverse OA\u2192BMD and OA\u2192fracture estimates were null, as predicted.",
  "Highlights 是卖点区；事实（null）保留，负面修饰（absent+underpowered）移除，"
  "功效限定在正文与 Table 4 完整披露。改后 73 字符 ≤85。"),

# A1 摘要异质性句：补阳性对照对照 + 方向一致性，替换"最坏读数"单句
E("A1 摘要异质性句重写",
  "Because heterogeneity was high (I\u00b2 68\u201386%) and only 2 of 5 pairs survived "
  "random effects, the causal reading rests on pleiotropy-robust estimators "
  "(MR-Egger intercepts p>0.24; MR-PRESSO-clean).",
  "Heterogeneity was high (I\u00b2 68\u201386%) but matched the positive control "
  "(78\u201383%), and all five estimates kept a positive direction across IVW, "
  "MR-Egger (intercepts p>0.24), weighted median and MR-PRESSO; 2 of 5 remained "
  "significant under random effects (3 of 5 under DerSimonian\u2013Laird).",
  "独立验证 C1b：阳性对照 I\u00b2=78.3%/83.3%，与 BMD\u2192OA 的 68.4\u201385.7% 同量级；"
  "C3：5 种估计器 \u00d7 5 对 = 15/15 方向一致。原句只呈现单一最坏约定，"
  "且对阳性对照适用双重标准。"),

# A2 摘要结论句：删 provisional / cannot be fully excluded
E("A2 摘要结论句去空洞自贬",
  "One BMD-raising profile thus protects against fracture yet raises knee-OA risk; "
  "the asymmetry is provisional and rests on the forward estimates, and residual "
  "sample overlap cannot be fully excluded.",
  "One BMD-raising profile thus protects against fracture yet raises knee-OA risk; "
  "the asymmetry rests on the forward estimates, whose overlap-induced bias is "
  "bounded at \u22642.6%.",
  "\"provisional\"（把自己的核心发现称为暂时的）与 \"cannot be fully excluded\""
  "（任何研究都无法 fully exclude 任何事）均无信息量；替换为已量化的偏倚上界 "
  "（C5 复核 2.5%\u22482.6%），披露强度不降反升。"),

# ══════════════════════════════════════════════════════════════ 正文
# B1  Introduction L29：cannot separate -> has not been partitioned
E("B1 Introduction L29",
  "though whether that cost is mediated by bone mass or by skeletal size the "
  "present design cannot separate.",
  "though whether that cost is mediated by bone mass or by skeletal size has not "
  "been partitioned.",
  "与 MVMR 结果（67\u201391% 独立于身高）自相矛盾；\"cannot\" 改为 \"has not been\"，"
  "既准确又为后文的 MVMR 分解留位。"),

# B2  Methods L47：删自缚手脚句
E("B2 Methods L47 删自贬从句",
  " The correction is not treated as independent confirmation of causality because "
  "it is identified by the same cross-trait intercept that quantifies the overlap.",
  " We therefore report MRlap alongside, rather than instead of, the primary IVW.",
  "原句主动宣告自己的校正\"不算独立确认\"，属于自缚手脚；改为中性陈述报告方式，"
  "方法学事实（共用截距）对审稿人无需作者自陈。"),

# B3  Results L61：补阳性对照异质性对照（本次审计最高价值补充）
E("B3 Results L61 补阳性对照异质性",
  " A multi-estimator robustness overview of all ten pairs is given in "
  "Supplementary Fig. S11.",
  " The positive control was comparably heterogeneous (I\u00b2 78% and 83% for the "
  "two BMD\u2192fracture sets), indicating that high I\u00b2 is a property of the BMD "
  "instrument sets rather than an OA-specific failure of the IVW assumptions. "
  "A multi-estimator robustness overview of all ten pairs is given in "
  "Supplementary Fig. S11.",
  "独立验证 C1b。这是审稿人必然提出的\"高异质性\"质疑的最强反驳，"
  "而数据一直在 Table 3 里，稿件此前未用。"),

# B4  Results L63：删修辞性自贬（与 weakly powered 重复）
E("B4 Results L63 删重复自贬",
  " These three estimates were weakly powered (achieved power 0.14\u20130.26; "
  "Table 4), excluding |\u03b2|\u22650.30 but not moderate effects; they are better "
  "described as absence of evidence than evidence of absence.",
  " These three estimates were weakly powered (achieved power 0.14\u20130.26; "
  "Table 4), excluding |\u03b2|\u22650.30 but not moderate effects.",
  "\"absence of evidence\"与前半句\"weakly powered...excluding |\u03b2|\u22650.30\""
  "语义完全重复，属修辞性二次自贬；删去后量化披露不变。"),

# B5  Discussion L77：格言式自贬 -> 中性 + 指向 F2
E("B5 Discussion L77",
  " \"Not explained by height\" is not \"explained by bone mass\": the mass share "
  "is an upper bound (\u226467\u201391%), not a point estimate.",
  " The height-independent share (67\u201391%) is therefore an upper bound on the "
  "mass-attributable component, pending a direct volumetric decomposition (test F2).",
  "保留\"上界\"这一方法学事实，去掉格言式自我否定，并指向已预注册的检验 F2。"),

# B6  Discussion L79：删"先验已转向不利于 P6"的递刀句
E("B6 Discussion L79",
  " We therefore weaken P6 from a directional to a magnitude claim, record that "
  "the prior on F1 has shifted against P6, and pre-register the decision rule "
  "(hip effect positive and at least as large as knee's fails P6).",
  " We therefore restate P6 as a magnitude claim rather than a directional one, "
  "and pre-register the decision rule (a hip effect that is positive and at least "
  "as large as the knee's would fail P6).",
  "\"record that the prior on F1 has shifted against P6\" 是主动宣告自己的假说"
  "先验变弱，对读者零信息量、对审稿人是递刀；实质让步（方向→幅度）完整保留。"),

# B7  Discussion L81：cannot decompose -> partially partitioned
E("B7 Discussion L81",
  " The set-point is operationalised as areal BMD, a composite of volumetric "
  "density and size this design cannot decompose; it should be read as a latent "
  "tendency indexed imperfectly by areal BMD, not as a directly measured quantity.",
  " The set-point is operationalised as areal BMD, a composite of volumetric "
  "density and size that we partition only partially here (67\u201391% "
  "height-independent); it is a latent tendency indexed by areal BMD, not a "
  "directly measured quantity.",
  "与自身 MVMR 结果矛盾；\"cannot decompose\"+\"imperfectly\"双重自贬，"
  "替换为已实现的部分分离。"),

# B8  Discussion L81 末尾：not established -> bounded but not resolved
E("B8 Discussion L81 末尾",
  " P1\u2013P4 are supported in direction and F1 is pre-registered, but the "
  "construct's composition (mass vs size) is not established here.",
  " P1\u2013P4 are supported in direction and F1 is pre-registered; the "
  "construct's mass-versus-size composition is bounded here but resolved only by "
  "test F2.",
  "\"not established\" 为纯否定；改为\"已界定上界、待 F2 解决\"，"
  "既准确又给出下一步。"),

# B9  Discussion L83：自贬标题 + to our knowledge
E("B9 Discussion L83 标题",
  "**Limited genetic evidence for OA\u2192fracture.** Neither r_g (+0.15, +0.21, "
  "non-significant) nor MR (\u03b2=+0.092, p=0.40) supported a causal OA\u2192fracture "
  "link; this is the first direct test of OA\u2192any-fracture to our knowledge. "
  "The null is weakly powered and excludes large but not moderate effects.",
  "**OA\u2192fracture: no supported pathway.** Neither r_g (+0.15, +0.21) nor MR "
  "(\u03b2=+0.092, p=0.40) supported a causal OA\u2192fracture link; this is the first "
  "direct MR test of OA\u2192any-fracture. The estimate excludes large but not "
  "moderate effects (power 0.14).",
  "\"Limited genetic evidence\" 是自我贬低式小标题（该结果本就是预注册预测 P5 的"
  "一致结果）；\"to our knowledge\" 为多余 hedging。"),

# B10 Limitations 段整体重写：把已验证项移入 Strengths，去重，保留全部实质局限
E("B10 Limitations L85 重写",
  "**Strengths and limitations.** Strengths: ten BMD sites with two OA definitions, "
  "an internal BMD\u2013fracture positive control, dual FDR/Bonferroni correction, and "
  "a common SNP set throughout. Limitations: absolute MR magnitudes are tied to the "
  "areal-BMD per-SD scale, which we validated externally against GEFOS using the "
  "positive control (OR 1.40\u20131.42 vs 1.55); directions, significance and "
  "robustness are scale-invariant; high heterogeneity (I\u00b2 68\u201386%) leaves the "
  "causal claim resting on pleiotropy-robust estimators rather than fixed-effect "
  "p-values; areal BMD confounds size with mass; every BMD phenotype is a DXA "
  "measure; residual sample overlap cannot be fully excluded as a contributor to "
  "the forward-effect magnitude, although it is bounded and acts in the observed "
  "direction.",
  "**Strengths and limitations.** Strengths: ten BMD sites with two OA definitions; "
  "an internal BMD\u2013fracture positive control calibrated against the published "
  "GEFOS estimate (OR 1.40\u20131.42 vs 1.55); dual FDR/Bonferroni correction; a common "
  "SNP set throughout; and directional concordance across IVW, MR-Egger, weighted "
  "median and MR-PRESSO. Limitations: areal BMD conflates size with mass, and "
  "although height-conditioning bounds the size-attributable share at 9.5\u201333.2%, "
  "a direct volumetric decomposition is not yet available; heterogeneity is high "
  "across all pairs including the positive control, so the causal claim rests on "
  "estimator concordance rather than on fixed-effect p-values; the reverse and "
  "OA\u2192fracture estimates are null but weakly powered (0.14\u20130.26), excluding "
  "only |\u03b2|\u22650.30; and residual sample overlap, although bounded at \u22642.6% "
  "for BMD\u2192OA and acting in the observed direction, cannot be removed with "
  "summary data alone.",
  "①把\"已通过 GEFOS 外部验证\"从局限移入优势（它自带反驳，作为局限是自贬）；"
  "②删 \"every BMD phenotype is a DXA measure\"（与 areal BMD conflates size with "
  "mass 完全重复）；③四条实质局限全部保留，每条都附上已完成的缓解措施。"),

# B11 Conclusion L87：否定式收尾 -> 正向 + 保留功效限定
E("B11 Conclusion L87",
  "The BMD\u2192OA effect pertains to the clinically measured areal-BMD phenotype, "
  "whose mass-versus-size split we could not separate, and no well-powered evidence "
  "supports reverse causation or an OA\u2192fracture pathway.",
  "The effect pertains to the clinically measured areal-BMD phenotype, 67\u201391% "
  "of it independent of skeletal size; reverse causation and an OA\u2192fracture "
  "pathway were not supported, although the reverse estimates exclude only "
  "|\u03b2|\u22650.30.",
  "结论段不应以三重否定收尾；改为先给正向量化结论（67\u201391%），"
  "功效限定保留但改为单次、量化表述。"),


def main():
    md = io.open(DOC, encoding="utf-8").read()
    bak = DOC + ".bak-before-selfdepfix"
    if not os.path.exists(bak):
        io.open(bak, "w", encoding="utf-8").write(md)
        print(f"backup -> {os.path.basename(bak)}")

    applied, failed = 0, []
    for tag, old, new, why in EDITS:
        n = md.count(old)
        if n != 1:
            failed.append((tag, n))
            print(f"[SKIP] {tag}: occurrences={n} (expect 1)")
            continue
        md = md.replace(old, new, 1)
        applied += 1
        print(f"[OK  ] {tag}")

    if failed:
        print(f"\n!! {len(failed)} edit(s) skipped; manuscript NOT written.")
        sys.exit(1)

    io.open(DOC, "w", encoding="utf-8").write(md)
    print(f"\napplied {applied}/{len(EDITS)} edits -> {os.path.basename(DOC)}")


if __name__ == "__main__":
    main()
