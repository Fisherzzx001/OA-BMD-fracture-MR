# -*- coding: utf-8 -*-
"""
130 · 第十二轮整改（正文文本项）
------------------------------------------------------------------
落地 docs/115_第十二轮预审报告_OAC.md 中的文本类条目：
  P0-1  L83 误挂 [35]
  P0-5  正文侧 BMI 通路披露（Methods + Limitations）
  P1-1  正文报出 BMD→OA 的 OR（Results；Table 3 由 131 脚本重建）
  P1-2  L95 结论按量级降调
  P1-4  5/10 位点入 Limitations
  P1-5  性别/性别报告（Methods + Limitations）
  P1-8  [23] 无效引用
  P1-9  Steiger 不显著的过度推论
  P1-10 F6 措辞统一为"未检出异质性"并明写判别力有限
  P2-1  引用 [7] Bowden 2016（weighted median）
  P2-2  L303 F1–F5 → F1–F7
  P2-3  P6 改写的前瞻性条款
  P2-4  L83 "introduced below" 措辞
  P2-5  标题口径（13 phenotypes 为 LDSC 全集）
  P2-7  SOS 降格为启发式 + 补 vBMD 形式的定量内容
  P2-8  [33] 改为"综述中回顾"式引用
  P2-6  Acknowledgements
  P2-10 生成式 AI 使用声明

每条替换前断言 old 唯一命中；结束后自检残留。
"""
import io
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MD = os.path.join(ROOT, "docs", "08_投稿稿_OA_BMD_fracture.md")

t = io.open(MD, encoding="utf-8").read()
orig = t
L = t.split("\n")
n_edits = 0


def rep(old, new, label, must_line_startswith=None):
    """唯一命中替换；可限定所在行前缀（行号会因前序替换漂移）。"""
    global t, n_edits
    hits = [i for i, ln in enumerate(t.split("\n")) if old in ln]
    if must_line_startswith is not None:
        hits = [i for i in hits if t.split("\n")[i].startswith(must_line_startswith)]
    assert len(hits) == 1, f"[{label}] 命中 {len(hits)} 次（应为 1）: {old[:70]!r}"
    t = t.replace(old, new, 1)
    n_edits += 1
    print(f"  OK  {label}")


print("=" * 88)
print("130 · 第十二轮整改：正文文本项")
print("=" * 88)

# ---------------------------------------------------------------- P2-5 标题
rep("Genetically higher bone mineral density raises knee osteoarthritis risk yet protects "
    "against fracture: bidirectional Mendelian randomisation of 13 skeletal phenotypes",
    "Genetically higher bone mineral density raises knee osteoarthritis risk yet protects "
    "against fracture: genetic correlation and bidirectional Mendelian randomisation across "
    "13 skeletal phenotypes",
    "P2-5 标题口径限定为 13 表型（LDSC 全集）")

# ---------------------------------------------------------------- P1-8 / P2-8 引言引用
rep("Osteoarthritis (OA) and bone mineral density (BMD) are the two most prevalent skeletal "
    "traits in ageing populations [23], and the nature",
    "Osteoarthritis (OA) is among the most prevalent joint disorders of ageing [23], and areal "
    "bone mineral density (BMD) — the clinical measure of skeletal strength — is measured in "
    "the same populations; the nature",
    "P1-8 [23] 不再支撑未被其论证的命题")
rep("Observational studies have linked higher areal BMD to increased OA prevalence and "
    "incidence [33],",
    "Observational studies, reviewed in [33], have linked higher areal BMD to increased OA "
    "prevalence and incidence,",
    "P2-8 [33] 改为述评式引用（不再以社论支撑实证命题）")

# ---------------------------------------------------------------- P1-5 性别（Methods）
rep("This secondary analysis required no new ethical approval.",
    "This secondary analysis required no new ethical approval. All four source GWAS contribute "
    "sex-combined primary estimates and none of the deposited statistics is stratified by sex, "
    "so every estimate in this study is sex-averaged and effect modification by sex cannot be "
    "tested here — a material limitation given that both BMD and OA are sex-dimorphic (§4).",
    "P1-5 数据来源补性别构成/性别分层不可行",
    must_line_startswith="**Data sources.**")

# ---------------------------------------------------------------- P2-1 引 [7]
rep("weighted median and leave-one-out, with Cochran Q",
    "weighted median [7] and leave-one-out, with Cochran Q",
    "P2-1 引用 weighted median 方法学原文 [7]")

# ---------------------------------------------------------------- P0-5 BMI（Methods）
rep("a grid (n_h ∈ {0,100,300,1,000}) and both EUR and EUR_noUKB height sources are in the "
    "Supplementary Material.",
    "a grid (n_h ∈ {0,100,300,1,000}) and both EUR and EUR_noUKB height sources are in the "
    "Supplementary Material. Height was chosen over body-mass index (BMI) deliberately: BMI is "
    "itself a major knee-OA risk factor and part of the adiposity pathway through which a "
    "BMD-raising process could raise knee-OA risk, whereas height indexes skeletal size without "
    "adiposity. The consequence is stated explicitly rather than absorbed silently: this "
    "multivariable analysis does not condition on BMI, so the height-independent residual may "
    "still contain a BMI-mediated component, and the size-attributable share is bounded from "
    "below only (§4). A BMI-conditioned multivariable analysis is the natural extension and "
    "requires BMI summary statistics that are not part of this deposit.",
    "P0-5 Methods 明示 BMI 未条件化及其后果")

# ---------------------------------------------------------------- P1-1 正文报 OR
rep("Direction was consistent across all five. MR-Egger intercepts were all p>0.24",
    "Direction was consistent across all five. On the odds-ratio scale the effect is small at "
    "the individual level — OR 1.017–1.028 per SD higher BMD across the five sites, with "
    "per-site 95% confidence intervals in Table 3 — against OR 1.40 (1.20–1.63) per SD in the "
    "largest one-sample estimate [22], a 12–20-fold difference on the log-odds scale adjudicated "
    "in §4. MR-Egger intercepts were all p>0.24",
    "P1-1 Results 报出 BMD→OA 的 OR 与量级")

# ---------------------------------------------------------------- P1-9 Steiger
rep("and was non-significant for the three reverse pairs (|z|<0.1, p≈0.48–0.49), consistent "
    "with the null and underpowered reverse MR and with BMD, not OA, driving the association.",
    "and was non-significant for the three reverse pairs (|z|<0.1, p≈0.48–0.49). That "
    "non-result does not discriminate direction: with 69 instruments and a null reverse "
    "estimate the test lacks the power to separate the forward from the reverse model, and it "
    "is not read here as evidence that BMD rather than OA drives the association.",
    "P1-9 Steiger 不显著不再被读作支持性证据")

# ---------------------------------------------------------------- P0-1 / P2-4 / P2-3
rep("Zhou et al. [5] reported higher femoral-neck BMD protected against hip OA [35],",
    "Zhou et al. [5] reported higher femoral-neck BMD protected against hip OA,",
    "P0-1 删除误挂的 [35]（Faber 为髋形状研究）")
rep("P6 and the falsification rule F1 are introduced below and in Table S11.",
    "P6 and the falsification rule F1 are specified in Table S11 and adjudicated in §4.",
    "P2-4 P6/F1 指针措辞（G16 门禁同步更新）")
rep("The directional form was the pre-registered prediction; its restatement here deviates from "
    "that registered form and is prompted by the contradictory hip evidence.",
    "The directional form was the pre-registered prediction; its restatement here deviates from "
    "that registered form and is prompted by the contradictory hip evidence. The restated, "
    "magnitude form is itself a new prediction: its decision rule scores prospectively collected "
    "hip estimates and is not applied retrospectively to the published estimates adjudicated "
    "above.",
    "P2-3 P6 改写的前瞻性条款")

# ---------------------------------------------------------------- P2-7 SOS 定量内容
rep("The set-point is operationalised as areal BMD — a composite of volumetric density and "
    "size — and is a latent tendency rather than a measured quantity.",
    "The set-point is operationalised as areal BMD — a composite of volumetric density and "
    "size — and is a latent tendency rather than a measured quantity; it is used here as a "
    "heuristic that organises the predictions below, not as a fitted model. Its quantitative "
    "content is this: if the causal quantity is volumetric density, then an areal-BMD estimate "
    "is diluted by skeletal size, so instruments for volumetric BMD (or areal BMD conditioned "
    "on bone area) should yield a BMD→OA effect at least as large as the areal one — the "
    "quantitative form of F2, which the present areal data cannot test.",
    "P2-7 SOS 明确为启发式 + 补 vBMD 形式的定量预测")

# ---------------------------------------------------------------- P1-10 F6 措辞
rep("the five site-specific BMD→knee-OA estimates are statistically indistinguishable "
    "(Q=6.92 on 4 df, p=0.14, with the covariance induced by shared instruments handled "
    "exactly; τ²=0), whereas a site-specific rival predicts heterogeneity the data do not show.",
    "heterogeneity among the five site-specific BMD→knee-OA estimates was not detected "
    "(Q=6.92 on 4 df, p=0.14, τ²=0, with the covariance induced by shared instruments handled "
    "exactly). F6 is therefore not rejected, but it is a weak test: with five estimates and "
    "site-mean r_g of 0.69–0.96, a single systemic set-point and correlated site-specific "
    "processes both predict homogeneity, so F6 discriminates little even before its limited "
    "power to detect heterogeneity is considered.",
    "P1-10 F6 措辞统一为未检出 + 明写判别力有限")

# ---------------------------------------------------------------- P1-4 / P1-5 / P0-5 局限性
rep("Areal BMD conflates size with mass [18,32,37];",
    "Five of the ten prespecified BMD sites entered MR because strict clumping completed for "
    "them in the compute environment; the five deferred sites showed comparable genome-wide "
    "signal, so the analysed subset is a partition determined by compute availability rather "
    "than by biology and is not random with respect to that availability (Table S12). Areal BMD "
    "conflates size with mass [18,32,37];",
    "P1-4 5/10 位点入 Limitations")
rep("so the 67–91% height-independent residual is not established as bone mass (resolved only "
    "by test F2).",
    "so the 67–91% height-independent residual is not established as bone mass (resolved only "
    "by test F2). Conditioning was on height, not BMI, so a BMI-mediated component of the "
    "BMD→OA effect is not excluded and that residual is a lower bound on the non-size share "
    "rather than an estimate of bone mass. All estimates are sex-averaged and effect "
    "modification by sex is untested, because the deposited statistics are not sex-stratified.",
    "P0-5 + P1-5 Limitations 补 BMI 与性别")

# ---------------------------------------------------------------- P1-2 结论降调
rep("Trials of treatments that raise areal BMD should record incident OA at the joint of "
    "interest, not fracture alone.",
    "The joint-endpoint effect is small in individual terms (OR 1.017–1.028 per SD higher BMD, "
    "95% CI 1.007–1.038 across sites) but attaches to a ubiquitous and modifiable exposure; "
    "these estimates do not support changing prescribing practice, and excluding a clinically "
    "meaningful joint cost requires trials or pharmacogenetic designs powered for incident OA "
    "recorded at the joint of interest alongside fracture.",
    "P1-2 结论按效应量降调")

# ---------------------------------------------------------------- P2-2 附件清单
rep("- **Table S11.** Systemic Osteogenic Set-point (SOS) theory: predictions (P1–P6) and "
    "falsifying tests (F1–F5).",
    "- **Table S11.** Systemic Osteogenic Set-point (SOS) theory: predictions (P1–P6) and "
    "falsifying tests (F1–F7).",
    "P2-2 附件清单 F1–F5 → F1–F7")

# ---------------------------------------------------------------- P2-6 / P2-10 声明
_kw = "**Keywords.**"
assert t.count(_kw) == 1
_concl = "\n**Conclusion.**"
assert t.count(_concl) == 1
_add = ("\n**Acknowledgements.** The authors thank the investigators and participants of the "
        "source GWAS consortia for making summary statistics publicly available. No funding "
        "was received for this work.\n\n**Declaration of generative AI in scientific writing.** "
        "An AI writing assistant was used during preparation for language editing, formatting "
        "and internal-consistency checking; all analyses, results, interpretations and "
        "conclusions are the authors' own and were verified against the deposited code and "
        "outputs. The authors take full responsibility for the content of the article.\n")
_marker = "\n## Author contributions"
if _marker in t:
    t = t.replace(_marker, _add + _marker, 1)
    n_edits += 1
    print("  OK  P2-6 Acknowledgements + P2-10 生成式 AI 声明（置于 Author contributions 前）")
else:
    print("  !!  未找到 '## Author contributions'，声明未插入")
    sys.exit(1)

# ---------------------------------------------------------------- 自检
io.open(MD, "w", encoding="utf-8").write(t)
print()
print("---- 自检：禁止残留 ----")
_bad = []
for s, lab in [("protected against hip OA [35]", "P0-1 [35] 误挂"),
               ("introduced below and in Table S11", "P2-4 旧措辞"),
               ("falsifying tests (F1–F5)", "P2-2 F1–F5"),
               ("statistically indistinguishable", "P1-10 indistinguishable"),
               ("driving the association", "P1-9 Steiger 旧措辞")]:
    if s in t:
        _bad.append(lab)
assert not _bad, f"残留: {_bad}"
print("  无残留 OK")
assert "1.017–1.028" in t, "OR 未写入"
assert "BMI" in t, "BMI 披露未写入"
assert "sex-averaged" in t, "性别披露未写入"
assert "[7]" in t and "weighted median [7]" in t, "Bowden 引用未写入"
print(f"  必需项齐备 OK；本脚本共完成 {n_edits} 处修改")
print(f"  词数变化：{len([w for w in orig.split() if re.search(r'[0-9A-Za-z]', w)])} → "
      f"{len([w for w in t.split() if re.search(r'[0-9A-Za-z]', w)])}")
