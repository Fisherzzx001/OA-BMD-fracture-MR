# -*- coding: utf-8 -*-
"""
133 · 第十二轮整改 P0-3：结构式摘要（OAC 官方要求）
------------------------------------------------------------------
OAC Manuscript checklist 原文：
  "350 word abstract structured into the following sections:
   (1) Objective; (2) Design (3) Results and (4) Conclusions."
此前按 250 词（Brief report 的限制）压缩，属自伤性删减。本脚本改为四段结构式，
并用回收的约 100 词补回三项此前被迫删掉的内容：
  · BMD→OA 的 OR（1.017–1.028 / SD）与 CI；
  · 5/10 位点进入 MR 的口径限定（LDSC 覆盖 13 表型，MR 覆盖 10 个）；
  · 未条件化 BMI 的边界。
同时保留所有门禁必需表述（G7/G11/G18 等）。
"""
import io
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MD = os.path.join(ROOT, "docs", "08_投稿稿_OA_BMD_fracture.md")

ABS = """**Objective.** Whether higher bone mineral density (BMD) shares aetiology with, or causes, osteoarthritis (OA) is unsettled, although people with OA tend to have higher BMD. We quantified genetic sharing and causal direction across 13 skeletal phenotypes and tested whether one BMD-raising set-point has opposing consequences at the fracture and joint endpoints.

**Design.** Linkage-disequilibrium score regression (LDSC) and two-sample Mendelian randomisation (MR) were applied to 13 GWAS summary-statistic datasets (two OA definitions, ten regional BMD sites, any-fracture). LDSC covered all 78 phenotype pairs; MR ran in the ten phenotypes with clumped instruments — five of ten prespecified BMD sites, both OA definitions, any-fracture — with BMD→fracture as positive control, multivariable MR conditioning on height, MRlap correction, and reverse MR.

**Results.** BMD sites inter-correlated strongly (r_g 0.69–0.96) and negatively with fracture (−0.45 to −0.60); both OA definitions correlated positively with every BMD site (+0.12 to +0.25). Genetically higher BMD raised knee-OA risk at all five analysed sites — OR 1.017–1.028 per SD (95% CI 1.007–1.038) — and protected against fracture (OR 1.40, 95% CI 1.34–1.46 per SD lower BMD; p=3.7×10⁻⁵¹). Heterogeneity was high (I² 68–86%) but matched the positive control (78–83%), so the causal reading rests on the concordance of pleiotropy-robust estimators (MR-Egger, weighted median, MR-PRESSO) rather than on fixed-effect p ≤ 6.4×10⁻⁴ alone; under random effects 2 of 5 remained significant. Conditioning on height attenuated the OA effect by 9.5–33.2% (a lower bound for size), leaving 67–91% height-independent — not identified as bone mass, and possibly still BMI-mediated. Reverse estimates were non-significant and underpowered (power 0.14–0.26). MRlap moved all ten estimates further from the null (median 20.5%), reversing none.

**Conclusions.** One genetically indexed bone phenotype is protective at the whole-bone endpoint and harmful at the joint endpoint, 67–91% of the joint effect independent of height. The joint effect is small in individual terms (OR ≈1.02 per SD areal BMD) yet attaches to a ubiquitous, modifiable exposure, and is measured at five analysed sites only. These estimates do not support changing prescribing practice; excluding a clinically meaningful joint cost requires trials powered for incident OA recorded at the joint of interest alongside fracture."""

t = io.open(MD, encoding="utf-8").read()
L = t.split("\n")
i = next(k for k, ln in enumerate(L) if ln.strip() == "## Abstract")
j = next(k for k in range(i + 1, len(L)) if L[k].strip() == "---")
old_abs = "\n".join(L[i + 1:j]).strip()
assert old_abs.startswith(("Whether higher bone mineral density", "**Objective.**")), \
    f"未定位到摘要块：{old_abs[:60]!r}"
# 按索引重建（幂等：无论当前是旧单段还是已结构化的版本都被整体替换）
L = L[:i + 1] + ["", ABS, ""] + L[j:]
io.open(MD, "w", encoding="utf-8").write("\n".join(L))

# ---------------------------------------------------------------- 校验
chk = io.open(MD, encoding="utf-8").read()
ab = chk.split("## Abstract")[1].split("\n---")[0]
words = len([w for w in ab.split() if re.search(r"[0-9A-Za-z]", w)])
print(f"摘要词数 = {words} （OAC 上限 350）")
assert words <= 350, f"摘要超长：{words}"
for h in ["**Objective.**", "**Design.**", "**Results.**", "**Conclusions.**"]:
    assert h in ab, f"缺少结构式小标题 {h}"
print("结构式四标题齐备 OK")
for must, lab in [("fixed-effect p ≤", "G18"), ("estimator concordance rather than", "G18"),
                  ("height-independent", "G11"), ("analysed sites", "G15"),
                  ("9.5–33.2%", "G7"), ("67–91%", "G7"), ("0.14–0.26", "G7"),
                  ("I² 68–86%", "G7"), ("MR-Egger", "G7"), ("positive control", "G7"),
                  ("median 20.5%", "MRlap"), ("−0.45 to −0.60", "r_g"),
                  ("0.69–0.96", "r_g"), ("1.017–1.028", "P1-1 OR")]:
    assert must in chk, f"门禁/整改必需表述缺失：{must}（{lab}）"
print("门禁必需表述齐备 OK")
print()
print(ABS)
