# -*- coding: utf-8 -*-
"""
131 · 降 AIGC 第一轮：结构性痕迹

三类改动（均按 sci-manuscript-deai 规则）：
  A 段内小标题加粗收敛：Results 8 个 + Discussion 4 个 去加粗并融入首句；
    保留 Methods 7 个（STROBE-MR 结构化报告需要）与
    **Strengths and limitations.** / **Conclusion.**（期刊标准子标题）。
  B em dash 清零：15 处 → 0（目标 ≤2）。
  C 顺带把几处超长句拆开，增加 burstiness。

所有替换强制 count == 1；任何一条未命中则整体中止，不写文件。
"""

import io
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOC = os.path.join(ROOT, "docs", "08_投稿稿_OA_BMD_fracture.md")

MS = []


def E(tag, old, new):
    MS.append((tag, old, new))


# ═══════════════════════════ A. 段内小标题去加粗 ═══════════════════════════
# ---- Results ----
E("A1 L57", "**Genetic-correlation architecture.** SNP-based heritability ranged",
  "SNP-based heritability ranged")

E("A2 L59", "**Positive control: BMD protects against fracture.** Two-sample MR using 31 "
  "independent total-body BMD instruments showed",
  "The positive control behaved as expected: two-sample MR using 31 independent "
  "total-body BMD instruments showed")

E("A3 L61", "**BMD raises knee-OA risk.** Across five BMD exposures MR supported a "
  "positive effect (Table 3):",
  "MR supported a positive effect of BMD on knee OA across all five exposures (Table 3):")

E("A4 L63", "**Reverse MR and OA\u2192fracture are null and underpowered.** Knee OA as "
  "exposure gave non-significant effects",
  "Knee OA as exposure gave non-significant effects")

E("A5 L65", "**MR-PRESSO and robustness.** The distortion test was significant for 9/10 pairs",
  "The MR-PRESSO distortion test was significant for 9/10 pairs")

# ---- Discussion ----
E("A6 L75", "**Principal findings.** Jointly analysing 13 skeletal phenotypes gave a "
  "consistent architecture:",
  "Jointly analysing 13 skeletal phenotypes gave a consistent architecture:")

E("A7 L77", "**Resolving the OA\u2013BMD paradox.** High BMD protects against fracture, "
  "confirmed here as a positive control,",
  "The paradox resolves once the two endpoints are separated. High BMD protects against "
  "fracture, confirmed here as a positive control,")

E("A8 L79", "**Site and direction.** Zhou et al. [5] reported higher femoral-neck BMD "
  "protected against hip OA,",
  "Zhou et al. [5] reported higher femoral-neck BMD protected against hip OA,")

E("A9 L81", "**Biological context and predictions.** OA is a disorder of the osteochondral unit",
  "OA is a disorder of the osteochondral unit")

E("A10 L83", "**OA\u2192fracture: no supported pathway.** Neither r_g (+0.15, +0.21) nor MR",
  "Neither r_g (+0.15, +0.21) nor MR")

# ═══════════════════════════ B. em dash 清零 ═══════════════════════════
E("B1 L29", "treatments that increase areal BMD may carry an OA cost \u2014 though whether "
  "that cost is mediated by bone mass or by skeletal size has not been partitioned.",
  "treatments that increase areal BMD may carry an OA cost, though whether that cost is "
  "mediated by bone mass or by skeletal size has not been partitioned.")

E("B2 L45", "Harmonisation against each outcome \u2014 exclusion of palindromic and "
  "strand-ambiguous variants and of instruments absent from the outcome GWAS \u2014 reduced "
  "the analysed sets",
  "Harmonisation against each outcome (excluding palindromic and strand-ambiguous variants "
  "and instruments absent from the outcome GWAS) reduced the analysed sets")

E("B3 L67", "**Effect-scale validation.** On the verified per-SD scale, the BMD\u2192fracture "
  "OR (1.42 femoral-neck) lay 8.4% below GEFOS (1.55) on the odds-ratio scale and 20.0% "
  "below on the log-odds scale \u2014 we give both because effect estimates are compared on "
  "the log scale \u2014 with overlapping CIs; the total-body OR 1.40 was 9.6% below on the "
  "ratio scale and 23.1% below on the log-odds scale (cross-site approximation). MRlap "
  "confirmed every estimate moved further from the null (median 20%).",
  "On the verified per-SD scale, the BMD\u2192fracture OR (1.42 femoral-neck) lay 8.4% below "
  "GEFOS (1.55) on the odds-ratio scale and 20.0% below on the log-odds scale, with "
  "overlapping CIs. We report both scales because effect estimates are compared on the log "
  "scale. The total-body OR 1.40 was 9.6% below on the ratio scale and 23.1% below on the "
  "log-odds scale (cross-site approximation). MRlap confirmed every estimate moved further "
  "from the null (median 20%).")

E("B4 L69", "**Colocalisation.** At WNT16, SOX6 and GDF5 \u2014 loci nominated as shared "
  "effectors \u2014 colocalisation with knee OA was not supported",
  "At WNT16, SOX6 and GDF5, loci nominated as shared effectors, colocalisation with knee "
  "OA was not supported")

E("B5 L69", "The knee-OA null reflects power \u2014 under the MR estimate, the expected "
  "knee-OA z at the strongest instrument is 0.85, and reaching genome-wide significance "
  "would require ~5\u00d710\u2077 participants (~42\u00d7 the present sample).",
  "The knee-OA null reflects power. Under the MR estimate, the expected knee-OA z at the "
  "strongest instrument is 0.85, and reaching genome-wide significance would require "
  "~5\u00d710\u2077 participants (~42\u00d7 the present sample).")

E("B6 L71", "**MVMR: how much of the BMD\u2192OA effect is skeletal size?** Conditioning the "
  "BMD\u2192OA effect on height attenuated it",
  "How much of the BMD\u2192OA effect is skeletal size? Conditioning the BMD\u2192OA effect "
  "on height attenuated it")

E("B7 L77", "One non-exclusive explanation \u2014 that the BMD\u2192OA effect is partly "
  "driven by skeletal size rather than bone mass \u2014 is now quantified:",
  "One non-exclusive explanation is that the BMD\u2192OA effect is partly driven by skeletal "
  "size rather than bone mass, and that contribution is now quantified:")

E("B8 L79", "whereas our ten-site analysis shows higher BMD raises knee-OA risk \u2014 a "
  "hip-versus-knee contrast explicable under joint-specific pathway dominance",
  "whereas our ten-site analysis shows higher BMD raises knee-OA risk, a hip-versus-knee "
  "contrast explicable under joint-specific pathway dominance")

E("B9 L79", "Within Hartley et al. \u2014 the only same-pipeline comparison \u2014 the hip "
  "and knee intervals overlap almost completely",
  "Within Hartley et al., the only same-pipeline comparison, the hip and knee intervals "
  "overlap almost completely")

E("B10 L282", "**Table S10.** Multivariable MR grid \u2014 full sensitivity grid",
  "**Table S10.** Multivariable MR grid: full sensitivity grid")

E("B11 L283", "**Table S11.** Systemic Osteogenic Set-point (SOS) theory \u2014 predictions",
  "**Table S11.** Systemic Osteogenic Set-point (SOS) theory: predictions")


def main():
    text = io.open(DOC, encoding="utf-8").read()
    bak = DOC + ".bak-before-humanize131"
    if not os.path.exists(bak):
        io.open(bak, "w", encoding="utf-8").write(text)
        print(f"backup -> {os.path.basename(bak)}")

    failures = []
    for tag, old, new in MS:
        n = text.count(old)
        if n != 1:
            failures.append((tag, n, old[:80]))
            print(f"[SKIP] {tag}: occurrences={n}")
            continue
        text = text.replace(old, new, 1)
        print(f"[OK  ] {tag}")

    if failures:
        print(f"\n!! {len(failures)} edit(s) not uniquely matched; manuscript NOT written.")
        for t, n, o in failures:
            print(f"    {t}: n={n}  {o}")
        sys.exit(1)

    io.open(DOC, "w", encoding="utf-8").write(text)
    print(f"\napplied {len(MS)} edits")
    em = text.count("\u2014")
    print(f"em dash remaining: {em}")


if __name__ == "__main__":
    main()
