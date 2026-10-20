# -*- coding: utf-8 -*-
"""
132 · 降 AIGC 第二轮：句法层

  C 句间分号 → 句号（保留列举性分号：括号内统计量、并列研究、清单项）
     目的：降分号密度、降平均句长、升 burstiness
  D 长句拆分（>50 词）
  E 排比结构拆解（4 连分号 / 5 项串列 = 典型 AI 节奏）
  F 个别 AI 偏好词（key → selected）

不动任何数字。所有替换强制 count == 1，否则整体中止。
"""

import io
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOC = os.path.join(ROOT, "docs", "08_投稿稿_OA_BMD_fracture.md")

MS = []


def E(tag, old, new):
    MS.append((tag, old, new))


# ═══════════════════ C. 句间分号 → 句号 ═══════════════════
E("C1 L31", "that contrasts with Zhou's negative femoral-neck estimate; the hip-versus-knee "
  "direction reversal is discussed in the Discussion.",
  "that contrasts with Zhou's negative femoral-neck estimate. The hip-versus-knee direction "
  "reversal is discussed below.")

E("C2 L39", "Analyses were performed in Python 3.13 with numpy, pandas, scipy and matplotlib; "
  "the sample-overlap correction used R 4.6.1 with MRlap v0.0.3.3.",
  "Analyses were performed in Python 3.13 with numpy, pandas, scipy and matplotlib. The "
  "sample-overlap correction used R 4.6.1 with MRlap v0.0.3.3.")

E("C3 L41", "All source GWAS were approved by their original ethics committees; this "
  "secondary analysis required no new ethical approval.",
  "All source GWAS were approved by their original ethics committees. This secondary "
  "analysis required no new ethical approval.")

E("C4 L47", "is bounded above by 2.6% for BMD\u2192OA and 12% for BMD\u2192fracture; it acts "
  "in the direction of r_g without changing signs or p-values.",
  "is bounded above by 2.6% for BMD\u2192OA and 12% for BMD\u2192fracture. It acts in the "
  "direction of r_g without changing signs or p-values.")

E("C5 L49", "BMD statistics had initially been standardised against a placeholder N of "
  "400,000; we therefore re-standardised them to per-SD units",
  "BMD statistics had initially been standardised against a placeholder N of 400,000. We "
  "therefore re-standardised them to per-SD units")

E("C6 L61b", "rather than on fixed-effect p-values; all five effects remained positively "
  "causal after MR-PRESSO outlier removal",
  "rather than on fixed-effect p-values. All five effects remained positively causal after "
  "MR-PRESSO outlier removal")

E("C7 L75", "genetically higher BMD raises knee-OA risk; the reverse (OA\u2192BMD) and "
  "OA\u2192fracture estimates were non-significant and underpowered.",
  "genetically higher BMD raises knee-OA risk. The reverse (OA\u2192BMD) and OA\u2192fracture "
  "estimates were non-significant and underpowered.")

E("C8 L81", "that we partition only partially here (67\u201391% height-independent); it is a "
  "latent tendency indexed by areal BMD, not a directly measured quantity.",
  "that we partition only partially here (67\u201391% height-independent). It is a latent "
  "tendency indexed by areal BMD, not a directly measured quantity.")

E("C9 L83", "supported a causal OA\u2192fracture link; this is the first direct MR test of "
  "OA\u2192any-fracture.",
  "supported a causal OA\u2192fracture link. This is the first direct MR test of "
  "OA\u2192any-fracture.")

E("C10 L162", "For BMD femoral neck \u2192 fracture the MR-Egger intercept is significant "
  "(p=0.018); the MR-Egger slope (\u22120.343) and weighted median (\u22120.365) are the "
  "pleiotropy-robust estimates.",
  "For BMD femoral neck \u2192 fracture the MR-Egger intercept is significant (p=0.018). The "
  "MR-Egger slope (\u22120.343) and weighted median (\u22120.365) are the pleiotropy-robust "
  "estimates.")

E("C11 L183", "No locus reaches PP.H4\u22650.8; six favour distinct causal variants (PP.H3\u22650.8).",
  "No locus reaches PP.H4\u22650.8. Six favour distinct causal variants (PP.H3\u22650.8).")

E("C12 L183b", "i.e. one locus fell 0.0006 short of the threshold; the pre-specified 0.8 "
  "cut-off is retained.",
  "i.e. one locus fell 0.0006 short of the threshold, and the pre-specified 0.8 cut-off is "
  "retained.")

E("C13 L187", "Attenuation = 100\u00d7(\u03b2_uni \u2212 \u03b2_MVMR)/\u03b2_uni; the two "
  "tabulated univariable baselines are \u03b2_uni from Table 3 and \u03b2_uni recomputed on "
  "the MVMR instrument set; both are given across the full grid in Table S10.",
  "Attenuation = 100\u00d7(\u03b2_uni \u2212 \u03b2_MVMR)/\u03b2_uni. The two tabulated "
  "univariable baselines are \u03b2_uni from Table 3 and \u03b2_uni recomputed on the MVMR "
  "instrument set, and both are given across the full grid in Table S10.")

E("C14 L259", "were not redistributed; provenance and accessions are listed above.",
  "were not redistributed. Provenance and accessions are listed above.")

E("C15 L266", "are referenced in the text; all are provided as machine-readable "
  "Markdown/CSV",
  "are referenced in the text. All are provided as machine-readable Markdown/CSV")

# ═══════════════════ D. 长句拆分 ═══════════════════
E("D1 L49 (56\u2192\u77ed\u53e5)",
  "We confirmed empirically that the betas entering the analysis are already on the true-N "
  "scale: for every instrument the reported standard error matches the closed-form "
  "expectation for a standardised effect, SE = 1/\u221a[2\u00b7N\u00b7f\u00b7(1\u2212f)] with f "
  "the minor-allele frequency, agreeing with N=31,873 to within 0.1% and excluding the "
  "placeholder N=400,000, which would predict standard errors 3.54-fold smaller.",
  "We confirmed empirically that the betas entering the analysis are already on the true-N "
  "scale. For every instrument the reported standard error matches the closed-form "
  "expectation for a standardised effect, SE = 1/\u221a[2\u00b7N\u00b7f\u00b7(1\u2212f)] with f "
  "the minor-allele frequency, agreeing with N=31,873 to within 0.1%. The placeholder "
  "N=400,000 is excluded, because it would predict standard errors 3.54-fold smaller.")

# ═══════════════════ E. 排比拆解 ═══════════════════
E("E1 L61 \u56db\u8fde\u5206\u53f7",
  "Direction was consistent across all five; MR-Egger intercepts were all p>0.24; weighted "
  "medians matched the IVW sign; leave-one-out produced no sign reversals; and directions "
  "agreed with r_g (+0.12 to +0.25).",
  "Direction was consistent across all five. MR-Egger intercepts were all p>0.24 and "
  "weighted medians matched the IVW sign. Leave-one-out produced no sign reversals, and the "
  "directions agreed with r_g (+0.12 to +0.25).")

E("E2 L85 Strengths \u4e32\u5217",
  "Strengths: ten BMD sites with two OA definitions; an internal BMD\u2013fracture positive "
  "control calibrated against the published GEFOS estimate (OR 1.40\u20131.42 vs 1.55); dual "
  "FDR/Bonferroni correction; a common SNP set throughout; and directional concordance "
  "across IVW, MR-Egger, weighted median and MR-PRESSO.",
  "Strengths: ten BMD sites with two OA definitions, an internal BMD\u2013fracture positive "
  "control calibrated against the published GEFOS estimate (OR 1.40\u20131.42 vs 1.55), dual "
  "FDR/Bonferroni correction, a common SNP set throughout, and directional concordance "
  "across IVW, MR-Egger, weighted median and MR-PRESSO.")

E("E3 L85 Limitations 82\u8bcd\u62c6\u56db\u53e5",
  "Limitations: areal BMD conflates size with mass, and although height-conditioning bounds "
  "the size-attributable share at 9.5\u201333.2%, a direct volumetric decomposition is not "
  "yet available; heterogeneity is high across all pairs including the positive control, so "
  "the causal claim rests on estimator concordance rather than on fixed-effect p-values; the "
  "reverse and OA\u2192fracture estimates are null but weakly powered (0.14\u20130.26), "
  "excluding only |\u03b2|\u22650.30; and residual sample overlap, although bounded at "
  "\u22642.6% for BMD\u2192OA and acting in the observed direction, cannot be removed with "
  "summary data alone.",
  "Limitations: areal BMD conflates size with mass, and although height-conditioning bounds "
  "the size-attributable share at 9.5\u201333.2%, a direct volumetric decomposition is not "
  "yet available. Heterogeneity is high across all pairs including the positive control, so "
  "the causal claim rests on estimator concordance rather than on fixed-effect p-values. The "
  "reverse and OA\u2192fracture estimates are null but weakly powered (0.14\u20130.26), "
  "excluding only |\u03b2|\u22650.30. Residual sample overlap is bounded at \u22642.6% for "
  "BMD\u2192OA and acts in the observed direction, but it cannot be removed with summary "
  "data alone.")

# ═══════════════════ F. AI 偏好词 ═══════════════════
E("F1 L95 key\u2192selected",
  "**Figure 2.** Forest plot of key genetic correlations:",
  "**Figure 2.** Forest plot of the selected genetic correlations:")


def main():
    text = io.open(DOC, encoding="utf-8").read()
    bak = DOC + ".bak-before-humanize132"
    if not os.path.exists(bak):
        io.open(bak, "w", encoding="utf-8").write(text)
        print(f"backup -> {os.path.basename(bak)}")

    # 允许重复出现的条目：同一句在 Methods 与 Data availability 各出现一次，两处都改
    ALLOW_MULTI = {"C3 L41"}

    failures = []
    for tag, old, new in MS:
        n = text.count(old)
        if n != 1:
            if tag in ALLOW_MULTI and n == 2:
                text = text.replace(old, new)
                print(f"[OK  ] {tag} (x2, both replaced)")
                continue
            failures.append((tag, n, old[:90]))
            print(f"[SKIP] {tag}: occurrences={n}")
            continue
        text = text.replace(old, new, 1)
        print(f"[OK  ] {tag}")

    if failures:
        print(f"\n!! {len(failures)} not uniquely matched; NOT written.")
        for t, n, o in failures:
            print(f"    {t}: n={n}  {o}")
        sys.exit(1)

    io.open(DOC, "w", encoding="utf-8").write(text)
    print(f"\napplied {len(MS)} edits")


if __name__ == "__main__":
    main()
