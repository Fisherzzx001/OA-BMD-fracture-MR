# -*- coding: utf-8 -*-
"""
116 · 把新图 S9–S12 接入稿件：图例、补充材料清单、正文引用
每处 text.count(old)==1 断言。
"""
import io
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DOC = os.path.normpath(os.path.join(HERE, "..", "docs", "08_投稿稿_OA_BMD_fracture.md"))

with io.open(DOC, encoding="utf-8") as f:
    t = f.read()

edits = [
    # ---- 1. Methods：工具筛选引用 S9
    ("47 of the 116 knee-OA instruments (41%) were lost at this "
     "step, an attrition we report explicitly.",
     "47 of the 116 knee-OA instruments (41%) were lost at this "
     "step, an attrition we report explicitly (Supplementary Fig. S9)."),

    # ---- 2. Results · 稳健性段引用 S11
    ("We therefore base the causal reading on the concordance of pleiotropy-robust "
     "estimators rather than on fixed-effect p-values; all five effects remained "
     "positively causal after MR-PRESSO outlier removal, with four of five strengthening.",
     "We therefore base the causal reading on the concordance of pleiotropy-robust "
     "estimators rather than on fixed-effect p-values; all five effects remained "
     "positively causal after MR-PRESSO outlier removal, with four of five strengthening. "
     "A multi-estimator robustness overview of all ten pairs is given in "
     "Supplementary Fig. S11."),

    # ---- 3. Results · MVMR 段引用 S10
    ("(13.5\u201333.2% against the primary Table 3 baselines and 9.5\u201321.2% against "
     "the MVMR-set baselines; primary n_h=100/EUR; Supplementary Tables S9 and S10)",
     "(13.5\u201333.2% against the primary Table 3 baselines and 9.5\u201321.2% against "
     "the MVMR-set baselines; primary n_h=100/EUR; Supplementary Tables S9 and S10 "
     "and Supplementary Fig. S10)"),

    # ---- 4. Discussion · 文献对照引用 S12
    ("Funck-Brentano et al. [36] likewise reported positive effects at both joints "
     "(knee OR 1.20, 1.05\u20131.37; hip OR 1.23, 1.12\u20131.36).",
     "Funck-Brentano et al. [36] likewise reported positive effects at both joints "
     "(knee OR 1.20, 1.05\u20131.37; hip OR 1.23, 1.12\u20131.36); a quantitative "
     "comparison on the native scales is given in Supplementary Fig. S12."),

    # ---- 5. 图例：补 S9–S12
    ("**S7.** Pleiotropy and heterogeneity diagnostics. **S8.** Colocalisation of "
     "BMD with knee OA and its power. **Graphical Abstract.**",
     "**S7.** Pleiotropy and heterogeneity diagnostics. **S8.** Colocalisation of "
     "BMD with knee OA and its power. **S9.** Instrument selection and attrition flow "
     "(STROBE-MR), tallied from the deposited clump and harmonised instrument tables. "
     "**S10.** Multivariable MR decomposition of the BMD\u2192OA effect conditional on "
     "height: (a) univariable vs MVMR estimates; (b) attenuation percentages with the "
     "BMD\u2192fracture positive control. **S11.** Multi-estimator robustness grid for "
     "the ten MR pairs (IVW fixed/random, MR-Egger, weighted median, PRESSO "
     "outlier-removed, MRlap-corrected, MVMR). **S12.** Quantitative comparison with "
     "published MR estimates (Hartley 2022; Funck-Brentano 2018; Zhou 2025) on their "
     "native scales. **Graphical Abstract.**"),

    # ---- 6. Supplementary Material 清单：补 4 项 + 范围 S1–S12
    ("- **Figure S8.** Colocalisation of BMD with knee OA and its power.\n",
     "- **Figure S8.** Colocalisation of BMD with knee OA and its power.\n"
     "- **Figure S9.** Instrument selection and attrition flow (STROBE-MR).\n"
     "- **Figure S10.** Multivariable MR decomposition of the BMD\u2192OA effect "
     "conditional on height.\n"
     "- **Figure S11.** Multi-estimator robustness grid for the ten MR pairs.\n"
     "- **Figure S12.** Quantitative comparison with published MR estimates.\n"),
    ("Supplementary figures S1\u2013S8 and the Graphical Abstract are provided as "
     "high-resolution PNG",
     "Supplementary figures S1\u2013S12 and the Graphical Abstract are provided as "
     "high-resolution PNG"),
]

for i, (old, new) in enumerate(edits, 1):
    c = t.count(old)
    assert c == 1, f"edit {i}: expected 1 occurrence, found {c}\n---\n{old[:120]}"
    t = t.replace(old, new)
    print(f"[ok] edit {i}")

with io.open(DOC, "w", encoding="utf-8") as f:
    f.write(t)
print("\nWrote", DOC)
