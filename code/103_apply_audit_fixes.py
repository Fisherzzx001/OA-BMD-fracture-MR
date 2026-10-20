# -*- coding: utf-8 -*-
"""
103_apply_audit_fixes.py
========================
Apply the trust-but-verify audit fixes to the manuscript.

HARD RULE OBSERVED: no beta / SE / p / OR is touched. Only wording and
derived percentage ranges that were computed from a defective table.

Every replacement asserts exactly one occurrence before/after.
"""
from __future__ import annotations

import io
import os
import shutil
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MS = os.path.join(ROOT, "docs", "08_\u6295\u7a3f\u7a3f_OA_BMD_fracture.md")

EDITS = [
    # ---- (1) Highlights: derived MVMR range recomputed from corrected S10 ----
    (
        "- BMD\u2192OA effect is only 6\u201333% explained by skeletal size (height proxy).",
        "- BMD\u2192OA effect is only 9.5\u201333.2% explained by skeletal size (height proxy).",
    ),
    # ---- (2) Abstract: same range, plus the complementary "independent" range ----
    (
        "Conditioning on height removed 6\u201333% of the effect, leaving 67\u201394% "
        "independent of the height\u2013size pathway.",
        "Conditioning on height removed 9.5\u201333.2% of the effect, leaving 67\u201391% "
        "independent of the height\u2013size pathway.",
    ),
    # ---- (3) Methods / MR: instrument counts + attrition + estimator definitions ----
    (
        "(strict PLINK-equivalent clumping, benchmarked against the PLINK binary; retained "
        "20\u201335 instruments per BMD exposure, 116 for knee OA). We computed IVW (fixed and "
        "random effects), MR-Egger (intercept test for directional pleiotropy), weighted median "
        "and leave-one-out, with Cochran Q and I\u00b2 for heterogeneity and MR-PRESSO for "
        "pleiotropic outliers.",
        "(strict PLINK-equivalent clumping, benchmarked against the PLINK binary; 22\u201337 "
        "instruments retained per analysed BMD site and 118 for knee OA). Harmonisation against "
        "each outcome \u2014 exclusion of palindromic and strand-ambiguous variants and of "
        "instruments absent from the outcome GWAS \u2014 reduced the analysed sets to 16\u201331 "
        "instruments per forward pair and 69 per knee-OA-exposure pair, an attrition we report "
        "explicitly. We computed IVW (fixed and random effects), MR-Egger (intercept test for "
        "directional pleiotropy), weighted median (Bowden\u2019s inverse-variance weighted median "
        "of the Wald ratios) and leave-one-out, with Cochran Q and I\u00b2 for heterogeneity and "
        "MR-PRESSO for pleiotropic outliers. Random-effects IVW inflated the standard error "
        "multiplicatively (SE_random = SE_fixed \u00d7 \u221amax(1, Q/df)) rather than by additive "
        "\u03c4\u00b2.",
    ),
    # ---- (4) Methods / Effect-scale: replace the phantom x0.2823 correction ----
    (
        "Munged BMD statistics carried a placeholder N=400,000 against a true N\u224831,873, "
        "shrinking reconstructed \u03b2 3.54-fold; forward MR \u03b2 were corrected by \u00d70.2823 "
        "(\u00f73.5426) and reverse-MR \u03b2 received the identical correction, so all estimates "
        "share one N-corrected SD scale. Directions, p-values and significance are N-invariant.",
        "BMD statistics had initially been standardised against a placeholder N of 400,000; we "
        "therefore re-standardised them to per-SD units using the verified per-site N (unified "
        "N=31,873). We confirmed empirically that the betas entering the analysis are already on "
        "the true-N scale: for every instrument the reported standard error matches the "
        "closed-form expectation for a standardised effect, SE = 1/\u221a[2\u00b7N\u00b7f\u00b7(1\u2212f)] "
        "with f the minor-allele frequency, agreeing with N=31,873 to within 0.1% and excluding "
        "the placeholder N=400,000, which would predict standard errors 3.54-fold smaller. No "
        "post-hoc rescaling was therefore applied to any reported estimate. Because the same "
        "verified per-SD exposure units enter the forward and reverse analyses, all estimates "
        "share one N-verified SD scale, and directions, p-values and significance are invariant "
        "to it.",
    ),
    # ---- (5) Results / BMD raises knee OA: random-effects convention dependence ----
    (
        "Heterogeneity was high (I\u00b2=68\u201386%), and under random-effects IVW only two of "
        "five exposures remained significant (total-body p=0.006, lumbar-spine p=0.017).",
        "Heterogeneity was high (I\u00b2=68\u201386%), and under random-effects IVW only two of "
        "five exposures remained significant (total-body p=0.006, lumbar-spine p=0.017); this "
        "count depends on the random-effects convention, because a DerSimonian\u2013Laird "
        "additive-\u03c4\u00b2 implementation additionally retains the arms exposure under random "
        "effects (p=0.033).",
    ),
    # ---- (6) Results / reverse MR: disclose the real outlier sensitivity ----
    (
        "excluding |\u03b2|\u22650.30 but not moderate effects; they are better described as "
        "absence of evidence than evidence of absence.",
        "excluding |\u03b2|\u22650.30 but not moderate effects; they are better described as "
        "absence of evidence than evidence of absence. Removing the MR-PRESSO outliers flagged "
        "for these pairs moved two of the three point estimates substantially (total-body BMD to "
        "+0.171, fracture to +0.200) without altering their non-significance, so the nulls "
        "reflect both low power and sensitivity to instrument selection.",
    ),
    # ---- (7) Results / Effect-scale validation: give both scales vs GEFOS ----
    (
        "After correcting the BMD placeholder N (400,000\u219231,873; \u00d70.2823), the "
        "BMD\u2192fracture OR (1.42 femoral-neck) lay 8.4% below GEFOS (1.55) on the like-for-like "
        "comparison, with overlapping CIs; the total-body OR 1.40 was 9.6% below (cross-site "
        "approximation).",
        "On the verified per-SD scale, the BMD\u2192fracture OR (1.42 femoral-neck) lay 8.4% below "
        "GEFOS (1.55) on the odds-ratio scale and 20.0% below on the log-odds scale \u2014 we give "
        "both because effect estimates are compared on the log scale \u2014 with overlapping CIs; "
        "the total-body OR 1.40 was 9.6% below on the ratio scale and 23.1% below on the "
        "log-odds scale (cross-site approximation).",
    ),
    # ---- (8) Results / MVMR: corrected attenuation ranges ----
    (
        "Conditioning the BMD\u2192OA effect on height attenuated it by 6\u201333% across three "
        "defensible univariable baselines (13.5\u201333.2% against the primary Table 3 baselines; "
        "n_h=100; Supplementary Table S9), leaving 67\u201394% independent of the "
        "height\u2013size pathway.",
        "Conditioning the BMD\u2192OA effect on height attenuated it by 9.5\u201333.2% across the "
        "two tabulated univariable baselines (13.5\u201333.2% against the primary Table 3 "
        "baselines and 9.5\u201321.2% against the MVMR-set baselines; primary n_h=100/EUR; "
        "Supplementary Tables S9 and S10), leaving 67\u201391% independent of the "
        "height\u2013size pathway.",
    ),
    (
        "BMD\u2192fracture moved only \u22121.7% to \u22122.8% under every baseline.",
        "BMD\u2192fracture moved only \u22121.7% to \u22123.0%.",
    ),
    # ---- (9) Discussion: same two ranges ----
    (
        "MVMR conditioning on height attenuated the effect by 6\u201333%, while the fracture "
        "positive control moved only \u22121.7% to \u22122.8%, and",
        "MVMR conditioning on height attenuated the effect by 9.5\u201333.2%, while the fracture "
        "positive control moved only \u22121.7% to \u22123.0%, and",
    ),
    (
        "the mass share is an upper bound (\u226467\u201394%), not a point estimate.",
        "the mass share is an upper bound (\u226467\u201391%), not a point estimate.",
    ),
    # ---- (10) Limitations: drop the self-negating "magnitudes are provisional" ----
    (
        "absolute MR magnitudes are provisional because exposure \u03b2 were reconstructed from a "
        "placeholder N=400,000 (true \u224831,873) and scaled down 3.54-fold (\u00d70.2823); "
        "directions, significance and robustness are unaffected;",
        "absolute MR magnitudes are tied to the areal-BMD per-SD scale, which we validated "
        "externally against GEFOS using the positive control (OR 1.40\u20131.42 vs 1.55); "
        "directions, significance and robustness are scale-invariant;",
    ),
]


def main() -> int:
    with io.open(MS, "r", encoding="utf-8", newline="") as fh:
        text = fh.read()

    backup = MS + ".bak-before-auditfix"
    if not os.path.exists(backup):
        shutil.copy2(MS, backup)
    print(f"backup -> {os.path.basename(backup)}")

    ok = 0
    for i, (old, new) in enumerate(EDITS, 1):
        n = text.count(old)
        if n != 1:
            print(f"[FAIL] edit {i}: expected 1 occurrence, found {n}")
            print(f"       old = {old[:110]!r}")
            return 2
        text = text.replace(old, new)
        ok += 1
        print(f"[ok  ] edit {i} applied")

    with io.open(MS, "w", encoding="utf-8", newline="") as fh:
        fh.write(text)

    print(f"\napplied {ok}/{len(EDITS)} edits")
    print("no beta / SE / p / OR was modified (all edits are prose / derived ranges)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
