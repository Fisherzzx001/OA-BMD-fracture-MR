#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
446_n3n7_remediation.py -- N3/N4/N5/N6/N7 (MAJOR) second-pass remediation.

Idempotent. Fixes the residual statistical-definition gaps left by code/444:

  N3  PP.H4>=0.8 stated as a posterior rule, not error-rate control (Methods);
      the paired-Wilcoxon p of the 9/31 vs 0/31 split labelled unadjusted /
      exploratory (it is a frequentist test outside both multiplicity families).
  N4  "directional concordance" removed from Strengths; legs WM CI crossing zero
      propagated to Results; Abstract no longer says the causal reading *rests on*
      pleiotropy-robust estimators.
  N5  PRESSO "removing outliers did not reverse any direction" no longer offered
      as robustness; S11/S20 figure captions demoted from "robustness grid" to
      "comparison grid (sign != magnitude)".
  N6  MDE (=2.8 SE, a design-sensitivity quantity) no longer used as an exclusion
      bound; the exclusion limit is the 95% CI upper bound 0.0897/0.0937/0.1047.
      Post-hoc (observed) power flagged as circular; fixed-effect SE caveat;
      forward/reverse beta units explicitly incommensurable; Steiger's
      "minimum detectable effect" no longer confused with an effect-size MDE.
  N7  eight-of-ten negative lambda disclosed in Methods; quantitative
      adjudication (Pearson r=-0.23 between lambda and correction size) showing
      the away-from-null shift is driven by weak-instrument/winner's-curse
      correction, not by the lambda sign.

Hard constraints: no reported beta/SE/p/OR/CI changed; no statistics re-run;
body word count must not increase (3993/4000 md, 3998/4000 docx); the 310 lock
string "achieved power is computed post-hoc" is preserved verbatim.

Usage: python3 code/446_n3n7_remediation.py
"""
import os
import re
import sys
import shutil

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "docs", "08_投稿稿_OA_BMD_fracture.md")
MIRROR = os.path.join(ROOT, "Manuscript_Ready", "01_Manuscript.md")
DOCXS = [
    os.path.join(ROOT, "submission_package", "Manuscript.docx"),
    os.path.join(ROOT, "submission_package", "Submission_MainText.docx"),
]
SUPP_PAIRS = [
    (os.path.join(ROOT, "docs", "Supplementary", "Table_S17_sensitivity.md"),
     os.path.join(ROOT, "Manuscript_Ready", "Supplementary", "Table_S17_sensitivity.md")),
    (os.path.join(ROOT, "docs", "Supplementary", "Table_S14_ancestry.md"),
     os.path.join(ROOT, "Manuscript_Ready", "Supplementary", "Table_S14_ancestry.md")),
]

# --------------------------------------------------------------- body edits
# (name, old, new). Net word delta is NEGATIVE (computed below).
BODY = [
 ("N4-abstract",
  "so the causal reading rests on pleiotropy-robust estimators",
  "so pleiotropy-robust estimators bound the causal reading"),
 ("N3-L59-posterior-rule",
  "**Colocalisation.** Bayesian colocalisation (coloc.abf; Giambartolomei et al. 2014 [22]) tested the three a priori loci (WNT16, SOX6, GDF5) and all 31 BMD instrument loci (±500 kb; PP.H4≥0.8), with a BMD–fracture control.",
  "**Colocalisation.** Bayesian coloc.abf [22] tested WNT16, SOX6, GDF5 and 31 BMD instrument loci (±500 kb), with a BMD–fracture control; PP.H4≥0.8 is a posterior rule, not error-rate control."),
 ("N3-L79-unadjusted-p",
  "(median 0.153; paired Wilcoxon one-sided p=0.0043, exact two-sided p=0.0086)",
  "(median 0.153; unadjusted exploratory paired Wilcoxon one-sided p=0.0043, exact two-sided p=0.0086)"),
 ("N4-L69-legs-CI",
  "weighted medians matched the IVW sign but were attenuated (~one-third at arms/legs; S17e)",
  "weighted medians matched the IVW sign but shrank ~one-third at arms/legs (legs CI spans zero; S17e)"),
 ("N4/N5-L99-strength",
  "dual FDR/Bonferroni correction and directional concordance across IVW, MR-Egger, weighted median and MR-PRESSO, and cross-site homogeneity",
  "dual FDR/Bonferroni correction, sign-only estimator agreement with attenuated magnitudes (S17e), and cross-site homogeneity"),
 ("N5-L75-presso",
  "and identified 0–3 outliers per pair; removing outliers did not reverse any direction, and strict-clumped estimates reproduced the conservative 1-Mb clumping.",
  "and identified 0–3 outliers per pair; removing them preserved signs, not robustness against systematic pleiotropy, and strict-clumped estimates reproduced the 1-Mb clumping."),
 ("N6-L71-exclusion-bound",
  "these were weakly powered (0.14–0.26; Table 3) and excluded only |β|≥0.09 in SD of BMD per unit OA liability.",
  "these were weakly powered (0.14–0.26; Table 3) and exclude only |β| above their 95% CI upper bounds (0.090–0.105; Table 3 note)."),
 ("N6-L73-steiger-not-MDE",
  "with 69 instruments and a null reverse estimate the test is underpowered (minimum detectable effect exceeds the null), not evidence BMD drives OA.",
  "with 69 instruments and a null reverse estimate Steiger cannot rank exposure–outcome variance explained (not an effect-size MDE), not evidence BMD drives OA."),
 ("N6-L101-conclusion-bound",
  "the reverse estimates exclude only |β|≥0.09; the positive knee-OA × fracture r_g is genetic correlation, not causation",
  "the reverse estimates exclude only |β| above 0.090–0.105; the positive knee-OA × fracture r_g is genetic correlation, not causation"),
 ("N7-L55-negative-lambda",
  "participant overlap is accommodated through its internally estimated cross-trait LDSC intercept λ, with no overlap fraction supplied (configuration and default-setting sensitivity in Table S14).",
  "participant overlap enters through its cross-trait LDSC intercept λ (eight of ten negative, |λ|≤0.037; S14), with no overlap fraction supplied."),
]

# ---------------------------------------------- non-counted blocks (after ## Figure legends)
TABLE3_NOTE = (
 "MDE = (1.96 + 0.84) × SE is the minimum |β| detectable at 80% power "
 "(two-sided α=0.05); achieved power is computed post-hoc at the observed "
 "point estimate and should be read alongside the MDE, not as prospective "
 "power. |β|/SE is N-scale-invariant. K = strict instrument count. Achieved "
 "power refers to fixed-effect IVW; under the primary random-effects IVW only "
 "two of five forward BMD→knee-OA estimates remain significant (Table S17-b; §3)."
)
TABLE3_ADD = (
 " **MDE, exclusion limits and post-hoc power are three different quantities "
 "(N6 adjudication; Table S17-f/g).** MDE = (1.96+0.84) × SE is a "
 "*design-sensitivity* quantity: were the true effect equal to the MDE, a "
 "two-sided α=0.05 test would reject with 80% probability. It is not the bound "
 "on what a null result excludes. That bound is the 95% CI upper limit "
 "|β| + 1.96 SE, which for the three reverse pairs is 0.0897 (OA→fracture), "
 "0.0937 (OA→BMD_total) and 0.1047 (OA→BMD_femoral_neck) SD of BMD per unit OA "
 "liability — all at or above the MDE (0.0896), and 16.9% above it at the "
 "femoral-neck site. Quoting 0.09 as the exclusion limit therefore understates "
 "the moderate effects still compatible with the data. Achieved (observed) power "
 "is Φ(|z|−1.96)+Φ(−|z|−1.96) with z the *observed* |β|/SE; it is a "
 "deterministic function of the same p value and is necessarily low whenever a "
 "result is non-significant, so it cannot independently establish that a null "
 "was underpowered (the observed-power fallacy). It is reported for transparency "
 "only. Both MDE and achieved power use fixed-effect IVW SE; under the primary "
 "random-effects convention the SE, and hence the MDE, is larger. Forward β "
 "(log-odds of OA per SD BMD) and reverse β (SD of BMD per unit OA liability) "
 "have different units, so their MDEs are not commensurable and must not be read "
 "on one scale; the forward–reverse power asymmetry is driven by exposure "
 "variance explained (R²_exposure 0.034–0.063 forward versus 0.0003 reverse, "
 "reflecting h² 0.26–0.39 versus 0.035), not by instrument count."
)

S11_OLD = ("**S11.** Multi-estimator robustness grid for the ten MR pairs "
           "(IVW fixed/random, MR-Egger, weighted median, PRESSO outlier-removed, "
           "MRlap-corrected, MVMR).")
S11_NEW = ("**S11.** Multi-estimator comparison grid for the ten MR pairs "
           "(IVW fixed/random, MR-Egger, weighted median, PRESSO outlier-removed, "
           "MRlap-corrected, MVMR); sign concordance is not magnitude concordance "
           "(Table S17-e/g).")

S20_OLD = ("All ten pairs retain a consistent direction across every estimator; "
           "the three reverse/null pairs (C) stay null and heterogeneous under all five.")
S20_NEW = ("All ten pairs retain a consistent direction across every estimator; "
           "the three reverse/null pairs (C) stay null and heterogeneous under all "
           "five. Sign agreement is descriptive, not robustness evidence: weighted-"
           "median magnitudes fall to ~one-third of IVW at arms/legs with the legs "
           "bootstrap interval crossing zero, and MR-PRESSO removes detected "
           "outliers only, not systematic pleiotropy (9/10 distortion tests "
           "significant). See Table S17-e/g.")

S14NOTE_ANCHOR = (
 "*Sensitivity to the default.* All ten |λ| are small and of mixed sign, so the "
 "internally estimated overlap signal is minor relative to the reported effects. "
 "Because MRlap's default does not take a user-specified overlap fraction, we "
 "assess sensitivity against the independently tabulated counts: the Burgess "
 "bounds above (≤2.6–8.4% overlap; biases ≤0.7–4.5%) bound the away-from-null "
 "component, and every corrected estimate retained its direction (10/10; median "
 "+20.5% away from the null; `code/output/mrlap_direction_verification.txt`), "
 "so no headline conclusion depends on the default."
)
S14NOTE_ADD = (
 " *Adjudication of the negative intercepts (N7).* λ is a cross-trait LDSC "
 "intercept, not a sample-overlap fraction, and its sign does not map "
 "monotonically onto overlap-bias direction. Eight of ten λ are negative yet all "
 "ten corrected estimates moved away from the null; across the ten pairs λ and "
 "the relative correction are essentially uncorrelated (Pearson r = −0.23): the "
 "most negative λ (−0.0361, BMD_total→fracture) shows a +15.3% shift, while a "
 "near-zero λ (+0.0010, BMD_lumbar→OA_knee) shows +21.8%. The away-from-null "
 "shift is therefore attributable to MRlap's weak-instrument and winner's-curse "
 "correction (§2), not to the sign of λ. Negative values are compatible with "
 "sampling variation around a near-zero intercept and with assortative-mating / "
 "dynastic effects; they neither demonstrate negative genetic sharing nor "
 "license the claim that no positive overlap-driven inflation exists. Overlap "
 "bias is bounded instead by the independently tabulated UKB counts and Burgess "
 "bounds above, which do not depend on λ."
)

INDEX_OLD = ("- **Table S17.** Sensitivity analyses requested at review: (a) uncertainty "
             "of the effect-scale reconstruction constant K; (b) significance counts "
             "under three random-effects conventions; (c) change in the IVW estimate "
             "on removal of MR-PRESSO outliers for all ten pairs; (d) colocalisation "
             "counts as a function of the PP.H4 threshold.")
INDEX_NEW = ("- **Table S17.** Sensitivity analyses requested at review: (a) uncertainty "
             "of the effect-scale reconstruction constant K; (b) significance counts "
             "under three random-effects conventions; (c) change in the IVW estimate "
             "on removal of MR-PRESSO outliers for all ten pairs; (d) colocalisation "
             "counts as a function of the PP.H4 threshold; (e) pleiotropy-robust "
             "estimator attenuation and bootstrap-CI asymmetry; (f) reverse-MR power "
             "and the MDE / exclusion-limit distinction; (g) posterior thresholds, "
             "sign concordance and PRESSO outlier removal as non-evidence.")

# ------------------------------------------------------------ supplement edits
S17E_OLD = (
 "The weighted-median bootstrap intervals are visibly right-skewed (point "
 "estimates near the lower bound, e.g. arms +0.026 against 0.019–0.217): Bowden "
 "weights (β_exposure/SE_outcome)² are highly uneven with 16–31 instruments, so "
 "a few strong instruments dominate and the percentile bootstrap distribution "
 "has a long right tail. The asymmetry is a small-instrument-number property of "
 "the estimator, not evidence of a second effect direction; we therefore report "
 "the intervals as computed and do not force symmetry."
)
S17E_ADD = (
 " **Inferential consequence (N4).** A percentile-bootstrap interval is not "
 "required to be centred on its point estimate, so a point estimate near the "
 "lower bound is expected here rather than a defect of the computation; a "
 "symmetric Wald interval built around the weighted median would be misleading. "
 "The substantive consequence is that the legs interval (−0.003 to +0.132) "
 "includes zero while IVW is significant (p=2.3×10⁻⁵): at that site the most "
 "pleiotropy-robust single estimator does not support a non-null effect. Sign "
 "agreement across estimators therefore establishes only that no estimate "
 "reverses direction; it establishes neither stability of magnitude nor "
 "stability of significance, and is not offered as robustness evidence "
 "(Fig. S11/S20 are comparison grids)."
)

S17F_OLD = (
 "The reverse MDE (minimum detectable |β|=0.089–0.090 in Table 3) is computed "
 "from these instruments' actual SE (≈0.032); the quoted exclusion bound "
 "|β|≥0.09 SD therefore corresponds to this MDE and is not an arbitrary "
 "cut-off. Forward and reverse estimates are on the same as-deposited scale, so "
 "the power asymmetry is biological (h²), not a scale artefact."
)
S17F_NEW = (
 "The reverse MDE (0.089–0.090 in Table 3) is 2.8 × the fixed-effect IVW SE "
 "(≈0.032) and is a **design-sensitivity** quantity, not an exclusion bound: "
 "the exclusion limit is the 95% CI upper bound |β| + 1.96 SE = 0.0897 / 0.0937 "
 "/ 0.1047 for OA→fracture / OA→BMD_total / OA→BMD_femoral_neck, all at or "
 "above the MDE and 16.9% above it at the femoral-neck site. Achieved power is "
 "post-hoc (observed) power, a deterministic function of the observed p value "
 "that cannot independently explain a null; see the Table 3 note for the full "
 "adjudication. Forward and reverse β are **not** on one scale — forward β is "
 "log-odds of OA per SD BMD, reverse β is SD of BMD per unit OA liability "
 "(the reverse exposure is a binary liability, not a 1-SD continuous trait) — "
 "so forward and reverse MDEs are not commensurable. The power asymmetry is "
 "driven by exposure variance explained (low OA h²), not by instrument count "
 "and not by a scale difference."
)

S17G_ANCHOR = "### S17-f. Why reverse MR has more instruments yet lower power"
S17G_BLOCK = (
 "### S17-g. Three statistical distinctions behind the N3/N4/N5/N6 remediation\n\n"
 "**(i) A Bayesian posterior threshold is not a multiplicity control (N3).** "
 "PP.H4≥0.8 is a locus-level decision rule on a posterior probability under "
 "fixed priors (p1=p2=1×10⁻⁴, p12=5×10⁻⁶); it controls no frequentist FWER or "
 "FDR across the 31 loci. Multiplicity control in this study is confined to the "
 "two analysis-defined families in §5 (the 78 LDSC pairs and the ten primary MR "
 "pairs). Accordingly the 9/31-versus-0/31 endpoint split and its paired "
 "Wilcoxon p value (one-sided 0.0043, exact two-sided 0.0086) are unadjusted, "
 "exploratory descriptive summaries, not multiplicity-adjusted inference, and "
 "are labelled as such in §2/§3. No colocalisation claim is a confirmatory "
 "test.\n\n"
 "**(ii) Sign concordance is not robustness evidence (N4/N5).** Weighted-median "
 "point estimates shrink to ~one-third of IVW at arms/legs, the legs bootstrap "
 "interval crosses zero, and MR-PRESSO removes only the 0–3 detected outliers "
 "per pair while 9/10 global distortion tests remain significant. Outlier "
 "removal therefore addresses detectable outlying variants, not systematic "
 "(balanced) horizontal pleiotropy, and is not offered as robustness evidence. "
 "The multi-estimator grids (Figs S11/S20) show direction agreement only.\n\n"
 "**(iii) MDE, exclusion limit and post-hoc power are distinct (N6).** MDE is "
 "design sensitivity at 80% power; the exclusion limit of a null result is the "
 "95% CI upper bound; achieved power is post-hoc and circular as an explanation "
 "of non-significance. These are quantified in the Table 3 note and §S17-f. "
 "No reported β/SE/p/OR/CI was changed by this remediation.\n\n"
)

S14_SUPP_OLD = (
 "A negative cross-trait intercept (BMD_total→fracture, λ=−0.0361) runs "
 "opposite to the positive bias expected from sample overlap; in cross-trait "
 "LDSC such negative values arise from sampling variation around a near-zero "
 "intercept and from genetic nurture/assortative mating, not from negative "
 "genetic sharing. Because the sign is mixed and |λ| small, no positive "
 "overlap-driven inflation needs to be subtracted here, and the negative value "
 "does not contradict the conclusion that residual overlap is minor."
)
S14_SUPP_NEW = (
 "A negative cross-trait intercept (BMD_total→fracture, λ=−0.0361) runs "
 "opposite to the positive bias expected from sample overlap; in cross-trait "
 "LDSC such negative values are compatible with sampling variation around a "
 "near-zero intercept and with assortative-mating / dynastic (genetic nurture) "
 "effects, not with negative genetic sharing. **Correction of an earlier "
 "over-reading (N7):** λ is an intercept, not a sample-overlap fraction, so "
 "its sign does not license the inference that no positive overlap-driven "
 "inflation exists. Eight of ten λ are negative yet all ten corrected estimates "
 "moved away from the null, and λ is essentially uncorrelated with the relative "
 "correction across the ten pairs (Pearson r = −0.23; the most negative λ, "
 "−0.0361, shows only a +15.3% shift, while a near-zero λ, +0.0010, shows "
 "+21.8%). The away-from-null displacement is therefore attributable to MRlap's "
 "weak-instrument and winner's-curse correction rather than to the λ sign. "
 "Overlap-bias magnitude is bounded by the independently tabulated UKB counts "
 "and the Burgess bounds below, which do not depend on λ; the negative values "
 "do not contradict the conclusion that residual overlap is minor."
)


def words(s):
    return len([x for x in s.split() if re.search(r"[^\W_]", x, re.U)])


def edit_text(path, triples, label):
    txt = open(path, encoding="utf-8").read()
    orig = txt
    out = []
    for name, old, new in triples:
        if new in txt and old not in txt:
            out.append(f"  {name}: SKIP (already applied)")
            continue
        n = txt.count(old)
        if n == 0:
            out.append(f"  {name}: MISSING anchor")
            continue
        txt = txt.replace(old, new, 1)
        out.append(f"  {name}: APPLIED ({n} match)")
    if txt != orig:
        open(path, "w", encoding="utf-8", newline="\n").write(txt)
    print(f"[{label}] {os.path.relpath(path, ROOT)}")
    print("\n".join(out))
    return txt


def edit_docx(path, triples, label):
    from docx import Document
    d = Document(path)
    out = []
    changed = False
    for name, old, new in triples:
        hits = 0
        for p in d.paragraphs:
            full = "".join(r.text for r in p.runs)
            if old in full and new not in full:
                merged = full.replace(old, new, 1)
                if p.runs:
                    p.runs[0].text = merged
                    for r in p.runs[1:]:
                        r.text = ""
                else:
                    p.add_run(merged)
                hits += 1
                changed = True
        if hits:
            out.append(f"  {name}: APPLIED ({hits} para)")
        elif new in "".join(p.text for p in d.paragraphs):
            out.append(f"  {name}: SKIP (already applied)")
        else:
            out.append(f"  {name}: MISSING anchor")
    if changed:
        d.save(path)
    print(f"[{label}] {os.path.relpath(path, ROOT)}")
    print("\n".join(out))


def main():
    os.chdir(ROOT)
    src0 = open(SRC, encoding="utf-8").read()

    # --- word-budget guard (must not increase) ---
    delta = sum(words(n) - words(o) for _, o, n in BODY)
    print(f"[budget] body word delta = {delta:+d}  (3993 -> {3993+delta} md, "
          f"3998 -> {3998+delta} docx)")
    if delta > 0:
        raise SystemExit("ABORT: body edits would increase word count")

    # --- 1. source manuscript: body + non-counted blocks ---
    extra = [
        ("N6-Table3-note", TABLE3_NOTE, TABLE3_NOTE + TABLE3_ADD),
        ("N4/N5-S11-caption", S11_OLD, S11_NEW),
        ("N4/N5-S20-caption", S20_OLD, S20_NEW),
        ("N7-S14-note", S14NOTE_ANCHOR, S14NOTE_ANCHOR + S14NOTE_ADD),
        ("S17-index-e/f/g", INDEX_OLD, INDEX_NEW),
    ]
    edit_text(SRC, BODY + extra, "docs/08")

    # --- 2. byte mirror ---
    shutil.copyfile(SRC, MIRROR)
    print(f"[mirror] byte copy -> {os.path.relpath(MIRROR, ROOT)}")

    # --- 3. both docx (run-level merge) ---
    for p in DOCXS:
        edit_docx(p, BODY + extra, "docx")

    # --- 4. supplementary files (two copies each) ---
    s17 = [
        ("N4-S17e-consequence", S17E_OLD, S17E_OLD + S17E_ADD),
        ("N6-S17f-scale-correction", S17F_OLD, S17F_NEW),
    ]
    s14 = [("N7-S14-adjudication", S14_SUPP_OLD, S14_SUPP_NEW)]
    for a, b in SUPP_PAIRS:
        for path, triples in ((a, s17 if "S17" in a else s14),
                              (b, s17 if "S17" in b else s14)):
            edit_text(path, triples, "supp")
        # append S17-g once, guarded by its heading
        for path in ((a, b) if "S17" in a else ()):
            txt = open(path, encoding="utf-8").read()
            if "### S17-g." not in txt:
                txt = txt.rstrip("\n") + "\n\n" + S17G_BLOCK
                open(path, "w", encoding="utf-8", newline="\n").write(txt)
                print(f"[supp] S17-g APPENDED -> {os.path.relpath(path, ROOT)}")
            else:
                print(f"[supp] S17-g SKIP -> {os.path.relpath(path, ROOT)}")

    # --- 5. consistency checks ---
    a = open(SUPP_PAIRS[0][0], encoding="utf-8").read()
    b = open(SUPP_PAIRS[0][1], encoding="utf-8").read()
    print("[check] S17 dual copies identical:", a == b)
    a = open(SUPP_PAIRS[1][0], encoding="utf-8").read()
    b = open(SUPP_PAIRS[1][1], encoding="utf-8").read()
    print("[check] S14 dual copies identical:", a == b)
    print("[check] mirror == source:", open(SRC, "rb").read() == open(MIRROR, "rb").read())
    print("[check] 310 lock string preserved:",
          "achieved power is computed post-hoc" in open(SRC, encoding="utf-8").read())
    print("DONE")


if __name__ == "__main__":
    main()
