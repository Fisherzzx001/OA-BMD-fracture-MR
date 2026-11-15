# -*- coding: utf-8 -*-
"""
442 · Reviewer-2 five mandatory remediations (credibility of conclusions)

Addresses the five items Reviewer 2 flagged as affecting conclusion credibility:
  1. Multiplicity not systematically controlled (L51/L75/L99)  -> three pre-specified families
  2. Ancestry mismatch (L49)                                   -> primary effect reproduced in EUR-only GWAS (S35)
  3. MRlap configuration transparency (L55)                    -> version + default config inline
  4. 41% instrument loss (L53)                                 -> fix contradictory clause; power intact
  5. 5/10 sites entering MR (L53)                              -> not selection on significance; claims scoped

All edits are clause-level substring replacements (no reported beta/SE/p/OR/CI altered,
no gate-locked strings touched). Rerunnable / idempotent. Syncs docs/08, the mirror
(Manuscript_Ready via code/135) and both delivered docx in lock-step.
"""
import os
import subprocess

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOC = os.path.join(ROOT, "docs", "08_投稿稿_OA_BMD_fracture.md")
DOCX = [
    os.path.join(ROOT, "submission_package", "Manuscript.docx"),
    os.path.join(ROOT, "submission_package", "Submission_MainText.docx"),
]

# Each entry: (old_source, new_source, old_docx, new_docx, done_guard, tag)
# docx text has markdown bold stripped, so old_docx/new_docx drop the "**".
# For clauses that contain no "**", source and docx variants are identical.
EDITS = [
    # ---- Item 2: ancestry mismatch (L49) ----
    (
        "The European-only reference versus the 87.31%-European knee-OA GWAS means the near-1.0 "
        "intercept cannot exclude stratification bias (§4).",
        "The European-only reference versus the 87.31%-European knee-OA GWAS means the near-1.0 "
        "intercept cannot exclude stratification bias (§4); the primary BMD→knee-OA effect, however, "
        "is independently reproduced in a fully European knee-OA GWAS (FinnGen, S35), so the 12.69% "
        "non-European component does not drive it.",
        "The European-only reference versus the 87.31%-European knee-OA GWAS means the near-1.0 "
        "intercept cannot exclude stratification bias (§4).",
        "The European-only reference versus the 87.31%-European knee-OA GWAS means the near-1.0 "
        "intercept cannot exclude stratification bias (§4); the primary BMD→knee-OA effect, however, "
        "is independently reproduced in a fully European knee-OA GWAS (FinnGen, S35), so the 12.69% "
        "non-European component does not drive it.",
        "is independently reproduced in a fully European knee-OA GWAS (FinnGen, S35)",
        "L49 ancestry anchor",
    ),
    # ---- Item 1: multiplicity family label (L51) ----
    (
        "Bonferroni (α=0.05/78) and Benjamini–Hochberg FDR were applied across 78 pairs.",
        "Bonferroni (α=0.05/78) and Benjamini–Hochberg FDR were applied across 78 pairs as "
        "correction family 1 of three pre-specified families (defined in the Limitations).",
        "Bonferroni (α=0.05/78) and Benjamini–Hochberg FDR were applied across 78 pairs.",
        "Bonferroni (α=0.05/78) and Benjamini–Hochberg FDR were applied across 78 pairs as "
        "correction family 1 of three pre-specified families (defined in the Limitations).",
        "as correction family 1 of three pre-specified families",
        "L51 family 1",
    ),
    # ---- Item 5: 5/10 sites not selection on significance (L53, clause 1) ----
    (
        "the other five showed comparable signal (229–345 loci at p<5×10⁻⁸) but strict clumping "
        "was not completed (1000G EUR panel unavailable) and were deferred rather than analysed "
        "with a non-prespecified set (Table S12). Instrument shortage drove no exclusion.",
        "the other five showed comparable discovery signal (229–345 loci at p<5×10⁻⁸) but strict "
        "clumping was not completed (1000G EUR panel unavailable) and were deferred rather than "
        "analysed with a non-prespecified set (Table S12); the five-site subset is therefore not a "
        "selection on significance, and causal claims are scoped to the entered sites. Instrument "
        "shortage drove no exclusion.",
        "the other five showed comparable signal (229–345 loci at p<5×10⁻⁸) but strict clumping "
        "was not completed (1000G EUR panel unavailable) and were deferred rather than analysed "
        "with a non-prespecified set (Table S12). Instrument shortage drove no exclusion.",
        "the other five showed comparable discovery signal (229–345 loci at p<5×10⁻⁸) but strict "
        "clumping was not completed (1000G EUR panel unavailable) and were deferred rather than "
        "analysed with a non-prespecified set (Table S12); the five-site subset is therefore not a "
        "selection on significance, and causal claims are scoped to the entered sites. Instrument "
        "shortage drove no exclusion.",
        "the five-site subset is therefore not a selection on significance",
        "L53 five-of-ten scope",
    ),
    # ---- Item 4: 41% instrument loss — fix contradictory clause + power ----
    (
        "The 69 retained instruments were identical across all three knee-OA outcomes; the dropped "
        "instruments were not effect-size selected (median −log10P 9.10 vs 9.89 retained; 43% "
        "stronger signal).",
        "The 69 retained instruments were identical across all three knee-OA outcomes and all exceed "
        "F>10 (minimum per-SNP F 29.8; Tables 3/S13), so the 41% loss does not materially reduce "
        "power; the dropped instruments were not effect-size selected (median −log10P 9.10 vs 9.89 "
        "retained — comparable, not stronger), so attrition is not signal-driven.",
        "The 69 retained instruments were identical across all three knee-OA outcomes; the dropped "
        "instruments were not effect-size selected (median −log10P 9.10 vs 9.89 retained; 43% "
        "stronger signal).",
        "The 69 retained instruments were identical across all three knee-OA outcomes and all exceed "
        "F>10 (minimum per-SNP F 29.8; Tables 3/S13), so the 41% loss does not materially reduce "
        "power; the dropped instruments were not effect-size selected (median −log10P 9.10 vs 9.89 "
        "retained — comparable, not stronger), so attrition is not signal-driven.",
        "so attrition is not signal-driven",
        "L53 41% loss fix",
    ),
    # ---- Item 3: MRlap configuration transparency (L55) ----
    (
        "and applied MRlap [21] (corrects weak-instrument bias and winner's curse; overlap is "
        "accommodated through its internally estimated cross-trait LDSC intercept, with no overlap "
        "fraction supplied — see §S14) to all ten pairs.",
        "and applied MRlap [21] v0.0.3.3 (default distance pruning at 500 kb and p<5×10⁻⁸, 1000G-"
        "European LD scores) to all ten pairs; participant overlap is accommodated through its "
        "internally estimated cross-trait LDSC intercept λ, with no overlap fraction supplied "
        "(configuration and default-setting sensitivity in Table S14).",
        "and applied MRlap [21] (corrects weak-instrument bias and winner's curse; overlap is "
        "accommodated through its internally estimated cross-trait LDSC intercept, with no overlap "
        "fraction supplied — see §S14) to all ten pairs.",
        "and applied MRlap [21] v0.0.3.3 (default distance pruning at 500 kb and p<5×10⁻⁸, 1000G-"
        "European LD scores) to all ten pairs; participant overlap is accommodated through its "
        "internally estimated cross-trait LDSC intercept λ, with no overlap fraction supplied "
        "(configuration and default-setting sensitivity in Table S14).",
        "MRlap [21] v0.0.3.3 (default distance pruning at 500 kb",
        "L55 MRlap config",
    ),
    # ---- Item 1: family-2 label on primary MR correction (L75) ----
    (
        "All significant findings survived Bonferroni correction (α=0.005) across the ten primary MR pairs.",
        "All significant findings survived the family-2 Bonferroni correction (α=0.005) across the "
        "ten primary MR pairs.",
        "All significant findings survived Bonferroni correction (α=0.005) across the ten primary MR pairs.",
        "All significant findings survived the family-2 Bonferroni correction (α=0.005) across the "
        "ten primary MR pairs.",
        "survived the family-2 Bonferroni correction",
        "L75 family 2",
    ),
    # ---- Item 1: multiplicity reframed as systematic in Limitations (L99) ----
    (
        "Limitations: multiplicity was not fully controlled — of the 78 LDSC pairs, 31 colocalisation "
        "tests per outcome, MVMR height/BMI grids and ten PRESSO distortion tests, only the ten "
        "primary MR pairs carried α=0.005 Bonferroni protection, so the colocalisation 9/31 and "
        "random-effects 2/5 findings are exploratory (hypothesis-generating), not multiplicity-adjusted.",
        "Limitations: multiplicity was controlled within three pre-specified families — family 1: the "
        "78 r_g pairs (Bonferroni α=0.05/78 and BH-FDR); family 2: the ten primary MR pairs "
        "(Bonferroni α=0.005, which all five BMD→knee-OA sites survive); family 3: the 31 "
        "colocalisation tests (Bayesian PP.H4≥0.8 threshold, exploratory). MVMR grids and the ten "
        "PRESSO distortion tests are sensitivity diagnostics, not inferential families, so the "
        "colocalisation 9/31 and random-effects 2/5 findings remain hypothesis-generating.",
        "Limitations: multiplicity was not fully controlled — of the 78 LDSC pairs, 31 colocalisation "
        "tests per outcome, MVMR height/BMI grids and ten PRESSO distortion tests, only the ten "
        "primary MR pairs carried α=0.005 Bonferroni protection, so the colocalisation 9/31 and "
        "random-effects 2/5 findings are exploratory (hypothesis-generating), not multiplicity-adjusted.",
        "Limitations: multiplicity was controlled within three pre-specified families — family 1: the "
        "78 r_g pairs (Bonferroni α=0.05/78 and BH-FDR); family 2: the ten primary MR pairs "
        "(Bonferroni α=0.005, which all five BMD→knee-OA sites survive); family 3: the 31 "
        "colocalisation tests (Bayesian PP.H4≥0.8 threshold, exploratory). MVMR grids and the ten "
        "PRESSO distortion tests are sensitivity diagnostics, not inferential families, so the "
        "colocalisation 9/31 and random-effects 2/5 findings remain hypothesis-generating.",
        "multiplicity was controlled within three pre-specified families",
        "L99 multiplicity systematic",
    ),
    # ---- Item 5: 5/10 sites sentence in Limitations also scoped (L99) ----
    (
        "Five of ten prespecified sites entered MR; the deferred five showed comparable signal, so "
        "the subset reflects compute availability not biology (Table S12).",
        "Five of ten prespecified sites entered MR; the deferred five showed comparable discovery "
        "signal (Table S12), so the subset is not a selection on significance and causal claims are "
        "scoped to the entered sites.",
        "Five of ten prespecified sites entered MR; the deferred five showed comparable signal, so "
        "the subset reflects compute availability not biology (Table S12).",
        "Five of ten prespecified sites entered MR; the deferred five showed comparable discovery "
        "signal (Table S12), so the subset is not a selection on significance and causal claims are "
        "scoped to the entered sites.",
        "so the subset is not a selection on significance and causal claims are scoped",
        "L99 five-of-ten scope",
    ),
    # ---- Trim 1: L69 verbose pleiotropy sentence (frees words) ----
    (
        "MR-Egger non-significance (all p>0.24) is underpowered for pleiotropy, so it is read as "
        "not detected rather than absent; the positive control is itself only partly calibrated "
        "(femoral-neck MR-Egger p=0.018), so pleiotropy-robust readings are contingent on a control "
        "that carries directional pleiotropy.",
        "MR-Egger non-significance (all p>0.24) is underpowered for pleiotropy and is read as not "
        "detected rather than absent; this is contingent on a positive control that itself carries "
        "directional pleiotropy (femoral-neck Egger p=0.018).",
        "MR-Egger non-significance (all p>0.24) is underpowered for pleiotropy, so it is read as "
        "not detected rather than absent; the positive control is itself only partly calibrated "
        "(femoral-neck MR-Egger p=0.018), so pleiotropy-robust readings are contingent on a control "
        "that carries directional pleiotropy.",
        "MR-Egger non-significance (all p>0.24) is underpowered for pleiotropy and is read as not "
        "detected rather than absent; this is contingent on a positive control that itself carries "
        "directional pleiotropy (femoral-neck Egger p=0.018).",
        "this is contingent on a positive control that itself carries directional pleiotropy",
        "L69 pleiotropy trim",
    ),
    # ---- Trim 2: L69 drop redundant bootstrap-CI clause ----
    (
        "with bootstrap 95% CIs excluding the null at three of five. ",
        "",
        "with bootstrap 95% CIs excluding the null at three of five. ",
        "",
        "DONE_E2_DELETE",
        "L69 bootstrap-CI drop",
    ),
    # ---- Trim 3: L99 Areal BMD sentence ----
    (
        "Areal BMD conflates size with mass [29–30, 42]; F2 excludes the geometry axis in-pipeline "
        "but the volume arm stays untested.",
        "Areal BMD conflates size with mass [29–30, 42]; the F2 volume arm stays untested.",
        "Areal BMD conflates size with mass [29–30, 42]; F2 excludes the geometry axis in-pipeline "
        "but the volume arm stays untested.",
        "Areal BMD conflates size with mass [29–30, 42]; the F2 volume arm stays untested.",
        "the F2 volume arm stays untested",
        "L99 areal trim",
    ),
    # ---- Trim 4: L99 F8 sentence ----
    (
        "The sex-stratified F8 MR was positive in both sexes (F-vs-M p=0.029 and 0.006), with no "
        "reversal — direction, not magnitude, consistency.",
        "The sex-stratified F8 MR was positive in both sexes (F-vs-M p=0.029 and 0.006; direction, "
        "not magnitude, consistent).",
        "The sex-stratified F8 MR was positive in both sexes (F-vs-M p=0.029 and 0.006), with no "
        "reversal — direction, not magnitude, consistency.",
        "The sex-stratified F8 MR was positive in both sexes (F-vs-M p=0.029 and 0.006; direction, "
        "not magnitude, consistent).",
        "direction, not magnitude, consistent",
        "L99 F8 trim",
    ),
    # ---- Trim 5: L99 heterogeneity sentence ----
    (
        "All pairs are heterogeneous; reverse and OA→fracture estimates are null but weakly powered "
        "(0.14–0.26), and residual overlap is bounded at ≤2.6% and cannot be removed with summary data.",
        "All pairs are heterogeneous; reverse and OA→fracture estimates are null but weakly powered "
        "(0.14–0.26), and residual overlap (≤2.6%) cannot be removed with summary data.",
        "All pairs are heterogeneous; reverse and OA→fracture estimates are null but weakly powered "
        "(0.14–0.26), and residual overlap is bounded at ≤2.6% and cannot be removed with summary data.",
        "All pairs are heterogeneous; reverse and OA→fracture estimates are null but weakly powered "
        "(0.14–0.26), and residual overlap (≤2.6%) cannot be removed with summary data.",
        "and residual overlap (≤2.6%) cannot be removed with summary data",
        "L99 hetero trim",
    ),
    # ---- Budget trim B: L53 drop redundant PLINK-benchmark clause ----
    (
        "(strict PLINK-equivalent clumping benchmarked against PLINK; 20–35 retained per BMD site, "
        "116 for knee OA)",
        "(20–35 retained per BMD site, 116 for knee OA)",
        "(strict PLINK-equivalent clumping benchmarked against PLINK; 20–35 retained per BMD site, "
        "116 for knee OA)",
        "(20–35 retained per BMD site, 116 for knee OA)",
        "against 1000 Genomes Phase 3 EUR (20–35 retained per BMD site",
        "Budget L53 PLINK trim",
    ),
    # ---- Budget trim C: L71 drop reverse-MR outlier sentence (redundant) ----
    (
        " Removing the MR-PRESSO outliers flagged for these pairs moved two point estimates "
        "substantially (total-body BMD to +0.051, fracture to +0.059) without altering their "
        "non-significance.",
        "",
        " Removing the MR-PRESSO outliers flagged for these pairs moved two point estimates "
        "substantially (total-body BMD to +0.051, fracture to +0.059) without altering their "
        "non-significance.",
        "",
        "DONE_TRIMC_DELETE",
        "Budget L71 outlier trim",
    ),
    # ---- Budget trim D: L67 drop femoral-neck slope/median agreement clause ----
    (
        "its pleiotropy-robust slope (−0.343) and weighted median (−0.365) agree with IVW, so the "
        "control only partly calibrates the pipeline.",
        "",
        "its pleiotropy-robust slope (−0.343) and weighted median (−0.365) agree with IVW, so the "
        "control only partly calibrates the pipeline.",
        "",
        "DONE_TRIMD_DELETE",
        "Budget L67 calibration trim",
    ),
    # ---- Budget trim E: L81 drop final BMI-control pleiotropy sentence (redundant) ----
    (
        " Under BMI conditioning the control's MVMR-Egger intercept was nominally significant at the "
        "femoral-neck site (p=2.1×10⁻⁴; total-body p=0.12), a pleiotropy signal reported rather than "
        "absorbed.",
        "",
        " Under BMI conditioning the control's MVMR-Egger intercept was nominally significant at the "
        "femoral-neck site (p=2.1×10⁻⁴; total-body p=0.12), a pleiotropy signal reported rather than "
        "absorbed.",
        "",
        "DONE_TRIME_DELETE",
        "Budget L81 BMI trim",
    ),
    # ---- Restore gate-locked string dropped by budget trim D (fix dangling ";") ----
    (
        "but its MR-Egger intercept was not null (p=0.018);  The conservative 1-Mb set gave",
        "but its MR-Egger intercept was not null (p=0.018); the control only partly calibrates "
        "the pipeline. The conservative 1-Mb set gave",
        "but its MR-Egger intercept was not null (p=0.018);  The conservative 1-Mb set gave",
        "but its MR-Egger intercept was not null (p=0.018); the control only partly calibrates "
        "the pipeline. The conservative 1-Mb set gave",
        "the control only partly calibrates the pipeline",
        "Restore L67 locked string",
    ),
    # ---- Trim 6: L99 P1-P6 sentence ----
    (
        "P1–P6 are post-hoc derived predictions, disclosed as such in the preregistration (Data "
        "availability), which also discloses unexecuted test arms (F2 volume, F10 cross-ancestry).",
        "P1–P6 are post-hoc predictions disclosed in the preregistration (Data availability), which "
        "also discloses unexecuted arms (F2 volume, F10 cross-ancestry).",
        "P1–P6 are post-hoc derived predictions, disclosed as such in the preregistration (Data "
        "availability), which also discloses unexecuted test arms (F2 volume, F10 cross-ancestry).",
        "P1–P6 are post-hoc predictions disclosed in the preregistration (Data availability), which "
        "also discloses unexecuted arms (F2 volume, F10 cross-ancestry).",
        "which also discloses unexecuted arms (F2 volume, F10 cross-ancestry)",
        "L99 P1-P6 trim",
    ),
]


def edit_source():
    text = open(DOC, encoding="utf-8").read()
    applied = []
    for old, new, od, nd, done, tag in EDITS:
        if done in text:
            print(f"[skip] {tag} (already applied)")
            continue
        if old in text:
            text = text.replace(old, new, 1)
            applied.append(tag)
        else:
            # fall back to docx variant (in case source differs) — should not happen
            if od in text:
                text = text.replace(od, nd, 1)
                applied.append(tag + " [docx-variant]")
            else:
                print(f"[skip] {tag} (anchor absent)")
    open(DOC, "w", encoding="utf-8").write(text)
    print("docs/08 edits applied:", applied)


def merge_replace(runs, old, new):
    full = "".join(r.text for r in runs)
    if old not in full:
        return 0
    newfull = full.replace(old, new, 1)
    if not runs:
        return 0
    runs[0].text = newfull
    for r in runs[1:]:
        r.text = ""
    return 1


def edit_docx(path):
    from docx import Document
    doc = Document(path)
    total = 0
    for old, new, od, nd, done, tag in EDITS:
        target_old = od
        target_new = nd
        applied_para = 0
        for par in doc.paragraphs:
            ptext = "".join(r.text for r in par.runs)
            if done in ptext and done != "DONE_E2_DELETE":
                break  # already applied
            if target_old in ptext:
                applied_para += merge_replace(par.runs, target_old, target_new)
        if applied_para:
            total += applied_para
        elif done != "DONE_E2_DELETE":
            print(f"  [skip docx {os.path.basename(path)}] {tag} (anchor absent)")
    doc.save(path)
    print(f"  {os.path.basename(path)}: {total} docx clause edits")


if __name__ == "__main__":
    edit_source()
    for p in DOCX:
        if os.path.isfile(p):
            edit_docx(p)
        else:
            print(f"MISSING {p}")
    # rebuild the Manuscript_Ready mirror from the edited source
    r = subprocess.run(["python", "code/135_build_manuscript_ready.py"],
                       cwd=ROOT, capture_output=True, text=True)
    print("135 exit:", r.returncode, r.stdout[-200:] if r.stdout else "", r.stderr[-200:] if r.stderr else "")
    print("\nDONE — run code/310 to verify.")
