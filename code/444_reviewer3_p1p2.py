# -*- coding: utf-8 -*-
"""Reviewer-3 adversarial review (docs/151) — P1 + P2 remediations.

P1:
  N3  family 3 (Bayesian PP.H4) is NOT a multiplicity-controlled family: recast
      the "three pre-specified families" claim as two analysis-defined frequentist
      families + exploratory Bayesian thresholding; drop the misleading
      "pre-specified" time-ordering (HARKing risk) -> "analysis-defined".
  N4  weighted medians are attenuated (about one-third at arms/legs) with
      right-skewed bootstrap CIs: state this in the body and give the full
      attenuation table + CI-asymmetry explanation in Table S17e.
  N5  MR-PRESSO: stop using "outlier removal left direction unchanged" as a
      robustness claim while also saying PRESSO misses non-outlier pleiotropy.
  N8  FinnGen: S35 reframed as within-ancestry direction consistency, not an
      independent replication; shared instruments / residual overlap persist.
  N9  F2 conditional MVMR (conditional F<4): "direction-indicative" has no
      statistical basis -> "non-informative under weak conditional instruments".

P2:
  N6  reverse-MR power paradox: explained in Table S17f (OA h2 0.035 -> reverse
      instrument R2_exposure 0.0003 vs forward 0.034-0.063).
  N7  negative cross-trait LDSC lambda: explained in Table S14.
  N10 +/-20% absolute-scale uncertainty attached to headline OR in Conclusion.
  N11 1,980 knee-OA "loci" were NOT clumped: honestly restated as LD-correlated
      SNPs, so the 94.4% sign-concordance is descriptive, not independent loci.

Idempotent. Body edits go to docs/08 + both docx (run-level merge); supplement
edits go to the markdown supplement files (docx strips the supplement block;
mirror rebuilt separately via code/135). No reported beta/SE/p/OR/CI is changed.
Gate-locked strings (FinnGen "replicated", "only partly calibrates", Highlights
"2/5 sites", etc.) are deliberately untouched.
"""
from pathlib import Path
from docx import Document

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "docs/08_投稿稿_OA_BMD_fracture.md"
S14 = ROOT / "docs/Supplementary/Table_S14_ancestry.md"
S17 = ROOT / "docs/Supplementary/Table_S17_sensitivity.md"
S35 = ROOT / "docs/Supplementary/Table_S35_eur_only_mr.md"
DOCX = [ROOT / "submission_package/Manuscript.docx",
        ROOT / "submission_package/Submission_MainText.docx"]

# ---------------- BODY (md + docx) ----------------
BODY = [
    # N11 — honest restatement of 1,980 loci as LD-correlated SNPs
    ("N11-1980",
     "across 1,980 genome-wide-significant knee-OA loci, 94.4% had concordant signed Z",
     "across 1,980 genome-wide-significant knee-OA SNPs (LD-correlated, not independent loci), 94.4% had concordant signed Z"),

    # N4 — weighted-median attenuation stated in body; keep word cost low
    ("N4-WM",
     "MR-Egger intercepts [19] were all p>0.24 and weighted medians matched the IVW sign at all five sites, MR-Egger non-significance",
     "MR-Egger intercepts [19] were null and weighted medians matched the IVW sign but were attenuated (~one-third at arms/legs; S17e); MR-Egger non-significance"),

    # N4/N5 — trim redundant hedge to reclaim words
    ("N4-trim",
     "is underpowered for pleiotropy and is read as not detected rather than absent",
     "is underpowered for pleiotropy and means undetected, not absent"),

    # N5 — PRESSO outlier removal no longer sold as robustness evidence
    ("N5-PRESSO",
     "All five effects remained positively causal after MR-PRESSO outlier removal (Table S17).",
     "Outlier removal left signs unchanged (Table S17), not robustness evidence given the distortion below."),

    # N3 — L51 two analysis-defined families
    ("N3-L51",
     "as correction family 1 of three pre-specified families (defined in the Limitations)",
     "as correction family 1 of two analysis-defined multiplicity families (defined in the Limitations)"),

    # N3 — L99 family 3 recast as exploratory Bayesian thresholding
    ("N3-L99",
     "multiplicity was controlled within three pre-specified families — family 1: the 78 r_g pairs (Bonferroni α=0.05/78 and BH-FDR); family 2: the ten primary MR pairs (Bonferroni α=0.005; all five sites survive under fixed effects, none under random effects); family 3: the 31 colocalisation tests (Bayesian PP.H4≥0.8 threshold, exploratory). MVMR grids and the ten PRESSO distortion tests are sensitivity diagnostics, not inferential families, so the colocalisation 9/31 and random-effects 2/5 findings remain hypothesis-generating.",
     "multiplicity was controlled within two analysis-defined families — family 1: the 78 r_g pairs (Bonferroni α=0.05/78 and BH-FDR); family 2: the ten primary MR pairs (Bonferroni α=0.005; all five sites survive under fixed effects, none under random effects). Colocalisation (31 tests) applied exploratory Bayesian posterior thresholding (PP.H4≥0.8), not frequentist error control, and MVMR grids and the ten PRESSO distortion tests are sensitivity diagnostics; the colocalisation 9/31 and random-effects 2/5 findings therefore remain hypothesis-generating."),

    # N9 — F2 conditional MVMR: direction-indicative -> non-informative
    ("N9-F2",
     "with conditional F2 analysis direction-indicative not conclusive (conditional F below 4 in all but one mode; Table S28, Part B)",
     "with conditional F2 analysis non-informative under weak conditional instruments (conditional F below 4 in all but one mode; Table S28, Part B)"),

    # N10 — attach scale uncertainty to the headline OR in the Conclusion
    ("N10-scale",
     "95% CI 1.025–1.135 across sites) but attaches to a ubiquitous, modifiable exposure",
     "95% CI 1.025–1.135 across sites, ±~20% absolute-scale) but attaches to a ubiquitous, modifiable exposure"),
]

# ---------------- SUPPLEMENT (md only) ----------------
S14_OLD = "which ranged from −0.0361 to +0.0153 (|λ| ≤ 0.037) across the ten pairs."
S14_NEW = ("which ranged from −0.0361 to +0.0153 (|λ| ≤ 0.037) across the ten pairs. "
           "A negative cross-trait intercept (BMD_total→fracture, λ=−0.0361) runs opposite to the "
           "positive bias expected from sample overlap; in cross-trait LDSC such negative values arise "
           "from sampling variation around a near-zero intercept and from genetic nurture/assortative "
           "mating, not from negative genetic sharing. Because the sign is mixed and |λ| small, no "
           "positive overlap-driven inflation needs to be subtracted here, and the negative value does "
           "not contradict the conclusion that residual overlap is minor.")

S17_ANCHOR = "The five highest PP.H4 values against knee OA are rs11228269 (0.7994), rs9594759 (0.6271), rs12600549 (0.5530), rs9293511 (0.4518), rs6040061 (0.2941). Against the fracture endpoint 9/31 loci reach PP.H4 ≥0.8, so the endpoint-split contrast reported in the manuscript is preserved at every threshold from 0.5 to 0.9: relaxing the cut-off to 0.7 adds exactly one locus and leaves the interpretation unchanged."
S17_ADD = """

### S17-e. Pleiotropy-robust estimators: magnitude attenuation and bootstrap-CI asymmetry

The body reports that weighted medians share the IVW sign at all five sites; the full pattern is a material attenuation at some sites. From Table 2:

| Pair | IVW β | Weighted median β | WM/IVW | Weighted-median 95% CI |
| --- | ---: | ---: | ---: | --- |
| total → knee | +0.078 | +0.065 | 0.83 | +0.027 to +0.177 |
| lumbar → knee | +0.095 | +0.115 | 1.21 | +0.014 to +0.189 |
| arms → knee | +0.075 | +0.026 | 0.35 | +0.019 to +0.217 |
| legs → knee | +0.060 | +0.026 | 0.43 | −0.003 to +0.132 |
| femoral neck → knee | +0.058 | +0.052 | 0.90 | −0.000 to +0.109 |

At arms and legs the most pleiotropy-robust single estimator is roughly one-third of the IVW point estimate, and the legs interval crosses zero. This is disclosed rather than hidden: the headline OR 1.060–1.099 is an IVW quantity, and a reader who prefers the more conservative robust estimator should read the effect as smaller and less certain at those two sites. The weighted-median bootstrap intervals are visibly right-skewed (point estimates near the lower bound, e.g. arms +0.026 against 0.019–0.217): Bowden weights (β_exposure/SE_outcome)² are highly uneven with 16–31 instruments, so a few strong instruments dominate and the percentile bootstrap distribution has a long right tail. The asymmetry is a small-instrument-number property of the estimator, not evidence of a second effect direction; we therefore report the intervals as computed and do not force symmetry.

### S17-f. Why reverse MR has more instruments yet lower power

The reverse analysis uses 69 OA instruments versus 16–31 per forward BMD site, yet its achieved power is 0.14–0.26 against 0.93–1.00 forward. Instrument count is not the driver of power: the exposure variance explained is. From `instrument_strength_steiger.csv`, forward BMD→OA pairs have R²_exposure 0.034–0.063, whereas the OA→BMD pairs have R²_exposure≈0.0003, because OA liability has low SNP heritability (h²=0.035) versus 0.26–0.39 for BMD. With near-zero exposure R² the non-centrality parameter stays small regardless of instrument number. The reverse MDE (minimum detectable |β|=0.089–0.090 in Table 3) is computed from these instruments' actual SE (≈0.032); the quoted exclusion bound |β|≥0.09 SD therefore corresponds to this MDE and is not an arbitrary cut-off. Forward and reverse estimates are on the same as-deposited scale, so the power asymmetry is biological (h²), not a scale artefact.
"""

S35_OLD = "Because FinnGen knee OA and the GO knee-OA GWAS are independent European samples, this table is a within-ancestry consistency check; cross-ancestry replication remains unexecuted (§S32, F10)."
S35_NEW = ("Because FinnGen knee OA and the GO knee-OA GWAS are independent European samples, this table "
           "is a within-ancestry direction-consistency check, not an independent replication of the causal "
           "claim: the BMD exposure instruments are the same Qian 2025 set used throughout, so the "
           "shared-instrument and residual UK Biobank overlap limitations (Table S14) persist; "
           "cross-ancestry replication remains unexecuted (§S32, F10).")

SUPP = [
    (S14, "N7-lambda", S14_OLD, S14_NEW),
    (S17, "S17e/f", S17_ANCHOR, S17_ANCHOR + S17_ADD),
    (S35, "N8-FinnGen", S35_OLD, S35_NEW),
]


def edit_text(path, triples):
    t = path.read_text(encoding="utf-8")
    out = []
    for name, old, new in triples:
        if new in t:
            out.append(f"  {name}: SKIP (already applied)")
        elif old in t:
            t = t.replace(old, new, 1)
            out.append(f"  {name}: APPLIED")
        else:
            out.append(f"  {name}: MISSING anchor")
    path.write_text(t, encoding="utf-8")
    return out


def edit_docx(path, triples):
    d = Document(path)
    out = []
    changed = False
    for name, old, new in triples:
        hits = 0
        for par in d.paragraphs:
            runs = par.runs
            full = "".join(r.text for r in runs)
            if old in full:
                merged = full.replace(old, new, 1)
                if runs:
                    runs[0].text = merged
                    for r in runs[1:]:
                        r.text = ""
                hits += 1
                changed = True
        if hits:
            out.append(f"  {name}: APPLIED ({hits})")
        else:
            out.append(f"  {name}: MISSING/elsewhere")
    if changed:
        d.save(path)
    return out


def main():
    print("== docs/08 body ==")
    for x in edit_text(SRC, BODY):
        print(x)
    print("== supplements ==")
    for path, name, old, new in SUPP:
        for x in edit_text(path, [(name, old, new)]):
            print(f"[{path.name}] {x}")
    for p in DOCX:
        print(f"== {p.name} ==")
        for x in edit_docx(p, BODY):
            print(x)


if __name__ == "__main__":
    main()
