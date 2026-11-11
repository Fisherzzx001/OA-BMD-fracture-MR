# -*- coding: utf-8 -*-
"""
code/337 - Word-neutral compression of docs/08 IMRaD prose to OAC <=4000.

Hard guardrail: every accepted replacement must leave the multiset of
numeric-bearing tokens in the IMRaD body EXACTLY unchanged (same as the
code/336 snapshot). This guarantees no reported beta/SE/p/OR/CI/r_g/N,
reference number, table/figure number or p-value is touched. Replacements that
would alter a numeric token are skipped and reported.

Run:  python3 code/337_compress_wordneutral.py            # dry run (report)
      python3 code/337_compress_wordneutral.py --apply     # write docs/08
"""
import re, sys
from collections import Counter

F = "docs/08_投稿稿_OA_BMD_fracture.md"

REPL = [
    # ---- Introduction ----
    ("Osteoarthritis (OA) is a prevalent disorder of ageing [1]; areal bone mineral density (BMD), the clinical measure of skeletal strength, is measured in the same populations.",
     "Osteoarthritis (OA) and areal bone mineral density (BMD, the clinical measure of skeletal strength) co-occur in ageing populations [1]."),
    ("Observational studies, reviewed in [3], link higher areal BMD to increased OA prevalence and incidence, although osteoporosis and OA were long regarded as opposing conditions [4].",
     "Observational studies [3] link higher areal BMD to increased OA prevalence, though the two were long viewed as opposing [4]."),
    ("The direction matters clinically: if high BMD predisposes to OA, joint endpoints deserve monitoring in bone-anabolic trials, though whether such a cost is mediated by bone mass or skeletal size is itself unresolved.",
     "Clinically, if high BMD predisposes to OA, joint endpoints warrant monitoring in bone-anabolic trials; whether the cost reflects bone mass or skeletal size is unresolved."),
    ("Applied to OA and BMD, scoping and methodological reviews have given partial answers [7, 8]:",
     "Reviews give partial answers [7, 8]:"),
    ("and a 2025 two-sample MR reported total BMD as a significant positive mediator of knee-OA risk (mediation proportion 3.59%) [11].",
     "and a 2025 MR reported total BMD as a positive knee-OA mediator (3.59% mediation) [11]."),
    ("Our BMD→knee-OA analysis, by contrast, finds a positive effect contrasting with Zhou's negative femoral-neck estimate (discussed below).",
     "Our analysis instead finds a positive BMD→knee-OA estimate (discussed below)."),
    ("Our study extends that work in four ways: (i) ten regional BMD sites rather than one composite, so site specificity is tested directly (five of ten entered MR; §2.3); (ii) conditioning on adult height and separately on BMI, separating size from adiposity; (iii) an internal BMD→fracture positive control; and (iv) a full 13×13 LDSC matrix rather than pairwise estimates.",
     "We extend this in four ways: (i) ten regional sites, not one composite, testing site specificity directly (five of ten entered MR; §2.3); (ii) conditioning separately on height and BMI to separate size from adiposity; (iii) an internal BMD→fracture positive control; (iv) a full 13×13 LDSC matrix."),
    ("Three questions remain open. First, most analyses included at most two BMD sites, so site specificity is unknown despite OA being joint-localised and BMD site-regulated. Second, OA has usually been represented by a single definition, leaving their bone relationship uncertain. Third, the OA–BMD–fracture triangle has not been analysed jointly; modelling all three tests the pipeline and whether the OA–BMD association carries fracture relevance.",
     "Three questions remain. First, most analyses included at most two BMD sites, so site specificity is unknown despite OA being joint-localised and BMD site-regulated. Second, OA was usually represented by a single definition. Third, the OA–BMD–fracture triangle was not analysed jointly; modelling all three tests the pipeline and whether the OA–BMD association carries fracture relevance."),
    # ---- Methods ----
    ("All four sources contribute sex-combined primary estimates, so every estimate is sex-averaged; a sex-stratified sensitivity analysis (F8) confirmed BMD→OA positivity in both sexes (§4).",
     "All four contribute sex-combined estimates, so all are sex-averaged; a sex-stratified sensitivity (F8) confirmed BMD→OA positivity in both sexes (§4)."),
    ("All LDSC used 1000 Genomes Phase 3 European LD scores; ancestry composition is in Table S14. The European-only reference against the 87.31%-European knee-OA GWAS means the near-1.0 intercept does not exclude population-stratification bias (primary limitation, §4).",
     "All LDSC used 1000 Genomes Phase 3 European LD scores; ancestry is in Table S14. The European-only reference versus the 87.31%-European knee-OA GWAS means the near-1.0 intercept cannot exclude stratification bias (§4)."),
    ("Genetic correlation (r_g) was estimated by LDSC [5] in Python (Bulik-Sullivan) on 690,150 HapMap3 SNPs with 1000 Genomes Phase 3 European LD scores, using delete-one-block jackknife (100 blocks) for SEs.",
     "Genetic correlation (r_g) was estimated by LDSC [5] on 690,150 HapMap3 SNPs with 1000 Genomes Phase 3 European LD scores (delete-one-block jackknife, 100 blocks)."),
    ("Both Bonferroni (α=0.05/78) and Benjamini–Hochberg FDR were applied across 78 phenotype pairs.",
     "Bonferroni (α=0.05/78) and Benjamini–Hochberg FDR were applied across 78 pairs."),
    ("A spectral decomposition of the resulting 13×13 r_g matrix was performed as an exploratory probe of latent structure (§3, Fig. S22).",
     "Spectral decomposition of the 13×13 r_g matrix probed latent structure (§3, Fig. S22)."),
    # NB: STROBE-MR item 5 gate (code/182) requires the literals
    # "Valid instruments must satisfy three assumptions" / "relevance" /
    # "independence from exposure" / "exclusion restriction". Keep them all;
    # compress only the trailing "We address each below" and the dash clause.
    ("**MR.** Valid instruments must satisfy three assumptions: relevance; independence from exposure–outcome confounders; and exclusion restriction — no effect on the outcome except through the exposure. We address each below.",
     "**MR.** Valid instruments must satisfy three assumptions — relevance, independence from exposure–outcome confounders, and exclusion restriction (no effect on the outcome except via the exposure)."),
    ("Clumping completed for five sites (arms, femoral neck, lumbar spine, legs, total body), each ≥15 independent instruments; the remaining five showed comparable genome-wide signal (229–345 loci at p<5×10⁻⁸) but strict clumping was not completed (1000G EUR panel unavailable) and they were deferred rather than analysed with a non-prespecified set (Table S12).",
     "Clumping completed for five sites (arms, femoral neck, lumbar spine, legs, total body; each ≥15 instruments); the other five showed comparable signal (229–345 loci at p<5×10⁻⁸) but strict clumping was not completed (1000G EUR panel unavailable) and were deferred rather than analysed with a non-prespecified set (Table S12)."),
    ("Harmonisation against the outcome GWAS reduced the analysed sets to 16–31 instruments per forward pair and 69 per knee-OA-exposure pair; 47 of 116 knee-OA instruments (41%) were lost (Fig. S9), an attrition we report.",
     "Harmonisation reduced the analysed sets to 16–31 instruments per forward pair and 69 for knee-OA-as-exposure; 47 of 116 knee-OA instruments (41%) were lost (Fig. S9)."),
    ("Under the Burgess approximation [18], the IVW bias equals the overlap fraction × r_g and acts in the direction of r_g (pairwise UKB overlaps and bias bounds in Table S14); because r_g shares the sign of the MR effect in all ten pairs (as does the SNP-effect correlation ρ_GE, Fig. S17), it acts away from the null.",
     "Under the Burgess approximation [18], IVW bias equals the overlap fraction × r_g and, because r_g shares the MR effect's sign in all ten pairs (as does ρ_GE, Fig. S17), acts away from the null (pairwise overlaps and bounds in Table S14)."),
    ("The joint correction moved every estimate further from the null (median +20.5%; observed-versus-corrected difference p<0.05 for 9/10 pairs), showing towards-null components dominate the away-from-null overlap component, reversing no direction or significance.",
     "The joint correction moved every estimate from the null (median +20.5%; observed-versus-corrected difference p<0.05 for 9/10), showing towards-null components dominate; no direction or significance reversed."),
    ("Primary IVW estimates are conservative for the combined bias, not overlap alone; the internally estimated per-pair cross-trait intercepts and the default-setting sensitivity are in Table S14.",
     "Primary IVW estimates are conservative for the combined bias, not overlap alone (per-pair intercepts and default-setting sensitivity in Table S14)."),
    ("Worst-case Burgess biases from the tabulated overlaps are ≤0.7% of the BMD→knee-OA estimate and ≤4.5% of the positive-control estimate (Table S14); MRlap retained all ten directions (median +20.5% away from the null), so the positive joint-endpoint effect survives overlap correction.",
     "Worst-case Burgess bias is ≤0.7% of the BMD→knee-OA estimate and ≤4.5% of the positive control (Table S14); MRlap retained all ten directions (median +20.5%), so the joint-endpoint effect survives overlap correction."),
    ("The Z ratio is unity by construction and cannot identify K; it is a same-source check only.",
     "The Z ratio is unity by construction, a same-source check that cannot identify K."),
    ("Each correction multiplies β and SE by one constant, preserving the per-SNP z; for the femoral-neck, lumbar and arms BMD→OA pairs it also removed one instrument each (26→25, 22→21, 17→16), so their pooled z and p use the final set.",
     "Each correction multiplies β and SE by one constant, preserving per-SNP z; femoral-neck, lumbar and arms BMD→OA pairs each lost one instrument (26→25, 22→21, 17→16), so pooled z and p use the final set."),
    ("**MVMR.** Multivariable MR conditioned the BMD→OA effect on adult height (Yengo et al. 2022, N=5.0M) [23] as a skeletal-size proxy, using the MVMR-IVW estimator [24] and MVMR-Egger [25], judged by the Sanderson–Smith conditional F (Fc≥10).",
     "**MVMR.** Multivariable MR conditioned BMD→OA on adult height (Yengo 2022, N=5.0M) [23] as a size proxy, using MVMR-IVW [24] and MVMR-Egger [25], judged by Sanderson–Smith conditional F (Fc≥10)."),
    ("The primary configuration requested up to n_h=100 height instruments (49–50 survived harmonisation); a grid (n_h ∈ {0,100,300,1,000}) and both EUR and EUR_noUKB sources are in the Supplementary Material.",
     "The primary configuration used n_h=100 height instruments (49–50 survived harmonisation); a grid (n_h ∈ {0,100,300,1,000}) and EUR/EUR_noUKB sources are in Supplementary."),
    ("Height was the primary size proxy: BMI is a major knee-OA risk factor and part of the adiposity pathway through which a BMD-raising process could raise knee-OA risk, whereas height indexes skeletal size alone.",
     "Height was the primary size proxy: BMI is itself a knee-OA risk factor in the adiposity pathway through which BMD could raise OA risk, whereas height indexes skeletal size alone."),
    ("A second MVMR conditioned the same effect on BMI (GIANT–UK Biobank, Pulit et al. 2019, N≈806,800) [26] with the same estimators, n_bmi=100 primary (61–63 instruments survived) and a grid (n_bmi ∈ {100,300,800}); conditional F statistics and estimates are in Supplementary Table S16.",
     "A second MVMR conditioned the effect on BMI (GIANT–UK Biobank, Pulit 2019, N≈806,800) [26], n_bmi=100 primary (61–63 survived), grid (n_bmi ∈ {100,300,800}); conditional F and estimates are in Table S16."),
    ("Height and BMI were conditioned in separate analyses, so each proxy is adjusted one at a time rather than jointly.",
     "Height and BMI were conditioned separately, each adjusted one at a time."),
    # ---- Results ----
    ("Because LDSC uses a European-only reference while the knee-OA GWAS is multi-ancestry (Table S14), the near-1.0 intercept does not exclude population-stratification or ancestry-mismatch bias (§4).",
     "Because LDSC uses a European-only reference while knee-OA GWAS is multi-ancestry (Table S14), the near-1.0 intercept cannot exclude stratification/ancestry-mismatch bias (§4)."),
    ("The two definitions also agreed at the variant level: across 1,980 genome-wide-significant knee-OA loci, 94.4% had concordant signed Z with self-reported OA.",
     "The definitions agreed at variant level: across 1,980 genome-wide-significant knee-OA loci, 94.4% had concordant signed Z with self-reported OA."),
    ("A spectral decomposition resolved two orthogonal axes — a bone-mass/fracture axis (PC1 66.7% of variance; all ten BMD sites loading +0.30, fracture −0.18, both OA definitions ≈0) and an OA axis (PC2 14.1%; OA_knee +0.69, OA_selfrep +0.68, BMD ≈0) — so BMD and OA occupy separate genetic dimensions (Fig. S22; exploratory, OA h² 0.035).",
     "Spectral decomposition resolved two orthogonal axes — a bone-mass/fracture axis (PC1 66.7%; ten BMD sites loading +0.30, fracture −0.18, OA ≈0) and an OA axis (PC2 14.1%; OA_knee +0.69, OA_selfrep +0.68, BMD ≈0) — so BMD and OA occupy separate genetic dimensions (Fig. S22; exploratory, OA h² 0.035)."),
    ("but its MR-Egger intercept was not null (p=0.018), a directional-pleiotropy signal we report rather than absorb; its pleiotropy-robust slope (−0.343) and weighted median (−0.365) agree with IVW, so the control calibrates the pipeline despite its own non-clean diagnostic.",
     "but its MR-Egger intercept was not null (p=0.018), a directional-pleiotropy signal reported rather than absorbed; its pleiotropy-robust slope (−0.343) and weighted median (−0.365) agree with IVW, so the control calibrates the pipeline despite its non-clean diagnostic."),
    ("with 16–31 instruments the intercept test has limited power, so a non-significant intercept excludes detectable unbalanced pleiotropy rather than pleiotropy itself.",
     "with 16–31 instruments the intercept test has limited power, so a non-significant intercept excludes detectable — not all — unbalanced pleiotropy."),
    ("indicating high I² is a property of the BMD instrument sets rather than an OA-specific IVW failure.",
     "indicating high I² reflects the BMD instrument sets, not an OA-specific IVW failure."),
    ("That non-result does not discriminate direction: with 69 instruments and a null reverse estimate the test lacks power to separate the models and is not read as evidence BMD rather than OA drives the association.",
     "That non-result does not discriminate direction: with 69 instruments and a null reverse estimate, the test lacks power and is not read as evidence for BMD over OA as the driver."),
    ("How much of the BMD→OA effect is skeletal size? Conditioning the BMD→OA effect on height attenuated it by 9.5–33.2% across the two tabulated univariable baselines (13.5–33.2% against the primary Table 2 baselines and 9.5–21.2% against the MVMR-set baselines; primary n_h=100/EUR; Supplementary Tables S9 and S10 and Supplementary Fig. S10), leaving 67–91% independent of height.",
     "How much of the BMD→OA effect is skeletal size? Height conditioning attenuated it 9.5–33.2% across two univariable baselines (13.5–33.2% vs primary Table 2; 9.5–21.2% vs MVMR-set; primary n_h=100/EUR; Tables S9/S10, Fig. S10), leaving 67–91% height-independent."),
    ("Because height indexes stature but not bone width or geometry [28], this residual is not identified as bone mass; test F2 addresses the geometry component (§4, Table S28).",
     "Because height indexes stature not bone width/geometry [28], this residual is not identified as bone mass; F2 addresses geometry (§4, Table S28)."),
    ("The positive control confirmed that the established BMD→fracture effect survives the same conditioning: BMD→fracture moved only −1.7% to −3.0% across the same two baselines (−2.8% against the MVMR-set baseline), the control being insensitive to baseline choice.",
     "The positive control survived the same conditioning: BMD→fracture moved only −1.7% to −3.0% across the same baselines (−2.8% vs MVMR-set), insensitive to baseline choice."),
    ("Because conditioning on height captured a genuine, direct OA-raising height effect (θ_H=+0.145 to +0.165, p≈10⁻¹⁶ to 10⁻¹⁹), the attenuation reflects size confounding, not dilution by a null covariate.",
     "Height conditioning captured a genuine direct OA-raising effect (θ_H=+0.145 to +0.165, p≈10⁻¹⁶ to 10⁻¹⁹), so attenuation reflects size confounding, not dilution by a null covariate."),
    # ---- Discussion ----
    ("The opposing-endpoint pattern is old [2, 4]; what is new is its demonstration for one instrument set in one calibrated pipeline.",
     "The opposing-endpoint pattern is old [2, 4]; what is new is its demonstration in one calibrated pipeline."),
    ("so the hip evidence is not on one measurement scale; we adjudicate rather than cite selectively, and F1 demands a positive, knee-like hip effect.",
     "so the hip evidence is not on one scale; we adjudicate rather than cite selectively, and F1 demands a positive, knee-like hip effect."),
    ("This is a difference of scale and design, not biology: Hartley instrumented heel-ultrasound eBMD (bone size), we DXA areal BMD (area-normalised); Hartley's one-sample 2SLS (overlap 100%) contrasts with ours bounded at ≤2.6%, MRlap moving every estimate a median 20.5% further from the null.",
     "This is a difference of scale and design, not biology: Hartley instrumented heel-ultrasound eBMD (bone size), we DXA areal BMD; Hartley's one-sample 2SLS (overlap 100%) contrasts with ours bounded at ≤2.6%, MRlap moving every estimate a median 20.5% from the null."),
    ("Conditioning the OA→fracture effect on BMD_total in a multivariable MR (69 OA instruments; Table S31) left it unchanged (β=+0.034, 95% CI −0.029 to 0.097), confirming BMD is not a confounder of the null and that the positive r_g is non-causal; Egger-MVMR showed no horizontal pleiotropy (intercept p=0.20).",
     "Conditioning OA→fracture on BMD_total in multivariable MR (69 OA instruments; Table S31) left it unchanged (β=+0.034, 95% CI −0.029 to 0.097), confirming BMD is not a confounder and the positive r_g non-causal; Egger-MVMR showed no horizontal pleiotropy (p=0.20)."),
    ("Higher BMD in OA has not meant fewer fractures, the residual association attributed to falls [41].",
     "Higher BMD in OA has not meant fewer fractures; the residual association is attributed to falls [41]."),
]


def body_span(text):
    ai = text.find("**Objective.**")
    ci = text.find("**Conclusions.**", ai)
    tail = text[ci:]
    em = re.search(r"\n---", tail)
    ae = ci + (em.start() if em else len(tail))
    leg = text.find("## Figure legends")
    return ae, leg


def prose_lines(body):
    return "\n".join(l for l in body.splitlines() if not l.strip().startswith("|"))


def bwords(body):
    return len([t for t in prose_lines(body).split()
                if re.search(r"[^\W_]", t, re.UNICODE)])


def bnum(body):
    # "Numbers unchanged" = multiset of every numeric literal (integer, decimal,
    # thousands-separated, scientific mantissa/exponent digits). This is the real
    # hard constraint (no reported beta/SE/p/OR/CI/r_g/N/ref-no changes) and is
    # insensitive to adjacent punctuation, e.g. `[1];` vs `[1].` or `(100` vs `100`.
    return Counter(re.findall(r"\d+(?:[.,]\d+)*", prose_lines(body)))


def bpunct_tokens(body):
    # Secondary, informational: whitespace tokens containing a digit (punctuation-
    # sensitive). Differences here are reported but do NOT block a replacement,
    # provided the numeric-literal multiset (bnum) is identical.
    return Counter(t for t in prose_lines(body).split() if re.search(r"\d", t))


def main():
    apply = "--apply" in sys.argv
    text = open(F, encoding="utf-8").read()
    ae, leg = body_span(text)
    body = text[ae:leg]
    before_w = bwords(body)
    before_n = bnum(body)

    newbody = body
    applied, misses, conflicts = 0, [], []
    for old, new in REPL:
        if old not in newbody:
            misses.append(old[:60])
            continue
        cand = newbody.replace(old, new, 1)
        if bnum(cand) != before_n:
            diff = (bnum(cand) - before_n, before_n - bnum(cand))
            conflicts.append((old[:55], diff))
            continue
        newbody = cand
        applied += 1

    after_w = bwords(newbody)
    after_n = bnum(newbody)
    punct_before = bpunct_tokens(body)
    punct_after = bpunct_tokens(newbody)
    punct_changed = sum((punct_before - punct_after).values())
    print(f"=== code/337 word-neutral compression {'[APPLY]' if apply else '[DRY RUN]'} ===")
    print(f"replacements applied : {applied}/{len(REPL)}")
    print(f"missing (no match)   : {len(misses)}")
    for m in misses:
        print("   MISS:", m)
    print(f"skipped (numeric-token change): {len(conflicts)}")
    for o, d in conflicts:
        print("   CONFLICT:", o, "->", d)
    print(f"body words: {before_w} -> {after_w}  (delta {after_w - before_w:+d})")
    print(f"numeric-literal multiset identical: {before_n == after_n}  "
          f"({sum(before_n.values())} numeric literals immutable)")
    print(f"punctuation-adjacent token changes (numbers equal): {punct_changed}")
    if before_n != after_n:
        print("!! GUARD FAIL - refusing to write")
        sys.exit(1)
    if apply:
        open(F, "w", encoding="utf-8").write(text[:ae] + newbody + text[leg:])
        print("WROTE", F)
    else:
        need = after_w - 4000
        print(f"DRY RUN only. Over 4000 by {need:+d} words after these edits.")


if __name__ == "__main__":
    main()
