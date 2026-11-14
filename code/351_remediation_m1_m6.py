#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
code/351 — M1–M6 remediation edits for docs/08 (对抗性续轮审稿 docs/149).

Convention (matches code/338 / code/349):
  - DRY RUN by default. Pass --apply to write.
  - Backs up docs/08 to .workbuddy/backups/ before any write.
  - Each replacement must match EXACTLY ONCE in the source; if not,
    it is reported and SKIPPED (no partial corruption).
  - Word delta is estimated (whitespace token count) for budgeting
    against the 4000-word body gate (code/148).
"""
import sys, os, shutil, datetime

SRC = "docs/08_投稿稿_OA_BMD_fracture.md"
BACKUP_DIR = ".workbuddy/backups"

def wc(s):
    return len(s.split())

# (label, old, new)
REPLS = [
    # ---------- M1 : effect-scale reconstruction (L57) ----------
    ("M1_L57_placeholder",
     "the deposited BMD \u03b2 and SE are on a placeholder N=400,000 scale; restoring the verified N=31,873 requires multiplying both by \u221a(400,000/31,873)=3.543 (per-SNP z invariant).",
     "the deposited BMD \u03b2 and SE are on a placeholder N=400,000 scale (all 691,715 rows carry N=400,000; a deposition convention, not the analysed N); restoring the verified N=31,873 requires multiplying both by \u221a(400,000/31,873)=3.543 (per-SNP z invariant). The absolute scale is anchored to the fracture positive control, not variant-level-validated: the reconstructed BMD\u2192fracture OR (1.40) matches the independent GEFOS estimate (1.55) within a 23.1% log-odds gap and excludes a 3.5\u00d7 over-inflation counterfactual (OR 3.30), and all five BMD\u2192OA ORs survive \u00b120% scale perturbation (Table S34)."),

    # ---------- M1 : GEFOS honesty (L77) ----------
    ("M1_L77_gefsos",
     "On the verified per-SD scale, the BMD\u2192fracture OR (1.42 femoral-neck) lay 8.4% below GEFOS (1.55) on the odds-ratio scale and 20.0% below on the log-odds scale, with overlapping CIs (we compare on the log scale). The total-body OR 1.40 was 9.6% below on the ratio and 23.1% on the log-odds scale (cross-site approximation).",
     "On the verified per-SD scale, the BMD\u2192fracture OR (1.42 femoral-neck) lay 8.4% below GEFOS (1.55) on the odds-ratio scale and 20.0% below on the log-odds scale, with overlapping CIs (we compare on the log scale). The total-body OR 1.40 was 9.6% below on the ratio and 23.1% on the log-odds scale (cross-site approximation). The 20\u201323% log-odds gap is not sampling error alone (distinct BMD sub-phenotypes and munge pipelines) and bounds, rather than certifies, the absolute scale; a variant-level cross-check against an independent BMD GWAS remains the definitive validation and is outstanding (Table S34)."),

    # ---------- M2 : title (L1) ----------
    ("M2_title",
     "Genetically higher BMD raises knee osteoarthritis risk yet protects against fracture on pleiotropy-robust evidence",
     "Genetically higher BMD raises knee osteoarthritis risk yet protects against fracture on pleiotropy-aware Mendelian-randomisation evidence"),

    # ---------- M2 : Highlight L13 (coloc endpoint-specific) ----------
    ("M2_HL13",
     "- BMD variants colocalise with fracture (9/31), not knee OA (0/31): endpoint-specific.",
     "- BMD variants colocalise with fracture (9/31) but not with knee OA (0/31) \u2014 a split driven by colocalisation power, not established endpoint-specificity (the strongest instrument gives expected knee-OA z\u22480.85, needing ~5\u00d710\u2077 participants for genome-wide significance)."),

    # ---------- M2 : Highlight L15 (set-point) ----------
    ("M2_HL15",
     "- A systemic bone-forming set-point is a testable hypothesis consistent with data.",
     "- A systemic bone-forming set-point is a hypothesis whose directional prediction is consistent with data, but whose core claim \u2014 one variant set fixing opposite endpoint signs \u2014 is untested (F2 volume arm, F10)."),

    # ---------- M2/M6 : L85 opening ----------
    ("M2_L85",
     "on pleiotropy-robust evidence, genetically higher BMD raises knee-OA risk",
     "on pleiotropy-aware evidence, genetically higher BMD raises knee-OA risk"),

    # ---------- M2/M6 : L93 SOS kernel (two sub-edits) ----------
    ("M2_L93a",
     "SOS is thus a falsifiable causal assertion, not a restatement of the pattern.",
     "SOS is thus a falsifiable hypothesis, not a demonstrated mechanism."),
    ("M2_L93b",
     "and distinct variant sets or concordant signs would refute it (Fig. S23).",
     "and distinct variant sets or concordant signs would refute it (Fig. S23); the colocalisation split (9/31 fracture vs 0/31 knee-OA) is a power phenomenon, not sign-opposition evidence (expected knee-OA z\u22480.85; genome-wide significance needs ~5\u00d710\u2077 participants), and the set-point core is untested (F2 volume arm, F10)."),

    # ---------- M3 : L55 MRlap honesty ----------
    ("M3_L55",
     "The joint correction moved every estimate from the null (median +20.5%; observed-versus-corrected difference p<0.05 for 9/10), showing towards-null components dominate; no direction or significance reversed.",
     "The joint correction moved every estimate from the null (median +20.5%; observed-versus-corrected difference p<0.05 for 9/10); MRlap corrects weak-instrument bias and winner's curse but not horizontal pleiotropy, so this shift also accommodates residual directional pleiotropy \u2014 'pleiotropy-robust' is read as 'pleiotropy-aware' (the femoral-neck control's Egger p=0.018 shows directional pleiotropy is present). No direction or significance reversed."),

    # ---------- M3 : L69 Egger ----------
    ("M3_L69",
     "with 16\u201331 instruments the intercept test has limited power, excluding only detectable unbalanced pleiotropy.",
     "with 16\u201331 instruments the intercept test has limited power and excludes only detectable unbalanced pleiotropy, not pleiotropy per se."),

    # ---------- M3/M4 : L75 PRESSO + multiplicity pointer ----------
    ("M3_L75",
     "All significant findings survived Bonferroni correction (\u03b1=0.005). The IVW estimate is therefore not read as a single clean causal parameter (Table 2, Fig. S11, S20).",
     "All significant findings survived Bonferroni correction (\u03b1=0.005) across the ten primary MR pairs. The IVW estimate is therefore not read as a single clean causal parameter (Table 2, Fig. S11, S20); the 9/10 significant PRESSO distortion tests indicate widespread horizontal pleiotropy whose non-outlier component is not removed by outlier exclusion, and only the ten primary MR pairs carried \u03b1=0.005 protection (see multiplicity note in Limitations)."),

    # ---------- M4 : L99 limitations multiplicity ----------
    ("M4_L99",
     "Limitations: the European-only BMD instruments and LDSC reference against an 87.31%-European knee-OA outcome mean forward estimates generalise most directly to European-ancestry populations (Table S14).",
     "Limitations: the European-only BMD instruments and LDSC reference against an 87.31%-European knee-OA outcome mean forward estimates generalise most directly to European-ancestry populations (Table S14). Multiplicity was not fully controlled: of the 78 LDSC pairs, 31 colocalisation tests per outcome, MVMR height/BMI grids and ten PRESSO distortion tests, only the ten primary MR pairs carried \u03b1=0.005 Bonferroni protection, so the colocalisation 9/31 and random-effects 2/5 findings are descriptive rather than multiplicity-adjusted."),

    # ---------- M5 : L61 MVMR collider caveat ----------
    ("M5_L61",
     "Height and BMI were conditioned separately, each adjusted one at a time.",
     "Height and BMI were conditioned separately, each adjusted one at a time. Because height and BMI are themselves skeletal-system traits, conditioning on them can act as collider conditioning on the BMD\u2192OA path and induce bias; we therefore read MVMR attenuation as a sensitivity bound on size mediation, not a deconfounded effect."),

    # ---------- M5 : L81 residual not bone mass (strengthen) ----------
    ("M5_L81",
     "Because height indexes stature not bone width/geometry [28], this residual is not identified as bone mass; F2 probes geometry (\u00a74, Table S28).",
     "Because height indexes stature not bone width/geometry [28], this residual is not identified as bone mass; it may reflect residual size, geometry or other skeletal pathways rather than a clean non-size component, and F2 probes geometry (\u00a74, Table S28)."),

    # ---------- M5 : Abstract L27 soften ----------
    ("M5_abs27",
     "Conditioning on height attenuated the OA effect by 9.5\u201333.2% (a lower bound for size); conditioning separately on BMI left it positive and significant at all five sites (\u22123.7% to 46.0% attenuation), so the 67\u201391% height-independent residual is neither height- nor BMI-mediated \u2014 not identified as bone mass.",
     "Conditioning on height attenuated the OA effect by 9.5\u201333.2% (a lower bound for size); conditioning separately on BMI left it positive and significant at all five sites (\u22123.7% to 46.0% attenuation), so 67\u201391% of the joint effect is independent of height conditioning and not BMI-mediated in MVMR, though not resolved into bone mass versus other skeletal pathways."),

    # ---------- M2 : L101 Conclusion ----------
    ("M2_L101",
     "**Conclusion.** Genetically higher BMD protects against fracture and, on the pleiotropy-robust evidence, raises knee-OA risk (2/5 sites significant under random-effects IVW).",
     "**Conclusion.** Genetically higher BMD protects against fracture and, on the pleiotropy-aware evidence, raises knee-OA risk (2/5 sites significant under random-effects IVW)."),

    # ---------- M6 : S18 legend over-claim ----------
    ("M6_S18",
     "confirming the instrument set shares causal variants with bone but not joint",
     "consistent with the instrument set sharing causal variants with bone but not joint under current power"),

    # ---------- TOC : add S34 entry (outside body count) ----------
    ("TOC_S34",
     "**Table S33.** Section-2 review remediation audit: overlap-bias denominators and MRlap positive-control transfer (P1-2), EUR-only ancestry sensitivity (P1-3), K bootstrap and dual-scale sensitivity (P1-4), and reverse-MR instrument-set attrition (P2-2).",
     "**Table S33.** Section-2 review remediation audit: overlap-bias denominators and MRlap positive-control transfer (P1-2), EUR-only ancestry sensitivity (P1-3), K bootstrap and dual-scale sensitivity (P1-4), and reverse-MR instrument-set attrition (P2-2).\n- **Table S34.** Effect-scale reconstruction independent validation (M1): placeholder-N file signature (691,715 rows at N=400,000), GEFOS positive-control concordance (reconstructed OR 1.40 vs 1.55; 23.1% log-odds gap; 3.5\u00d7 over-inflation counterfactual ruled out), and \u00b120% scale-sensitivity of all five BMD\u2192OA MR estimates (code/350; Fig. M1 in Figures_v18_Suppl/)."),

    # ---------- line 636 summary: extend table range ----------
    ("SUM_636",
     "supplementary tables S8\u2013S17, S25\u2013S32 and the STROBE-MR checklist",
     "supplementary tables S8\u2013S17, S25\u2013S34 and the STROBE-MR checklist"),

    # ================= COMPRESSION CUTS (fund word budget) =================
    # C1 : L91 trim
    ("C1_L91",
     "This is a difference of scale and design, not biology: Hartley instrumented heel-ultrasound eBMD (bone size), we DXA areal BMD; Hartley's one-sample 2SLS (overlap 100%) contrasts with ours bounded at \u22642.6%, MRlap moving every estimate a median 20.5% from the null. Directions concord, magnitudes unresolved.",
     "This is a difference of scale and design, not biology (Hartley's heel-ultrasound eBMD and one-sample 2SLS versus our DXA areal BMD and \u22642.6% overlap). Directions concord, magnitudes unresolved."),

    # C2 : L89 trim
    ("C2_L89",
     "In the only same-pipeline comparison [12], hip and knee intervals overlap (1.05\u20131.57 vs 1.20\u20131.63), the hip log-odds 27% below the knee's (0.247 vs 0.336) \u2014 positivity, not equality. Our cross-site hip estimate (same study/instruments) resolves F1: the same 31 instruments give BMD on hip OA IVW \u03b2=+0.097 (OR 1.102, 1.068\u20131.136 per SD), positive and knee-like (hip minus knee z=0.94, p=0.35; Table S27), restoring P6's pre-registered direction.",
     "In the only same-pipeline comparison [12], hip and knee intervals overlap (1.05\u20131.57 vs 1.20\u20131.63); our cross-site hip estimate (same 31 instruments) gives BMD\u2192hip-OA IVW \u03b2=+0.097 (OR 1.102, 1.068\u20131.136), positive and knee-like (hip\u2212knee z=0.94, p=0.35; Table S27), restoring P6's pre-registered direction."),

    # C3 : S30 point (vi) trim
    ("C3_S30vi",
     " (vi) The formation-side SOS reading gains endotype and mechanistic support. Genetically distinct hypertrophic (osteophyte-rich) OA endotypes [46] corroborate an osteophyte-dominant joint pathway, and the RSPO3\u2013NFATC1 subchondral-bone switch [47] \u2014 whereby genetic liability to low bone mass protects against OA by increasing subchondral compliance \u2014 supplies, by reverse inference, the mechanistic anchor for a high-bone set-point as OA-promoting. A companion five-layer OA\u2013eBMD\u2013fracture funnel [48] (the authors' own unpublished deposit) reports a reverse OA\u2194bone signal under eBMD instruments; this bounds rather than refutes the constitutional-asymmetry reading of P5, which here concerns DXA areal BMD and fracture under low power (0.14\u20130.26).",
     " (vi) The formation-side SOS reading gains endotype support (hypertrophic osteophyte-rich OA endotypes [46]; the RSPO3\u2013NFATC1 subchondral-bone switch [47]); a companion OA\u2013eBMD\u2013fracture funnel [48] reports a reverse OA\u2194bone signal under eBMD instruments, bounding rather than refuting the low-power P5 reading."),

    # C4 : S28 Part A interpretation trim
    ("C4_S28",
     "**Interpretation of Part A.** The pooled geometry arm is an order of magnitude smaller than the BMD arm and carries no consistent direction (hip and knee pooled estimates have opposite signs), so the BMD->OA effect is not accounted for by the hip-geometry axes tested. Individual modes, however, are heterogeneous (I2 75-95%) and several show significant MR-Egger intercepts, so the per-mode estimates are pleiotropy-affected and are not read causally; the single-variable design is informative rather than a clean falsification, and it is reported as such. The conditional BMD-plus-hip-shape MVMR in Part B now supplies the formal test and reaches the same qualitative conclusion for knee OA (and is inconclusive for hip OA). The literal volumetric-instrument form of F2 (areal BMD replaced by a volume- or area-conditioned BMD instrument) remains untested because no fit-for-purpose public common-SNP summary-statistics file for a volumetric BMD trait was found in the GWAS Catalog; the manuscript's height conditioning (9.5-33.2% attenuation) and BMI conditioning already cover the stature and adiposity axes.",
     "**Interpretation of Part A.** The pooled geometry arm is an order of magnitude smaller than the BMD arm and direction-inconsistent (hip and knee pooled estimates have opposite signs), so hip geometry does not account for the BMD->OA effect. Individual modes are heterogeneous (I2 75\u201395%) with significant MR-Egger intercepts, so per-mode estimates are pleiotropy-affected and not read causally; the single-variable design is informative, not a clean falsification. The conditional MVMR in Part B formalises the test (knee OA positive; hip OA inconclusive). The volumetric-instrument form of F2 remains untested (no public volumetric-BMD summary statistics)."),

    # C5 : S25 intro trim
    ("C5_S25",
     "Tissue specificity of the 13-phenotype genetic architecture was tested by a competitive gene-set enrichment with permutation calibration (MAGMA-style statistic; code/342): each GWAS SNP was assigned to the nearest gene body within \u00b150 kb; gene-level association = \u2212log10 of the minimum p across assigned SNPs; each GO set's mean was compared with background, and significance was calibrated by 10,000 permutations matched on gene length and chromosome (the previous one-sided Mann\u2013Whitney p was optimistic because the min-p statistic depends on gene length/SNP density and sets overlap). Positive enrichment \u0394 means the set's genes are more associated than background. Because cartilage_development and chondrocyte_differentiation overlap (Jaccard 0.44) and the broad skeletal-system-development set is 67% shared with cartilage sets, we additionally tested a direct cartilage-minus-bone contrast with overlapping genes removed. Bone-side sets: bone development, ossification. Mixed skeletal set (not bone-specific): skeletal-system development \u2014 28/42 genes shared with cartilage sets. Cartilage-side sets: cartilage development, chondrocyte differentiation, extracellular-matrix organisation, skeletal-joint morphogenesis. Source GWAS: OA knee, total-body BMD, any-fracture (sumstats in OA_FLAGSHIP/07_reference/ldsc; code/171, code/342).",
     "Tissue specificity was tested by competitive gene-set enrichment with permutation calibration (code/342): SNPs were assigned to the nearest gene (\u00b150 kb); gene-level association = \u2212log10 of the minimum p across assigned SNPs; each GO set's mean was compared with background, calibrated by 10,000 permutations matched on gene length and chromosome. Because cartilage_development and chondrocyte_differentiation overlap (Jaccard 0.44) and skeletal-system-development is 67% shared with cartilage sets, a direct cartilage-minus-bone contrast (overlapping genes removed) was also tested. Bone-side sets: bone development, ossification; mixed: skeletal-system development (28/42 shared with cartilage); cartilage-side: cartilage development, chondrocyte differentiation, extracellular-matrix organisation, skeletal-joint morphogenesis. Sources: OA knee, total-body BMD, any-fracture."),

    # C6 : Abstract Design trim
    ("C6_abs25",
     "**Design.** LDSC and two-sample MR were applied to 13 GWAS summary-statistic datasets (two OA definitions, ten regional BMD sites, any-fracture). LDSC covered all 78 phenotype pairs; MR ran in eight phenotypes with clumped instruments \u2014 five of ten prespecified BMD sites, both OA definitions and any-fracture \u2014 using BMD\u2192fracture as the positive control. Further configurations conditioned on height and separately on BMI, applied MRlap sample-overlap correction, and reverse MR.",
     "**Design.** LDSC and two-sample MR were applied to 13 GWAS summary-statistic datasets. LDSC covered all 78 phenotype pairs; MR ran in eight phenotypes with clumped instruments (five of ten BMD sites, both OA definitions, any-fracture), using BMD\u2192fracture as the positive control. Further configurations conditioned on height and BMI separately, applied MRlap overlap correction, and reverse MR."),

    # D1 : S30 intro trim
    ("D1_S30",
     "Genetic correlations (r_g) between total-body DXA BMD (GCST90568455; verified N = 31,986) and the four site-stratified OA outcomes \u2014 hip, knee, spine and hand (Hatzikotoulas et al. 2025, GWAS Catalog GCST90566798/6800/6801/6803, GRCh37) \u2014 were estimated with LD Score regression (ldsc v2.0.1, Bulik-Sullivan et al. 2015 [5]) using European-ancestry 1000 Genomes Phase 3 HapMap3 LD Scores (eur_w_ld_chr_1\u201322; 1,290,028 regression SNPs; the HLA region excluded). All traits were munged against the same w_hm3.snplist (1,217,311 SNPs), so the regression SNP set is identical across traits (685,394\u2013685,477 SNPs per BMD\u2013OA pair). Because r_g = b_gencov / \u221a(b_h2,1 \u00b7 b_h2,2) is a ratio of LD Score regression slopes, a global rescaling of N cancels exactly; the verified N is therefore required only before quoting h\u00b2, not r_g.",
     "Genetic correlations (r_g) between total-body DXA BMD (GCST90568455; verified N = 31,986) and the four site-stratified OA outcomes (hip, knee, spine, hand; Hatzikotoulas et al. 2025, GCST90566798/6800/6801/6803, GRCh37) were estimated by LD Score regression (ldsc v2.0.1 [5]) on European-ancestry 1000 Genomes Phase 3 HapMap3 LD Scores (eur_w_ld_chr_1\u201322; 1,290,028 regression SNPs; HLA excluded). All traits were munged against the same w_hm3.snplist (1,217,311 SNPs), so the regression SNP set is identical across traits. Because r_g is a ratio of LD Score regression slopes, a global rescaling of N cancels exactly; the verified N is required only before quoting h\u00b2, not r_g."),

    # D2 : S29 notes trim
    ("D2_S29",
     "*Notes.* n = 31 BMD_total instruments per site. The spine and hand (osteophyte-dominant) ORs are numerically larger than the knee (cartilage-dominant) OR by \u22484\u20134.5% on the odds scale (spine\u2212knee log-OR difference = +0.043; hand\u2212knee = +0.039). All three MR-Egger intercepts are non-significant, so the gradient is not driven by directional pleiotropy; confidence intervals overlap broadly and the differences are not individually significant.",
     "*Notes.* n = 31 BMD_total instruments per site; spine and hand (osteophyte-dominant) ORs exceed the knee (cartilage-dominant) OR by \u22484\u20134.5% on the odds scale (spine\u2212knee log-OR +0.043; hand\u2212knee +0.039). All three MR-Egger intercepts are non-significant (gradient not driven by directional pleiotropy), but CIs overlap and the differences are not individually significant."),

    # D4 : S27 verification trim
    ("D4_S27",
     "**Verification against the published knee-OA row.** Re-deriving the knee-OA estimate from the raw GWAS Catalog file reproduces the manuscript's Table 2 values: IVW beta +0.078 (this table +0.078), OR 1.081 (1.055-1.108), I2 81% (80.8%) and weighted median +0.065 (+0.027 to +0.177) (+0.066, +0.027 to +0.174). The MR-Egger intercept is non-significant in both (this table p = 0.63; manuscript p = 0.810).",
     "**Verification against the published knee-OA row.** Re-deriving the knee-OA estimate from the raw GWAS Catalog file reproduces Table 2: IVW beta +0.078 (OR 1.081, 1.055\u20131.108), I2 81% (80.8%), weighted median +0.065 (+0.027 to +0.177); MR-Egger intercept non-significant (this table p = 0.63; manuscript p = 0.810)."),

    # E1 : S32 status trim
    ("E1_S32",
     "**Status.** The registered F10 cross-ancestry replication is **not executed**: Hatzikotoulas 2025 released site-stratified (hip/knee/spine/hand) but no ancestry-stratified OA summary statistics, and the previously staged `GCST90566811` is the X-chromosome knee-OA GWAS (0/31 autosomal BMD instruments matched), not a European-restricted autosomal outcome. This table reports a narrower, executable question: is the pooled BMD\u2192knee-OA estimate driven by the variants whose summary statistics were contributed by the most ancestry-diverse set of cohorts?",
     "**Status.** The registered F10 cross-ancestry replication is **not executed**: Hatzikotoulas 2025 released site-stratified but no ancestry-stratified OA summary statistics, and the previously staged `GCST90566811` is the X-chromosome knee-OA GWAS (0/31 autosomal BMD instruments matched), not a European-restricted autosomal outcome. This table asks a narrower question: is the pooled BMD\u2192knee-OA estimate driven by variants from the most ancestry-diverse cohorts?"),
]

def main():
    if not os.path.exists(SRC):
        print("ERROR: source not found:", SRC); sys.exit(2)
    with open(SRC, encoding="utf-8") as f:
        text = f.read()

    total_delta = 0
    print("=== DRY/RUN preview ===")
    for label, old, new in REPLS:
        n = text.count(old)
        d = wc(new) - wc(old)
        total_delta += d
        flag = "OK" if n == 1 else ("FAIL(n=0)" if n == 0 else f"WARN(n={n})")
        print(f"[{flag}] {label}: delta_words={d:+d} (old={wc(old)} new={wc(new)})")
    print(f"--- projected net word delta in edited regions: {total_delta:+d} ---")

    if "--apply" not in sys.argv:
        print("\nDRY RUN: no file written. Re-run with --apply to commit.")
        sys.exit(0)

    # validate all match exactly once before writing
    bad = [(l, text.count(o)) for l, o, n in REPLS if text.count(o) != 1]
    if bad:
        print("ABORT: non-unique matches:", bad); sys.exit(3)

    os.makedirs(BACKUP_DIR, exist_ok=True)
    ts = datetime.datetime.now().strftime("%Y-%m-%d_%H%M%S")
    bak = f"{BACKUP_DIR}/docs08_pre351_{ts}.md"
    shutil.copy(SRC, bak)
    print("backup ->", bak)

    out = text
    for label, old, new in REPLS:
        out = out.replace(old, new)
    with open(SRC, "w", encoding="utf-8") as f:
        f.write(out)
    print("APPLIED", len(REPLS), "replacements.")

if __name__ == "__main__":
    main()
