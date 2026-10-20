# -*- coding: utf-8 -*-
import re
SRC = "docs/08_投稿稿_OA_BMD_fracture.md"
lines = open(SRC, encoding="utf-8").read().split("\n")

def words(s):
    return len(re.findall(r"[A-Za-z][A-Za-z0-9'\-]*", s))

def replace_line(anchor, new_line):
    hits = [i for i,l in enumerate(lines) if anchor in l]
    assert len(hits)==1, f"anchor not unique ({len(hits)}): {anchor!r}"
    i = hits[0]
    assert anchor in lines[i]
    lines[i] = new_line

def replace_sub(anchor_line, old, new):
    hits = [i for i,l in enumerate(lines) if anchor_line in l]
    assert len(hits)==1
    i = hits[0]
    assert old in lines[i], f"substring not found in line {i}: {old!r}"
    lines[i] = lines[i].replace(old, new, 1)

# --- L93: SOS paragraph (whole-line condense) ---
replace_line(
"OA is a disorder of the osteochondral unit",
"OA is a disorder of the osteochondral unit [27–30], in which subchondral bone remodelling feeds cartilage degradation. A systemic, polygenic bone-forming set-point readable as areal BMD would protect at the whole-bone endpoint and raise OA risk at the joint endpoint — the Systemic Osteogenic Set-point (SOS [34]), operationalised as areal BMD and used here as a heuristic, not a fitted model. If the causal quantity is volumetric density, an areal-BMD estimate is diluted by skeletal size, so volumetric-BMD (or area-conditioned) instruments should yield a BMD→OA effect at least as large as the areal one — the quantitative form of F2, which the present areal data cannot test. Height captures neither bone width nor geometry [35], so the height-independent 67–91% residual is not established as bone mass (resolved only by F2). The reading generates six predictions (P1–P6) and seven falsifying tests (F1–F7; Table S11); P1 and P4 are supported (P4 subject to the mass-versus-size split that only F2 resolves) and P6 is restated as a magnitude claim, while P2 is under-powered and P3 non-discriminating. Two tests are executed here for the first time: F6 (cross-site readout homogeneity) found no heterogeneity among the five site-specific BMD→knee-OA estimates (Q=6.92, 4 df, p=0.14, τ²=0) and F7 (joint-locus exclusion) found GDF5 not genome-wide-significant for BMD at any site, and excluding WNT16 or SOX6 changes no estimate's direction or significance."
)

# --- L95: boundaries paragraph (whole-line condense) ---
replace_line(
"Two boundaries limit what these predictions can discriminate.",
"Two boundaries limit what these predictions can discriminate. First, the tests executed here separate the set-point architecture from site-specific, local-mechanical and joint-locus rivals, but not from a generic systemic rival in which any BMD-raising process also acts at the joint; that separation requires F2–F5 (volumetric decomposition, tissue enrichment, eQTL colocalisation, osteophyte-stratified outcome), unavailable with the present data. Second, P2's site-ordering clause is under-powered (arms minus femoral neck +0.005, one-sided p=0.18) and P3 is non-discriminating (a null colocalisation is predicted under any diffuse mechanism); one pairwise contrast (lumbar spine versus legs, p=0.022) reaches nominal significance — chance-compatible under Bonferroni but not excluded. Among the pre-registered predictions the genuinely discriminating one is P4, resolved only by F2. The lowest-cost, highest-yield test is F5, pre-registered here: re-run with the OA outcome stratified into osteophyte-predominant versus joint-space-narrowing-predominant subtypes and test whether the BMD→OA estimate concentrates in the former, against the same fracture positive control; F5 needs an OA GWAS with radiographic subtypes and is pre-registered rather than performed, so we claim no result from it."
)

# --- L91: Hartley magnitude paragraph (whole-line condense) ---
replace_line(
"The directional agreement with Hartley et al.",
"The directional agreement with Hartley et al. [22] is matched by a magnitude disagreement: their knee estimate is OR 1.40 (1.20–1.63) per SD higher eBMD — 0.336 on the log-odds scale — against 0.017–0.028 here, a 12–20-fold difference. Three features make this a difference of scale and design rather than measured biology. First, the exposures differ: Hartley instrumented heel-ultrasound eBMD, scaling with bone size, whereas we instrumented DXA areal BMD, area-normalised, so the two per-SD effects are not denominated in a common unit. Second, the designs differ: Hartley's estimate is one-sample two-stage least squares on allele scores in the same individuals as the outcome (overlap 100% by construction), whereas ours is bounded at ≤2.6% and an MRlap correction moved every estimate a median 20.5% further from the null. Third, our pipeline is calibrated: the BMD→fracture positive control returns OR 1.40–1.42 against GEFOS 1.55, so a systematic deflation would have to be specific to the OA outcome. Inverting Hartley's reverse estimate (β 4.16 SD per doubling of knee-OA risk, implying OR ≈1.18 per SD higher eBMD) leaves a residual gap of roughly five-fold. We read the directions as concordant and the magnitudes as unresolved, and record the disagreement as a limitation."
)

# --- L89: Zhou/Hartley/Funck paragraph (whole-line condense) ---
replace_line(
"Zhou et al. [5] reported higher femoral-neck BMD",
"Zhou et al. [5] reported higher femoral-neck BMD protected against hip OA, whereas our ten-site analysis shows higher BMD raises knee-OA risk — a hip-versus-knee contrast cautioning against generalising any single BMD→OA estimate across joints. Hartley et al. [22] reported a positive BMD effect on both hip (OR 1.28, 1.05–1.57) and knee OA (OR 1.40, 1.20–1.63); Funck-Brentano et al. [36] likewise reported a positive overall effect (all-OA OR 1.14, 1.06–1.22), significant at both knee and hip; a quantitative comparison on the native scales is in Fig. S12. These contradict the directional form of prediction P6, motivated by Zhou's protective femoral-neck finding; P6 and falsification rule F1 are in Table S11 and adjudicated in §4. We adjudicate rather than cite selectively. First, the hip evidence disagrees with itself (protective [5], harmful [22], [36]). Second, decisively for F1, the criterion is conjunctive: the hip effect must be positive AND of knee-like magnitude. Within Hartley et al., the only same-pipeline comparison, the hip and knee intervals overlap almost completely (1.05–1.57 vs 1.20–1.63) and on the log-odds scale the hip estimate lies 27% below the knee's (0.247 vs 0.336); Hartley establishes positivity, not equality. We therefore restate P6 as a magnitude claim and pre-register the decision rule (a hip effect positive and at least as large as the knee's would fail P6); the directional form was pre-registered, its restatement here deviates from that form and is prompted by the contradictory hip evidence, and the restated form scores only prospectively collected hip estimates."
)

# --- L69: MR results paragraph (whole-line condense) ---
replace_line(
"MR supported a positive effect of BMD on knee OA",
"MR supported a positive effect of BMD on knee OA across all five exposures (Table 2): total-body β=+0.023 (p=3.9×10⁻¹⁰), lumbar-spine +0.028 (p=7.7×10⁻⁹), arms +0.022 (p=1.7×10⁻⁶), legs +0.018 (p=2.3×10⁻⁵), femoral-neck +0.017 (p=6.4×10⁻⁴) — OR 1.017–1.028 per SD higher BMD, against OR 1.40 (1.20–1.63) per SD in the largest one-sample estimate [22], a 12–20-fold log-odds difference (adjudicated §4). MR-Egger intercepts were all p>0.24 and weighted medians matched the IVW sign at all five sites, with bootstrap 95% CIs excluding the null at three of five; with 16–31 instruments the intercept test is under-powered, so a null intercept excludes unbalanced pleiotropy of detectable size rather than pleiotropy itself. Leave-one-out produced no sign reversals and directions agreed with r_g (+0.12 to +0.25). Heterogeneity was high (I²=68–86%); under random-effects IVW only two of five exposures remained significant (total-body p=0.006, lumbar-spine p=0.017; a second random-effects convention additionally retains arms, p=0.033), so the count depends on the random-effects convention and we base the causal reading on estimator concordance rather than fixed-effect p-values. All five effects remained positively causal after MR-PRESSO outlier removal, with four of five strengthening. The positive control was comparably heterogeneous (I² 78% and 83%), indicating high I² is a property of the BMD instrument sets rather than an OA-specific IVW failure; a multi-estimator robustness overview of all ten pairs is in Supplementary Fig. S11."
)

# --- L99: Strengths and limitations (whole-line condense) ---
replace_line(
"**Strengths and limitations.** Strengths:",
"**Strengths and limitations.** Strengths: ten BMD sites with two OA definitions, an internal BMD–fracture positive control calibrated against GEFOS (OR 1.40–1.42 vs 1.55), dual FDR/Bonferroni correction and directional concordance across IVW, MR-Egger, weighted median and MR-PRESSO. Limitations: the BMD instruments are European-only whereas the knee-OA outcome is 87.31% multi-ancestry (Table S14) and LDSC uses a European-only reference, so the forward estimates generalise most directly to European-ancestry populations (the near-1.0 intercept does not exclude ancestry-mismatch bias). Five of ten prespecified BMD sites entered MR because strict clumping completed for them; the five deferred sites showed comparable genome-wide signal, so the analysed subset reflects compute availability rather than biology (Table S12). Areal BMD conflates size with mass [18,32,37]; height-conditioning bounds the size share at 9.5–33.2% (a lower bound), leaving the 67–91% height-independent residual unresolved into bone mass versus geometry (test F2), and because conditioning was on height not BMI a BMI-mediated component is not excluded [35]. All estimates are sex-averaged and effect modification by sex is untested (deposited statistics are not sex-stratified). Heterogeneity is high across all pairs including the positive control, so the causal claim rests on estimator concordance rather than fixed-effect p-values; the reverse and OA→fracture estimates are null but weakly powered (0.14–0.26), and residual overlap is bounded at ≤2.6% for BMD→OA but cannot be removed with summary data. The tests executed here (F6, F7) separate the set-point from site-specific and joint-locus rivals but not from a generic systemic-pleiotropy rival, requiring the external-data tests F2–F5."
)

# --- L87: drop the redundant height-independent sentence ---
replace_sub(
"The opposing-endpoint pattern is an old clinical observation",
" The height-independent residual (67–91%) is not identified as bone mass; its mass-versus-geometry split is resolved only by test F2.",
""
)

open(SRC, "w", encoding="utf-8").write("\n".join(lines))

# ---- re-measure ----
t = "\n".join(lines)
ai = t.find("**Objective.**")
ci = t.find("**Conclusions.**", ai)
tail = t[ci:]
em = re.search(r"\n---", tail)
abstract_end = ci + (em.start() if em else len(tail))
leg = t.find("## Figure legends")
body = t[abstract_end:leg]
body_lines = [l for l in body.splitlines() if not l.strip().startswith("|")]
print("MAIN BODY words:", words("\n".join(body_lines)))
print("abstract words:", words(t[ai:abstract_end]))

# gatekeeper phrase checks
checks = {
    "height-independent": t.count("height-independent"),
    "prespecified canonical DXA regions": t.count("prespecified canonical DXA regions"),
    "a site enters MR only if": t.count("a site enters MR only if"),
    "9.5–33.2%": t.count("9.5–33.2%"),
    "67–91%": t.count("67–91%"),
    "20.5%": t.count("20.5%"),
    "0.799": t.count("0.799"),
    "87.31%": t.count("87.31%"),
    "12–20-fold": t.count("12–20-fold"),
    "OR 1.017–1.028": t.count("OR 1.017–1.028"),
}
print("GATEKEEPERS:", checks)
