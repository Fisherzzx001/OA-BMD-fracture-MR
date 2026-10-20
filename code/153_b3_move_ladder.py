# -*- coding: utf-8 -*-
"""B3: move SOS prediction-falsification ladder forward into main text (§4).
Replaces scattered L93/L95 ladder mentions with one enumerated prose subsection,
and trims L91/L99 redundancy, keeping main-body words <=4000 and display items <=8.
All gatekeeper phrases (height-independent, 9.5-33.2%, 67-91%, 20.5%, 87.31%,
12-20-fold, OR 1.017-1.028) are preserved.
"""
import re, sys
F = "docs/08_投稿稿_OA_BMD_fracture.md"
t = open(F, encoding="utf-8").read()

def replace_once(old, new, label):
    global t
    n = t.count(old)
    assert n == 1, f"[{label}] expected 1 match, found {n}"
    t = t.replace(old, new, 1)
    print(f"  [{label}] replaced ({len(old.split())} -> {len(new.split())} words approx)")

# ---- OP_A: trim L93 trailing ladder-status sentence to a pointer ----
old_a = ("The reading generates six predictions (P1\u2013P6) and seven falsifying tests "
         "(F1\u2013F7; Table S11); P1 and P4 are supported (P4 subject to the mass-versus-size "
         "split that only F2 resolves) and P6 is restated as a magnitude claim, while P2 is "
         "under-powered and P3 non-discriminating. Two tests are executed here for the first "
         "time: F6 (cross-site readout homogeneity) found no heterogeneity among the five "
         "site-specific BMD\u2192knee-OA estimates (Q=6.92, 4 df, p=0.14, \u03c4\u00b2=0) and F7 "
         "(joint-locus exclusion) found GDF5 not genome-wide-significant for BMD at any site, "
         "and excluding WNT16 or SOX6 changes no estimate's direction or significance.")
new_a = ("The reading generates six predictions (P1\u2013P6) and seven falsifying tests "
         "(F1\u2013F7); the full ladder, with designs and current judgement, is set out below "
         "and in Table S11.")
replace_once(old_a, new_a, "OP_A")

# ---- OP_B: replace L95 with consolidated enumerated ladder subsection ----
old_b = ("Two boundaries limit what these predictions can discriminate. First, the tests "
         "executed here separate the set-point architecture from site-specific, local-mechanical "
         "and joint-locus rivals, but not from a generic systemic rival in which any BMD-raising "
         "process also acts at the joint; that separation requires F2\u2013F5 (volumetric "
         "decomposition, tissue enrichment, eQTL colocalisation, osteophyte-stratified outcome), "
         "unavailable with the present data. Second, P2's site-ordering clause is under-powered "
         "(arms minus femoral neck +0.005, one-sided p=0.18) and P3 is non-discriminating (a null "
         "colocalisation is predicted under any diffuse mechanism); one pairwise contrast (lumbar "
         "spine versus legs, p=0.022) reaches nominal significance \u2014 chance-compatible under "
         "Bonferroni but not excluded. Among the pre-registered predictions the genuinely "
         "discriminating one is P4, resolved only by F2. The lowest-cost, highest-yield test is "
         "F5, pre-registered here: re-run with the OA outcome stratified into osteophyte-predominant "
         "versus joint-space-narrowing-predominant subtypes and test whether the BMD\u2192OA estimate "
         "concentrates in the former, against the same fracture positive control; F5 needs an OA "
         "GWAS with radiographic subtypes and is pre-registered rather than performed, so we claim "
         "no result from it.")
new_b = ("**Pre-registered SOS predictions and falsifying tests.** *Predictions (derived).* "
         "**P1** opposing endpoints \u2014 supported (positive-control OR 1.40 per SD lower BMD; "
         "BMD\u2192knee-OA 5/5 positive). **P2** systemic not local-mechanical \u2014 systemic clause "
         "supported (non-weight-bearing arms \u03b2 = +0.022, p = 1.7\u00d710\u207b\u2076), ordering "
         "clause under-powered (arms \u2212 femoral neck +0.005, p = 0.18). **P3** diffuse polygenic "
         "\u2014 0/31 loci colocalise (PP.H4 \u2265 0.8); direction-consistent, non-discriminating. "
         "**P4** size a separable partial component \u2014 height attenuation 9.5\u201333.2% (fracture "
         "control \u22121.7% to \u22123.0%); 67\u201391% residual unresolved into mass versus geometry (F2). "
         "**P5** constitutional asymmetry \u2014 reverse null (power 0.14\u20130.26), consistent but not "
         "evidence. **P6** site-pathway dominance \u2014 weakened to a magnitude claim (hip evidence "
         "contradictory). *Falsifying tests.* **F1** hip-OA replication; **F2** volumetric/size "
         "decomposition; **F3** bone-versus-cartilage tissue enrichment (LDSC-SEG); **F4** eQTL "
         "colocalisation; **F5** osteophyte-stratified outcome \u2014 F1\u2013F5 prospective, requiring "
         "external data, F5 pre-registered here as the immediate next step and not performed. "
         "**F6** cross-site homogeneity (executed: Q = 6.92, 4 df, p = 0.14, \u03c4\u00b2 = 0) and "
         "**F7** joint-locus exclusion (executed: GDF5 not a BMD instrument at any site; excluding "
         "WNT16 or SOX6 changes no estimate's direction or significance) both passed. P4 is the "
         "genuinely discriminating prediction, resolved only by F2; F1\u2013F5 are required to "
         "separate the set-point from a generic systemic-pleiotropy rival.")
replace_once(old_b, new_b, "OP_B")

# ---- OP_C: trim L91 (Hartley magnitude disagreement) ----
old_c = ("The directional agreement with Hartley et al. [22] is matched by a magnitude "
         "disagreement: their knee estimate is OR 1.40 (1.20\u20131.63) per SD higher eBMD \u2014 "
         "0.336 on the log-odds scale \u2014 against 0.017\u20130.028 here, a 12\u201320-fold difference. "
         "Three features make this a difference of scale and design rather than measured biology. "
         "First, the exposures differ: Hartley instrumented heel-ultrasound eBMD, scaling with bone "
         "size, whereas we instrumented DXA areal BMD, area-normalised, so the two per-SD effects are "
         "not denominated in a common unit. Second, the designs differ: Hartley's estimate is "
         "one-sample two-stage least squares on allele scores in the same individuals as the outcome "
         "(overlap 100% by construction), whereas ours is bounded at \u22642.6% and an MRlap correction "
         "moved every estimate a median 20.5% further from the null. Third, our pipeline is calibrated: "
         "the BMD\u2192fracture positive control returns OR 1.40\u20131.42 against GEFOS 1.55, so a "
         "systematic deflation would have to be specific to the OA outcome. Inverting Hartley's "
         "reverse estimate (\u03b2 4.16 SD per doubling of knee-OA risk, implying OR \u22481.18 per SD "
         "higher eBMD) leaves a residual gap of roughly five-fold. We read the directions as concordant "
         "and the magnitudes as unresolved, and record the disagreement as a limitation.")
new_c = ("The directional agreement with Hartley et al. [22] is matched by a magnitude disagreement: "
         "their knee estimate is OR 1.40 (1.20\u20131.63) per SD higher eBMD against 0.017\u20130.028 here, "
         "a 12\u201320-fold difference. Three features make this a difference of scale and design rather "
         "than measured biology. First, the exposures differ: Hartley instrumented heel-ultrasound eBMD, "
         "scaling with bone size, whereas we instrumented DXA areal BMD, area-normalised. Second, the "
         "designs differ: Hartley's estimate is one-sample two-stage least squares (overlap 100%) whereas "
         "ours is bounded at \u22642.6% and an MRlap correction moved every estimate a median 20.5% further "
         "from the null. Third, our pipeline is calibrated: the BMD\u2192fracture positive control returns "
         "OR 1.40\u20131.42 against GEFOS 1.55. We read the directions as concordant and the magnitudes as "
         "unresolved, and record the disagreement as a limitation.")
replace_once(old_c, new_c, "OP_C")

# ---- OP_D: delete duplicate boundary sentence at end of L99 ----
old_d = ("The tests executed here (F6, F7) separate the set-point from site-specific and "
         "joint-locus rivals but not from a generic systemic-pleiotropy rival, requiring the "
         "external-data tests F2\u2013F5.")
new_d = ""
replace_once(old_d, new_d, "OP_D")

open(F, "w", encoding="utf-8").write(t)

# ---- re-measure ----
def words(s): return len(re.findall(r"[A-Za-z][A-Za-z0-9'\-]*", s))
ai = t.find("**Objective.**"); ci = t.find("**Conclusions.**", ai)
tail = t[ci:]; em = re.search(r"\n---", tail)
abstract_end = ci + (em.start() if em else len(tail))
leg = t.find("## Figure legends"); body = t[abstract_end:leg]
bl = [l for l in body.splitlines() if not l.strip().startswith("|")]
print("MAIN BODY words after B3:", words("\n".join(bl)))
print("display items (Figure/Table in main):",
      sorted(set(re.findall(r"\bFigure \d+\b", t[:leg]))),
      "tables:", sorted(set(re.findall(r"\n## Table \d+\.", t[:leg]))))
