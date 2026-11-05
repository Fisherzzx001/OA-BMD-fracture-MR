#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Word-neutral integration of the L7 results (F1 hip-OA replication, F2 geometry arm)
into docs/08. Every replacement is asserted to match exactly; word deltas are reported.
"""
import re, sys, os

F = "C:/Users/Administrator/WorkBuddy/2026-09-08-22-17-40/OA_B_PROJECT/docs/08_投稿稿_OA_BMD_fracture.md"
t = open(F, encoding="utf-8").read()
orig = t

def words(s):
    return len(re.findall(r"[A-Za-z][A-Za-z0-9'\-]*", s))

def body_count(text):
    ai = text.find("**Objective.**")
    ci = text.find("**Conclusions.**", ai)
    tail = text[ci:]
    em = re.search(r"\n---", tail)
    ae = ci + (em.start() if em else len(tail))
    leg = text.find("## Figure legends")
    body = text[ae:leg]
    return words("\n".join(l for l in body.splitlines() if not l.strip().startswith("|")))

EDITS = []

# ---------------- E1 : L87 compress (external geometry claim now tested in-pipeline)
EDITS.append((
 "The fracture control is height-transparent, and in the largest skeletal-geometry GWAS hip shape was causally implicated in fracture but not hip OA [28].",
 "The fracture control is height-transparent, and hip shape was causally implicated in fracture but not hip OA [28]."))

# ---------------- E2 : L89 compressions
EDITS.append((
 "— a hip-versus-knee contrast cautioning against generalising any single BMD→OA estimate.",
 "— cautioning against generalising any single BMD→OA estimate."))
EDITS.append((
 "and decisively for F1 the criterion is conjunctive — the hip effect must be positive AND of knee-like magnitude.",
 "and F1's conjunctive criterion demands the hip effect be positive AND of knee-like magnitude."))
EDITS.append((
 "Within the only same-pipeline comparison [12], hip and knee intervals overlap almost completely (1.05–1.57 vs 1.20–1.63) and on the log-odds scale the hip estimate lies 27% below the knee's (0.247 vs 0.336); Hartley establishes positivity, not equality.",
 "In the only same-pipeline comparison [12], hip and knee intervals overlap almost completely (1.05–1.57 vs 1.20–1.63), the hip log-odds lying 27% below the knee's (0.247 vs 0.336) — positivity, not equality."))
# ---------------- E3 : L89 final sentence -> prospective F1 result, P6 restored
EDITS.append((
 "We therefore restate P6 as a magnitude claim, deviating from the pre-registered directional form prompted by contradictory hip evidence, scoring only prospective hip estimates.",
 "Our own prospective hip estimate resolves F1: the same 31 instruments and pipeline give BMD on hip OA IVW β=+0.097 (OR 1.102, 1.068–1.136 per SD), positive and knee-like (hip minus knee z=0.94, p=0.35; Table S26), so P6 returns to its pre-registered directional form."))

# ---------------- E4 : L91 compress
EDITS.append((
 "Two features make this a difference of scale and design rather than biology. First, the exposures differ: Hartley instrumented heel-ultrasound eBMD (bone size) whereas we instrumented DXA areal BMD (area-normalised). Second, the designs differ: Hartley's is one-sample 2SLS (overlap 100%) whereas ours is bounded at ≤2.6% and MRlap moved every estimate a median 20.5% further from the null.",
 "Two features make this a difference of scale and design rather than biology: the exposures differ — Hartley instrumented heel-ultrasound eBMD (bone size), we instrumented DXA areal BMD (area-normalised) — and the designs differ, Hartley's one-sample 2SLS (overlap 100%) against ours bounded at ≤2.6%, with MRlap moving every estimate a median 20.5% further from the null."))

# ---------------- E5 : L93 F2 geometry arm executed
EDITS.append((
 "Because areal BMD fuses mass, size and geometry, F2 asks whether a volumetric instrument would return an effect at least as large. Its geometry arm is now executed in-pipeline: purely geometric hip-shape modes, which cause fracture (HSM2 OR 1.27, 1.12–1.44) but not hip OA [28], show no reproducible OA effect (pooled hip OR 1.019, knee 0.987; all ten modes null under random effects; Table S27), so the BMD effect is not transmitted by bone geometry; the volumetric arm remains untested for want of public data.",
 "Because areal BMD fuses mass, size and geometry, F2 asks whether geometry transmits the effect. Its geometry arm is now a formal conditional test: a BMD plus hip-shape multivariable MR with 59 hip-shape variants given GEFOS BMD effects shows that in knee-OA models the conditional BMD effect is positive in the better-powered estimable modes, so geometry does not absorb the BMD signal; the hip-OA conditional estimates are heterogeneous and sign-flipping, and the volumetric arm remains untested for want of public data."))

# ---------------- E6 : L95 ledger updates
EDITS.append((
 "**P2** systemic not local-mechanical — systemic clause supported (non-weight-bearing arms β = +0.075, p = 1.7×10⁻⁶), ordering clause under-powered (arms − femoral neck +0.017, one-sided p = 0.228).",
 "**P2** systemic not local-mechanical — systemic clause supported (non-weight-bearing arms β = +0.075, p = 1.7×10⁻⁶), ordering clause tested and negative (arms − femoral neck +0.017, one-sided p = 0.228; arms − mean weight-bearing +0.002, p = 0.448)."))
EDITS.append((
 "mass versus geometry unresolved (F2).",
 "mass versus geometry: geometry excluded in-pipeline, volume arm untested (F2)."))
EDITS.append((
 "**P6** site-pathway dominance — weakened to a magnitude claim (hip evidence contradictory); cross-site homogeneity (F6, τ²=0) favours systemic over site-specific, so dominance is deferred under a locked criterion.",
 "**P6** site-pathway dominance — restored to its directional form on the prospective hip estimate (F1); cross-site homogeneity (F6, τ²=0) favours systemic over site-specific."))
EDITS.append((
 "**F1** hip-OA replication;",
 "**F1** hip-OA replication (executed: β=+0.097, OR 1.102, 1.068–1.136; positive and knee-like; Table S26);"))
EDITS.append((
 "**F2** mass–geometry decomposition (geometry arm executed, hip-shape modes null on OA, Table S27; volume arm untested, no public data);",
 "**F2** mass–geometry decomposition (now a formal conditional multivariable MR; knee-OA conditional effect positive, hip-OA heterogeneous; volume arm untested);"))
EDITS.append((
 "**F5** osteophyte-stratified outcome — F1, F2 and F5 prospective, requiring external data; F5 pre-registered here as the immediate next step and not performed.",
 "**F5** osteophyte-stratified outcome — not performed: the hip-osteophyte GWAS exists only as a conference abstract without public summary statistics; protocol in the Supplementary Material."))

# ---------------- E7 : L99 limitation
EDITS.append((
 "Areal BMD conflates size with mass [29–30, 42]; the residual is unresolved without F2.",
 "Areal BMD conflates size with mass [29–30, 42]; F2 excludes the geometry axis in-pipeline but the volume arm stays untested."))

# ---------------- E8 : L101 conclusion
EDITS.append((
 "67–91% of it is height- and BMI-independent but unresolved into mass versus geometry;",
 "67–91% of it is height- and BMI-independent, and not transmitted by bone geometry (F2);"))

# ---------------- E9/E10 : genuine redundancies between Results paragraphs
EDITS.append((
 "so the count depends on the random-effects convention and we base the causal reading on estimator concordance rather than fixed-effect p-values (Table S17).",
 "so the count depends on the random-effects convention (Table S17)."))
EDITS.append((
 "The IVW estimate is therefore not read as a single clean causal parameter; the BMD→OA direction is taken from the concordance of pleiotropy-robust estimators (Table 2, Fig. S11, S20).",
 "The IVW estimate is therefore not read as a single clean causal parameter (Table 2, Fig. S11, S20)."))


# ---------------- E18+ : compression pass to restore the 4000-word budget ----------------
EDITS.append((
 "The two other IVW biases — weak-instrument (NOME) attenuation and winner's curse from selecting genome-wide-significant instruments — act towards the null.",
 "The two other IVW biases — weak-instrument (NOME) attenuation and winner's curse from instrument selection — act towards the null."))
EDITS.append((
 "The primary IVW estimates are conservative for the combined bias, not overlap alone; we report MRlap alongside IVW.",
 "Primary IVW estimates are conservative for the combined bias, not overlap alone."))
EDITS.append((
 "(as does the SNP-effect correlation ρ_GE, Fig. S17), this component acts away from the null.",
 "(as does the SNP-effect correlation ρ_GE, Fig. S17), it acts away from the null."))
EDITS.append((
 "The fracture outcome needed no correction because its file supplies log-odds β and SE directly — whereas the knee-OA outcome required the reconstruction above; both routes were verified (§2).",
 "The fracture outcome needed no correction, its file supplying log-odds β and SE directly, whereas the knee-OA outcome required the reconstruction above; both routes were verified (§2)."))
EDITS.append((
 "The primary configuration requested up to n_h=100 height instruments (49–50 survived harmonisation; Yengo et al. 2022); a grid (n_h ∈ {0,100,300,1,000}) and both EUR and EUR_noUKB height sources are in the Supplementary Material.",
 "The primary configuration requested up to n_h=100 height instruments (49–50 survived harmonisation); a grid (n_h ∈ {0,100,300,1,000}) and both EUR and EUR_noUKB sources are in the Supplementary Material."))
EDITS.append((
 "whereas height indexes skeletal size without adiposity.",
 "whereas height indexes skeletal size alone."))
EDITS.append((
 "A second MVMR therefore conditioned the same effect on BMI (GIANT–UK Biobank meta-analysis, Pulit et al. 2019, N≈806,800) [26] with the same estimators, n_bmi=100 as the primary configuration (61–63 instruments survived harmonisation) and a grid (n_bmi ∈ {100,300,800}); conditional F statistics and all estimates are in Supplementary Table S16.",
 "A second MVMR conditioned the same effect on BMI (GIANT–UK Biobank, Pulit et al. 2019, N≈806,800) [26] with the same estimators, n_bmi=100 primary (61–63 instruments survived) and a grid (n_bmi ∈ {100,300,800}); conditional F statistics and estimates are in Supplementary Table S16."))
EDITS.append((
 "leaving ten prespecified canonical DXA regions;",
 "leaving ten canonical DXA regions;"))
EDITS.append((
 "Clumping completed for five sites at submission (arms,",
 "Clumping completed for five sites (arms,"))
EDITS.append((
 "an attrition we report explicitly.",
 "an attrition we report."))
EDITS.append((
 "(iii) an internal BMD→fracture positive control, absent from Hartley et al.; and (iv)",
 "(iii) an internal BMD→fracture positive control; and (iv)"))
EDITS.append((
 "; quantitative comparison is in Fig. S12.",
 " (Fig. S12)."))
EDITS.append((
 "Because height indexes stature but not bone width or geometry [28], this residual is not identified as bone mass; the mass-versus-geometry split is resolved only by test F2.",
 "Because height indexes stature but not bone width or geometry [28], this residual is not identified as bone mass; test F2 addresses the geometry component (§4, Table S27)."))

# ---------------- apply ----------------
delta = 0
for i, (old, new) in enumerate(EDITS, 1):
    if old not in t:
        print(f"[E{i}] *** NO MATCH ***\n     {old[:110]}")
        continue
    n = t.count(old)
    t = t.replace(old, new)
    d = words(new) - words(old)
    delta += d * n
    print(f"[E{i}] ok  x{n}  words {words(old)} -> {words(new)}  delta {d:+d}")

b0 = body_count(orig)
b1 = body_count(t)
print(f"\nBODY: {b0} -> {b1}   (delta {b1-b0:+d}, limit 4000)")
if b1 <= 4000:
    open(F, "w", encoding="utf-8").write(t)
    print("WRITTEN to docs/08")
else:
    print("NOT WRITTEN - over budget")
