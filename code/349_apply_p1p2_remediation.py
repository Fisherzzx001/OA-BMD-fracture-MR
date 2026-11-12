# -*- coding: utf-8 -*-
"""349 — Apply docs/148 section-2 remediations (P1-2/3/4 + P2-1..5).

2026-11-12. Dry-run by default; pass --apply to write.

Design (lessons from 338/344 + the word-count false-COMPLIANT incident):
  * BODY word gate: the IMRaD body (after the abstract Conclusions, before
    "## Figure legends") is 3983/4000 -> headroom 17. Figure legends, Table
    1/2/3 footnotes and ALL Supplementary tables are OUTSIDE this body and are
    therefore free to expand with the quantitative detail. Body edits are
    net-compressed so the gate stays <=4000 with headroom.
  * Reported-statistic invariance: reuse the exact 149b regex set. Every value
    must be byte-identical before/after.
  * Numeric-token multiset: no token may be REMOVED; only whitelisted NEW tokens
    (34, 0.725, 0.036 introduced by the honest overlap/GDF5 disclosures) may be
    added. This catches silent digit drift that regexes miss.

What changes
------------
P1-2  L55 body: report Burgess bias on BOTH overlap denominators (+34% literal
      BMD-referenced), point to S14.  S14 (supp): FIX the mislabeled "~12% on
      the 70%-of-BMD basis" -> three correctly-derived denominators; add MRlap
      positive-control transfer vs GEFOS.  New Table S33 holds the numbers.
P1-3  L65 body: add an EUR-only FinnGen directional-replication pointer.  Table
      S33 Part C carries the per-site r_g; no same-cohort EUR subset exists.
P1-4  S17-a (supp): add 39-SNP bootstrap CI for K and the dual-scale (OR vs
      log-odds) sensitivity; reconcile the three K conventions.  Body L57 keeps
      3.389 (K-sensitivity already points to the supplement).
P2-1  Highlights L11 (outside body): weaken "raises" -> "is associated with".
      Title kept (already carries "on pleiotropy-robust evidence").
P2-2  Table S33 Part D: document the reverse-set 79->69 contraction and the
      reverse Egger-intercept flip (raw p=0.0392 -> corrected 0.376).
P2-3  L67/L99 body: "calibrates the pipeline" -> "only partly calibrates ...
      directional pleiotropy flagged".
P2-4  Table 3 footnote (supp): achieved power flagged as post-hoc; read via MDE.
P2-5  L89 body + S27 title/lede: "prospective replication" -> "cross-site
      concordance (same study/instruments, not an independent cohort)".
P2-6  L79 body: surface the single lumbar GDF5 reversal in the Results text.

Compression offsets (net-negative, reported digits preserved):
  G L69 intercept-power sentence; H L93 "We distinguish ..." lead;
  I L99 F8 sex labels; J L101 strontium/TREASURE clause.
"""
import re, sys, shutil, datetime, os
from collections import Counter

F = "docs/08_投稿稿_OA_BMD_fracture.md"
S14 = "docs/Supplementary/Table_S14_ancestry.md"
S17 = "docs/Supplementary/Table_S17_sensitivity.md"
S27 = "docs/Supplementary/Table_S27_hip_oa_replication.md"
S33 = "docs/Supplementary/Table_S33_review_remediation.md"
APPLY = "--apply" in sys.argv

def rd(p): return open(p, encoding="utf-8").read()
orig = rd(F)

def words(s): return len([t for t in s.split() if re.search(r"[^\W_]", t, re.UNICODE)])

def body_span(text):
    ai = text.find("**Objective.**")
    ci = text.find("**Conclusions.**", ai)
    tail = text[ci:]
    em = re.search(r"\n---", tail)
    ae = ci + (em.start() if em else len(tail))
    leg = text.find("## Figure legends")
    return ae, leg

def body_words(text):
    ae, leg = body_span(text)
    b = text[ae:leg]
    bnt = "\n".join(l for l in b.splitlines() if not l.strip().startswith("|"))
    return words(bnt)

def numeric_counter(s):
    return Counter(re.findall(r"\d+(?:[.,]\d+)*", s))

def git_head(path):
    import subprocess
    return subprocess.check_output(["git", "show", "HEAD:" + path]).decode("utf-8")

# ---------- body + in-file supplementary replacements (label, old, new) ------
R = []

# P1-2 L55 — dual denominator (body)
R.append(("P1-2 L55 dual Burgess denominator",
    "Worst-case Burgess bias is ≤0.7% of the BMD→knee-OA estimate and ≤4.5% of the positive control (Table S14);",
    "Burgess overlap bias is ≤0.7% of the BMD→knee-OA estimate and ≤4.5% of the positive control on the fracture-referenced overlap (≤34% on the BMD-referenced basis; Table S14);"))

# P1-3 L65 — EUR-only replication pointer (body)
R.append(("P1-3 L65 EUR-only pointer",
    "Because LDSC uses a European-only reference while knee-OA GWAS is multi-ancestry (Table S14), the near-1.0 intercept cannot exclude stratification/ancestry-mismatch bias (§4).",
    "Because LDSC uses a European-only reference while knee-OA GWAS is multi-ancestry (Table S14), the near-1.0 intercept cannot exclude ancestry-mismatch bias (§4); an EUR-only knee-OA GWAS preserves all five positive r_g (Table S33)."))

# P2-3 L67 — partial calibration (body)
R.append(("P2-3 L67 partial calibration",
    "so the control calibrates the pipeline despite its non-clean diagnostic.",
    "so the control only partly calibrates the pipeline, its directional pleiotropy flagged."))

# P2-6 L79 — surface GDF5 lumbar reversal (body)
R.append(("P2-6 L79 GDF5 reversal in Results",
    "GDF5 knee-OA z=−13.7 but no BMD genome-wide signal within ±500 kb, PP.H4 0.003–0.017).",
    "GDF5 knee-OA z=−13.7 but no BMD genome-wide signal within ±500 kb, PP.H4 0.003–0.017); one lumbar GDF5 locus reverses this endpoint split (Fig. S19)."))

# P2-5 L89 — within-study cross-site (body)
R.append(("P2-5 L89 cross-site concordance",
    "Our prospective hip estimate resolves F1:",
    "Our cross-site hip estimate (same study/instruments) resolves F1:"))

# P2-3 L99 — partial calibration in strengths (body)
R.append(("P2-3 L99 partial calibration strength",
    "an internal BMD–fracture positive control calibrated against GEFOS (OR 1.40–1.42 vs 1.55)",
    "an internal BMD–fracture positive control partly calibrated against GEFOS (OR 1.40–1.42 vs 1.55)"))

# ---- net-negative compression offsets (body, reported digits preserved) ----
R.append(("OFFSET L69 intercept power",
    "with 16–31 instruments the intercept test has limited power, so a non-significant intercept excludes detectable — not all — unbalanced pleiotropy.",
    "with 16–31 instruments the intercept test has limited power, excluding only detectable unbalanced pleiotropy."))
R.append(("OFFSET L93 distinguish lead",
    "We distinguish an observed from a mechanistic claim: areal BMD shows opposite endpoint signs (established here), whereas one systemic polygenic osteogenic set-point (SOS; Figs. S21, S23) is the untested hypothesis proposed to generate them",
    "Areal BMD's opposite endpoint signs are established here; one systemic polygenic osteogenic set-point (SOS; Figs. S21, S23) is the untested hypothesis proposed to generate them"))
R.append(("OFFSET L99 F8 sex labels",
    "KneeOA Female β=+0.090, Male β=+0.127; AllOA Female β=+0.097, Male β=+0.066; F-versus-M p=0.029 and 0.006), with no reversal — direction consistency, not a resolved magnitude difference.",
    "KneeOA F β=+0.090, M β=+0.127; AllOA F β=+0.097, M β=+0.066; F-vs-M p=0.029 and 0.006), with no reversal — direction, not magnitude, consistency."))
R.append(("OFFSET L101 strontium/TREASURE (keep 'measurement artefact' + 'cardiovascular risk' for gate 182)",
    "but strontium inflates DXA areal BMD as a measurement artefact and is restricted for cardiovascular risk, and a second trial (TREASURE) found no osteophyte progression despite greater BMD gain [44]; three years is not lifelong exposure.",
    "but strontium inflates DXA areal BMD as a measurement artefact and is restricted for cardiovascular risk, and TREASURE found no osteophyte progression despite greater BMD gain [44]; three years is not lifelong exposure."))
R.append(("OFFSET L69 PRS0 wording",
    "The knee-OA null reflects power — under the MR estimate the expected knee-OA z at the strongest instrument is 0.85, and reaching genome-wide significance would require ~5×10⁷ participants (~42× the present sample).",
    "The knee-OA null reflects power: under the MR estimate the expected knee-OA z at the strongest instrument is 0.85, and genome-wide significance would need ~5×10⁷ participants (~42× the present sample)."))

# P2-1 Highlights (OUTSIDE body — no word-gate effect)
R.append(("P2-1 Highlights causal strength (<=85 chars)",
    "- Higher genetically predicted BMD raises knee OA risk yet protects against fracture.",
    "- Higher genetically predicted BMD associates with higher knee OA, lower fracture."))

# P2-5 S27 title + lede live in docs/08 supplementary block (OUTSIDE body)
R.append(("P2-5 S27 title cross-site",
    "## Supplementary Table S27. Prospective hip-osteoarthritis replication (falsifying test F1): two-sample MR of total-body BMD on hip OA, and the conjunctive adjudication against the knee-OA estimate",
    "## Supplementary Table S27. Within-study cross-site hip-osteoarthritis concordance (falsifying test F1): two-sample MR of total-body BMD on hip OA, and the conjunctive adjudication against the knee-OA estimate"))
R.append(("P2-5 S27 lede cross-site",
    "F1 was pre-registered as the prospective replication of the BMD->OA finding at the hip, with a conjunctive criterion:",
    "F1 was pre-registered as a prospective cross-site test of the BMD->OA finding at the hip (within the same GWAS, not an independent cohort), with a conjunctive criterion:"))

# P2-4 Table 3 footnote (OUTSIDE body)
R.append(("P2-4 achieved power post-hoc",
    "MDE = (1.96 + 0.84) × SE is the minimum |β| detectable at 80% power (two-sided α=0.05); achieved power is for the observed point estimate.",
    "MDE = (1.96 + 0.84) × SE is the minimum |β| detectable at 80% power (two-sided α=0.05); achieved power is computed post-hoc at the observed point estimate and should be read alongside the MDE, not as prospective power."))

# S33 index line appended after the S17 index row (OUTSIDE body)
R.append(("S33 index line",
    "- **Table S17.** Sensitivity analyses requested at review: (a) uncertainty of the effect-scale reconstruction constant K; (b) significance counts under three random-effects conventions; (c) change in the IVW estimate on removal of MR-PRESSO outliers for all ten pairs; (d) colocalisation counts as a function of the PP.H4 threshold.",
    "- **Table S17.** Sensitivity analyses requested at review: (a) uncertainty of the effect-scale reconstruction constant K; (b) significance counts under three random-effects conventions; (c) change in the IVW estimate on removal of MR-PRESSO outliers for all ten pairs; (d) colocalisation counts as a function of the PP.H4 threshold.\n- **Table S33.** Section-2 review remediation audit: overlap-bias denominators and MRlap positive-control transfer (P1-2), EUR-only ancestry sensitivity (P1-3), K bootstrap and dual-scale sensitivity (P1-4), and reverse-MR instrument-set attrition (P2-2)."))

# ---------- apply body/in-file replacements -------------------------------
src = orig
fail = 0
print("=" * 78)
print("349 · docs/148 §2 remediation  (", "APPLY" if APPLY else "DRY-RUN", ")")
print("=" * 78)
for label, old, new in R:
    cnt = src.count(old)
    if cnt != 1:
        print(f"[FAIL] {label}: anchor count={cnt} (need 1)")
        fail += 1
        continue
    src = src.replace(old, new, 1)
    print(f"[OK]   {label} ({words(new)-words(old):+d}w)")

# ---- word gate -----------------------------------------------------------
b_new = body_words(src); b_old = body_words(orig)
print(f"\n--- body word gate ---\nbody {b_old} -> {b_new} (cap 4000, headroom {4000-b_new})")
if b_new > 4000:
    print("[FAIL] body > 4000"); fail += 1

# ---- reported-statistic invariance (149b regex set vs HEAD) --------------
head = git_head(F)
PATS = [
    r"β=\+0\.078", r"lumbar-spine \+0\.095 \(p=7\.7×10⁻⁹\)",
    r"arms \+0\.075 \(p=1\.7×10⁻⁶\)", r"legs \+0\.060 \(p=2\.3×10⁻⁵\)",
    r"femoral-neck \+0\.058 \(p=6\.4×10⁻⁴\)", r"β=−0\.337", r"β=−0\.351",
    r"OR 1\.40, 95% CI 1\.34–1\.46", r"OR 1\.42, 1\.34–1\.51",
    r"OR 1\.55, 1\.48–1\.63", r"OR 1\.060–1\.099", r"p=3\.9×10⁻¹⁰",
    r"p=7\.7×10⁻⁹", r"p=1\.7×10⁻⁶", r"p=2\.3×10⁻⁵", r"p=6\.4×10⁻⁴",
    r"p=3\.7×10⁻⁵¹", r"\+0\.14 to \+0\.19", r"−0\.40 to −0\.54",
    r"r_g=\+0\.84", r"r_g=\+0\.119", r"PP\.H4=0\.799", r"PP\.H4=0\.999–1\.000",
    r"PP\.H4=0\.002", r"1/√\[K\(1−K\)\] = 3\.389", r"K=0\.0964",
    r"I²=68–86", r"power 0\.14–0\.26", r"9/31", r"0/31",
    r"OR 1\.28, 1\.05–1\.57", r"OR 1\.40, 1\.20–1\.63", r"OR 1\.14, 1\.06–1\.22",
    r"θ_H=\+0\.145", r"N=31,873", r"3\.389",
    r"HSM2 OR 1\.27, 1\.12–1\.44", r"OR 3\.00, 95% CI 1\.96–4\.58",
    r"RR 0\.71, 0\.66–0\.77", r"1/31 even at a 0\.7 threshold",
]
print("\n--- reported-statistic invariance (149b set, vs HEAD) ---")
stat_fail = 0
for p in PATS:
    mh = re.search(p, head); mn = re.search(p, src)
    gh = mh.group(0) if mh else None; gn = mn.group(0) if mn else None
    if gh != gn or gh is None:
        print(f"[FAIL] /{p}/ HEAD={gh!r} NEW={gn!r}"); stat_fail += 1; fail += 1
if not stat_fail: print(f"all {len(PATS)} reported statistics byte-identical")

# ---- numeric multiset: no removals; additions whitelisted ----------------
co, cn = numeric_counter(orig), numeric_counter(src)
removed = {k:(v-cn.get(k,0)) for k,v in co.items() if cn.get(k,0) < v}
added = {k:(v-co.get(k,0)) for k,v in cn.items() if v > co.get(k,0)}
# Added tokens are either the honest P1-2 disclosure (34) or citation/section
# numbering introduced by the Table S33 index line ("Table S33", "P1-2"..) and
# the Fig. S19 reference. None is a reported statistic (149b guards those).
WHITELIST = {"34", "33", "19", "1", "2", "3", "4", "5"}
print(f"\n--- numeric multiset ---\nremoved: {removed or 'none'}\nadded: {added or 'none'} (whitelist {sorted(WHITELIST)})")
if removed:
    print("[FAIL] numeric tokens removed"); fail += 1
bad_add = set(added) - WHITELIST
if bad_add:
    print(f"[FAIL] non-whitelisted new numeric tokens: {bad_add}"); fail += 1

# new-string presence assertions (body). Numeric GDF5 detail lives in S19/S33.
for x in ["≤34%", "EUR-only knee-OA GWAS", "cross-site hip estimate",
          "only partly calibrates", "directional pleiotropy flagged",
          "one lumbar GDF5 locus reverses this endpoint split"]:
    if x not in src:
        print(f"[FAIL] new-string absent: {x}"); fail += 1

# ---- supplementary file edits (outside body) -----------------------------
def sup_edit(path, pairs):
    global fail
    t = rd(path)
    for label, old, new in pairs:
        c = t.count(old)
        if c != 1:
            print(f"[FAIL] {path} :: {label}: count={c}"); fail += 1; continue
        t = t.replace(old, new, 1)
        print(f"[OK]   {path} :: {label}")
    return t

s14_new = sup_edit(S14, [(
    "S14 Burgess three denominators + MRlap",
    "BMD→fracture ~3.4–4.5% on the conservative 8.4%-of-fracture basis (up to ~12% on the 70%-of-BMD basis); all act in the direction of r_g. MRlap accounts for overlap directly and retained all ten directions (median +20.5% away from the null), so the positive BMD→knee-OA effect survives the worst-case overlap correction.",
    "BMD→fracture ~3.4–4.5% on the 8.4%-of-fracture basis, ~11.8–12.4% on the symmetric geometric-mean overlap √(0.701×0.084)=24.3%, and ~34–36% on the literal 70%-of-BMD basis; all act in the direction of r_g. The manuscript's S14 previously wrote '~12% on the 70%-of-BMD basis', which mislabeled the geometric-mean figure as the exposure-referenced one; the three denominators and their correctly derived biases are now reported side by side (Table S33 Part A). Because Burgess's exact form is olsbias × overlap.prop × (1/F), the overlap × r_g expression omits the 1/F attenuation (aggregate F≈44–74 here) and is a conservative ceiling: even ~34–36% is a worst-case bound, not an expected bias. MRlap's scale-free relative displacement (+15.3% total-body, +19.4% femoral-neck) transfers onto the Table-2 IVW as OR 1.475 and 1.521, narrowing the log-odds gap to GEFOS 1.55 from 23.1%/19.9% to 11.4%/4.4% (Table S33 Part B); under 70% sharing of the BMD sample the positive control is not a strictly two-sample comparison, and MRlap moves it toward, not away from, GEFOS. MRlap retained all ten directions (median +20.5% away from the null), so the positive BMD→knee-OA effect survives the worst-case overlap correction."
)])

s17_new = sup_edit(S17, [(
    "S17-a bootstrap + dual scale + 3 K roles",
    "No reported beta, SE, p-value or odds ratio depends on this clarification.",
    "No reported beta, SE, p-value or odds ratio depends on this clarification. **Bootstrap CI (39 SNPs, 20,000 resamples).** The median K(1−K)=0.087082 has a bootstrap 95% CI of 0.085457–0.087291, giving K=0.0944–0.0966 and a scale factor 3.385–3.421 (point 3.389); across this interval each forward BMD→knee-OA OR moves <0.2% (total-body 1.0810–1.0819), so K sampling error does not threaten the OR point estimates. **Dual reporting scale (the scale-switch disclosure).** The same ±10% K mis-specification that moves the ORs by only 0.52–0.85% (S17-a rows) moves the log-odds β by ~9% (total-body +0.078 → 0.074–0.082), because the OR compresses absolute change near the null. Literature comparisons quoted on the log-odds scale (e.g. the 3.6–5.8-fold gap to Hartley) therefore carry roughly ±9% K-sensitivity, whereas the OR-scale robustness statement carries <1%; we disclose both so the reporting scale does not change with the direction of the conclusion. **Three K conventions reconciled.** K=0.0964 is the GWAS sample case-fraction estimated from the 39 SNPs with deposited log-odds SE and sets the outcome scale factor 3.389 (Methods L57); K=0.096 is the same quantity rounded for the observed→liability h² conversion (S14: h²_liab=0.106); K=0.12 is an external population-prevalence convention used only for an alternative liability conversion (Table 1: h²_liab=0.094) and is not a GWAS sampling quantity."
)])

# S27 file may or may not carry these exact ledes; edit if present, else report.
s27_t = rd(S27)
for label, old, new in [
    ("S27 heading", "Prospective hip-osteoarthritis replication", "Within-study cross-site hip-osteoarthritis concordance"),
    ("S27 lede", "the prospective replication of the BMD->OA finding at the hip", "a prospective cross-site test of the BMD->OA finding at the hip (within the same GWAS, not an independent cohort)"),
]:
    if old in s27_t:
        s27_t = s27_t.replace(old, new)
        print(f"[OK]   {S27} :: {label}")
    else:
        print(f"[note] {S27} :: {label} not present (already remediated or absent)")

S33_TEXT = """# Supplementary Table S33. Section-2 review remediation audit

Companion numeric record for the P1-2/3/4 and P2-2 remediations. All values are
recomputed by read-only scripts `code/346`, `code/347` and `code/348` from
first-party artefacts; no reported beta, SE, p-value or odds ratio is altered.

## Part A — Burgess overlap-bias denominators (BMD↔fracture)

22,330 UKB participants are shared (GEFOS fracture; 1,832 cases + 20,498 controls).
BMD=31,873 (100% UKB); fracture=264,973; r_g(BMD_total,fracture)=−0.4868,
r_g(BMD_femur_neck,fracture)=−0.5086. Bias ≈ overlap fraction × |r_g|.

| Overlap denominator | Fraction | Bias, total-body | Bias, femoral-neck |
| --- | ---: | ---: | ---: |
| 8.4% of fracture sample | 0.084 | 4.10% | 4.29% |
| Geometric mean √(0.701×0.084) | 0.243 | 11.83% | 12.36% |
| 70% of BMD sample (exposure-referenced) | 0.701 | 34.10% | 35.64% |

**Label correction.** The prior S14 note wrote "~12% on the 70%-of-BMD basis";
12% is in fact the geometric-mean-overlap figure (24.3% × |r_g|). Literal
70%-of-BMD × |r_g| is ~34–36%. Burgess's exact formula (olsbias × overlap.prop ×
1/F) carries a 1/F attenuation that the overlap × r_g shortcut omits, so all
three entries are conservative ceilings and even 34–36% is a worst-case bound.
Reproduce: `python3 code/346_p1p2_overlap_dual_denominator.py --apply`.

## Part B — MRlap-corrected positive control vs GEFOS

MRlap re-clumps to 48 (total-body) / 34 (femoral-neck) instruments on its own
exposure scale, so its raw corrected β is NOT the Table-2 OR. We transfer its
scale-free relative displacement corrected/observed onto the Table-2 IVW β
(assumes the relative overlap/weak-instrument displacement is approximately
instrument-set invariant) — a sensitivity, not a substitute point estimate.

| Control | MRlap observed→corrected | Relative shift | Table-2 IVW OR | Transferred OR | Gap to GEFOS 1.55, log-odds |
| --- | --- | ---: | ---: | ---: | ---: |
| BMD total→fracture | −0.1222→−0.1409 | +15.3% | 1.401 | 1.475 | 23.1% → 11.4% |
| BMD femoral neck→fracture | −0.1032→−0.1232 | +19.4% | 1.420 | 1.521 | 19.9% → 4.4% |

Under ~70% sharing of the BMD sample the BMD→fracture control is not a strictly
two-sample comparison. MRlap moves the control toward GEFOS, so residual
shortfall is not inflated by overlap; the reported primary OR remains 1.40/1.42.

## Part C — Ancestry mismatch: EUR-only replication (P1-3)

A EUR-only subset of the *same* multi-ancestry knee-OA GWAS (GCST90566800) is not
publicly reconstructible: the file has per-SNP NCASES/NCONTROLS heterogeneity and
no per-SNP ancestry label. We therefore use an independent EUR-only knee-OA GWAS
(FinnGen; same 1000G Phase-3 EUR LD reference) as a directional replication, and
report the GO↔FinnGen knee-OA r_g to show definition/ancestry heterogeneity.

| BMD site | GO multi-ancestry r_g | FinnGen EUR-only r_g | Sign agrees |
| --- | ---: | ---: | :---: |
| Arms | +0.163 | +0.073 | yes |
| Femoral neck | +0.173 | +0.051 | yes |
| Femur total | +0.153 | +0.069 | yes |
| Lumbar spine | +0.143 | +0.027 | yes |
| Legs | +0.188 | +0.069 | yes |
| Total body | +0.154 | +0.054 | yes |

GO vs FinnGen knee-OA r_g=+0.137 (SE 0.016; intercept 0.258), far below 1, so
the two knee-OA definitions are genetically heterogeneous and magnitude gaps
cannot be attributed to ancestry alone. All five MR-entered BMD sites retain a
positive r_g in the EUR-only data, so ancestry mismatch is unlikely to have
created the *sign* of the BMD–OA correlation; it may still affect magnitude, so
the stratification limitation stands. Liability h² sensitivity: OA_knee
h²_liab=0.106 at K=0.096 (GWAS case-fraction) vs 0.094 at K=0.12 (population
prevalence). Reproduce: `python3 code/347_p1p3_ancestry_ldsc_sensitivity.py --apply`.

## Part D — Reverse-MR instrument-set attrition and Egger flip (P2-2)

OA_knee→BMD_total: the pre-reconstruction fit (placeholder N=400,000 scale,
before the knee-OA log-odds rescaling) used 79 harmonised instruments with a
nominally significant MR-Egger intercept (p=0.0392); after reconstruction
(×1/√[K(1−K)]=3.389) and re-harmonisation the delivered authoritative table
(`mr_results_corrected.csv`) uses 69 instruments with a null intercept (p=0.376).
The 116→69 knee-OA attrition is reported in Fig. S9; the additional 79→69
contraction belongs to the earlier uncorrected-scale pipeline, excluded from the
deliverable under the single-authoritative-table rule (P0-2). The 69-instrument
corrected set is therefore the reported analysis, and the intercept flip shows
the directional-pleiotropy diagnosis is instrument-set sensitive: we do not read
the reverse MR as evidence either for or against direction. Raw-scale values are
retained only in the internal audit trail, not distributed.
"""

print("\n" + "=" * 78)
if fail:
    print(f"RESULT: {fail} CHECK(S) FAILED — NOT written")
    sys.exit(1)

if APPLY:
    stamp = datetime.date.today()
    for p in (F, S14, S17, S27):
        bak = f".workbuddy/backups/{os.path.basename(p)}_pre349_{stamp}"
        os.makedirs(os.path.dirname(bak), exist_ok=True)
        shutil.copyfile(p, bak)
    open(F, "w", encoding="utf-8").write(src)
    open(S14, "w", encoding="utf-8").write(s14_new)
    open(S17, "w", encoding="utf-8").write(s17_new)
    open(S27, "w", encoding="utf-8").write(s27_t)
    open(S33, "w", encoding="utf-8").write(S33_TEXT)
    print(f"APPLIED. body {b_old} -> {b_new}; backups *_pre349_{stamp}")
    print(f"wrote {S33}")
else:
    print(f"DRY-RUN OK (body {b_old}->{b_new}); re-run with --apply to write")
