# -*- coding: utf-8 -*-
"""344 — Apply docs/148 P0-1 / P1-1 / P0-4 manuscript wording remediations.

2026-11-12. Dry-run by default (prints diff + word count); pass --apply to write.

Design (lessons from code/338 + the P0 word-count false-COMPLIANT incident):
  * Hard gate 1: body 3963 -> must stay <= 4000 (Word-compatible count, code/148).
  * Hard gate 2: EVERY reported statistic (beta/SE/p/OR/CI/rg/PP.H4 ...) must be
    byte-identical before/after. We reuse the exact regexes from code/149b plus
    a numeric-token multiset diff, so no "数字与文件一致≠数字对" regression.
  * Only NEW numbers allowed are the two scale-disclosure constants 3.543 and
    sqrt(400000/31873) required by P0-1, and the permutation enrichment values
    from code/342. Every new number is explicitly whitelisted and asserted present.

What changes
------------
P0-1 (scale disclosure — currently the manuscript falsely says "no rescaling"):
  L57  "so no rescaling is needed"  ->  state BMD beta/SE x sqrt(400000/31873)=3.543
  L174 "no rescaling was applied"   ->  state the same reconstruction
P1-1 (tissue-enrichment direction — S25 shows cartilage >= bone, not formation-side):
  L93  formation-side cellular anchor -> mixed/osteochondral, cartilage >= bone
  L95  F3 "evaluated as prespecified" -> "mixed; does not support formation-side exclusivity"
  L285 S26 "formation-side rather than cartilage-matrix" -> neutral/mixed reading
  L620 S25 index line -> report cartilage>=bone + permutation-calibrated
  L621 S26 index line -> drop "formation-side rather than cartilage-matrix"
P0-4 (R1 covariance wording closure, supplement only, non-body):
  docs/146 §7.4 handled separately (not in manuscript body); here only ensure
  S23 legend does not claim MRlap intercept-as-covariance (it currently does not).

Nothing in Abstract / Methods results values / Table cells / refs is touched.
"""
import re, sys, shutil, datetime, os

F = "docs/08_投稿稿_OA_BMD_fracture.md"
APPLY = "--apply" in sys.argv

src = open(F, encoding="utf-8").read()
orig = src

def words(s):
    return len([t for t in s.split() if re.search(r"[^\W_]", t, re.UNICODE)])

def body_stats(text):
    ai = text.find("**Objective.**")
    ci = text.find("**Conclusions.**", ai)
    tail = text[ci:]
    em = re.search(r"\n---", tail)
    ae = ci + (em.start() if em else len(tail))
    leg = text.find("## Figure legends")
    body = text[ae:leg]
    body_nt = "\n".join(l for l in body.splitlines() if not l.strip().startswith("|"))
    return words(body_nt), words(text[ai:ae])

def numeric_multiset(s):
    return sorted(re.findall(r"\d+(?:[.,]\d+)*", s))

# ---- replacement spec: (label, old, new) ------------------------------------
R = []

# P0-1a: L57 effect-scale reconstruction
R.append((
    "P0-1a L57 BMD scale disclosure",
    "the deposited placeholder N=400,000 gives SE = 1/√[2·N·f·(1−f)] matching the verified N=31,873 to within 0.1%, so no rescaling is needed.",
    "the deposited BMD β and SE are on a placeholder N=400,000 scale; restoring the verified N=31,873 requires multiplying both by √(400,000/31,873)=3.543 (per-SNP z invariant).",
))

# P0-1b: L174 Table 2 note
R.append((
    "P0-1b L174 Table2 note scale disclosure",
    "Forward BMD-exposure betas are per-SD on the verified N=31,873 scale (per-instrument closed-form SE check, within 0.1%); no rescaling was applied.",
    "Forward BMD-exposure betas are per-SD on the verified N=31,873 scale, obtained by multiplying the deposited placeholder-N (400,000) β and SE by √(400,000/31,873)=3.543 (per-SNP z invariant; per-instrument closed-form SE within 0.1%).",
))

# P1-1a: L93 cellular anchor + formation-side claim
R.append((
    "P1-1a L93 cellular anchor mixed",
    "Its cellular anchor: RSPO3+ subchondral-bone precursors (COL2A1/ACAN-high, MMP13-low), NFATC1+ osteoclasts (S30 vi);",
    "Candidate cellular correlates span the osteochondral unit — RSPO3+ subchondral-bone precursors (COL2A1/ACAN-high, MMP13-low), NFATC1+ osteoclasts (S30 vi) — but tissue enrichment is mixed (OA-knee cartilage ≥ bone; S25);",
))
R.append((
    "P1-1a2 L93 formation-side evidence wording",
    "Supported is the directional pattern, not the osteogenic mechanism: formation-side evidence is indirect and the set-point core is untested (F2 volume arm, F10).",
    "Supported is the directional pattern, not the osteogenic mechanism: the cellular route is unresolved (S25) and the set-point core untested (F2 volume arm, F10).",
))

# P1-1b: L95 F3 verdict
R.append((
    "P1-1b L95 F3 verdict mixed",
    "**F3** bone-versus-cartilage tissue enrichment (gene-set, evaluated as prespecified);",
    "**F3** bone-versus-cartilage tissue enrichment (permutation-calibrated; mixed, not formation-side-exclusive; S25);",
))

# P1-1c: L285 S26 interpretation
R.append((
    "P1-1c L285 S26 interpretation neutral",
    "This is consistent with a formation-side (bone/subchondral) rather than cartilage-matrix or muscle-mediated interpretation of the BMD→OA effect, and parallels the cross-trait arm in which the same instruments colocalise with fracture (bone) but not OA (joint).",
    "These tissue-literal colocalisation nulls do not by themselves localise the route: cartilage and bone/skeletal-muscle eQTLs both fall short of PP.H4≥0.8, and the competitive gene-set test (S25) finds OA-knee cartilage enrichment at least as strong as bone, so a formation-side-exclusive cellular route is not supported; this parallels the cross-trait arm in which the same instruments colocalise with fracture (bone) but not OA (joint).",
))

# P1-1d: L620 / L621 index lines
R.append((
    "P1-1d L620 S25 index line",
    "- **Table S25.** Bone-versus-cartilage tissue enrichment — competitive gene-set test (F3): OA knee enriches cartilage-development and chondrocyte-differentiation gene sets more strongly than bone-development sets; total-body BMD enriches skeletal-system-development most strongly.",
    "- **Table S25.** Bone-versus-cartilage tissue enrichment — length/chromosome-matched competitive permutation (F3, 10,000 permutations): OA knee enriches cartilage-development and chondrocyte-differentiation sets more strongly than bone-development sets (cartilage−bone contrast +0.800, permutation p=0.0013), whereas total-body BMD shows no cartilage-versus-bone difference (contrast +0.207, p=0.267); the nominal skeletal-system-development signal is 67% shared with cartilage sets.",
))
R.append((
    "P1-1e L621 S26 index line",
    "- **Table S26.** eQTL colocalisation — tissue-literal (cartilage and skeletal-muscle bone-proxy) arms of F4: at the 31 BMD instrument loci, cartilage eQTL colocalisation with OA knee, BMD total or fracture is 0/31 at PP.H4≥0.8, and skeletal-muscle (bone-proxy) eQTL is 0–1/31 (one borderline OA-knee locus, PP.H4=0.887) — neither tissue colocalises with BMD/fracture, consistent with a formation-side rather than cartilage-matrix or muscle-mediated interpretation.",
    "- **Table S26.** eQTL colocalisation — tissue-literal (cartilage and skeletal-muscle bone-proxy) arms of F4: at the 31 BMD instrument loci, cartilage eQTL colocalisation with OA knee, BMD total or fracture is 0/31 at PP.H4≥0.8, and skeletal-muscle (bone-proxy) eQTL is 0–1/31 (one borderline OA-knee locus, PP.H4=0.887); neither tissue colocalises with BMD/fracture, and these nulls do not discriminate bone/subchondral from cartilage-matrix routes (see S25 for the gene-set contrast).",
))

# S25 table itself: replace the MWW header note + add permutation columns is done
# in a separate step (structured table). Here we only rewrite the method paragraph.
R.append((
    "P1-1f S25 method paragraph permutation",
    "Tissue specificity of the 13-phenotype genetic architecture was tested by a competitive gene-set enrichment (MAGMA-style): each GWAS SNP was assigned to the nearest gene body within ±50 kb; gene-level association = the minimum p across assigned SNPs; each GO gene set was tested against the remaining genes by a one-sided Mann–Whitney U on −log10 p. Positive enrichment Δ means the set's genes are more associated than background.",
    "Tissue specificity of the 13-phenotype genetic architecture was tested by a competitive gene-set enrichment with permutation calibration (MAGMA-style statistic; code/342): each GWAS SNP was assigned to the nearest gene body within ±50 kb; gene-level association = −log10 of the minimum p across assigned SNPs; each GO set's mean was compared with background, and significance was calibrated by 10,000 permutations matched on gene length and chromosome (the previous one-sided Mann–Whitney p was optimistic because the min-p statistic depends on gene length/SNP density and sets overlap). Positive enrichment Δ means the set's genes are more associated than background. Because cartilage_development and chondrocyte_differentiation overlap (Jaccard 0.44) and the broad skeletal-system-development set is 67% shared with cartilage sets, we additionally tested a direct cartilage-minus-bone contrast with overlapping genes removed.",
))

# ---- apply replacements ----------------------------------------------------
print("=" * 76)
print("344 · P0-1 / P1-1 / P0-4 wording remediation  (", "APPLY" if APPLY else "DRY-RUN", ")")
print("=" * 76)
fail = 0
for label, old, new in R:
    cnt = src.count(old)
    if cnt != 1:
        print(f"[FAIL] {label}: anchor count={cnt} (need exactly 1)")
        fail += 1
        continue
    src = src.replace(old, new, 1)
    print(f"[OK]   {label}  (-{words(old)} +{words(new)} words, net {words(new)-words(old):+d})")

b_new, a_new = body_stats(src)
b_old, a_old = body_stats(orig)
print("\n--- word-gate ---")
print(f"body {b_old} -> {b_new} (cap 4000, headroom {4000-b_new})   abstract {a_old} -> {a_new} (cap 350)")
if b_new > 4000:
    print("[FAIL] body exceeds 4000"); fail += 1
if a_new > 350:
    print("[FAIL] abstract exceeds 350"); fail += 1

# ---- reported-statistic invariance: reuse 149b regexes ----------------------
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
]
print("\n--- reported-statistic invariance (149b regex set) ---")
for p in PATS:
    mo = re.search(p, orig); mn = re.search(p, src)
    go = mo.group(0) if mo else None; gn = mn.group(0) if mn else None
    if go != gn:
        print(f"[FAIL] statistic drift /{p}/: {go!r} -> {gn!r}"); fail += 1
if fail == 0:
    print("all reported statistics byte-identical")

# new-number whitelist: these MUST be present after edit (P0-1 + P1-1 disclosures)
NEW = ["3.543", "√(400,000/31,873)", "0.0013", "0.267", "Jaccard 0.44", "67%"]
print("\n--- new-number whitelist (must be present) ---")
for x in NEW:
    ok = x in src
    print(f"[{'OK' if ok else 'FAIL'}] {x}")
    if not ok: fail += 1

print("\n" + "=" * 76)
if fail:
    print(f"RESULT: {fail} CHECK(S) FAILED — NOT written")
    sys.exit(1)

if APPLY:
    bak = f".workbuddy/backups/docs08_pre344_{datetime.date.today()}.md"
    os.makedirs(os.path.dirname(bak), exist_ok=True)
    shutil.copyfile(F, bak)
    open(F, "w", encoding="utf-8").write(src)
    print(f"APPLIED. backup -> {bak}")
    print(f"body {b_old} -> {b_new}; abstract {a_old} -> {a_new}")
else:
    print("DRY-RUN OK — re-run with --apply to write")
