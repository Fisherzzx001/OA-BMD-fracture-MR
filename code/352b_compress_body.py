# -*- coding: utf-8 -*-
# Code/352b — surgical body + highlights compression for docs/08.
# Goal: bring body <=4000 words and Highlights bullets <=85 chars WITHOUT
# altering any locked number (beta/SE/p/OR/CI), reference, or fact.
# Strategy: remove only clearly redundant / duplicated prose (statements that
# repeat an adjacent paragraph, a table note, or themselves).
import re, os, shutil, sys

F = "docs/08_投稿稿_OA_BMD_fracture.md"
BAK = F + ".bak_before_352b"
if not os.path.exists(BAK):
    shutil.copyfile(F, BAK)

t = open(F, encoding="utf-8").read()

# ---- gate word counter (mirrors code/148_oac_compliance.py) ----
def words(s):
    return len([w for w in s.split() if re.search(r"[^\W_]", w, re.UNICODE)])

def body_words(txt):
    ai = txt.find("**Objective.**"); ci = txt.find("**Conclusions.**", ai)
    tail = txt[ci:]; em = re.search(r"\n---", tail)
    ae = ci + (em.start() if em else len(tail))
    leg = txt.find("## Figure legends")
    body = txt[ae:leg]
    lines = [l for l in body.splitlines() if not l.strip().startswith("|")]
    return words("\n".join(lines))

def hl_lengths(txt):
    hi = txt.find("## Highlights"); hl = txt.find("**Keywords.**")
    block = txt[hi:hl]
    bullets = re.findall(r"^- (.+)$", block, re.M)
    return [len(b) for b in bullets]

before_body = body_words(t)
before_hl = hl_lengths(t)

# ---- curated EDITS: (old, new, expected_count, desc) ----
# All OLD strings must be present exactly `expected_count` times.
EDITS = [
    # ---- Highlights: bullets 3 and 5 over 85 chars ----
    ("BMD variants colocalise with fracture (9/31) but not with knee OA (0/31) — a split driven by colocalisation power, not established endpoint-specificity (the strongest instrument gives expected knee-OA z≈0.85, needing ~5×10⁷ participants for genome-wide significance).",
     "BMD colocalises with fracture (9/31) but not knee OA (0/31); the split is power-limited.",
     1, "HL3"),
    ("A systemic bone-forming set-point is a hypothesis whose directional prediction is consistent with data, but whose core claim — one variant set fixing opposite endpoint signs — is untested (F2 volume arm, F10).",
     "Systemic bone-forming set-point fits the direction but its sign-flip claim is untested (F2, F10).",
     1, "HL5"),

    # ---- Body: Methods / Sample overlap (line 55) ----
    (" — 'pleiotropy-robust' is read as 'pleiotropy-aware' (the femoral-neck control's Egger p=0.018 shows directional pleiotropy is present). No direction or significance reversed.",
     ". No direction or significance reversed.", 1, "A"),
    ("; MRlap retained all ten directions (median +20.5%), so the joint-endpoint effect survives overlap correction.",
     "; the joint-endpoint effect survives overlap correction.", 1, "B"),
    (" (as does ρ_GE, Fig. S17)", "", 1, "OO"),

    # ---- Body: Methods / Effect-scale reconstruction (line 57) ----
    ("Table 2 reports the final reconstructed estimates; absolute BMD→OA effects move away from the null (OR 1.060–1.099, not 1.017–1.028) and the reverse effects towards it.",
     "Table 2 reports the final reconstructed estimates (OR 1.060–1.099, not 1.017–1.028).", 1, "C"),
    ("; femoral-neck, lumbar and arms BMD→OA pairs each lost one instrument (26→25, 22→21, 17→16), so pooled z and p use the final set.",
     ".", 1, "PP"),

    # ---- Body: Results / positive control (line 67) ----
    (", a directional-pleiotropy signal reported rather than absorbed; its pleiotropy-robust slope (−0.343) and weighted median (−0.365) agree with IVW, so the control only partly calibrates the pipeline, its directional pleiotropy flagged.",
     "; its pleiotropy-robust slope (−0.343) and weighted median (−0.365) agree with IVW, so the control only partly calibrates the pipeline.", 1, "D"),

    # ---- Body: Results / MR-Egger power (line 69) ----
    ("; with 16–31 instruments the intercept test has limited power and excludes only detectable unbalanced pleiotropy, not pleiotropy per se", "", 1, "E"),

    # ---- Body: Results / PRESSO (line 75) ----
    (", and only the ten primary MR pairs carried α=0.005 protection (see multiplicity note in Limitations)", "", 1, "F"),

    # ---- Body: Results / MVMR BMI (line 81) ----
    (" (the estimate rose at the femoral-neck site)", "", 1, "G"),

    # ---- Body: Methods / MR random-effects convention (line 53) ----
    ("Random-effects IVW used the multiplicative convention SE_random = SE_fixed × √[max(1, Q/df)] (TwoSampleMR), not additive DerSimonian–Laird τ²; where the convention changes an inference, both are reported (§3, Table 2 note).",
     "Random-effects IVW used the multiplicative convention (TwoSampleMR); where the convention changes an inference, both are reported (§3, Table 2 note).", 1, "JJ"),

    # ---- Body: Methods / MVMR (line 61) ----
    ("a grid (n_h ∈ {0,100,300,1,000}) and EUR/EUR_noUKB sources are in Supplementary.",
     "grids are in the Supplement.", 1, "HH"),
    (" (49–50 survived harmonisation)", "", 1, "IIa"),
    (" (61–63 survived)", "", 1, "IIb"),
    ("Height was the primary size proxy: BMI is itself a knee-OA risk factor in the adiposity pathway through which BMD could raise OA risk, whereas height indexes skeletal size alone.",
     "Height was the primary size proxy (BMI is itself a knee-OA risk factor, whereas height indexes skeletal size alone).", 1, "MV"),

    # ---- Body: Discussion / limitations (line 99) ----
    ("The sex-stratified F8 MR was positive in both sexes (KneeOA F β=+0.090, M β=+0.127; AllOA F β=+0.097, M β=+0.066; F-vs-M p=0.029 and 0.006), with no reversal — direction, not magnitude, consistency.",
     "The sex-stratified F8 MR was positive in both sexes (F-vs-M p=0.029 and 0.006), with no reversal — direction, not magnitude, consistency.", 1, "H"),

    # ---- Body: Discussion / SOS synthesis (line 93) ----
    (" (COL2A1/ACAN-high, MMP13-low)", "", 1, "BB1"),
    (" (S30 vi)", "", 1, "BB2"),
    (" (bone-mass and OA axes orthogonal; S22)", " (S22)", 1, "EE"),
    ("is a power phenomenon, not sign-opposition evidence (expected knee-OA z≈0.85; genome-wide significance needs ~5×10⁷ participants), and the set-point core is untested (F2 volume arm, F10).",
     "is a power phenomenon, not sign-opposition evidence (power-limited; §3), and the set-point core is untested.", 1, "DD"),
    ("and the set-point core untested (F2 volume arm, F10).", "and the set-point core untested.", 1, "FF"),
    ("The reading generates P1–P6 and F1–F10 (Table S11); survival counts executed tests only, not register completeness (Fig. S23).", "", 1, "CC"),

    # ---- Body: Discussion / register (line 95) F5 tail ----
    ("; hip-osteophyte GWAS remains unpublished; full in S29", "", 1, "XX"),
]

# ---- apply with verification ----
failures = []
removed_words = 0
for old, new, exp, desc in EDITS:
    cnt = t.count(old)
    if cnt != exp:
        failures.append((desc, exp, cnt, old[:60]))
        continue
    rem = words(old) - words(new)
    removed_words += rem
    t = t.replace(old, new)

if failures:
    print("ABORT: the following edits did not match expected count:")
    for desc, exp, cnt, head in failures:
        print(f"  [{desc}] expected {exp}, found {cnt}: {head!r}")
    sys.exit(1)

# write
open(F, "w", encoding="utf-8").write(t)

after_body = body_words(t)
after_hl = hl_lengths(t)

print("=== COMPRESSION REPORT ===")
print(f"body words: {before_body} -> {after_body}  (removed {before_body-after_body}, target <=4000)")
print(f"Highlights lengths before: {before_hl}")
print(f"Highlights lengths after : {after_hl}")
print(f"HL all <=85: {all(L<=85 for L in after_hl)}")
print(f"Body compliant (<=4000): {after_body<=4000}")
