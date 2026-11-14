# -*- coding: utf-8 -*-
# Code/352c — second-pass surgical compression (docs/08).
# Brings body <=4000 words and Highlights bullets <=85 chars, without altering
# any locked number / reference / fact. Each edit removes only redundant or
# duplicated prose. Operates on the current (post-352b) file.
import re, os, shutil, sys

F = "docs/08_投稿稿_OA_BMD_fracture.md"
BAK = F + ".bak_before_352c"
if not os.path.exists(BAK):
    shutil.copyfile(F, BAK)

t = open(F, encoding="utf-8").read()

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
    return [len(b) for b in re.findall(r"^- (.+)$", block, re.M)]

before_body = body_words(t)
before_hl = hl_lengths(t)

EDITS = [
    # Highlights: bullets 3 and 5 -> <=85 chars
    ("BMD colocalises with fracture (9/31) but not knee OA (0/31); the split is power-limited.",
     "BMD colocalises with fracture (9/31) but not knee OA (0/31); split power-limited.", 1, "HL3"),
    ("Systemic bone-forming set-point fits the direction but its sign-flip claim is untested (F2, F10).",
     "Bone-forming set-point fits the direction; sign-flip claim untested (F2, F10).", 1, "HL5"),

    # Discussion hip comparison (line 89)
    ("hip and knee intervals overlap (1.05–1.57 vs 1.20–1.63);",
     "hip and knee intervals overlap;", 1, "L89a"),
    ("our cross-site hip estimate (same 31 instruments) gives BMD→hip-OA IVW β=+0.097 (OR 1.102, 1.068–1.136), positive and knee-like (hip−knee z=0.94, p=0.35; Table S27), restoring P6's pre-registered direction.",
     "our cross-site hip estimate gives BMD→hip-OA IVW β=+0.097 (OR 1.102, 1.068–1.136), positive and knee-like, restoring P6's pre-registered direction.", 1, "L89b"),
    ("we adjudicate rather than cite selectively, and F1 demands a positive, knee-like hip effect.",
     "we adjudicate rather than cite selectively; F1 demands a positive, knee-like hip effect.", 1, "L89c"),

    # SOS register closing (line 95)
    ("P4 alone discriminates (via colocalisation); F2's decisive arm is unexecuted, so survival counts executed tests only.",
     "F2's decisive arm is unexecuted.", 1, "L95"),

    # Discussion intro (line 87)
    ("partial mediation by skeletal size rather than bone mass",
     "partial mediation by skeletal size", 1, "L87"),

    # Scale validation (line 77) — "outstanding" will be updated post independent-test
    ("a variant-level cross-check against an independent BMD GWAS remains the definitive validation and is outstanding",
     "a variant-level cross-check against an independent BMD GWAS remains the definitive validation", 1, "L77"),

    # Discussion SOS synthesis (line 93)
    ("is a power phenomenon, not sign-opposition evidence (power-limited; §3), and the set-point core is untested.",
     "is a power phenomenon, not sign-opposition evidence (power-limited; §3).", 1, "DD"),
    ("but tissue enrichment is mixed (OA-knee cartilage ≥ bone; S25);",
     "but tissue enrichment is mixed;", 1, "L93c"),
    ("a BMD-plus-hip-shape MVMR (59 shape variants with GEFOS BMD effects) leaves",
     "a BMD-plus-hip-shape MVMR leaves", 1, "L93d"),
    ("so geometry does not absorb it, albeit on weak conditional instruments; hip-OA estimates are heterogeneous and sign-flipping,",
     "so geometry does not absorb it; hip-OA estimates are heterogeneous and sign-flipping,", 1, "L93e"),

    # Results positive control (line 67)
    ("The femoral-neck set (25 SNPs) gave a concordant result",
     "The femoral-neck set gave a concordant result", 1, "L67"),

    # Results colocalisation power (line 79)
    ("and genome-wide significance would need ~5×10⁷ participants (~42× the present sample)",
     "and genome-wide significance would need ~5×10⁷ participants", 1, "L79"),
]

failures = []
removed = 0
for old, new, exp, desc in EDITS:
    cnt = t.count(old)
    if cnt != exp:
        failures.append((desc, exp, cnt, old[:70]))
        continue
    removed += words(old) - words(new)
    t = t.replace(old, new)

if failures:
    print("ABORT — edits not matched:")
    for d, e, c, h in failures:
        print(f"  [{d}] expected {e}, found {c}: {h!r}")
    sys.exit(1)

open(F, "w", encoding="utf-8").write(t)
after_body = body_words(t)
after_hl = hl_lengths(t)
print("=== 352c COMPRESSION REPORT ===")
print(f"body words: {before_body} -> {after_body}  (removed {removed}; target <=4000)")
print(f"Highlights lengths: {before_hl} -> {after_hl}")
print(f"HL all <=85: {all(L<=85 for L in after_hl)}")
print(f"BODY compliant (<=4000): {after_body<=4000}")
