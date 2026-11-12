# -*- coding: utf-8 -*-
"""
code/338 - DRAFT revisions from the sixth pre-review committee (2026-11-12).

Purpose: draft word-neutral, number-neutral edits addressing P0 x3 / P1 x3 / P2 box,
WITHOUT applying them to the single source of truth. Dry-run ONLY by default:
reports per-edit word delta and numeric-literal multiset delta, and enforces
  (a) total IMRaD word delta <= 0   (headroom is only 35 words: 3965/4000)
  (b) the multiset of ALL 673 numeric literals is EXACTLY unchanged
  (c) no reported beta/SE/p/OR/CI token is altered (subset of (b))
If you later want to apply, run with --apply (writes docs/08) — but that is a
separate, human-approved step and must still pass 148/149b/182/310/331/334/336.

Edits map to docs/146 (ledger) sections. Each item: (id, old, new).
"""
import re, sys
from collections import Counter

F = "docs/08_投稿稿_OA_BMD_fracture.md"

# ---- exact-string replacements (must match docs/08 verbatim) ----
REPL = [
    # ===== P1-7 (R3-3稻草人): "instead" -> consistent-with. Word-neutral-ish. =====
    ("P1-7",
     "and a 2025 MR reported total BMD as a positive knee-OA mediator (3.59% mediation) [11]. Our analysis instead finds a positive BMD→knee-OA estimate (discussed below).",
     "and a 2025 MR reported total BMD as a positive knee-OA mediator (3.59% mediation) [11], consistent with our positive BMD→knee-OA estimate (discussed below)."),

    # ===== P0-1 (R2-1 K循环): disclose SE-ratio identification + point to S17-a.
    #      Net cost paid by tightening the same paragraph elsewhere. =====
    ("P0-1a",
     "The Z ratio is unity by construction, a same-source check that cannot identify K.",
     "The Z ratio is unity by construction and cannot identify K; the log-odds SE does (K-sensitivity in the supplement; p-values K-invariant)."),
    ("P0-1b",
     "Two verified scale corrections were applied.",
     "Two scale corrections were applied."),

    # ===== P0-2 (R1-1概念置换): split observational vs mechanistic; drop the
    #      bisphosphonate inference that backfires. Net-neutral by trimming. =====
    ("P0-2a",
     "A systemic, polygenic bone-forming set-point readable as areal BMD would protect at the whole-bone endpoint and raise OA risk at the joint — the Systemic Osteogenic Set-point (SOS) (Figs. S21, S23), which we pose as a testable hypothesis distinct in kind from the bone-forming triad (osteophytes, enthesophytes and high bone mass; OR 3.00, 95% CI 1.96–4.58 [39]) and from the broader osteoporosis–osteoarthritis 'bone–cartilage paradox' [3, 4].",
     "We distinguish an observed from a mechanistic claim: areal BMD shows opposite endpoint signs (established here), whereas one systemic polygenic osteogenic set-point (SOS; Figs. S21, S23) is the untested hypothesis proposed to generate them — distinct in kind from the bone-forming triad (osteophytes, enthesophytes, high bone mass; OR 3.00, 95% CI 1.96–4.58 [39]) and the 'bone–cartilage paradox' [3, 4]."),
    ("P0-2b",
     "anti-resorptives also raise areal BMD, yet bisphosphonate users had a lower knee-replacement rate (RR 0.71, 0.66–0.77) [40], implicating formation rather than density.",
     "anti-resorptives also raise areal BMD yet bisphosphonate users had lower knee-replacement (RR 0.71, 0.66–0.77) [40], so areal BMD cannot by itself distinguish formation from resorption."),

    # ===== P0-3 (R2-2/R2-4部位外推+F2措辞): limit to five sites; soften F2. =====
    ("P0-3a",
     "The effect concerns the clinically measured areal-BMD phenotype at the five analysed sites; 67–91% of it is height- and BMI-independent, and not transmitted by bone geometry (F2), albeit on weak conditional instruments (conditional F below 4 in all but one mode), so this is direction-indicative rather than magnitude-precise (Table S28, Part B);",
     "The effect concerns the clinically measured areal-BMD phenotype and is established only at the five sites entering MR; 67–91% is height- and BMI-independent, with conditional F2 analysis direction-indicative not conclusive (conditional F below 4 in all but one mode; Table S28, Part B);"),
    ("P0-3b",
     "so geometry does not absorb the signal;",
     "so geometry does not absorb it, albeit on weak conditional instruments;"),

    # ===== P1-5 (R1-2/C证据两用): state orthogonality is neutral, not supportive. =====
    ("P1-5",
     "The discriminating kernel is not a label for the observed pattern; it is non-trivial because BMD, OA and fracture occupy distinct genetic dimensions (S22: the PC1 bone-mass/fracture axis is orthogonal to the PC2 OA axis), so a single set-point cannot be reduced to a generic bone–OA correlation.",
     "The discriminating kernel is not a label for the pattern; BMD, OA and fracture occupy distinct genetic dimensions (S22: PC1 bone-mass/fracture axis orthogonal to PC2 OA axis), but orthogonality alone is neutral about a single set-point — it neither establishes nor refutes one, consistent also with two coupled modules."),

    # ===== P1-4 (R4-E "9/10"计数口径): make explicit that only F4 discriminates
    #      and survival counts executed tests, not register completeness. =====
    ("P1-4",
     "P4 is the discriminating prediction; F2, its decisive test, remains unexecuted.",
     "P4 alone discriminates (via colocalisation); F2's decisive arm is unexecuted, so survival counts executed tests only."),

    # ===== P2箱 =====
    # P2-9 身高/BMI代理措辞 (L81): already careful; no numeric change, minor tighten
    ("P2-9a",
     "Because height indexes stature not bone width/geometry [28], this residual is not identified as bone mass; F2 addresses geometry (§4, Table S28).",
     "Because height indexes stature not bone width/geometry [28], this residual is not identified as bone mass; F2 probes geometry (§4, Table S28)."),
    # P2-10 "可干预 vs 不改处方"矛盾 (L29): resolve by explicit qualifier
    ("P2-10",
     "The joint effect is modest per SD (OR 1.06–1.10 areal BMD) yet attaches to a ubiquitous, modifiable exposure. These estimates do not support changing prescribing practice;",
     "The joint effect is modest per SD (OR 1.06–1.10 areal BMD) at a ubiquitous exposure. These estimates do not support changing prescribing practice;"),
    # P2-12 Steiger MDE (L73): disclose underpowered rather than bare "lacks power"
    ("P2-12",
     "That non-result does not discriminate direction: with 69 instruments and a null reverse estimate, the test lacks power and is not read as evidence for BMD over OA as the driver.",
     "That non-result does not discriminate direction: with 69 instruments and a null reverse estimate the test is underpowered (minimum detectable effect exceeds the null), not evidence BMD drives OA."),
]

# zeros / notes: items intentionally left as no-op are skipped in application.
NOOP = {id for id, old, new in REPL if old == new}


def imrad_segments(text):
    L = text.split("\n")
    s = next(i for i, x in enumerate(L) if x.startswith("## 1. Introduction"))
    e = next(i for i, x in enumerate(L) if x.startswith("## Figure legends"))
    return s, e, L


def wc(seg):
    return len([x for x in seg.split() if re.search(r"[^\W_]", x, re.UNICODE)])


def nums(seg):
    return Counter(re.findall(r"\d+(?:[.,]\d+)*", seg))


def main():
    apply = "--apply" in sys.argv
    text = open(F, encoding="utf-8").read()
    s, e, L = imrad_segments(text)
    body0 = "\n".join(x for x in L[s:e] if not x.strip().startswith("|"))
    w0, n0 = wc(body0), nums(body0)

    newtext = text
    applied, skipped = [], []
    for rid, old, new in REPL:
        if rid in NOOP:
            continue
        cnt = newtext.count(old)
        if cnt != 1:
            skipped.append((rid, f"match count={cnt}"))
            continue
        cand = newtext.replace(old, new, 1)
        # recompute body numeric multiset on candidate
        _, _, Lc = imrad_segments(cand)
        bodyc = "\n".join(x for x in Lc[s:e] if not x.strip().startswith("|"))
        nc = nums(bodyc)
        if nc != n0:
            skipped.append((rid, f"numeric multiset changed: missing={list((n0-nc).elements())[:6]} extra={list((nc-n0).elements())[:6]}"))
            continue
        dw = wc(bodyc) - w0
        # update baseline word count as we go so cumulative delta is tracked
        applied.append((rid, dw, wc(bodyc)))
        w0 = wc(bodyc)
        newtext = cand

    _, _, Le = imrad_segments(newtext)
    bodyf = "\n".join(x for x in Le[s:e] if not x.strip().startswith("|"))
    wf, nf = wc(bodyf), nums(bodyf)

    print("=" * 78)
    print(f"{'id':8} {'Δwords':>8} {'running total':>14}")
    print("-" * 78)
    run = 0
    for rid, dw, tot in applied:
        run += dw
        print(f"{rid:8} {dw:+8d} {run:+14d}")
    print("-" * 78)
    orig_total = wc("\n".join(x for x in L[s:e] if not x.strip().startswith("|")))
    print(f"IMRaD words: {orig_total} -> {wf}  (delta {wf-orig_total:+d}; headroom to 4000 = {4000-wf})")
    print(f"numeric-literal multiset identical: {nf == nums(chr(10).join(x for x in L[s:e] if not x.strip().startswith('|')))}")
    if skipped:
        print("\nSKIPPED (did not meet guardrails):")
        for rid, why in skipped:
            print(f"  {rid}: {why}")
    print("=" * 78)
    if nf != nums("\n".join(x for x in L[s:e] if not x.strip().startswith("|"))):
        print("FAIL: numeric literals changed"); return 2
    if wf > 4000:
        print("FAIL: over 4000 words"); return 3
    print("DRAFT GUARDRAILS PASS (word-neutral, number-neutral)")
    if apply:
        open(F, "w", encoding="utf-8").write(newtext)
        print("APPLIED to", F)
    else:
        print("(dry-run only; docs/08 NOT modified — use --apply after human approval)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
