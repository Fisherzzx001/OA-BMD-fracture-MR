#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
code/352 — main-text compression to keep body <= 4000 after M1-M6 additions.
Follows code/351 guardrails: dry-run default, --apply with backup,
each replacement must match EXACTLY ONCE.
Targets verbose / restated prose in Methods+Results+Discussion only
(the region code/148 counts as "body").
"""
import sys, os, shutil, datetime

SRC = "docs/08_投稿稿_OA_BMD_fracture.md"
BACKUP_DIR = ".workbuddy/backups"

def wc(s):
    import re
    return len([t for t in s.split() if re.search(r"[^\W_]", t, __import__("re").UNICODE)])

REPLS = [
    # L53 Methods: drop redundant knee-OA instrument attrition detail (in Fig S9)
    ("B_L53",
     "Harmonisation reduced the analysed sets to 16\u201331 instruments per forward pair and 69 for knee-OA-as-exposure; 47 of 116 knee-OA instruments (41%) were lost (Fig. S9). The 69 retained instruments were identical across all three knee-OA outcomes; the dropped instruments were not effect-size selected (median \u2212log10P 9.10 vs 9.89 retained; 43% stronger signal).",
     "Harmonisation left 16\u201331 instruments per forward pair and 69 for knee-OA-as-exposure, with no effect-size selection (Fig. S9)."),

    # L55 Methods: trim Burgess percentage tail (in Table S14)
    ("B_L55",
     "Primary IVW estimates are conservative for the combined bias, not overlap alone (per-pair intercepts and default-setting sensitivity in Table S14). Burgess overlap bias is \u22640.7% of the BMD\u2192knee-OA estimate and \u22644.5% of the positive control on the fracture-referenced overlap (\u226434% on the BMD-referenced basis; Table S14); MRlap retained all ten directions (median +20.5%), so the joint-endpoint effect survives overlap correction.",
     "Primary IVW estimates are conservative for the combined bias, not overlap alone (Table S14); MRlap retained all ten directions (median +20.5%), so the joint-endpoint effect survives overlap correction."),

    # L65 Results: delete variant-level agreement sentence (in Fig S24)
    ("B_L65a",
     "The definitions agreed at variant level: across 1,980 genome-wide-significant knee-OA loci, 94.4% had concordant signed Z with self-reported OA. ",
     ""),

    # L65 Results: delete 77/78 restatement (in Table S15)
    ("B_L65b",
     "77/78 pairs were significant at both BH-FDR and Bonferroni, the sole exception being self-reported-OA \u00d7 fracture. ",
     ""),

    # L67 Results: trim positive-control stability detail
    ("B_L67",
     "The estimate held under random effects (p=2.4\u00d710\u207b\u00b9\u00b2), the MR-Egger intercept was null (p=0.558), weighted median agreed (\u22120.364), and leave-one-out produced no sign reversal. ",
     "The estimate was stable (MR-Egger p=0.558; weighted median \u22120.364; leave-one-out no sign reversal). "),

    # L69 Results: trim random-effects convention detail (in Table S17)
    ("B_L69",
     "Heterogeneity was high (I\u00b2=68\u201386%); under random-effects IVW only two of five exposures remained significant (total-body p=0.006, lumbar-spine p=0.017; a second random-effects convention additionally retains arms, p=0.033), so the count depends on the random-effects convention (Table S17). All five effects remained positively causal after MR-PRESSO outlier removal, with four of five strengthening (Table S17).",
     "Heterogeneity was high (I\u00b2=68\u201386%); only two of five exposures remained significant under random-effects IVW (Table S17). All five remained positively causal after MR-PRESSO outlier removal, four of five strengthening (Table S17)."),

    # L81 Results: drop theta_H restatement (in Table S9)
    ("B_L81",
     " Height conditioning captured a genuine direct OA-raising effect (\u03b8_H=+0.145 to +0.165, p\u224810\u207b\u00b9\u2076 to 10\u207b\u00b9\u2079), so attenuation reflects size confounding, not dilution by a null covariate.",
     ""),

    # L89 Discussion: trim hip-estimate detail
    ("B_L89",
     "our cross-site hip estimate (same 31 instruments) gives BMD\u2192hip-OA IVW \u03b2=+0.097 (OR 1.102, 1.068\u20131.136), positive and knee-like (hip\u2212knee z=0.94, p=0.35; Table S27), restoring P6's pre-registered direction.",
     "our cross-site hip estimate (same 31 instruments) gives BMD\u2192hip-OA \u03b2=+0.097 (OR 1.102, 1.068\u20131.136), positive and knee-like (Table S27), restoring P6's pre-registered direction."),

    # L101 Conclusion: tighten prescribing-practice caveat
    ("B_L101",
     "these estimates do not support changing prescribing practice; one trial of a BMD-raising drug reported structural benefit [43], but strontium inflates DXA areal BMD as a measurement artefact and is restricted for cardiovascular risk, and TREASURE found no osteophyte progression despite greater BMD gain [44]; three years is not lifelong exposure.",
     "these estimates do not support changing prescribing practice; a BMD-raising trial reported structural benefit [43] but strontium inflates DXA areal BMD as an artefact and TREASURE found no osteophyte progression despite greater BMD gain [44]; three years is not lifelong exposure."),

    # Highlights: shorten to <=85 chars (code/148 rule)
    ("H_HL13",
     "- BMD variants colocalise with fracture (9/31) but not with knee OA (0/31) \u2014 a split driven by colocalisation power, not established endpoint-specificity (the strongest instrument gives expected knee-OA z\u22480.85, needing ~5\u00d710\u2077 participants for genome-wide significance).",
     "- BMD colocalises with fracture (9/31) but not knee-OA (0/31) \u2014 power, not endpoint-specific."),
    ("H_HL15",
     "- A systemic bone-forming set-point is a hypothesis whose directional prediction is consistent with data, but whose core claim \u2014 one variant set fixing opposite endpoint signs \u2014 is untested (F2 volume arm, F10).",
     "- A systemic bone-forming set-point is a hypothesis; its core claim is untested (F2,F10)."),
]

def main():
    with open(SRC, encoding="utf-8") as f:
        text = f.read()
    total = 0
    print("=== DRY/RUN preview (code/352) ===")
    for label, old, new in REPLS:
        n = text.count(old)
        d = wc(new) - wc(old)
        total += d
        print(f"[{'OK' if n==1 else ('FAIL' if n==0 else 'WARN'+str(n))}] {label}: delta={d:+d}")
    print(f"--- projected net body delta: {total:+d} ---")
    if "--apply" not in sys.argv:
        print("DRY RUN: no file written."); sys.exit(0)
    bad = [(l, text.count(o)) for l, o, n in REPLS if text.count(o) != 1]
    if bad:
        print("ABORT non-unique:", bad); sys.exit(3)
    os.makedirs(BACKUP_DIR, exist_ok=True)
    ts = datetime.datetime.now().strftime("%Y-%m-%d_%H%M%S")
    bak = f"{BACKUP_DIR}/docs08_pre352_{ts}.md"
    shutil.copy(SRC, bak)
    print("backup ->", bak)
    out = text
    for label, old, new in REPLS:
        out = out.replace(old, new)
    open(SRC, "w", encoding="utf-8").write(out)
    print("APPLIED", len(REPLS), "replacements.")

if __name__ == "__main__":
    main()
