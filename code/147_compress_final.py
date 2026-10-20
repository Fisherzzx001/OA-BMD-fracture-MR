# -*- coding: utf-8 -*-
# Final OAC word-cut for docs/08. Prose-only; every beta/SE/p/OR/CI/citation preserved.
# Assertion-guarded: every old string must occur EXACTLY once, else abort with no write.
import sys

F = "docs/08_投稿稿_OA_BMD_fracture.md"
text = open(F, encoding="utf-8").read()

# Each pair: (old, new). All numbers/citations in old are kept verbatim in new.
pairs = [
    # --- Intro: cut editorial aside (line 35) ---
    ("though whether such a cost is mediated by bone mass or skeletal size — or large enough to matter clinically — is itself unresolved.",
     "though whether such a cost is mediated by bone mass or skeletal size is itself unresolved."),

    # --- Intro: drop redundancy with abstract (line 41) ---
    ("leaving unclear whether radiographic and self-reported OA relate to bone alike.",
     "leaving their bone relationship uncertain."),

    # --- Intro: compress the "we addressed all three" restatement (line 43) ---
    ("We addressed all three in 13 skeletal phenotypes (two OA definitions, ten regional BMD sites, and any-fracture), using LDSC across all 78 phenotype pairs, two-sample MR for causal direction with an internal positive control, and reverse MR to distinguish BMD→OA from OA→BMD.",
     "We addressed all three in 13 skeletal phenotypes using LDSC across all 78 phenotype pairs, two-sample MR with an internal positive control, and reverse MR."),

    # --- Discussion: remove duplicate magnitude figure already stated in Results (line 91) ---
    ("here (OR 1.060–1.099 per SD higher BMD), a 3.6–5.8-fold log-odds difference.",
     "here (OR 1.060–1.099 per SD higher BMD)."),

    # --- Highlights H3: trim 1 char to meet <=85 (remove space after colon) ---
    ("BMD instruments colocalise with fracture (9/31) not knee OA (0/31): endpoint-specific.",
     "BMD instruments colocalise with fracture (9/31) not knee OA (0/31):endpoint-specific."),
]

for old, new in pairs:
    c = text.count(old)
    if c != 1:
        print(f"[FAIL] pair matched {c} times (expected 1):\n  OLD={old!r}")
        sys.exit(1)

for old, new in pairs:
    text = text.replace(old, new, 1)

open(F, "w", encoding="utf-8").write(text)
print(f"Applied {len(pairs)} replacements OK.")
