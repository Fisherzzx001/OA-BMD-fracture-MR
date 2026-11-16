#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
451_v28_release_prep.py -- bump release metadata to v2.8 BEFORE tagging.

This is the defence against the v2.6 incident: last time .zenodo.json still
said v2.6 when the v2.6 tag was minted, so Zenodo labelled the record v2.5
and a corrective v2.6.1 re-archive was required. Here both .zenodo.json and
the hard-coded tag/release prose in code/318 are advanced to v2.8, then
validated, before any push happens.

Read-only with respect to the manuscript. Idempotent.
"""
import os
import json

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ZEN = os.path.join(ROOT, ".zenodo.json")
REL = os.path.join(ROOT, "code", "318_github_ssh_release.py")

SHA = "3fd452b300ff03a65fa8097ba63d050f479cd9b29d80ffd9cd7ecaf1593164f4"

DESC_TAIL = (
    "v2.8 (this archive) supersedes v2.7 and adds a second, statistics-focused "
    "remediation pass plus a new numerical-consistency gate; it inherits all v2.7 "
    "content (descriptive multi-phenotype genetic-architecture plus triangulation "
    "reframing with the Systemic Osteogenic Set-point downgraded to a post-hoc "
    "interpretive hypothesis, five Reviewer-2 mandatory items, and the eighth-round "
    "P0 N1/N2 fixes with nine P1/P2 items). The v2.8 delta resolves five MAJOR "
    "statistical defects: (N3) PP.H4>=0.8 is restated as a Bayesian posterior "
    "decision rule, not frequentist multiplicity control, and the paired-Wilcoxon "
    "p is labelled unadjusted/exploratory; (N4) weighted-median estimates are "
    "disclosed as matching the IVW sign while shrinking to about 0.35-0.43 of its "
    "magnitude at arms/legs (legs bootstrap CI -0.003 to +0.132 crosses zero), so "
    "sign concordance is no longer offered as robustness evidence; (N5) MR-PRESSO "
    "outlier removal is explicitly not removal of systematic/balanced horizontal "
    "pleiotropy (9/10 distortion tests significant), and Figures S11/S20 are "
    "demoted from robustness grids to comparison grids; (N6) MDE = 2.8 x SE is "
    "restored to a design-sensitivity quantity, the exclusion limits are the 95% "
    "CI upper bounds 0.0897/0.0937/0.1047, and the observed-power fallacy, "
    "random-effects SE, forward/reverse beta incommensurability and "
    "Steiger-versus-effect-size-MDE confusion are corrected; (N7) the cross-trait "
    "LDSC intercept lambda is an intercept, not a sample-overlap fraction, eight "
    "of ten values are negative, and Pearson r = -0.23 shows the away-from-null "
    "shift is driven by weak-instrument/winner's-curse correction rather than the "
    "lambda sign. Two new verification scripts accompany this archive: code/448 "
    "(141 cross-carrier assertions) and code/450 (an eighth gate of 65 cross-table "
    "numerical-consistency assertions derived only from authoritative CSVs and "
    "supplementary tables, including five negative self-tests). Nine verification "
    "checks remain green (148 COMPLIANT: abstract 346/350, body 3992/4000, 50 refs, "
    "8 display items; 149b SAFE; 182 79/79; 310 99/99 NO REGRESSION; 331 31/31; "
    "334 18/18; 336 docx Word-like 3997/4000 COMPLIANT; code/448 141/141; "
    "code/450 65/65). No reported beta/SE/p/OR/CI changed relative to v2.7; no "
    "statistical analysis was re-run. Raw GWAS summary statistics are not "
    "redistributed; they are publicly published (cited in the manuscript) and "
    "re-derived by the pipeline."
)

NEW_NOTES = (
    "Submission v2.8 (2026-11-16). This is a substantive statistics-focused "
    "remediation archive, NOT a metadata-only re-archive: .zenodo.json version was "
    "bumped to v2.8 BEFORE tagging, so the minted record carries the correct "
    "version label (the earlier v2.6 -> v2.6.1 correction is not repeated). Two "
    "local commits since tag v2.7 (d745824, 4537d4c), plus this release-preparation "
    "commit, implement the N3-N7 second-pass remediation of five MAJOR statistical "
    "defects and add the eighth gate. All changes are confined to definitions, "
    "table captions/notes and interpretive self-consistency: no reported "
    "beta/SE/p/OR/CI changed and no statistical analysis was re-run. N3: PP.H4>=0.8 "
    "is a posterior rule, Wilcoxon p unadjusted/exploratory. N4: weighted-median "
    "magnitude shrinkage and legs CI crossing zero propagated to Abstract/Results/"
    "Strengths; sign-only agreement no longer robustness evidence. N5: MR-PRESSO "
    "outlier removal is not robustness against systematic pleiotropy; S11/S20 "
    "demoted to comparison grids. N6: MDE restored to design sensitivity; exclusion "
    "limits are 95% CI upper bounds 0.0897/0.0937/0.1047; observed-power fallacy, "
    "random-effects SE, beta-unit incommensurability and Steiger clarified. N7: "
    "lambda is an intercept, 8/10 negative; Pearson r -0.23 attributes the "
    "away-from-null shift to winner's-curse correction. New code/448 (141 "
    "cross-carrier checks) and code/450 (65 cross-table numerical-consistency "
    "checks with five negative self-tests). Nine checks green (148 COMPLIANT: "
    "abstract 346/350, body 3992/4000 / 149b SAFE / 182 79/79 / 310 99/99 / "
    "331 31/31 / 334 18/18 / 336 docx 3997 / 448 141/141 / 450 65/65). Manuscript "
    "SHA256 " + SHA + ". Concept DOI: https://doi.org/10.5281/zenodo.23400860 "
    "(always resolves to the most recent version; v2.7 archive 24248841 and earlier "
    "versions remain available)."
)

TAG_MSG = [
    "Submission v2.8 (2026-11-16): second-pass remediation of five MAJOR statistical",
    "defects (N3-N7) plus a new eighth gate for cross-table numerical consistency.",
    ".zenodo.json version was bumped to v2.8 BEFORE tagging so the Zenodo record",
    "carries the correct version label. N3: PP.H4>=0.8 is a posterior decision rule,",
    "not multiplicity control; Wilcoxon p unadjusted/exploratory. N4: weighted-median",
    "magnitude shrinkage at arms/legs (0.35/0.43) and legs CI crossing zero disclosed;",
    "sign concordance is not robustness evidence. N5: MR-PRESSO outlier removal is not",
    "removal of systematic pleiotropy; S11/S20 demoted to comparison grids. N6:",
    "MDE=2.8xSE is design sensitivity, not an exclusion bound; exclusion limits are",
    "95% CI upper bounds 0.0897/0.0937/0.1047; observed-power fallacy, random-effects",
    "SE, forward/reverse beta incommensurability and Steiger clarified. N7: lambda is a",
    "cross-trait LDSC intercept, not an overlap fraction, 8/10 negative; Pearson",
    "r -0.23 shows the away-from-null shift is winner's-curse driven. New code/448",
    "(141 cross-carrier assertions) and code/450 (65 numerical-consistency assertions",
    "with five negative self-tests). Two local commits since tag v2.7 (d745824,",
    "4537d4c) plus this release-preparation commit. Wording, table-definition and",
    "interpretive changes only: no reported beta/SE/p/OR/CI changed, no statistical",
    "analysis re-run. Nine checks green (148 346/3992 / 149b SAFE / 182 79/79 /",
    "310 99/99 / 331 31/31 / 334 18/18 / 336 docx 3997 / 448 141/141 / 450 65/65).",
    "Concept DOI: https://doi.org/10.5281/zenodo.23400860",
]

RELEASE_BODY = [
    "## Submission v2.8 (2026-11-16) - second-pass MAJOR statistics remediation + eighth gate",
    "",
    "This is a **substantive** archive, not a metadata-only re-archive. `.zenodo.json`",
    "was bumped to v2.8 **before** tagging, so the GitHub-Zenodo integration mints a",
    "record with the correct version label (the earlier v2.6 -> v2.6.1 correction is",
    "not repeated). Two local commits since tag v2.7 (`d745824`, `4537d4c`), plus this",
    "release-preparation commit.",
    "",
    "All changes are confined to **definitions, table captions/notes and interpretive",
    "self-consistency**: no reported beta/SE/p/OR/CI changed, and no statistical",
    "analysis was re-run.",
    "",
    "### Five MAJOR statistical defects (second pass)",
    "- **N3** - `PP.H4>=0.8` is a Bayesian posterior decision rule, **not** frequentist",
    "  multiplicity control; the paired-Wilcoxon p is labelled unadjusted/exploratory.",
    "- **N4** - Weighted median matches the IVW sign but shrinks to about 0.35-0.43 of",
    "  its magnitude at arms/legs (legs bootstrap CI -0.003 to +0.132 crosses zero);",
    "  sign concordance is no longer offered as robustness evidence.",
    "- **N5** - MR-PRESSO outlier removal does **not** remove systematic/balanced",
    "  horizontal pleiotropy (9/10 distortion tests significant); Figures S11/S20",
    "  demoted from robustness grids to comparison grids.",
    "- **N6** - MDE = 2.8 x SE restored to a **design-sensitivity** quantity; exclusion",
    "  limits are the 95% CI upper bounds **0.0897 / 0.0937 / 0.1047**; observed-power",
    "  fallacy, random-effects SE, forward/reverse beta incommensurability and",
    "  Steiger-versus-MDE confusion corrected.",
    "- **N7** - Cross-trait LDSC intercept lambda is an intercept, **not** a",
    "  sample-overlap fraction; 8/10 values negative; Pearson r = -0.23 attributes the",
    "  away-from-null shift to weak-instrument/winner's-curse correction rather than",
    "  the lambda sign.",
    "",
    "### New verification",
    "- `code/448`: 141 cross-carrier assertions (manuscript, mirror, two docx, dual",
    "  supplementary copies).",
    "- `code/450`: eighth gate, 65 cross-table numerical-consistency assertions derived",
    "  only from authoritative CSVs and supplementary tables, with five in-memory",
    "  negative self-tests proving each detector fires on corrupted data. This closes",
    "  the blind spot where the seven gates checked word counts and lock strings but",
    "  not cross-table arithmetic.",
    "",
    "Full audit trail: `docs/149` sections 8.10 / 8.10.6; review report `docs/151`.",
    "Audit trail is kept out of the manuscript.",
    "",
    "Nine checks green (148 COMPLIANT: abstract 346/350, body 3992/4000, 50 refs, 8",
    "display items / 149b SAFE / 182 79/79 / 310 99/99 NO REGRESSION / 331 31/31 /",
    "334 18/18 / 336 docx Word-like 3997 COMPLIANT / code/448 141/141 / code/450",
    "65/65). Concept DOI: https://doi.org/10.5281/zenodo.23400860 (resolves to the",
    "latest version; earlier archives including v2.7 / 24248841 remain available).",
]


def jstr(s):
    """JSON string literal without the surrounding quotes."""
    return json.dumps(s, ensure_ascii=False)[1:-1]


def prep_zenodo():
    raw = open(ZEN, encoding="utf-8").read()
    if '"version": "v2.8"' in raw and "v2.8 (this archive) supersedes" in raw:
        print("[zenodo] already v2.8 -- SKIP")
        return
    if '"version": "v2.7"' not in raw:
        raise SystemExit("[zenodo] unexpected version; expected v2.7")

    raw = raw.replace('"publication_date": "2026-11-15"',
                      '"publication_date": "2026-11-16"', 1)

    # description: keep everything up to the v2.7 narrative, replace the tail
    d0 = raw.index('"description": "') + len('"description": "')
    dcut = raw.index("v2.7 (this archive) supersedes", d0)
    dend = raw.index('",\n  "access_right"', dcut)
    preserved = raw[d0:dcut]
    raw = raw[:d0] + preserved + jstr(DESC_TAIL) + raw[dend:]

    # notes: replace the whole value
    n0 = raw.index('"notes": "') + len('"notes": "')
    nend = raw.index('"\n}', n0)
    raw = raw[:n0] + jstr(NEW_NOTES) + raw[nend:]

    raw = raw.replace('"version": "v2.7"', '"version": "v2.8"', 1)
    obj = json.loads(raw)  # hard validation
    open(ZEN, "w", encoding="utf-8", newline="\n").write(raw)
    print("[zenodo] version=%s pubdate=%s desc_len=%d notes_len=%d"
          % (obj["version"], obj["publication_date"],
             len(obj["description"]), len(obj["notes"])))


def replace_py_block(src, varname, lines, end_marker, inline_space=False):
    # The original blocks close their final string literal with ')' on the
    # SAME line (no standalone '        )' line), so delimit by the next
    # top-level statement instead.
    start = "        %s = (\n" % varname
    i = src.index(start)
    j = src.index(end_marker, i)
    content = list(lines)
    if inline_space:
        content = [x + " " for x in content[:-1]] + [content[-1]]
    # Python source needs full quoted literals; json.dumps is a valid
    # Python string literal for these ASCII/UTF-8 payloads.
    rendered = ["            " + json.dumps(x, ensure_ascii=False)
                for x in content]
    rendered[-1] += ")"
    block = start + "\n".join(rendered)
    return src[:i] + block + src[j:]


def prep_release_script():
    src = open(REL, encoding="utf-8").read()
    if 'RELEASE_TAG", "v2.8")' in src and TAG_MSG[0] in src:
        print("[code/318] already v2.8 -- SKIP")
        return
    if 'RELEASE_TAG", "v2.7")' not in src:
        raise SystemExit("[code/318] unexpected default TAG; expected v2.7")

    src = src.replace('RELEASE_TAG", "v2.7")', 'RELEASE_TAG", "v2.8")', 1)
    src = replace_py_block(
        src, "tag_msg", TAG_MSG,
        '\n        rc, out, err = git(["tag", "-a", TAG, "-m", tag_msg])',
        inline_space=True)
    src = replace_py_block(
        src, "bodytxt", RELEASE_BODY,
        '\n        rel_name = "submission-" + TAG')
    compile(src, REL, "exec")  # hard validation
    open(REL, "w", encoding="utf-8", newline="\n").write(src)
    print("[code/318] default TAG=v2.8, tag_msg/bodytxt rewritten, syntax OK")


if __name__ == "__main__":
    prep_zenodo()
    prep_release_script()
    print("DONE")
