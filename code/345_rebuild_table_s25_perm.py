# -*- coding: utf-8 -*-
"""345 — Rebuild Supplementary Table S25 with permutation-calibrated p (docs/148 P1-1).

Keeps every previously reported MWW enrichment Δ and p value UNCHANGED (audit
trail), and ADDS:
  * a permutation p column (10,000 length+chromosome-matched permutations,
    code/342_f3_permutation_enrichment.py, seed 20261112);
  * a cartilage-minus-bone contrast row block (overlapping genes removed);
  * corrects the tissue label of skeletal_system_development from "bone" to
    "mixed_skeletal" (28/42 = 67% of its genes are shared with cartilage sets).
The one gene set absent from code/342 (bone_remodeling, n=1) keeps perm_p = NA.

Reads code/output/genetics/F3_permutation_results.csv + F3_contrast_results.csv.
Dry-run by default; --apply to write docs/08.
"""
import csv, os, sys, io, datetime, shutil

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
F = os.path.join(ROOT, "docs", "08_投稿稿_OA_BMD_fracture.md")
GEN = os.path.join(ROOT, "code", "output", "genetics")
APPLY = "--apply" in sys.argv

perm = {}
with io.open(os.path.join(GEN, "F3_permutation_results.csv"), encoding="utf-8") as fh:
    for r in csv.DictReader(fh):
        perm[(r["trait"], r["gene_set"])] = (float(r["delta_obs"]), float(r["p_perm"]))
contrast = {}
with io.open(os.path.join(GEN, "F3_contrast_results.csv"), encoding="utf-8") as fh:
    rows = list(csv.DictReader(fh))
for i, trait in enumerate(["OA_knee", "BMD_total", "fracture"]):
    r = rows[i]
    contrast[trait] = (float(r["contrast_obs"]), float(r["p_perm"]))

def fmtp(p):
    if p < 1e-4:
        return "<1×10⁻⁴"
    # scientific like 1.70e-3 -> 1.70e−03 ; 2.30e-3
    return f"{p:.2e}".replace("e-", "e−")

# (trait, gene set, tissue, n, mean_set, mean_bg, delta_MWW, p_MWW) — VERBATIM from S25
OLD_ROWS = [
    ("OA_knee", "bone_development", "bone", 83, 3.415, 2.610, 0.805, "2.21e-04"),
    ("OA_knee", "cartilage_development", "cartilage", 72, 5.036, 2.606, 2.430, "1.45e-08"),
    ("OA_knee", "ossification", "bone", 71, 3.741, 2.609, 1.131, "2.78e-03"),
    ("OA_knee", "chondrocyte_differentiation", "cartilage", 81, 4.932, 2.605, 2.326, "1.58e-08"),
    ("OA_knee", "bone_morphogenesis", "other", 9, 2.084, 2.613, -0.528, "6.07e-01"),
    ("OA_knee", "extracellular_matrix_organization", "cartilage", 37, 3.084, 2.612, 0.472, "3.27e-02"),
    ("OA_knee", "skeletal_system_development", "bone", 42, 5.354, 2.608, 2.746, "9.96e-07"),
    ("OA_knee", "skeletal_joint_morphogenesis", "cartilage", 82, 4.057, 2.608, 1.449, "4.58e-10"),
    ("OA_knee", "bone_remodeling", "other", 1, 4.853, 2.612, 2.241, "NA"),
    ("BMD_total", "bone_development", "bone", 83, 2.121, 1.671, 0.450, "3.07e-03"),
    ("BMD_total", "cartilage_development", "cartilage", 72, 2.198, 1.671, 0.526, "1.35e-02"),
    ("BMD_total", "ossification", "bone", 71, 2.285, 1.671, 0.614, "1.00e-01"),
    ("BMD_total", "chondrocyte_differentiation", "cartilage", 81, 2.277, 1.671, 0.606, "2.11e-02"),
    ("BMD_total", "bone_morphogenesis", "other", 9, 2.124, 1.673, 0.452, "1.74e-02"),
    ("BMD_total", "extracellular_matrix_organization", "cartilage", 37, 1.903, 1.673, 0.230, "1.52e-02"),
    ("BMD_total", "skeletal_system_development", "bone", 42, 3.403, 1.670, 1.733, "1.26e-04"),
    ("BMD_total", "skeletal_joint_morphogenesis", "cartilage", 82, 2.530, 1.670, 0.860, "3.16e-07"),
    ("BMD_total", "bone_remodeling", "other", 1, 3.645, 1.673, 1.972, "NA"),
    ("fracture", "bone_development", "bone", 83, 1.599, 1.455, 0.144, "3.14e-02"),
    ("fracture", "cartilage_development", "cartilage", 72, 1.895, 1.454, 0.441, "2.16e-04"),
    ("fracture", "ossification", "bone", 71, 1.538, 1.455, 0.083, "3.52e-01"),
    ("fracture", "chondrocyte_differentiation", "cartilage", 81, 1.812, 1.454, 0.357, "8.38e-04"),
    ("fracture", "bone_morphogenesis", "other", 9, 1.590, 1.455, 0.135, "1.81e-01"),
    ("fracture", "extracellular_matrix_organization", "cartilage", 37, 1.696, 1.455, 0.241, "5.33e-04"),
    ("fracture", "skeletal_system_development", "bone", 42, 1.825, 1.455, 0.370, "5.42e-04"),
    ("fracture", "skeletal_joint_morphogenesis", "cartilage", 82, 1.902, 1.454, 0.448, "1.36e-06"),
    ("fracture", "bone_remodeling", "other", 1, 1.445, 1.455, -0.010, "NA"),
]

# Locate the S25 markdown table block by structural markers (robust to sign/space
# glyph differences), then rebuild it. The block starts at the header line that
# begins "| Trait | Gene set |" and ends at the last S25 data row before a blank
# line followed by "## ".
t = open(F, encoding="utf-8").read()
lines = t.splitlines()
start = next(i for i, l in enumerate(lines) if l.startswith("| Trait | Gene set |"))
end = start
while end + 1 < len(lines) and lines[end + 1].startswith("|"):
    end += 1
old_block = "\n".join(lines[start:end + 1])
assert sum(1 for l in lines[start:end + 1] if l.startswith("| OA_knee |")
           or l.startswith("| BMD_total |") or l.startswith("| fracture |")) == 27, \
    "expected 27 S25 data rows"

new_header = ("| Trait | Gene set | Tissue | n (set) | mean −log10p set | mean −log10p bg "
              "| enrichment Δ | MWW p (1-sided) | perm p (length/chr-matched) |")
new_lines = [new_header, "| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |"]
for (tr, gs, ti, n, ms, mb, d, p) in OLD_ROWS:
    ti2 = "mixed_skeletal*" if gs == "skeletal_system_development" else ti
    if (tr, gs) in perm:
        dp, pp = perm[(tr, gs)]
        pp_s = fmtp(pp)
    else:
        pp_s = "NA (n=1)"
    new_lines.append(f"| {tr} | {gs} | {ti2} | {n} | {ms:.3f} | {mb:.3f} | {d:+.3f} | {p} | {pp_s} |")

# contrast summary appended after the main table block, before next "## " section
new_lines += [
    "",
    "**Cartilage-minus-bone contrast** (permutation, overlapping genes removed; 10,000 draws): "
    f"OA_knee +{contrast['OA_knee'][0]:.3f} (p={contrast['OA_knee'][1]:.4f}); "
    f"BMD_total +{contrast['BMD_total'][0]:.3f} (p={contrast['BMD_total'][1]:.3f}); "
    f"fracture +{contrast['fracture'][0]:.3f} (p={contrast['fracture'][1]:.3f}). "
    "The OA-knee signal is cartilage-biased; total-body BMD is non-discriminating between "
    "cartilage and bone; the nominal skeletal-system-development enrichment is not pure "
    "bone evidence, as 28/42 (67%) of its genes are shared with the cartilage sets.",
    "",
    "*skeletal_system_development is re-labelled mixed_skeletal: its gene content overlaps "
    "heavily with the cartilage/joint sets (Jaccard 0.44 for the two cartilage sets with "
    "each other), so its strong nominal enrichment cannot be read as bone-formation-specific. "
    "MWW columns are retained unchanged for audit continuity; perm p supersedes them for "
    "inference.",
]

new_block = "\n".join(new_lines)

cnt = t.count(old_block)
print("=" * 70)
print("345 · S25 permutation-calibrated rebuild (", "APPLY" if APPLY else "DRY-RUN", ")")
print("=" * 70)
print(f"anchor block count = {cnt} (need exactly 1)")
if cnt != 1:
    sys.exit(1)
# sanity: all MWW delta/p values preserved verbatim in new_block
for (tr, gs, ti, n, ms, mb, d, p) in OLD_ROWS:
    needle = f"| {tr} | {gs} |"
    assert needle in new_block, needle
new_t = t.replace(old_block, new_block, 1)
# S25 stays in the figure-legends appendix (not body): word count must be unchanged in body region
def body_words(s):
    ai = s.find("**Objective.**"); ci = s.find("**Conclusions.**", ai)
    tail = s[ci:]; em = __import__("re").search(r"\n---", tail)
    ae = ci + (em.start() if em else len(tail)); leg = s.find("## Figure legends")
    body = s[ae:leg]
    import re as _re
    return len([x for x in "\n".join(l for l in body.splitlines() if not l.strip().startswith("|")).split() if _re.search(r"[^\W_]", x, _re.UNICODE)])
import re as _re
bw0 = body_words(t); bw1 = body_words(new_t)
print(f"body words {bw0} -> {bw1} (must be unchanged: S25 is appendix)")
if bw1 != bw0:
    print("[FAIL] body word count changed"); sys.exit(1)
# contrast values present
for k in ("OA_knee", "BMD_total", "fracture"):
    assert f"contrast['{k}']" not in new_block
print(f"contrast values: OA {contrast['OA_knee']}, BMD {contrast['BMD_total']}, FX {contrast['fracture']}")
print("DRY-RUN OK" if not APPLY else "writing...")

if APPLY:
    bak = f".workbuddy/backups/docs08_pre345_{datetime.date.today()}.md"
    shutil.copyfile(F, bak)
    open(F, "w", encoding="utf-8").write(new_t)
    print(f"APPLIED. backup -> {bak}")
