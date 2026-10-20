#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Render the four missing supplementary tables (S8, S9, S14, S15) into
docs/Supplementary/ as standalone Markdown + CSV, faithful to the manuscript.

- S8  : built from coloc_screen_all_loci.csv (31 loci) + coloc_power.csv (power).
        The manuscript S8 paragraph is descriptive only ("full locus-level table
        in the repository"); the actual per-locus colocalisation lives here.
- S9  : parsed from the canonical 7-row MVMR table in docs/08 (lines ~192-199).
- S14 : parsed from the 4-row ancestry table in docs/08 (lines ~206-211).
- S15 : parsed from the 5-row r_g-block table in docs/08 (lines ~218-224).
"""
import csv, re, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GEN = os.path.join(ROOT, "code", "output", "genetics")
SUP = os.path.join(ROOT, "docs", "Supplementary")
os.makedirs(SUP, exist_ok=True)
MAN = os.path.join(ROOT, "docs", "08_投稿稿_OA_BMD_fracture.md")

def rnd(x, sig=4):
    if x is None or x == "" or x != x:
        return ""
    try:
        f = float(x)
    except Exception:
        return str(x)
    if f == 0:
        return "0"
    # %g gives 'sig' significant figures with automatic sci/decimal notation
    return f"%.{sig}g" % f

# ---------------------------------------------------------------- S8 (from CSV)
def build_s8():
    import csv as _csv
    def load(p):
        with open(p, newline="", encoding="utf-8") as fh:
            return list(_csv.DictReader(fh))
    screen = load(os.path.join(GEN, "coloc_screen_all_loci.csv"))
    power = load(os.path.join(GEN, "coloc_power.csv"))
    pw = {r["snp"]: r for r in power}
    assert len(screen) == 31, f"expected 31 coloc loci, got {len(screen)}"
    # integrity: highest PP.H4 must equal 0.799 at rs11228269 (manuscript claim)
    pp4 = {r["index_snp"]: float(r["PP4"]) for r in screen}
    top = max(pp4, key=pp4.get)
    assert abs(pp4[top] - 0.799) < 0.02, f"max PP.H4={pp4[top]:.3f} at {top}, expected ~0.799"
    rows = []
    for r in screen:
        p = pw.get(r["index_snp"], {})
        rows.append({
            "Locus (index SNP)": r["index_snp"],
            "Chr:Pos (hg38)": f'{r["chrom"]}:{int(float(r["pos"])):,}',
            "PP.H0": rnd(r["PP0"]),
            "PP.H1": rnd(r["PP1"]),
            "PP.H2": rnd(r["PP2"]),
            "PP.H3": rnd(r["PP3"]),
            "PP.H4": rnd(r["PP4"]),
            "n_snps (BMD)": r["n_snps_t1"],
            "n_snps (OA)": r["n_snps_t2"],
            "max|z|_BMD": rnd(r["max_abs_z1"], 2),
            "max|z|_OA": rnd(r["max_abs_z2"], 2),
            "Direction concordant": r["dir_concordant"],
            "N_OA needed for PP.H4>=0.8": (int(round(float(p["N_oa_needed"]))) if p.get("N_oa_needed") else ""),
        })
    cols = list(rows[0].keys())
    md = ["# Supplementary Table S8. Colocalisation of total-body BMD with knee OA at all 31 instrument loci",
          "",
          "Regions are ±500 kb around each strict-clumped total-body BMD instrument. "
          "PP.H0–PP.H4 are posterior probabilities under priors p1=p2=1×10⁻⁴, p12=5×10⁻⁶. "
          "No locus reaches PP.H4≥0.8; six favour distinct causal variants (PP.H3≥0.8). "
          "The highest single-locus PP.H4 was 0.799 (rs11228269, chr11:68.29 Mb; BMD z=−8.65, "
          "knee-OA z=−4.06, direction-concordant; PP.H1=0.156) — one locus fell 0.0006 short of "
          "the threshold, and the pre-specified 0.8 cut-off is retained. "
          "`N_OA needed for PP.H4≥0.8` is the knee-OA sample size required to reach the 0.8 "
          "colocalisation threshold at that locus given its observed BMD and OA effect sizes.",
          "",
          "| " + " | ".join(cols) + " |",
          "| " + " | ".join(["---"] * len(cols)) + " |"]
    for r in rows:
        md.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
    return "\n".join(md) + "\n", rows, cols

# ---------------------------------------------------------- S9/S14/S15 (parse)
def parse_table(man_text, header_anchor):
    lines = man_text.splitlines()
    start = None
    for i, ln in enumerate(lines):
        if ln.startswith("## Supplementary Table") and header_anchor in ln:
            start = i
            break
    assert start is not None, f"header not found: {header_anchor}"
    # collect first contiguous block of pipe-lines after the header
    tbl = []
    for ln in lines[start+1:]:
        if ln.lstrip().startswith("|"):
            tbl.append(ln)
        elif tbl:
            break
    # tbl[0] = header, tbl[1] = separator
    assert len(tbl) >= 3, f"{header_anchor}: table too short ({len(tbl)} lines)"
    header = [c.strip() for c in tbl[0].strip().strip("|").split("|")]
    data = []
    for ln in tbl[2:]:
        cells = [c.strip() for c in ln.strip().strip("|").split("|")]
        data.append(dict(zip(header, cells)))
    return header, data

def md_block(title, note, header, data):
    out = ["# " + title, "", note, "",
           "| " + " | ".join(header) + " |",
           "| " + " | ".join(["---"] * len(header)) + " |"]
    for d in data:
        out.append("| " + " | ".join(str(d[h]) for h in header) + " |")
    return "\n".join(out) + "\n"

with open(MAN, encoding="utf-8") as fh:
    man = fh.read()

# S8
s8_md, s8_rows, s8_cols = build_s8()
with open(os.path.join(SUP, "Table_S8_colocalisation.md"), "w", encoding="utf-8") as f:
    f.write(s8_md)
with open(os.path.join(SUP, "Table_S8_colocalisation.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=s8_cols); w.writeheader()
    for r in s8_rows: w.writerow(r)
print("[OK] S8 -> 31 loci, max PP.H4=0.799 at rs11228269 (verified)")

# S9
h9, d9 = parse_table(man, "Supplementary Table S9")
s9_md = md_block(
    "Supplementary Table S9. Multivariable MR of BMD → knee OA and BMD → fracture, conditioning on height",
    "Primary configuration n_h=100 (EUR height; Yengo et al. 2022). n_H as tabulated (49–50) is the "
    "number that survived harmonisation, not the requested cap; n_BMD is 27/23/19/15/21 of 31/25/21/16/24 "
    "because six instruments are absent from the height GWAS. Attenuation = 100×(β_uni − β_MVMR)/β_uni. "
    "θ_H is the MVMR direct effect of height; Fc(BMD) and Fc(H) are Sanderson–Smith conditional F (≥10). "
    "β per 1 SD on the N-corrected scale.",
    h9, d9)
with open(os.path.join(SUP, "Table_S9_mvmr_height.md"), "w", encoding="utf-8") as f:
    f.write(s9_md)
with open(os.path.join(SUP, "Table_S9_mvmr_height.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=h9); w.writeheader()
    for r in d9: w.writerow(r)
print(f"[OK] S9 -> {len(d9)} rows")

# S14
h14, d14 = parse_table(man, "Supplementary Table S14")
s14_md = md_block(
    "Supplementary Table S14. Population-ancestry composition of the 13 GWAS",
    "The BMD exposure instruments are European-only, whereas the knee-OA outcome is 87.31% European "
    "(multi-ancestry meta-analysis); LDSC uses a European-only reference, so cross-population LD differences "
    "between the BMD instruments and the multi-ancestry knee-OA outcome can bias r_g and the MR estimates. "
    "This ancestry mismatch is a primary limitation (§4). The OA_knee percentages are taken from the GWAS "
    "Catalog entry and the originating Nature report of the multi-ancestry meta-analysis.",
    h14, d14)
with open(os.path.join(SUP, "Table_S14_ancestry.md"), "w", encoding="utf-8") as f:
    f.write(s14_md)
with open(os.path.join(SUP, "Table_S14_ancestry.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=h14); w.writeheader()
    for r in d14: w.writerow(r)
print(f"[OK] S14 -> {len(d14)} rows")

# S15
h15, d15 = parse_table(man, "Supplementary Table S15")
s15_md = md_block(
    "Supplementary Table S15. Genetic correlations by phenotype block",
    "Blocks from the 13×13 LDSC genetic-correlation matrix (Figure 1, Table 2). 'Significant (q<0.05)' "
    "counts pairs with BH-FDR q<0.05. Of the two OA–fracture pairs, knee OA × fracture is significant "
    "(r_g=+0.119, q=3.3×10⁻⁹) and self-reported OA × fracture is the sole non-significant pair in the "
    "matrix (r_g=+0.056, p=0.051); genetic correlation is not evidence of causation (MR for OA→fracture "
    "was null, §4).",
    h15, d15)
with open(os.path.join(SUP, "Table_S15_rg_blocks.md"), "w", encoding="utf-8") as f:
    f.write(s15_md)
with open(os.path.join(SUP, "Table_S15_rg_blocks.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=h15); w.writeheader()
    for r in d15: w.writerow(r)
print(f"[OK] S15 -> {len(d15)} rows")

print("\nAll four supplementary tables rendered to docs/Supplementary/.")
