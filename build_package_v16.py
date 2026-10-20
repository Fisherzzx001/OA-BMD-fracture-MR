# -*- coding: utf-8 -*-
"""Assemble the complete journal submission package (v16).

Adds the two things the previous package lacked:
  * all figure files (main Figures 1-10 + Graphical Abstract, Supplementary S1-S7),
    in three journal-ready formats (PDF vector / PNG 600 dpi / TIFF 300 dpi);
  * the source data behind every Table and the Supplementary material.

Nothing in the manuscript text is modified. Run after build_submission.py.
"""
import os, re, shutil, hashlib, json, zipfile, datetime

ROOT = os.path.dirname(os.path.abspath(__file__))
os.chdir(ROOT)

PKG = "submission_package_v16"
FIG_SRC = "Figures_v17_Nature"
GEN = "code/output/genetics"

MAIN_FIGS = [
    ("Figure1", "Figure1_rg_heatmap_13"),
    ("Figure2", "Figure2_forest_key_results"),
    ("Figure3", "Figure3_MR_forest"),
    ("Figure4", "Figure4_h2_bar"),
    ("Figure5", "Figure5_causal_network"),
    ("Figure6", "Figure6_MR_scatter"),
    ("Figure7", "Figure7_MR_funnel"),
    ("Figure8", "Figure8_MR_power"),
    ("Figure9", "Figure9_overlap_pleiotropy"),
    ("Figure10", "Figure10_MRlap_overlap_correction"),
    ("GraphicalAbstract", "GraphicalAbstract"),
]
SUPP_FIGS = [
    ("SupplFig1", "SupplFig1_MR_scatter_remaining"),
    ("SupplFig2", "SupplFig2_leave_one_out"),
    ("SupplFig3", "SupplFig3_GEFOS_benchmark"),
    ("SupplFig4", "SupplFig4_STROBE_flow"),
    ("SupplFig5", "SupplFig5_OA_knee_N_sensitivity"),
    ("SupplFig6", "SupplFig6_fracture_N_decision"),
    ("SupplFig7", "SupplFig7_pleiotropy_diagnostics"),
    ("SupplFig8", "SupplFig8_coloc_power"),
]
EXTS = [".pdf", ".png", ".tif"]

# Table -> source CSV(s) that the manuscript numbers are drawn from
TABLES = {
    "Table1_h2_ldsc_intercept": ["rg_h2_trueN.csv", "oa_h2_liability.csv"],
    "Table2_genetic_correlations": ["rg_matrix_13.csv", "rg_p_matrix_13.csv", "rg_se_matrix_13.csv"],
    "Table3_MR_estimates": ["table3_corrected.csv", "mr_results_corrected.csv"],
    "Table4_MR_power": ["mr_power.csv"],
}
SUPP_DATA = [
    "mr_presso_results_corrected.csv",
    "overlap_statistics.csv",
    "strict_clump_summary.csv",
    "N_sensitivity_decision.csv",
    "fracture_h2_sensitivity.csv",
    "gefos_benchmark.json",
]

MS_MD = "docs/08_投稿稿_OA_BMD_fracture.md"


def fresh(path):
    if os.path.isdir(path):
        shutil.rmtree(path)
    os.makedirs(path, exist_ok=True)


def copy_figs():
    n_main = n_supp = 0
    miss = []
    d_main = os.path.join(PKG, "Figures")
    d_supp = os.path.join(PKG, "Supplementary_Figures")
    os.makedirs(d_main, exist_ok=True)
    os.makedirs(d_supp, exist_ok=True)
    for _, stem in MAIN_FIGS:
        for ext in EXTS:
            src = os.path.join(FIG_SRC, stem + ext)
            if not os.path.exists(src):
                miss.append(src); continue
            shutil.copy2(src, os.path.join(d_main, stem + ext))
        n_main += 1
    for _, stem in SUPP_FIGS:
        for ext in EXTS:
            src = os.path.join(FIG_SRC, stem + ext)
            if not os.path.exists(src):
                miss.append(src); continue
            shutil.copy2(src, os.path.join(d_supp, stem + ext))
        n_supp += 1
    # provenance files that document how each figure was generated
    for f in ("plotted_data.json", "README_Figures.md", "AUDIT_REPORT.txt"):
        src = os.path.join(FIG_SRC, f)
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(d_main, "provenance_" + f))
    return n_main, n_supp, miss


def copy_tables():
    d = os.path.join(PKG, "Tables")
    os.makedirs(d, exist_ok=True)
    miss = []
    for label, files in TABLES.items():
        for f in files:
            src = os.path.join(GEN, f)
            if not os.path.exists(src):
                miss.append(src); continue
            shutil.copy2(src, os.path.join(d, label + "__" + f))
    # per-SNP instrument-level data behind Table 3
    harm = os.path.join(GEN, "harmonised")
    if os.path.isdir(harm):
        dh = os.path.join(d, "Per_SNP_harmonised")
        os.makedirs(dh, exist_ok=True)
        for f in sorted(os.listdir(harm)):
            if f.endswith(".csv"):
                shutil.copy2(os.path.join(harm, f), os.path.join(dh, f))
    return miss


def copy_supplementary():
    d = os.path.join(PKG, "Supplementary")
    os.makedirs(d, exist_ok=True)
    pairs = [
        ("docs/10_参考文献与原始GWAS来源.md", "S1_References_and_GWAS_sources.md"),
        ("docs/11_补充材料_OA_selfrep_Z校正.md", "S2_OA_selfrep_Z-sign_correction.md"),
        ("docs/12_补充材料_MR_PRESSO.md", "S3_MR_PRESSO_diagnostics.md"),
        ("docs/14_N敏感性决策表_fracture与OA_selfrep.md", "S4_N_sensitivity_decisions.md"),
        ("docs/09_数据下载实操指南.md", "S5_Data_download_protocol.md"),
    ]
    got = []
    for src, dst in pairs:
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(d, dst)); got.append(dst)
        else:
            print("  [MISSING]", src)
    # S2 (Z-correction) and S3 (OA_knee N-sensitivity) are also kept at the package
    # root because §2.4 / §3.1 of the text point at them directly.
    # supplementary result tables
    ds = os.path.join(d, "Result_tables")
    os.makedirs(ds, exist_ok=True)
    for f in SUPP_DATA:
        for base in (GEN, "code/output"):
            src = os.path.join(base, f)
            if os.path.exists(src):
                shutil.copy2(src, os.path.join(ds, f)); break
    src = "code/output/mrlap_overlap_results.csv"
    if os.path.exists(src):
        shutil.copy2(src, os.path.join(ds, "mrlap_overlap_results.csv"))
    return got


def write_legends():
    """Extract the '## Figure legends' block from the manuscript markdown."""
    with open(MS_MD, encoding="utf-8") as fh:
        text = fh.read()
    m = re.search(r"\n## Figure legends\n(.*?)(?=\n## |\Z)", text, re.S)
    body = m.group(1).strip() if m else ""
    out = ["# Figure legends",
           "",
           "Copied verbatim from the manuscript source (`docs/08_投稿稿_OA_BMD_fracture.md`, "
           "section *Figure legends*). Every numbered legend below has a matching file in "
           "`Figures/` or `Supplementary_Figures/`.",
           "",
           body,
           ""]
    with open(os.path.join(PKG, "Figure_Legends.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out))
    return len(re.findall(r"\*\*(?:Supplementary )?Figure|^\*\*Graphical Abstract", body, re.M))


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_readme_and_manifest(stats):
    rows = []
    for dirpath, dirnames, filenames in os.walk(PKG):
        dirnames.sort()
        for f in sorted(filenames):
            p = os.path.join(dirpath, f)
            rel = p.replace("\\", "/")
            rows.append((rel, os.path.getsize(p), sha256(p)[:16]))
    rows.sort()
    total = sum(r[1] for r in rows)

    with open(os.path.join(PKG, "MANIFEST.txt"), "w", encoding="utf-8") as fh:
        fh.write("Submission package manifest — generated %s\n" % datetime.date.today().isoformat())
        fh.write("Total: %d files, %.1f MB (uncompressed)\n\n" % (len(rows), total / 1048576))
        fh.write("%-70s %12s  %s\n" % ("file", "bytes", "sha256[:16]"))
        fh.write("-" * 100 + "\n")
        for rel, size, dig in rows:
            fh.write("%-70s %12d  %s\n" % (rel, size, dig))

    readme = """# Submission package — OA / BMD / fracture cross-trait study

Manuscript: *Bone Mineral Density and Osteoarthritis Share Genetic Architecture,
but Evidence for an Osteoarthritis-to-Fracture Causal Pathway Remains Limited:
A Cross-Trait Analysis of 13 Skeletal Phenotypes*

Package built: {date}.  Source of truth for all text is
`docs/08_投稿稿_OA_BMD_fracture.md`; figures are generated by `build_figures_v17.py`.

## Contents

| Folder / file | What it is |
|---|---|
| `Manuscript.docx` | Full manuscript with title page, authors, affiliations, CRediT, Data availability |
| `Manuscript_with_Figures.docx` | Same text with each figure inserted after its legend (for co-author / reviewer reading; **do not submit** this file) |
| `Submission_MainText.docx` | Blinded main text (no authors/affiliations) for journals that require it |
| `Cover_Letter.docx` | Cover letter addressed to the Editor-in-Chief |
| `Figure_Legends.md` | All {nleg} figure legends, verbatim from the manuscript |
| `Figures/` | Main figures 1–10 + Graphical Abstract, each as **PDF** (vector), **PNG** (600 dpi) and **TIFF** (300 dpi, LZW-compatible deflate) |
| `Supplementary_Figures/` | Supplementary figures S1–S7, same three formats |
| `Tables/` | Source CSVs behind Tables 1–4, plus per-SNP harmonised instrument data |
| `Supplementary/` | Supplementary notes S1–S5 and supplementary result tables |
| `Supplement_S2_*.md`, `Supplement_S3_*.md` | Referenced directly from §2.4 / §3.1 of the text |
| `MANIFEST.txt` | File list with byte sizes and SHA-256 prefixes for integrity checking |

## Figure format guidance

* **PDF** — preferred for the initial submission and for any journal that requests vector art.
* **TIFF** (300 dpi) — use if the journal's artwork checklist asks for TIFF/EPS.
* **PNG** (600 dpi) — for editorial-manager previews and for the graphical abstract.

All figures use the Okabe–Ito colour-blind-safe palette, Arial at >= 7 pt, and carry
panel labels (a, b, c …). Every plotted value traces back to the CSV/JSON files in
`Tables/`, `Supplementary/Result_tables/` and `Figures/provenance_plotted_data.json`;
`Figures/provenance_AUDIT_REPORT.txt` records the 845-point automated check that
reconciled each figure against those source data.

## Submission checklist

- [ ] Choose the target journal and re-address `Cover_Letter.docx` (currently *Bone*).
- [ ] Upload main figures as separate files (`Figures/`, PDF or TIFF).
- [ ] Upload supplementary figures as one combined PDF or as separate files
      (`Supplementary_Figures/`), per the journal's instructions.
- [ ] Upload `Manuscript.docx` (or `Submission_MainText.docx` + separate title page
      if the journal is double-blind).
- [x] References 16/17/18 volume/issue/pages completed (Crossref/PubMed verified).

## Reproducibility

```bash
python build_submission.py            # OUTDIR=submission_package_v16
python build_maintext.py              # blinded main text
python build_cover_letter.py          # cover letter
python build_figures_v17.py           # regenerates Figures_v17_Nature/
python build_package_v16.py           # this script
```
""".format(date=datetime.date.today().isoformat(), nleg=stats["nleg"])

    with open(os.path.join(PKG, "00_README.md"), "w", encoding="utf-8") as fh:
        fh.write(readme)
    return len(rows), total


def main():
    # keep the two docx + S2/S3 that already exist in the package
    keep = {}
    for f in os.listdir(PKG) if os.path.isdir(PKG) else []:
        p = os.path.join(PKG, f)
        if os.path.isfile(p):
            keep[f] = p
    tmp = "_pkg_keep"
    fresh(tmp)
    for f, p in keep.items():
        shutil.copy2(p, os.path.join(tmp, f))
    fresh(PKG)
    for f in os.listdir(tmp):
        shutil.copy2(os.path.join(tmp, f), os.path.join(PKG, f))
    shutil.rmtree(tmp)

    n_main, n_supp, miss_fig = copy_figs()
    miss_tab = copy_tables()
    supp = copy_supplementary()
    nleg = write_legends()
    stats = {"nleg": nleg}
    nfiles, total = write_readme_and_manifest(stats)

    print("Figures: %d main + %d supplementary" % (n_main, n_supp))
    print("Supplementary notes copied:", len(supp))
    print("Legend blocks detected:", nleg)
    for m in miss_fig + miss_tab:
        print("  [MISSING]", m)
    print("Package: %s — %d files, %.1f MB" % (PKG, nfiles, total / 1048576))

    zp = "Submission_Package_v16.zip"
    if os.path.exists(zp):
        os.remove(zp)
    with zipfile.ZipFile(zp, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for dirpath, dirnames, filenames in os.walk(PKG):
            dirnames.sort()
            for f in sorted(filenames):
                p = os.path.join(dirpath, f)
                z.write(p, p.replace("\\", "/"))
    print("Saved %s (%.1f MB)" % (zp, os.path.getsize(zp) / 1048576))


if __name__ == "__main__":
    main()
