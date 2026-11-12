# -*- coding: utf-8 -*-
"""Assemble a clean submission package zip from the current post-remediation artifacts.

Mirrors build_package_v16.py but targets the CURRENT 5-main-figure manuscript
(OAC display cap = 8) and bundles only artifacts that actually exist.

2026-11-12  (docs/148 P0-2 / P0-3 data-governance remediation)
---------------------------------------------------------------
P0-2: the repository previously shipped THREE mutually inconsistent
      mr_results*.csv (raw placeholder-N, corrected HEAD, pre-G6 backup).
      Only mr_results_corrected.csv is authoritative (Table 2 truth source;
      read by code/173,172,177,181,171c,178,14). We now ship ONLY the
      authoritative table and its matched PRESSO corrected table; the raw
      placeholder-N table and the _backup_pre_G6/ tree stay on disk but are
      EXCLUDED from the deliverable (documented in RESULTS_AUTHORITATIVE.md).
P0-3: the per-SNP harmonised instrument tables, from which Table 2's
      IVW/Egger/Q/I2/WM/LOO are derived (code/177, code/334), are NOW
      DELIVERED under Table_Source_Data/harmonised/ (STROBE-MR requirement).
"""
import os, zipfile, datetime, hashlib

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PKG = "Submission_Package_2026-10-17.zip"
OUT = os.path.join(ROOT, PKG)

MR = os.path.join(ROOT, "Manuscript_Ready")
SP = os.path.join(ROOT, "submission_package")
GEN = os.path.join(ROOT, "code", "output", "genetics")
HARM = os.path.join(GEN, "harmonised")

# P0-2: authoritative result tables ONLY (raw mr_results.csv and _backup_pre_G6
# deliberately excluded — see RESULTS_AUTHORITATIVE.md).
TABLE_CSVS = [
    "mr_results_corrected.csv", "mr_presso_results_corrected.csv",
    "rg_matrix_13.csv", "rg_h2_trueN.csv", "oa_h2_liability.csv",
    "rg_p_matrix_13.csv", "rg_se_matrix_13.csv", "table3_corrected.csv",
    "mr_power.csv", "overlap_statistics.csv", "strict_clump_summary.csv",
    "N_sensitivity_decision.csv", "fracture_h2_sensitivity.csv",
    "mvmr_height_main.csv", "mvmr_unify_baselines.csv",
    "coloc_screen_all_loci.csv", "plink_benchmark_summary.csv",
    "rg_blocks_summary.csv", "coloc_power.csv",
]

# Governance / provenance docs (independent OSF preregistration + K-a verification)
GOVERNANCE_DOCS = [
    ("docs/118_OSF预注册_OA_BMD_fracture.md", "Governance/OSF_preregistration_OA_BMD_fracture.md"),
    ("docs/119_K-a核验.md", "Governance/K-a_verification.md"),
    ("docs/RESULTS_AUTHORITATIVE.md", "Table_Source_Data/RESULTS_AUTHORITATIVE.md"),
]

# P0-3: per-SNP instrument tables backing every Table 2 cell (STROBE-MR).
HARMONISED_CSVS = [
    "BMD_total__fracture.csv", "BMD_femur_neck__fracture.csv",
    "BMD_total__OA_knee.csv", "BMD_lumbar__OA_knee.csv", "BMD_arms__OA_knee.csv",
    "BMD_leg__OA_knee.csv", "BMD_femur_neck__OA_knee.csv",
    "OA_knee__BMD_total.csv", "OA_knee__BMD_femur_neck.csv", "OA_knee__fracture.csv",
    "_loo.csv", "_validation.csv",
]


def sha12(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:12]


def add_tree(zf, real_dir, arc_root):
    n = 0
    for dp, _, fns in os.walk(real_dir):
        for fn in sorted(fns):
            rp = os.path.join(dp, fn)
            arc = os.path.join(arc_root, os.path.relpath(rp, real_dir)).replace("\\", "/")
            zf.write(rp, arc); n += 1
    return n


manifest_lines = []
with zipfile.ZipFile(OUT, "w", zipfile.ZIP_DEFLATED) as zf:
    n = 0
    # docx deliverables
    for fn in ("Manuscript.docx", "Cover_Letter.docx"):
        zf.write(os.path.join(SP, fn), fn); n += 1
    # manuscript_Ready mirror (text + figures + supplementary)
    n += add_tree(zf, os.path.join(MR, "Figures", "Main"), "Figures/Main")
    n += add_tree(zf, os.path.join(MR, "Figures", "Supplementary"), "Figures/Supplementary")
    n += add_tree(zf, os.path.join(MR, "Supplementary"), "Supplementary")
    for fn in ("01_Manuscript.md", "02_Cover_letter.md", "README.md"):
        zf.write(os.path.join(MR, fn), fn); n += 1
    # governance / provenance docs
    g = 0
    for real, arc in GOVERNANCE_DOCS:
        if os.path.exists(real):
            zf.write(real, arc); g += 1; n += 1
    # table source data — authoritative tables only (P0-2)
    m = 0
    for csv in TABLE_CSVS:
        rp = os.path.join(GEN, csv)
        if os.path.exists(rp):
            zf.write(rp, os.path.join("Table_Source_Data", csv)); m += 1
            manifest_lines.append(f"  {sha12(rp)}  Table_Source_Data/{csv}")
    # per-SNP harmonised instrument tables (P0-3)
    h = 0
    for csv in HARMONISED_CSVS:
        rp = os.path.join(HARM, csv)
        if os.path.exists(rp):
            zf.write(rp, os.path.join("Table_Source_Data/harmonised", csv)); h += 1
            manifest_lines.append(f"  {sha12(rp)}  Table_Source_Data/harmonised/{csv}")
    manifest = (
        f"Submission_Package_2026-10-17\n"
        f"Built: {datetime.datetime.now().isoformat(timespec='seconds')}\n"
        f"Manuscript main body words: 3963 (OAC cap 4000; code/148 Word-compatible count)\n"
        f"Abstract words: 339 (OAC cap 350, structured)\n"
        f"References: 50 (OAC cap 50)\n"
        f"Display items: Figure 1-5 + Table 1-3 = 8 (OAC cap 8)\n"
        f"--- data governance (docs/148 P0-2/P0-3) ---\n"
        f"AUTHORITATIVE MR table: mr_results_corrected.csv (Table 2 truth source)\n"
        f"AUTHORITATIVE PRESSO table: mr_presso_results_corrected.csv\n"
        f"EXCLUDED from deposit: mr_results.csv (raw placeholder-N), _backup_pre_G6/*\n"
        f"Per-SNP harmonised instrument tables delivered: {h} (Table 2 reproduction)\n"
        f"Reproduce Table 2: python3 code/334_data_level_trust_but_verify.py (18/18)\n"
        f"--- SHA256 (first 12) ---\n" + "\n".join(manifest_lines) + "\n"
        f"--- counts ---\n"
        f"Table source CSVs included: {m}\n"
        f"Harmonised per-SNP CSVs included: {h}\n"
        f"Governance docs included: {g}\n"
        f"Total zip entries: {n+m+h+1}\n"
    )
    zf.writestr("MANIFEST.txt", manifest)

print(f"[OK] wrote {PKG}")
print(manifest)
