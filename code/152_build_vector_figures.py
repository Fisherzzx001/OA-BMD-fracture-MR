#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Produce publication-grade VECTOR figure deliverables (PDF + TIFF) for the
17 supplementary figures (S1-S17), to satisfy the OAC/Elsevier vector-format
requirement. Output -> Manuscript_Ready/Figures/Supplementary/Vector/

Source mapping (consistent with the PNGs currently in Manuscript_Ready):
  S1-S8   : Figures_v17_Nature/SupplFig{1-8}.pdf  + .tif   (v17 originals)
  S9-S12  : Figures_v18_Suppl/SupplFig{9-12}.pdf  + .tif   (v18 originals)
  S13-S17 : Figures_v18_Suppl/SupplFig{13-17}.pdf (vector)
            TIFF for S13-S17 is NOT in the source tree -> rendered from the
            PDF at 600 DPI via PyMuPDF (high-res raster, journal-acceptable).

PDF is the true vector format and is provided for all 17 figures.
TIFF (raster, required by some submission systems) is provided for all 17:
12 from source files + 5 rendered from the vector PDFs.
"""
import os, shutil, importlib.util

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
V17 = os.path.join(ROOT, "Figures_v17_Nature")
V18 = os.path.join(ROOT, "Figures_v18_Suppl")
OUT = os.path.join(ROOT, "Manuscript_Ready", "Figures", "Supplementary", "Vector")
os.makedirs(OUT, exist_ok=True)

def have(p):
    return os.path.exists(p)

# descriptive base names (files carry subject suffixes)
NAMES = {
    1:  "SupplFig1_MR_scatter_remaining",
    2:  "SupplFig2_leave_one_out",
    3:  "SupplFig3_GEFOS_benchmark",
    4:  "SupplFig4_STROBE_flow",
    5:  "SupplFig5_OA_knee_N_sensitivity",
    6:  "SupplFig6_fracture_N_decision",
    7:  "SupplFig7_pleiotropy_diagnostics",
    8:  "SupplFig8_coloc_power",
    9:  "SupplFig9_instrument_flow",
    10: "SupplFig10_MVMR_decomposition",
    11: "SupplFig11_robustness_grid",
    12: "SupplFig12_literature_comparison",
    13: "SupplFig13_rg_forest",
    14: "SupplFig14_h2_bar",
    15: "SupplFig15_presso_funnel",
    16: "SupplFig16_mr_power",
    17: "SupplFig17_overlap_rhoGE",
}

pdf_n, tif_n, tif_gen = 0, 0, 0

# ---- S1-S8 : v17 pdf + tif
for i in range(1, 9):
    base = NAMES[i]
    sp = os.path.join(V17, base + ".pdf")
    st = os.path.join(V17, base + ".tif")
    assert have(sp) and have(st), f"missing v17 source for {base}"
    shutil.copy(sp, os.path.join(OUT, base + ".pdf"))
    shutil.copy(st, os.path.join(OUT, base + ".tif"))
    pdf_n += 1; tif_n += 1

# ---- S9-S12 : v18 pdf + tif
for i in range(9, 13):
    base = NAMES[i]
    sp = os.path.join(V18, base + ".pdf")
    st = os.path.join(V18, base + ".tif")
    assert have(sp) and have(st), f"missing v18 source for {base}"
    shutil.copy(sp, os.path.join(OUT, base + ".pdf"))
    shutil.copy(st, os.path.join(OUT, base + ".tif"))
    pdf_n += 1; tif_n += 1

# ---- S13-S17 : v18 pdf + rendered tif (600 DPI)
import pymupdf, io
from PIL import Image
DPI = 600
for i in range(13, 18):
    base = NAMES[i]
    sp = os.path.join(V18, base + ".pdf")
    assert have(sp), f"missing v18 pdf for {base}"
    shutil.copy(sp, os.path.join(OUT, base + ".pdf"))
    pdf_n += 1
    # render TIFF (LZW, lossless) from the vector PDF via PNG intermediate
    doc = pymupdf.open(sp)
    page = doc[0]
    pix = page.get_pixmap(dpi=DPI)
    png_bytes = pix.tobytes("png")
    doc.close()
    tif_path = os.path.join(OUT, base + ".tif")
    Image.open(io.BytesIO(png_bytes)).save(tif_path, "TIFF", compression="tiff_lzw")
    tif_gen += 1
    tif_n += 1
    print(f"  rendered {base}.tif @ {DPI} DPI -> {pix.width}x{pix.height}")

print(f"\n[OK] Vector figures written to: {OUT}")
print(f"  PDF  : {pdf_n} files (S1-S17, true vector)")
print(f"  TIFF : {tif_n} files ({tif_n - tif_gen} from source + {tif_gen} rendered from PDF)")
