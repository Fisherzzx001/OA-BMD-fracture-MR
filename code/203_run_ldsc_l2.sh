#!/bin/bash
# Compute 1000G phase-3 EUR baseline LD scores (.l2.ldscore.gz) for downstream
# cross-trait LDSC / MVMR. ldsc 2.0.1 (PyPI wheel, py3) installed --no-deps;
# runtime deps numpy/scipy/pandas/bitarray already in the custom python 3.13.
# NOTE: EUR.bim ships with CM (genetic distance) = 0, so we use --ld-wind-kb 1000
# (1 Mb physical window) instead of --ld-wind-cm; this also avoids the
# whole-chromosome LD Score mode (which OOMs on this box for 8.5M SNPs).
set -u
LDSC="C:/Users/Administrator/.workbuddy/binaries/python/versions/3.13.12/Scripts/ldsc.py"
OUTDIR="/g/REF/1000G_EUR/EUR_ldscore"
mkdir -p "$OUTDIR"
python3 "$LDSC" \
  --l2 \
  --bfile /g/REF/1000G_EUR/EUR \
  --ld-wind-kb 1000 \
  --out "$OUTDIR/EUR" \
  2>&1 | tee "$OUTDIR/ldsc_l2.log"
echo "EXIT_LDSC=${PIPESTATUS[0]}"
