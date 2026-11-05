#!/bin/bash
# One-command entry for the downstream munge + cross-trait LDSC rg pipeline.
# Default weights = standard Broad HapMap3 (--ref-ld-chr/--w-ld-chr eur_w_ld_chr_
# + munge --merge-alleles w_hm3.snplist), i.e. publication-grade configuration.
# Usage:
#   bash code/205_run_ldsc_munge_rg.sh --check      # verify prerequisites
#   bash code/205_run_ldsc_munge_rg.sh --smoke      # 200k-row end-to-end test
#   bash code/205_run_ldsc_munge_rg.sh              # full 4-site run (standard)
#   bash code/205_run_ldsc_munge_rg.sh --weights self   # exploratory contrast
# Everything else is passed straight through to the Python driver:
#   --sites knee spine   --no-h2   --skip-prep   --limit N   --fallback-M (self only)
set -u
cd "C:/Users/Administrator/WorkBuddy/2026-09-08-22-17-40/OA_B_PROJECT" || exit 2
exec python3 -u code/205_run_ldsc_munge_rg.py "$@"
