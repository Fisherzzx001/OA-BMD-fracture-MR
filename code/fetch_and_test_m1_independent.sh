#!/usr/bin/env bash
# fetch_and_test_m1_independent.sh
# ==================================
# Turnkey: download GEFOS total-body BMD (GCST005348; the exact "independent BMD
# GWAS" the reviewer asked for -- trait-matched to the Qian total-body BMD exposure)
# and run the M1 variant-level scale test (code/351b).
#
# WHY A CUSTOM DOWNLOADER: the GWAS Catalog / EBI HTTPS proxy truncates large or
# irregular responses. This script requests small (256KB) byte-ranges; because the
# server honours the Range `start`, every received byte is a valid prefix of its
# range, so advancing the offset by the ACTUAL bytes received keeps the file
# byte-contiguous even when individual requests are truncated. Short/zero-byte
# responses are retried with backoff; the whole thing resumes from the existing
# file size if restarted.
#
# USAGE (on a machine with open network access):
#   bash code/fetch_and_test_m1_independent.sh
#   # then read code/output/genetics/m1_audit/m1_independent_gefs_tbmd_test.json
#
# ALTERNATIVE (IEU ieu-a-980 femoral-neck BMD, GEFOS 2018, N~32.7k -- also
# independent and consistent with the manuscript's existing GEFOS positive control):
#   export OPENGWAS_JWT="<your token from https://api.opengwas.io/>"
#   ieugwasr::associations("ieu-a-980")   # in R, or the REST equivalent with Bearer token
#   then run code/351 --independent <that file>   # expects SNP,A1,A2,beta,se
set -u
cd "$(dirname "$0")/.." || exit 1

URL="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST005001-GCST006000/GCST005348/harmonised/29304378-GCST005348-EFO_0003923-build37.f.tsv.gz"
OUT="code/output/genetics/independent_bmd/GCST005348_harmonised_build37.f.tsv.gz"
TMP="code/output/genetics/independent_bmd/_dl_state"
mkdir -p "$(dirname "$OUT")" "$TMP"

CHUNK=5242880
TOTAL=468025799
MAX_RETRY=20
SLEEP_BASE=2

# Resume from actual file size (authoritative after any crash).
if [ -f "$OUT" ]; then offset=$(wc -c < "$OUT"); else offset=0; rm -f "$OUT"; fi
rm -f /tmp/_chunk.gz
echo "RESUME offset=$offset / $TOTAL" | tee -a "$TMP/log"

while [ "$offset" -lt "$TOTAL" ]; do
  end=$(( offset + CHUNK - 1 ))
  [ "$end" -ge "$TOTAL" ] && end=$((TOTAL-1))
  retry=0; got=0
  while [ "$retry" -lt "$MAX_RETRY" ]; do
    got=$(curl -sL -m 120 -r ${offset}-${end} "$URL" -o /tmp/_chunk.gz -w "%{size_download}")
    [ "$got" -gt 0 ] && break
    retry=$((retry+1)); sleep $((SLEEP_BASE * retry))
  done
  if [ "$got" -le 0 ]; then
    echo "FAILED at offset=$offset after $MAX_RETRY retries" | tee -a "$TMP/log"; exit 4
  fi
  cat /tmp/_chunk.gz >> "$OUT"
  offset=$(( offset + got ))
  if [ $(( offset % 10485760 )) -lt $CHUNK ]; then
    echo "$(date +%H:%M:%S) offset=$offset ($(( offset * 100 / TOTAL ))%)" | tee -a "$TMP/log"
  fi
done

echo "DOWNLOAD COMPLETE: $OUT" | tee -a "$TMP/log"
gzip -t "$OUT" && echo "GZIP OK" || { echo "GZIP CORRUPT"; exit 3; }

echo "=== RUNNING M1 independent scale test (code/351b) ==="
python3 code/351b_independent_scale_test.py
