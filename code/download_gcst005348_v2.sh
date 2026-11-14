#!/usr/bin/env bash
# Robust resumable downloader for GEFOS total-body BMD (GCST005348) harmonised stats.
# Problem: the EBI HTTPS proxy truncates large/irregular responses. Range requests
# honor `start`, so we advance `offset` by the ACTUAL bytes received each request,
# making the concatenation byte-contiguous even if individual requests are truncated.
set -u
URL="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST005001-GCST006000/GCST005348/harmonised/29304378-GCST005348-EFO_0003923-build37.f.tsv.gz"
OUT="code/output/genetics/independent_bmd/GCST005348_harmonised_build37.f.tsv.gz"
TMP="code/output/genetics/independent_bmd/_dl_state"
mkdir -p "$(dirname "$OUT")" "$TMP"

CHUNK=262144          # 256KB — well below the proxy's truncation threshold
MAX_TRUNC=44040191    # 42MB sanity cap (file is ~468MB; stop if something is wrong)
TOTAL=468025799
MAX_RETRY=12
SLEEP_BASE=2

# Resume from the actual size of OUT (standard resume logic: every successful
# append advances both the file and the offset, so the file is always authoritative).
if [ -f "$OUT" ]; then
  offset=$(wc -c < "$OUT")
else
  offset=0
  rm -f "$OUT"
fi
rm -f /tmp/_chunk.gz

echo "START offset=$offset TOTAL=$TOTAL" | tee -a "$TMP/log"

while [ "$offset" -lt "$TOTAL" ]; do
  end=$(( offset + CHUNK - 1 ))
  if [ "$end" -ge "$TOTAL" ]; then end=$((TOTAL-1)); fi
  retry=0
  got=0
  while [ "$retry" -lt "$MAX_RETRY" ]; do
    got=$(curl -sL -m 90 -r ${offset}-${end} "$URL" -o /tmp/_chunk.gz -w "%{size_download}")
    if [ "$got" -gt 0 ]; then break; fi
    retry=$((retry+1))
    sleep $((SLEEP_BASE * retry))
  done
  if [ "$got" -le 0 ]; then
    echo "FAILED at offset=$offset after $MAX_RETRY retries" | tee -a "$TMP/log"
    echo "$offset" > "$TMP/offset"
    exit 4
  fi
  # append
  cat /tmp/_chunk.gz >> "$OUT"
  offset=$(( offset + got ))
  echo "$offset" > "$TMP/offset"
  # light logging every ~10MB
  if [ $(( offset % 10485760 )) -lt $CHUNK ]; then
    echo "$(date +%H:%M:%S) offset=$offset ($(( offset * 100 / TOTAL ))%)" | tee -a "$TMP/log"
  fi
done

echo "DONE offset=$offset" | tee -a "$TMP/log"
gzip -t "$OUT" && echo "GZIP OK -> $OUT" || echo "GZIP CORRUPT"
rm -f "$TMP/offset"
