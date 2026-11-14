#!/usr/bin/env bash
# Parallel chunked download of GEFOS total-body BMD (GCST005348) harmonised summary stats.
# Rationale: EBI proxy throttles ~13KB/s per connection; parallelism multiplies throughput.
set -u
URL="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST005001-GCST006000/GCST005348/harmonised/29304378-GCST005348-EFO_0003923-build37.f.tsv.gz"
OUT="code/output/genetics/independent_bmd/GCST005348_harmonised_build37.f.tsv.gz"
TMP="code/output/genetics/independent_bmd/_chunks"
mkdir -p "$TMP"
TOTAL=468025799
N=30
CHUNK=$(( (TOTAL + N - 1) / N ))
echo "TOTAL=$TOTAL N=$N CHUNK=$CHUNK"

cleanup(){ rm -f "$TMP"/chunk_*.gz; }
trap cleanup EXIT

pids=()
for i in $(seq 0 $((N-1))); do
  start=$(( i * CHUNK ))
  end=$(( start + CHUNK - 1 ))
  if [ $end -ge $TOTAL ]; then end=$((TOTAL-1)); fi
  (
    curl -sL -m 600 --retry 6 --retry-delay 3 -r ${start}-${end} "$URL" -o "$TMP/chunk_$(printf '%02d' $i).gz"
  ) &
  pids+=($!)
done
echo "launched ${#pids[@]} parallel chunk downloads"
for p in "${pids[@]}"; do wait $p; done
echo "all chunks finished"

# verify and concatenate in order
cat "$TMP"/chunk_*.gz > "$OUT"
sz=$(wc -c < "$OUT")
echo "concatenated size=$sz (expected $TOTAL)"
if [ "$sz" -ne "$TOTAL" ]; then echo "SIZE MISMATCH"; exit 2; fi
gzip -t "$OUT" && echo "GZIP OK" || { echo "GZIP CORRUPT"; exit 3; }
echo "DOWNLOAD COMPLETE: $OUT"
