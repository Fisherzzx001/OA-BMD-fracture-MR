#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
199_download_gwas.py  --  robust parallel-Range downloader for GWAS Catalog .tsv.gz
Replicates the proven L7 method (parallel byte-range fetches from https://ftp.ebi.ac.uk,
which is NOT the throttled www.ebi.ac.uk host). Splits the file into N chunks, downloads
each range in parallel with retries, then concatenates in order. Output is a byte-exact
copy of the original gzip (range over a gzip just copies bytes).

Usage:
  python code/199_download_gwas.py <ACCESSION> <OUT_PATH> [N_CHUNKS]
Example:
  python code/199_download_gwas.py GCST90566801 G:/OA_B_PROJECT_cache/l7/GCST90566801.tsv.gz 24
"""
import os, sys, time, math, urllib.request, concurrent.futures

BASE = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90566001-GCST90567000"

def build_opener():
    op = urllib.request.build_opener(urllib.request.ProxyHandler({}))
    op.addheaders = [("User-Agent", "l8-dl/1.0")]
    return op

def get_size(op, url):
    # 0-0 range -> Content-Range: bytes 0-0/<TOTAL>
    for _ in range(5):
        try:
            req = urllib.request.Request(url, headers={"Range": "bytes=0-0"})
            with op.open(req, timeout=60) as r:
                cr = r.headers.get("Content-Range", "")
                if "/" in cr:
                    return int(cr.split("/")[-1])
        except Exception as e:
            time.sleep(3)
    raise RuntimeError("cannot determine size for " + url)

def fetch_chunk(op, url, s, e, out, retries=14):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Range": f"bytes={s}-{e}"})
            with op.open(req, timeout=180) as r, open(out, "wb") as f:
                while True:
                    buf = r.read(1 << 20)
                    if not buf:
                        break
                    f.write(buf)
            if os.path.getsize(out) == (e - s + 1):
                return True
        except Exception:
            time.sleep(2 + attempt * 2)
    return os.path.exists(out) and os.path.getsize(out) == (e - s + 1)

def main():
    if len(sys.argv) < 3:
        print("usage: 199_download_gwas.py <ACCESSION> <OUT_PATH> [N_CHUNKS]")
        sys.exit(2)
    acc = sys.argv[1]
    out = sys.argv[2]
    nch = int(sys.argv[3]) if len(sys.argv) > 3 else 24
    url = f"{BASE}/{acc}/{acc}.tsv.gz"
    tmpdir = out + ".chunks"
    os.makedirs(tmpdir, exist_ok=True)

    if os.path.exists(out):
        print(f"[skip] {out} already exists ({os.path.getsize(out)} bytes)")
        return

    op = build_opener()
    total = get_size(op, url)
    print(f"[size] {acc}: {total:,} bytes  ->  {out}", flush=True)

    # chunk boundaries (inclusive end)
    bounds = []
    step = total // nch
    for i in range(nch):
        s = i * step
        e = (s + step - 1) if i < nch - 1 else total - 1
        bounds.append((i, s, e))

    t0 = time.time()
    def work(b):
        i, s, e = b
        cpath = os.path.join(tmpdir, f"chunk_{i:03d}")
        ok = fetch_chunk(op, url, s, e, cpath)
        return i, ok, os.path.getsize(cpath) if os.path.exists(cpath) else 0

    results = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=nch) as ex:
        for i, ok, sz in ex.map(work, bounds):
            results[i] = (ok, sz)
            done = sum(1 for v in results.values() if v[0])
            print(f"  chunk {i:03d} ok={ok} {sz:,}B  ({done}/{nch}) {time.time()-t0:.0f}s",
                  flush=True)

    # retry any failed chunks singly
    failed = [b for b in bounds if not results[b[0]][0]]
    tries = 0
    while failed and tries < 6:
        tries += 1
        print(f"[retry round {tries}] {len(failed)} failed chunks", flush=True)
        nf = []
        for b in failed:
            i, s, e = b
            cpath = os.path.join(tmpdir, f"chunk_{i:03d}")
            if fetch_chunk(op, url, s, e, cpath):
                results[i] = (True, os.path.getsize(cpath))
            else:
                nf.append(b)
        failed = nf

    missing = [b for b in bounds if not results[b[0]][0]]
    if missing:
        print(f"[FAIL] {len(missing)} chunks still missing: {[b[0] for b in missing]}")
        sys.exit(1)

    # assemble in order
    assembled = 0
    with open(out, "wb") as fo:
        for i in range(nch):
            cpath = os.path.join(tmpdir, f"chunk_{i:03d}")
            with open(cpath, "rb") as fi:
                d = fi.read()
            fo.write(d)
            assembled += len(d)
    print(f"[assemble] {assembled:,} bytes  expected {total:,}  match={assembled==total}  "
          f"({time.time()-t0:.0f}s)", flush=True)
    if assembled != total:
        print("[FAIL] size mismatch after assembly")
        sys.exit(1)
    # cleanup chunk dir
    try:
        import shutil
        shutil.rmtree(tmpdir)
    except Exception:
        pass
    print(f"[done] {acc} -> {out}")

if __name__ == "__main__":
    main()
