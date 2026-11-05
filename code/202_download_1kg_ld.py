#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
202_download_1kg_ld.py  --  parallel-Range resume downloader for the MRC-IEU
1000G-EUR LD reference panel (1kg.v3.tgz, ~1.495 GB), stored at G:/REF/1000G_EUR.

Why this script exists:
  A single curl from http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz was found to be
  throttled to ~8.4 KB/s and then went stale (connection dropped, process alive
  but no socket, file frozen at 291.8 MB / 19.5%). fileserve honours Range
  requests (206 Partial Content), so we resume from the bytes already on disk and
  pull the REMAINING ~80% with N parallel connections (per-connection bandwidth
  shaping on these hosts -> parallelism is what buys throughput).

Behaviour:
  - Determines total size via a 0-0 Range HEAD.
  - If <out> already exists and is smaller than total, it is renamed to <out>.base
    (the bytes already downloaded are kept), and only [have, total) is split into
    N parallel chunks.
  - Assembles base + ordered chunks, verifies byte-exact size, then cleans up.
  - If <out> is missing/empty, downloads the whole file in parallel.

Usage:
  python code/202_download_1kg_ld.py [OUT_PATH] [N_CHUNKS]
Default OUT: G:/REF/1000G_EUR/1kg.v3.tgz ; N_CHUNKS: 16
"""
import os, sys, time, shutil, concurrent.futures, urllib.request

URL = "http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz"

def build_opener():
    op = urllib.request.build_opener(urllib.request.ProxyHandler({}))
    op.addheaders = [("User-Agent", "l8-dl/1.0")]
    return op

def get_size(op):
    for _ in range(6):
        try:
            req = urllib.request.Request(URL, headers={"Range": "bytes=0-0"})
            with op.open(req, timeout=60) as r:
                cr = r.headers.get("Content-Range", "")
                if "/" in cr:
                    return int(cr.split("/")[-1])
        except Exception:
            time.sleep(3)
    raise RuntimeError("cannot determine size for " + URL)

def fetch_chunk(op, s, e, out, retries=16):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(URL, headers={"Range": f"bytes={s}-{e}"})
            with op.open(req, timeout=300) as r, open(out, "wb") as f:
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
    out = sys.argv[1] if len(sys.argv) > 1 else "G:/REF/1000G_EUR/1kg.v3.tgz"
    nch = int(sys.argv[2]) if len(sys.argv) > 2 else 16
    op = build_opener()
    total = get_size(op)
    have = os.path.getsize(out) if os.path.exists(out) else 0
    print(f"[size] total={total:,} have={have:,} need={total-have:,} "
          f"({100.0*have/total:.1f}% already on disk)", flush=True)
    if have >= total:
        print("[skip] already complete"); return

    base_path = out + ".base"
    if have > 0 and not os.path.exists(base_path):
        os.rename(out, base_path)          # keep already-downloaded prefix
        print(f"[resume] keeping {have:,} bytes as base", flush=True)
    elif have == 0 and os.path.exists(out):
        os.remove(out)

    start = have
    tmpdir = out + ".chunks"
    os.makedirs(tmpdir, exist_ok=True)

    rem = total - start
    step = rem // nch
    bounds = []
    for i in range(nch):
        s = start + i * step
        e = (start + (i + 1) * step - 1) if i < nch - 1 else total - 1
        bounds.append((i, s, e))

    t0 = time.time()
    def work(b):
        i, s, e = b
        cpath = os.path.join(tmpdir, f"chunk_{i:03d}")
        ok = fetch_chunk(op, s, e, cpath)
        return i, ok, os.path.getsize(cpath) if os.path.exists(cpath) else 0

    results = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=nch) as ex:
        for i, ok, sz in ex.map(work, bounds):
            results[i] = (ok, sz)
            done = sum(1 for v in results.values() if v[0])
            print(f"  chunk {i:03d} ok={ok} {sz:,}B ({done}/{nch}) "
                  f"{time.time()-t0:.0f}s", flush=True)

    failed = [b for b in bounds if not results[b[0]][0]]
    tries = 0
    while failed and tries < 8:
        tries += 1
        print(f"[retry round {tries}] {len(failed)} failed chunks", flush=True)
        nf = []
        for b in failed:
            i, s, e = b
            cpath = os.path.join(tmpdir, f"chunk_{i:03d}")
            if fetch_chunk(op, s, e, cpath):
                results[i] = (True, os.path.getsize(cpath))
            else:
                nf.append(b)
        failed = nf

    missing = [b for b in bounds if not results[b[0]][0]]
    if missing:
        print(f"[FAIL] {len(missing)} chunks still missing: {[b[0] for b in missing]}")
        sys.exit(1)

    assembled = 0
    with open(out, "wb") as fo:
        if os.path.exists(base_path):
            with open(base_path, "rb") as fb:
                d = fb.read()
            fo.write(d)
            assembled += len(d)
        for i in range(nch):
            cpath = os.path.join(tmpdir, f"chunk_{i:03d}")
            with open(cpath, "rb") as fi:
                d = fi.read()
            fo.write(d)
            assembled += len(d)
    print(f"[assemble] {assembled:,} bytes  expected {total:,}  "
          f"match={assembled==total}  ({time.time()-t0:.0f}s)", flush=True)
    if assembled != total:
        print("[FAIL] size mismatch after assembly")
        sys.exit(1)
    try:
        if os.path.exists(base_path):
            os.remove(base_path)
        shutil.rmtree(tmpdir)
    except Exception:
        pass
    print(f"[done] -> {out}  ({total:,} bytes)")

if __name__ == "__main__":
    main()
