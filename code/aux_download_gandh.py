# -*- coding: utf-8 -*-
"""Resumable chunked downloader for GWAS Catalog harmonised sumstats.

Proxy throttles ftp.ebi.ac.uk to ~45 KB/s and kills long connections, so a
single streaming GET stalls. This downloads in small HTTP-Range chunks with
resume + retry, writing to a local cache file. Run in background.
"""
import os
import sys
import time
import socket
import urllib.request

URL = ("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
       "GCST90727001-GCST90728000/GCST90727231/harmonised/GCST90727231.h.tsv.gz")
DST = "G:/OA_B_PROJECT_cache/l7/GCST90727231.h.tsv.gz"
CHUNK = 2 * 1048576  # 2 MiB chunks
SOCK_TIMEOUT = 60
MAX_RETRY = 200


def main():
    socket.setdefaulttimeout(SOCK_TIMEOUT)
    os.makedirs(os.path.dirname(DST), exist_ok=True)
    # get total size
    total = None
    try:
        req = urllib.request.Request(URL, method="HEAD",
                                    headers={"User-Agent": "OA-BMD/f10"})
        with urllib.request.urlopen(req, timeout=SOCK_TIMEOUT) as r:
            total = int(r.headers.get("Content-Length"))
    except Exception as e:
        print("HEAD failed:", e, flush=True)
    if total is None:
        print("Cannot determine total size; abort.", flush=True)
        return
    print("total size:", total, flush=True)
    t0 = time.time()
    pos = os.path.getsize(DST) if os.path.exists(DST) else 0
    print("resuming from", pos, flush=True)
    retries = 0
    with open(DST, "ab" if pos else "wb") as f:
        if pos:
            f.seek(pos)
        while pos < total:
            end = min(pos + CHUNK - 1, total - 1)
            hdr = {"User-Agent": "OA-BMD/f10", "Range": f"bytes={pos}-{end}"}
            try:
                req = urllib.request.Request(URL, headers=hdr)
                r = urllib.request.urlopen(req, timeout=SOCK_TIMEOUT)
                data = r.read()
                if not data:
                    raise IOError("empty chunk")
                f.write(data)
                pos += len(data)
                retries = 0
                elapsed = time.time() - t0
                rate = (pos / 1048576) / elapsed if elapsed else 0
                print(f"  {pos/total*100:5.1f}%  {pos/1048576:.1f}MiB  "
                      f"{rate:.1f}MiB/s  t={elapsed:.0f}s", flush=True)
            except Exception as e:
                retries += 1
                if retries > MAX_RETRY:
                    print("TOO MANY RETRIES; abort at", pos, flush=True)
                    return
                print(f"  retry {retries} at {pos}: {type(e).__name__} "
                      f"{str(e)[:100]}", flush=True)
                time.sleep(2)
    print("DONE size", os.path.getsize(DST), flush=True)


if __name__ == "__main__":
    main()
