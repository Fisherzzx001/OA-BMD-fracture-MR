# -*- coding: utf-8 -*-
"""
F10 跨祖先敏感性 · MVP 非洲 OA（GCST90134287, n=35753）可行性核对

*** 结论（2026-11-10）：UNDERPOWERED / UNINFORMATIVE — 不作为结果报告 ***
------------------------------------------------------------------
实测 31 个欧洲 BMD 工具在该非洲祖先 OA GWAS 中仅 4–5 个可匹配；扩大
tabix 窗口到 ±512 kb 仍为 4 个，且 IVW 符号在窗口半径间翻转
（expand=0: β=−0.1205, P=0.374；expand=32: β=+0.0407, P=0.766）。
这既不是 build 偏移（坐标窗口已对上）也不是 rsID 更名（位置回退 0 救回），
而是该结局 GWAS 真实缺失这些常见欧洲变异。4–5 个工具的估计不可判读，
**不能**用来宣称方向稳健或反转，也**不能**关闭登记册 F10 或改动 S23 9/10。
本脚本与输出 JSON 仅作为「为什么 F10 仍不可执行」的阴性证据留痕，
不得写入稿件正文、不得作为跨祖先复制结果引用。

背景与诚实边界
--------------
登记册 F10 点名的是 Hatzikotoulas 2025 的 EAS 分层膝 OA，但该文件从未单独
发布（合并文件 GCST90566800 由 1.12M 欧裔主导）；纯 EAS OA（Liu TY 2025,
GCST90652199, n=221585）未发布 sumstats；MVP 的 EAS 分层仅 n=1704。
因此本脚本执行的是**可防御的跨祖先敏感性**（欧洲 BMD 工具 → 非欧 OA 结局），
**不是**登记册 F10 的 de novo EAS 复制。它检验「BMD→OA 方向在非欧人群是否反转」，
不能据此把 S23 的 9/10 改成 10/10。

为什么需要 tabix/BGZF 定点抓取
------------------------------
MVP 非洲 OA 的 harmonised 文件为 828 MB，而沙箱代理对 ftp.ebi.ac.uk 限速约
45 KB/s 并杀长连接，全量下载不可行。GWAS Catalog 同时提供 1.6 MB 的 .tbi
索引。本脚本：
  1. 解析 TBI（gzip 包着 TBI\\1），按 reg2bins 找 31 个工具位点的候选 bin；
  2. 把 chunk 虚拟偏移 (coffset<<16|uoffset) 转成 BGZF 块的实际字节区间；
  3. 合并区间后用 HTTP Range 只抓这些块（合计几 MB）；
  4. 逐块 gzip.decompress，按 rsID 精确匹配（build 无关）；
  5. 等位对齐后复用 code/13 的 IVW / Egger / LOO。

输入：G:/OA_B_PROJECT_cache/l7/GCST90134287.h.tsv.gz.tbi（1.6 MB，已下载）
      bmd_instruments.json（31 个锁定工具）
输出：code/output/genetics/L7_F10_cross_ancestry_MVP_AFR.json
"""
import os
import json
import gzip
import math
import struct
import socket
import urllib.request
import importlib.util
import numpy as np
from scipy.stats import norm

ROOT = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(ROOT)
CACHE = "G:/OA_B_PROJECT_cache/l7"
INST = os.path.join(CACHE, "bmd_instruments.json")
TBI = os.path.join(CACHE, "GCST90134287.h.tsv.gz.tbi")
BGZF_URL = ("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
            "GCST90134001-GCST90135000/GCST90134287/harmonised/"
            "GCST90134287.h.tsv.gz")
OUT = os.path.join(ROOT, "output", "genetics", "L7_F10_cross_ancestry_MVP_AFR.json")
socket.setdefaulttimeout(60)

_spec = importlib.util.spec_from_file_location(
    "mr13", os.path.join(ROOT, "13_mr_analysis.py"))
mr13 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mr13)


# ---------------- TBI ----------------
def reg2bins(beg, end):
    """tabix/htslib reg2bins, 0-based half-open."""
    end -= 1
    out = [0]
    for shift, base in ((26, 1), (23, 9), (20, 73), (17, 585), (14, 4681)):
        lo = base + (beg >> shift)
        hi = base + (end >> shift)
        out.extend(range(lo, hi + 1))
    return out


def parse_tbi(path):
    d = gzip.decompress(open(path, "rb").read())
    assert d[:4] == b"TBI\x01", "not a TBI index"
    off = 4
    n_ref, = struct.unpack_from("<i", d, off); off += 4
    fmt, = struct.unpack_from("<i", d, off); off += 4
    col_seq, col_beg, col_end, meta, skip, l_nm = struct.unpack_from(
        "<iiiiii", d, off); off += 24
    names = [x.decode() for x in d[off:off + l_nm].split(b"\x00") if x]
    off += l_nm
    refs = {}
    p = off
    for name in names:
        n_bin, = struct.unpack_from("<i", d, p); p += 4
        bins = {}
        for _ in range(n_bin):
            b, = struct.unpack_from("<I", d, p); p += 4
            nc, = struct.unpack_from("<i", d, p); p += 4
            chunks = []
            for _ in range(nc):
                a, z = struct.unpack_from("<QQ", d, p); p += 16
                chunks.append((a, z))
            bins[b] = chunks
        n_intv, = struct.unpack_from("<i", d, p); p += 4
        p += 8 * n_intv  # linear index, not needed here
        refs[name] = bins
    return refs, dict(n_ref=n_ref, format=fmt, col_seq=col_seq,
                      col_beg=col_beg, col_end=col_end, meta=meta, skip=skip,
                      names=names)


def chunks_for_loci(refs, loci, expand=0):
    """loci: list of (chrom str, bp int 1-based). -> list of (vo_beg, vo_end).

    Only the finest (16-kb) bins are used. Coarse bins (bin 0, level 1-4)
    carry chromosome-wide chunks spanning tens of MB, which would defeat the
    point of a targeted fetch under a throttled proxy. `expand` adds the
    neighbouring 16-kb windows on each side as a guard against a GRCh37/GRCh38
    coordinate offset between our instrument list and the harmonised outcome.
    """
    out = []
    missing = []
    for c, bp in loci:
        bins = refs.get(c) or refs.get("chr" + c)
        if not bins:
            missing.append((c, bp, "contig absent"))
            continue
        centre = (bp - 1) >> 14
        got = []
        for k in range(centre - expand, centre + expand + 1):
            got.extend(bins.get(4681 + k, []))
        if not got:
            # fall back to any *small* chunk among the coarse candidate bins
            for b in reg2bins(bp - 1, bp):
                for a, z in bins.get(b, []):
                    if (z >> 16) - (a >> 16) < (1 << 20):
                        got.append((a, z))
        if not got:
            missing.append((c, bp, "no chunk"))
        out.extend(got)
    if missing:
        print(f"  [warn] {len(missing)} loci without a fine-bin chunk: {missing[:5]}")
    return sorted(set(out))


# ---------------- BGZF ----------------
def bgzf_blocks(buf):
    """Yield decompressed payload of each BGZF member in buf."""
    p = 0
    while p + 18 <= len(buf):
        if buf[p:p + 4] != b"\x1f\x8b\x08\x04":
            break
        xlen, = struct.unpack_from("<H", buf, p + 10)
        q = p + 12
        bsize = None
        end_extra = p + 12 + xlen
        while q + 4 <= end_extra:
            si1, si2 = buf[q], buf[q + 1]
            slen, = struct.unpack_from("<H", buf, q + 2)
            if si1 == 66 and si2 == 67 and slen >= 2:
                bsize, = struct.unpack_from("<H", buf, q + 4)
            q += 4 + slen
        if bsize is None:
            break
        bend = p + bsize + 1
        if bend > len(buf):
            break
        try:
            yield gzip.decompress(buf[p:bend])
        except Exception:
            break
        p = bend


def http_range(url, start, end):
    hdr = {"User-Agent": "OA-BMD/f10", "Range": f"bytes={start}-{end}"}
    req = urllib.request.Request(url, headers=hdr)
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read()


def merge_intervals(iv, gap=1 << 20):
    """Merge (start,end) byte intervals separated by < gap."""
    if not iv:
        return []
    iv = sorted(iv)
    out = [list(iv[0])]
    for s, e in iv[1:]:
        if s <= out[-1][1] + gap:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def fetch_header():
    """Read the real header from the start of the BGZF file (one small range)."""
    buf = http_range(BGZF_URL, 0, 65536)
    text = b"".join(bgzf_blocks(buf)).decode("utf-8", "replace")
    return text.split("\n")[0].rstrip("\r").split("\t")


def fetch_rows(chunks, header, max_chunks=400):
    """Fetch BGZF blocks for chunks; return list of parsed dict rows.

    A targeted range usually starts *inside* a BGZF member (uoffset != 0), so
    its first line is a fragment. Rows whose column count differs from the real
    header are discarded; the real header is fetched separately from byte 0.
    """
    iv = []
    for vo_beg, vo_end in chunks[:max_chunks]:
        cs, ce = vo_beg >> 16, vo_end >> 16
        iv.append((cs, ce + 65536))
    merged = merge_intervals(iv)
    print(f"  {len(chunks)} chunks -> {len(merged)} merged byte ranges")
    rows = []
    n_bad = 0
    for i, (s, e) in enumerate(merged, 1):
        try:
            buf = http_range(BGZF_URL, s, e)
        except Exception as ex:
            print(f"  [warn] range {i} {s}-{e} failed: {type(ex).__name__}")
            continue
        text = b"".join(bgzf_blocks(buf)).decode("utf-8", "replace")
        for line in text.split("\n"):
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\r").split("\t")
            if len(parts) != len(header):
                n_bad += 1
                continue
            if parts[0] == "chromosome":  # header repeated mid-stream
                continue
            rows.append(dict(zip(header, parts)))
        print(f"  range {i}/{len(merged)} bytes={len(buf)} "
              f"decoded_rows={len(rows)}", flush=True)
    print(f"  discarded {n_bad} partial/fragment lines")
    return rows


def main():
    import sys
    expand = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    inst = json.load(open(INST, encoding="utf-8"))["inst"]
    print(f"[load] {len(inst)} locked BMD instruments (expand={expand} "
          f"x 16-kb windows each side)")
    refs, meta = parse_tbi(TBI)
    print(f"[tbi] n_ref={meta['n_ref']} names={meta['names'][:5]}... "
          f"col_beg={meta['col_beg']} skip={meta['skip']}")
    loci = [(str(s["chr"]), int(s["bp"])) for s in inst.values()]
    chunks = chunks_for_loci(refs, loci, expand=expand)
    print(f"[index] {len(chunks)} BGZF chunks cover the 31 instrument loci")
    header = fetch_header()
    print(f"[header] {len(header)} columns: {header}")
    rows = fetch_rows(chunks, header)
    print(f"[fetch] {len(rows)} rows decoded")

    wanted = set(inst.keys())
    matched = {}
    how = {}
    for r in rows:
        for key in ("rsid", "variant_id"):
            rs = r.get(key)
            if rs in wanted and rs not in matched:
                matched[rs] = r
                how[rs] = key
                break
    print(f"[match] {len(matched)}/{len(inst)} rsids matched by tabix window "
          f"(via rsid={sum(1 for v in how.values() if v=='rsid')}, "
          f"variant_id={sum(1 for v in how.values() if v=='variant_id')})")

    # ---- harmonise & MR ----
    recs = []
    for rs, s in inst.items():
        r = matched.get(rs)
        if not r:
            continue
        exp = str(s["exp_allele"]).lower()
        ea = str(r.get("effect_allele", "")).lower()
        oa = str(r.get("other_allele", "")).lower()
        try:
            bo = float(r["beta"]); so = float(r["standard_error"])
        except (KeyError, ValueError):
            continue
        if ea == exp:
            pass
        elif oa == exp:
            bo = -bo
        else:
            continue
        be, se_e = float(s["beta_e"]), float(s["se_e"])
        if be == 0 or so <= 0:
            continue
        recs.append(dict(rsid=rs, be=be, se_e=se_e, bo=bo, so=so))

    result = dict(
        design="F10 cross-ancestry SENSITIVITY (European BMD instruments -> "
               "non-European OA outcome); NOT the registered de novo EAS replication",
        outcome="MVP African-ancestry OA (GCST90134287, n=35753, "
                "McDonald 2022 Nat Genet)",
        exposure="31 European-only total-body DXA BMD instruments",
        n_instruments=len(inst), n_rsid_matched=len(matched),
        n_harmonised=len(recs),
        chunks=len(chunks), rows_decoded=len(rows),
    )
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    if len(recs) < 3:
        result["verdict"] = (f"INSUFFICIENT: only {len(recs)} harmonised SNPs; "
                             f"cannot run IVW/Egger")
        print(f"\n[VERDICT] {result['verdict']}")
        json.dump(result, open(OUT, "w", encoding="utf-8"),
                  ensure_ascii=False, indent=2)
        print(f"[write] {OUT}")
        return

    be = np.array([r["be"] for r in recs])
    se_e = np.array([r["se_e"] for r in recs])
    bo = np.array([r["bo"] for r in recs])
    so = np.array([r["so"] for r in recs])
    mr = mr13.run_mr(be, se_e, bo, so, snps=[r["rsid"] for r in recs])
    beta = mr["ivw"]; se = mr["ivw_se"]
    result.update(
        ivw_beta=float(beta), ivw_se=float(se), ivw_p=mr.get("ivw_p"),
        or_sd=math.exp(beta),
        ci_lo=math.exp(beta - 1.96 * se), ci_hi=math.exp(beta + 1.96 * se),
        egger_intercept=mr.get("egger_intercept"),
        egger_intercept_p=mr.get("egger_intercept_p"),
        Q=mr.get("Q"), Q_df=mr.get("Q_df"), Q_p=mr.get("Q_p"),
        loo_min=mr.get("loo_min"), loo_max=mr.get("loo_max"),
        loo_sign_flips=mr.get("loo_sign_flips"),
        snps=recs,
    )
    # Honest power gate: with only a handful of harmonised instruments the
    # IVW estimate is NOT interpretable. A null (or a sign that flips between
    # window radii) is not evidence of direction stability OR of reversal.
    n_ok = len(recs)
    p_ivw = mr.get("ivw_p")
    if n_ok < 10 or (p_ivw is not None and p_ivw > 0.05):
        result["verdict"] = (
            f"UNINFORMATIVE / UNDERPOWERED: only {n_ok}/31 European BMD "
            f"instruments are present in this African-ancestry OA GWAS "
            f"(IVW beta {beta:+.4f}, P={p_ivw:.2e}). The estimate is "
            "unstable (its sign flips between tabix window radii) and cannot "
            "support ANY claim about direction, stability or reversal. This "
            "analysis therefore does NOT close F10, does NOT change the S23 "
            "9/10 verdict, and must NOT be reported as a cross-ancestry "
            "replication. Recorded only as a negative feasibility check.")
    else:
        result["verdict"] = (
            f"IVW beta {beta:+.4f}, P={p_ivw:.2e} on {n_ok}/31 instruments; "
            "direction does not reverse, but this remains a cross-ancestry "
            "SENSITIVITY with European-ascertained instruments, NOT the "
            "registered de novo EAS replication; S23 9/10 is unchanged.")
    print(f"\n=== F10 cross-ancestry MR · MVP African OA ===")
    print(f"harmonised SNPs: {len(recs)}/{len(inst)}")
    print(f"IVW beta(SE)= {beta:.5f} ({se:.5f})  "
          f"OR/SD= {math.exp(beta):.4f} "
          f"({math.exp(beta-1.96*se):.4f}-{math.exp(beta+1.96*se):.4f})  "
          f"P= {mr.get('ivw_p'):.2e}")
    print(f"Egger intercept P= {mr.get('egger_intercept_p')}; "
          f"Q={mr.get('Q'):.2f} (df {mr.get('Q_df')}, P={mr.get('Q_p'):.3f})")
    print(f"LOO range [{mr.get('loo_min'):.4f}, {mr.get('loo_max'):.4f}], "
          f"sign flips={mr.get('loo_sign_flips')}")
    json.dump(result, open(OUT, "w", encoding="utf-8"),
              ensure_ascii=False, indent=2)
    print(f"\n[VERDICT] {result['verdict']}")
    print(f"[write] {OUT}")


if __name__ == "__main__":
    main()
