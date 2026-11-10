# -*- coding: utf-8 -*-
"""
F10 跨祖先复制 · 公开数据可得性探针（只读网络流，不落盘大文件）

STATUS: ABANDONED PROTOTYPE — 未执行完、未产出结果、未写入稿件。
------------------------------------------------------------------
2026-11-09 复核结论：现有非欧 / 祖先分层公开汇总统计均不满足登记册 F10 的
规格（见 `compliance_internal/public_data_recheck_2026-11-09.md`）：
非欧 BMD 只有 N=2633 的 heel QUS，祖先分层 OA 是 all-site / exome-wide，
用 31 个欧洲工具做符号一致性检验既功效不足，也不是 de novo 跨祖先复制，
故主动放弃该探针，避免把 portability check 误报为 replication。
可复跑的正式核对脚本为 `code/320_public_data_recheck.py`。本文件仅作
设计取舍留痕，**不产生任何已报数值**。

2026-11-10 追加：用户追问「F10 能否解决」后，**真正执行了一次**跨祖先 MR，
而非停留在检索。MVP 2022 祖先分层 OA（GCST90134287, 非洲 n=35753）用
tabix/BGZF 定点抓取实跑，结果仅 4–5/31 个欧洲工具可匹配、IVW 符号随窗口
半径翻转、P>0.37 → UNINFORMATIVE / UNDERPOWERED，不能关闭 F10、不改 S23
9/10、不得写入正文。详见
`compliance_internal/F10_cross_ancestry_execution_2026-11-10.md` 与
`code/323_f10_mvp_tabix_mr.py`。本探针的结论（登记册 F10 不可执行）因此
从「文献检索推断」升级为「实测验证」，但本文件仍不产生任何已报数值。

背景
----
稿件 Data availability 与预注册登记册目前记载：
  "F10's registered cross-ancestry replication not performed for want of
   public ancestry-stratified data"
本脚本对 2025-2026 年新释放的公开汇总统计做**流式**核查，判定该表述是否
仍然成立，以及能否真正执行一次跨祖先检验。

核查对象（GWAS Catalog / EBI FTP，全部为公开匿名汇总统计）
  1. Pan-UK Biobank (Karczewski 2025, Nat Genet, PMID 40968291)
     GCST90692328  South/Central Asian heel BMD, N=2633, GRCh37
  2. Genes & Health (Kim 2026, Nat Genet)
     GCST90727231  British South Asian OA excl. spine, N=44026, GRCh38,
                   **exome-wide sequencing**
  3. MVP + UKB (McDonald 2022, Nat Genet, PMID 36411363)
     GCST90134280  South Asian OA, N=7162, GRCh37 raw / GRCh38 harmonised
     GCST90134281  East Asian OA, N=1704
     GCST90134282  African OA, N=7188
     GCST90134287  African OA, N=35753

方法
----
按 **rsID** 匹配稿件锁定的 31 个欧洲 total-body BMD 工具（bmd_instruments.json），
逐文件统计：匹配数、等位基因可对齐数、beta 符号一致数、Pearson r。
不做任何 MR 估计，不改任何已报告 beta/SE/p/OR/CI。

诚实边界（必须随结果一并报告）
--------------------------------
* 暴露工具仍来自欧洲 total-body DXA BMD；非欧 BMD 可得的是 **heel QUS BMD**
  （不同测量部位/尺度），非欧 OA 可得的是 **全 OA 或 OA excl. spine**
  （非膝 OA 特异）。因此本探针最多支持「跨祖先可移植性 / 方向一致性」检验，
  **不构成**登记册 F10 所要求的 de novo 跨祖先复制（需非欧人群自建工具集）。
* G&H 为 exome-wide，常见调控区工具预期覆盖极低；该结果本身是有用的阴性证据。
"""
import gzip
import io
import json
import os
import sys
import time
import urllib.request

ROOT = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(ROOT)
CACHE = "G:/OA_B_PROJECT_cache/l7"
INST = os.path.join(CACHE, "bmd_instruments.json")
OUT = os.path.join(ROOT, "output", "genetics", "L7_F10_cross_ancestry_probe.json")

FTP = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
TARGETS = [
    # label, url, size_hint, build, ancestry, trait, technology, N
    ("PanUKB_SAS_heelBMD",
     FTP + "GCST90692001-GCST90693000/GCST90692328/harmonised/GCST90692328.h.tsv.gz",
     "484M", "GRCh37", "South/Central Asian", "heel BMD (UKB 3148)",
     "genome-wide array", 2633),
    ("GenesHealth_SAS_OA",
     FTP + "GCST90727001-GCST90728000/GCST90727231/harmonised/GCST90727231.h.tsv.gz",
     "47M", "GRCh38", "British South Asian", "OA excluding spine",
     "exome-wide sequencing", 44026),
]

HDR = {"User-Agent": "Mozilla/5.0 (F10 public-data availability probe)"}


def load_instruments():
    with open(INST, encoding="utf-8") as f:
        inst = json.load(f)["inst"]
    out = {}
    for rs, v in inst.items():
        out[rs] = {"chr": str(v["chr"]), "bp": int(v["bp"]),
                   "beta_eur": float(v["beta_e"]), "se_eur": float(v["se_e"]),
                   "eaf_eur": float(v["eaf"]), "exp_allele": v["exp_allele"]}
    return out


def stream_rows(url, max_seconds=900):
    """流式逐行解压远端 .gz，不落盘。yield 解码后的行。"""
    req = urllib.request.Request(url, headers=HDR)
    t0 = time.time()
    with urllib.request.urlopen(req, timeout=60) as resp:
        total = resp.headers.get("Content-Length")
        print(f"    Content-Length={total}", flush=True)
        gz = gzip.GzipFile(fileobj=resp)
        for raw in gz:
            if time.time() - t0 > max_seconds:
                raise TimeoutError(f"exceeded {max_seconds}s")
            yield raw.decode("utf-8", "replace")


def probe(label, url, build, anc, trait, tech, n, inst, size_hint):
    print(f"\n=== {label} ({size_hint}, {build}, N={n}) ===", flush=True)
    hit = {}
    columns = None
    nrows = 0
    t0 = time.time()
    for line in stream_rows(url):
        if nrows == 0:
            columns = line.rstrip("\n\r").split("\t")
            idx = {c: i for i, c in enumerate(columns)}
            need = ["beta", "standard_error", "rsid"]
            missing = [c for c in need if c not in idx]
            if missing:
                raise ValueError(f"{label}: missing columns {missing}; header={columns}")
            print(f"    columns={columns[:16]}", flush=True)
        nrows += 1
        if nrows == 1:
            continue
        f = line.rstrip("\n\r").split("\t")
        try:
            rs = f[idx["rsid"]]
        except IndexError:
            continue
        if rs in inst and rs not in hit:
            try:
                hit[rs] = {
                    "beta": float(f[idx["beta"]]),
                    "se": float(f[idx["standard_error"]]),
                    "chr": f[idx["chromosome"]] if "chromosome" in idx else None,
                    "bp": f[idx["base_pair_location"]] if "base_pair_location" in idx else None,
                    "ea": f[idx["effect_allele"]] if "effect_allele" in idx else None,
                    "oa": f[idx["other_allele"]] if "other_allele" in idx else None,
                    "eaf": f[idx["effect_allele_frequency"]] if "effect_allele_frequency" in idx else None,
                    "p": f[idx["p_value"]] if "p_value" in idx else None,
                }
            except (ValueError, IndexError):
                continue
        if len(hit) == len(inst):
            print("    all 31 instruments found -> early stop", flush=True)
            break
        if nrows % 5_000_000 == 0:
            print(f"    rows={nrows:,} hits={len(hit)} elapsed={time.time()-t0:.0f}s", flush=True)
    elapsed = time.time() - t0
    matched = len(hit)
    sign_conc = None
    r = None
    if matched >= 2:
        pairs = [(inst[rs]["beta_eur"], hit[rs]["beta"]) for rs in hit
                 if hit[rs]["se"] and hit[rs]["se"] > 0]
        sign_conc = sum(1 for a, b in pairs if (a > 0) == (b > 0))
        if len(pairs) >= 3:
            import math
            xs = [a for a, _ in pairs]; ys = [b for _, b in pairs]
            mx = sum(xs) / len(xs); my = sum(ys) / len(ys)
            num = sum((x - mx) * (y - my) for x, y in pairs)
            dx = math.sqrt(sum((x - mx) ** 2 for x in xs))
            dy = math.sqrt(sum((y - my) ** 2 for y in ys))
            r = num / (dx * dy) if dx > 0 and dy > 0 else None
    res = {
        "label": label, "accession": url.split("/")[-2] if "/harmonised/" in url else "",
        "url": url, "build": build, "ancestry": anc, "trait": trait,
        "technology": tech, "N": n, "rows_scanned": nrows,
        "elapsed_seconds": round(elapsed, 1),
        "instruments": len(inst), "matched": matched,
        "match_fraction": round(matched / len(inst), 4),
        "beta_sign_concordant": sign_conc, "beta_pearson_r": r,
        "matched_detail": hit,
    }
    print(f"    rows={nrows:,} matched={matched}/{len(inst)} "
          f"sign_conc={sign_conc} r={r} elapsed={elapsed:.0f}s", flush=True)
    return res


def main():
    inst = load_instruments()
    print(f"locked instruments: {len(inst)}", flush=True)
    results = []
    for label, url, size_hint, build, anc, trait, tech, n in TARGETS:
        try:
            results.append(probe(label, url, build, anc, trait, tech, n, inst, size_hint))
        except Exception as e:
            print(f"    ERROR {type(e).__name__}: {e}", flush=True)
            results.append({"label": label, "url": url, "error": f"{type(e).__name__}: {e}"})
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump({
            "generated": "2026-11-09",
            "purpose": "F10 cross-ancestry public-data availability probe",
            "instrument_source": "bmd_instruments.json (31 European total-body BMD instruments)",
            "honesty_boundary": (
                "Non-European exposure is heel QUS BMD and non-European outcomes are "
                "all-site OA / OA excluding spine; exposure instruments remain European. "
                "Any positive result is trans-ancestry portability / sign concordance, "
                "NOT the registered de novo cross-ancestry replication."),
            "results": results,
        }, f, ensure_ascii=False, indent=1)
    print(f"\nWROTE {OUT}", flush=True)


if __name__ == "__main__":
    main()
