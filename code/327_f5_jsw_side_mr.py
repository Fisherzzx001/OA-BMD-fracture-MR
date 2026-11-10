# -*- coding: utf-8 -*-
"""
F5「换推断层级」· JSN 侧：欧洲 31 BMD 工具 -> 髋最小关节间隙宽度 mJSW
============================================================================
结局 GCST90281365 (Faber 2023 EBioMedicine, PMID 37619450, n=50,745 EUR,
GRCh37)：hip mJSW 连续表型 (mm)，是软骨厚度的影像代理。注册检验 F5 要求
「髋骨赘为结局的 SNP-MR」，髋骨赘 GWAS 至今未公开；mJSW 是同一 Faber
结构内型框架的**关节间隙侧 (JSN side)**，是 F5「换推断层级」的可执行臂。

方向约定：mJSW 越大 = 关节间隙越宽 = 软骨越厚 = 越保护。故 BMD→mJSW
        β<0（高 BMD 使间隙变窄）与「高 BMD 促 OA」的 SOS 关节侧一致。

技术：复用 code/323 的 TBI/BGZF/HTTP-Range 定点抓取（265 MB BGZF 不全下），
      31 工具坐标 finest 16-kb bin -> Range -> 解压 -> rsid 匹配 -> 等位协调
      -> code/13 的 IVW/Egger/LOO。
诚实：这是 JSN 影像侧代理，不是注册的髋骨赘 SNP-MR；不改变 S23 9/10。
      若匹配 <10 或符号在窗口半径间翻转 -> UNINFORMATIVE，不作方向声称。

输入：G:/OA_B_PROJECT_cache/l7/GCST90281365.h.tsv.gz.tbi (1.5 MB, 已下载)
      bmd_instruments.json (31 锁定工具)
输出：code/output/genetics/L7_F5_JSN_mJSW_MR.json
"""
import os, json, math, gzip, struct, socket, urllib.request, importlib.util
import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
CACHE = "G:/OA_B_PROJECT_cache/l7"
INST = os.path.join(CACHE, "bmd_instruments.json")
TBI = os.path.join(CACHE, "GCST90281365.h.tsv.gz.tbi")
BGZF_URL = ("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
            "GCST90281001-GCST90282000/GCST90281365/harmonised/"
            "GCST90281365.h.tsv.gz")
OUT = os.path.join(ROOT, "output", "genetics", "L7_F5_JSN_mJSW_MR.json")
socket.setdefaulttimeout(60)

_spec = importlib.util.spec_from_file_location("mr13", os.path.join(ROOT, "13_mr_analysis.py"))
mr13 = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(mr13)
# reuse TBI/BGZF helpers from 323 without running its main
_s323 = importlib.util.spec_from_file_location("t323", os.path.join(ROOT, "323_f10_mvp_tabix_mr.py"))
t323 = importlib.util.module_from_spec(_s323); _s323.loader.exec_module(t323)


def run(expand):
    inst = json.load(open(INST, encoding="utf-8"))["inst"]
    refs, meta = t323.parse_tbi(TBI)
    loci = [(str(s["chr"]), int(s["bp"])) for s in inst.values()]
    chunks = t323.chunks_for_loci(refs, loci, expand=expand)
    # fetch header from the start of OUR BGZF (323.fetch_header hardcodes its own URL)
    buf = t323.http_range(BGZF_URL, 0, 65536)
    header = b"".join(t323.bgzf_blocks(buf)).decode("utf-8", "replace").split("\n")[0].rstrip("\r").split("\t")
    rows = []
    iv = []
    for vo_beg, vo_end in chunks[:400]:
        cs, ce = vo_beg >> 16, vo_end >> 16
        iv.append((cs, ce + 65536))
    merged = t323.merge_intervals(iv)
    print(f"[327] expand={expand}: {len(chunks)} chunks -> {len(merged)} ranges", flush=True)
    for s, e in merged:
        try:
            b = t323.http_range(BGZF_URL, s, e)
        except Exception as ex:
            print(f"  [warn] range {s}-{e}: {type(ex).__name__}"); continue
        text = b"".join(t323.bgzf_blocks(b)).decode("utf-8", "replace")
        for line in text.split("\n"):
            if not line or line.startswith("#") or line.startswith("chromosome"):
                continue
            parts = line.rstrip("\r").split("\t")
            if len(parts) != len(header): continue
            rows.append(dict(zip(header, parts)))
    wanted = set(inst.keys()); matched = {}
    for r in rows:
        for key in ("rsid", "rs_id", "variant_id"):
            rs = r.get(key)
            if rs in wanted and rs not in matched:
                matched[rs] = r; break
    recs = []
    for rs, s in inst.items():
        r = matched.get(rs)
        if not r: continue
        exp = str(s["exp_allele"]).lower()
        ea = str(r.get("effect_allele", "")).lower()
        oa = str(r.get("other_allele", "")).lower()
        try: bo = float(r["beta"]); so = float(r["standard_error"])
        except (KeyError, ValueError): continue
        if ea == exp: pass
        elif oa == exp: bo = -bo
        else: continue
        be, se_e = float(s["beta_e"]), float(s["se_e"])
        if be == 0 or so <= 0: continue
        recs.append(dict(rsid=rs, be=be, se_e=se_e, bo=bo, so=so))
    return dict(n_matched=len(matched), n_harmonised=len(recs), recs=recs,
                chunks=len(chunks), rows=len(rows))


def main():
    import sys
    expands = [int(x) for x in (sys.argv[1:] or ["0", "16", "32"])]
    runs = {}
    for ex in expands:
        r = run(ex)
        recs = r["recs"]
        entry = dict(n_matched=r["n_matched"], n_harmonised=r["n_harmonised"],
                     chunks=r["chunks"], rows_decoded=r["rows"])
        if len(recs) >= 3:
            be = np.array([x["be"] for x in recs]); se_e = np.array([x["se_e"] for x in recs])
            bo = np.array([x["bo"] for x in recs]); so = np.array([x["so"] for x in recs])
            mr = mr13.run_mr(be, se_e, bo, so, snps=[x["rsid"] for x in recs])
            beta = mr["ivw"]
            entry.update(ivw_beta=round(float(beta), 5), ivw_p=mr.get("ivw_p"),
                         or_per_mm=round(math.exp(beta), 5),
                         egger_intercept=mr.get("egger_intercept"),
                         egger_intercept_p=mr.get("egger_intercept_p"),
                         Q=mr.get("Q"), Q_df=mr.get("Q_df"), Q_p=mr.get("Q_p"),
                         loo_sign_flips=mr.get("loo_sign_flips"))
            entry["recs"] = recs
        runs[str(ex)] = entry
        print(f"[327] expand={ex}: matched {r['n_matched']}/31, harmonised {r['n_harmonised']}/31"
              + (f", IVW β={entry.get('ivw_beta')} P={entry.get('ivw_p'):.2e}" if "ivw_beta" in entry else "")
              + flush_str(), flush=True)
    # verdict across radii
    betas = [runs[str(e)].get("ivw_beta") for e in expands]
    nh = [runs[str(e)]["n_harmonised"] for e in expands]
    signs = [int(np.sign(b)) for b in betas if b is not None]
    stable = bool(signs) and (max(signs) == min(signs))
    # direction concordance with SOS joint arm: BMD raises OA == mJSW narrows (beta<0)
    sos_concordant = stable and signs[0] < 0
    if min(nh) < 10 or not stable:
        verdict = ("UNINFORMATIVE / UNDERPOWERED: <10 harmonised instruments or "
                   "IVW sign flips between tabix window radii; JSN-side direction "
                   "cannot be resolved. Not the registered hip-osteophyte SNP-MR; "
                   "S23 9/10 unchanged.")
    elif sos_concordant:
        verdict = (f"DIRECTION-CONCORDANT with SOS joint arm: higher BMD narrows "
                   f"hip joint space (mJSW β<0, stable across radii), the JSN-side "
                   f"imaging proxy of osteophyte-dominant formation. Still a JSN "
                   f"surrogate, NOT the registered hip-osteophyte SNP-MR; S23 9/10 unchanged.")
    else:
        verdict = ("DIRECTION-DISCORDANT: higher BMD widens hip joint space; JSN "
                   "proxy does not track the knee-OA risk arm. Reported honestly; "
                   "S23 9/10 unchanged.")
    out = dict(design="F5 JSN-side imaging proxy: 31 European BMD instruments -> hip mJSW",
               outcome="GCST90281365 hip minimum joint space width (Faber 2023, n=50745 EUR, mm)",
               note="mJSW larger = wider joint space = more cartilage = protective; NOT registered hip-osteophyte SNP-MR",
               runs=runs, verdict=verdict)
    json.dump(out, open(OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print(f"\n[VERDICT] {verdict}")
    print(f"[write] {OUT}")


def flush_str():
    return ""

if __name__ == "__main__":
    main()
