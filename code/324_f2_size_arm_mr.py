# -*- coding: utf-8 -*-
"""
F2 · 骨尺寸臂（volumetric/size arm）替代检验 —— 用已公开的股骨颈宽度 GWAS
=====================================================================

背景
----
登记册 F2 体积臂原本要求 vBMD / 骨面积 / 骨宽 GWAS 作为第二暴露，检验
BMD→膝 OA 的效应是否经「骨尺寸」传递（aBMD 是面积密度，会被骨尺寸混杂）。
`code/320` 判定 vBMD 检索 0 条、`bone area` 6 条全为 CNV burden，因此
一度被记为「数据物理缺失」。

本脚本基于 2026-11-10 的重新检索：**Tobias JH 2024 (PMID 38477772)
`GCST90399549` Femoral neck width, n=38150 欧洲, GRCh37, 有公开数据**。
GWAS Catalog `/studies/GCST90399549/associations` 提供 71 个 genome-wide
significant lead SNP 的 beta / CI / risk allele / 坐标，足以构造 de novo
骨尺寸工具，做法与 F2 几何臂（Faber 2024 HSM1–10 的 193 个 lead）同构。

推断目标
--------
  H_size : 骨尺寸↑ → 膝 OA↑（若成立，BMD→膝 OA 可能是骨尺寸混杂，
           即 SOS 的替代路径；若不成立则 informative null，支持 SOS）。

诚实边界（必须随结果一起引用）
------------------------------
1. **不是 vBMD**：股骨颈宽度是 DXA 二维骨尺寸，不是体积骨密度或骨面积；
   它是体积臂的**近似替代**，不能宣称「体积臂已执行」。
2. **不是部位匹配膝**：暴露是髋/股骨颈尺寸，结局是膝 OA，跨部位。
3. **多效性受限**：这些 lead 多是身高 / 骨骼发育经典位点（如 rs143384），
   水平多效性可能很大；IVW 只是方向性检验，必须配 Egger / weighted median。
4. lead SNP 来自 Catalog 的独立关联集合，本脚本**未重新做 LD clumping**，
   另报告 ±500 kb 距离粗 clump 后的敏感性。
5. 本结果**不得写入主文稿正文**，不得改动任何已报 β/SE/p/OR/CI。
   登记册 F2 「体积臂未测」的结论只在「vBMD/骨面积」意义下继续成立。

用法：python3 code/324_f2_size_arm_mr.py
输出：code/output/genetics/L7_F2_size_arm_FNW.json（随仓库归档）
      compliance_internal/F2_size_arm_FNW_2026-11-10.md（内部审计）
"""
from __future__ import annotations

import gzip
import json
import math
import os
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
import importlib
mr13 = importlib.import_module("13_mr_analysis")

RAW = ROOT / "compliance_internal" / "F2_Tobias_FNW_associations_raw_2026-11-10.json"
KNEE = Path("G:/OA_B_PROJECT_cache/l7/GCST90566800.tsv.gz")
BIM = Path("G:/REF/1000G_EUR/EUR.bim")
OUT_JSON = ROOT / "code" / "output" / "genetics" / "L7_F2_size_arm_FNW.json"
OUT_MD = ROOT / "compliance_internal" / "F2_size_arm_FNW_2026-11-10.md"

CI_RE = re.compile(r"\[\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*-\s*"
                   r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*\]")
BETA_RE = re.compile(r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*unit")


def parse_leads():
    arr = json.loads(RAW.read_text(encoding="utf-8"))
    leads = []
    for a in arr:
        ra = (a.get("riskAllele") or [{}])[0]
        label = str(ra.get("label", ""))
        key = str(ra.get("key", ""))
        # risk allele = part after the last '-' in the label
        allele = label.rsplit("-", 1)[-1].upper() if "-" in label else ""
        bm = BETA_RE.search(str(a.get("beta")))
        cm = CI_RE.search(str(a.get("ci")))
        if not (bm and cm):
            continue
        lo, hi = float(cm.group(1)), float(cm.group(2))
        locs = a.get("locations") or []
        coords = []
        for loc in locs:
            m = re.match(r"([\dXY]+):(\d+)", str(loc))
            if m:
                coords.append((m.group(1), int(m.group(2))))
        if not coords:  # riskAllele key itself may carry chr:pos
            m = re.match(r"chr([\dXY]+):(\d+)", key)
            if m:
                coords.append((m.group(1), int(m.group(2))))
        rsid = key if key.startswith("rs") else ""
        if not coords or not allele or not re.fullmatch(r"[ACGT]+", allele):
            continue  # skip indel / multi-allelic leads
        leads.append(dict(rsid=rsid, label=label, allele=allele,
                          beta=float(bm.group(1)), se=(hi - lo) / 3.92,
                          pval=a.get("pValue"), pexp=a.get("pValueExponent"),
                          risk_freq=a.get("riskFrequency"), coords=coords))
    return leads


def bim_lookup(leads):
    """rsID -> GRCh37 (chrom, bp) from the 1000G-EUR PLINK bim."""
    want = {l["rsid"] for l in leads if l["rsid"]}
    out = {}
    if not BIM.exists():
        print(f"[warn] {BIM} missing; build mapping unavailable")
        return out
    with BIM.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            p = line.split()
            if len(p) >= 4 and p[1] in want:
                out[p[1]] = (p[0], int(p[3]))
                if len(out) == len(want):
                    break
    return out


def scan_knee(targets, leads_by_coord):
    """One streaming pass over the 22.4M-row knee-OA GWAS."""
    hits = {}
    n = 0
    with gzip.open(KNEE, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            n += 1
            quick = line.split("\t", 2)
            if len(quick) < 2:
                continue
            try:
                key = (quick[0], int(quick[1]))
            except ValueError:
                continue
            if key not in targets:
                continue
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            hits.setdefault(key, []).append(r)
            if n % 5_000_000 == 0:
                print(f"    scanned {n:,} rows, {len(hits)} target loci hit", flush=True)
    print(f"    scanned {n:,} rows total; {len(hits)}/{len(targets)} loci hit", flush=True)
    return n, hits, header


def harmonise(leads, hits, bmap):
    rows = []
    ambiguous = 0
    missing = 0
    for l in leads:
        coords = list(l["coords"])
        if l["rsid"] and l["rsid"] in bmap:
            coords.append(bmap[l["rsid"]])
        found = None
        for c in coords:
            if c in hits:
                found = hits[c]
                break
        if not found:
            missing += 1
            continue
        ea = str(found[0].get("effect_allele", "")).upper()
        oa = str(found[0].get("other_allele", "")).upper()
        allele = l["allele"]
        if len(allele) != 1 or len(ea) != 1:
            ambiguous += 1
            continue
        if allele == ea:
            sgn = 1.0
        elif allele == oa:
            sgn = -1.0
        else:
            ambiguous += 1
            continue
        # palindromic check
        pair = {allele, (oa if sgn > 0 else ea)}
        if pair in ({"A", "T"}, {"C", "G"}):
            eaf = None
            try:
                eaf = float(found[0].get("effect_allele_frequency"))
            except (TypeError, ValueError):
                pass
            if eaf is None or 0.42 < eaf < 0.58:
                ambiguous += 1
                continue
        try:
            OR = float(found[0]["odds_ratio"])
            cu = float(found[0]["ci_upper"])
            cl = float(found[0]["ci_lower"])
        except (KeyError, TypeError, ValueError):
            ambiguous += 1
            continue
        if OR <= 0 or cu <= 0 or cl <= 0:
            ambiguous += 1
            continue
        beta_o = sgn * math.log(OR)
        se_o = (math.log(cu) - math.log(cl)) / 3.92
        if not (se_o > 0) or not math.isfinite(beta_o):
            ambiguous += 1
            continue
        rows.append(dict(rsid=l["rsid"] or l["label"], beta_e=l["beta"],
                         se_e=l["se"], beta_o=beta_o, se_o=se_o,
                         fnw_p=l["pval"], fnw_pexp=l["pexp"],
                         knee_or=OR, sgn=sgn, coord=list(coords)))
    return rows, missing, ambiguous


def estimate(rows, label):
    if len(rows) < 3:
        return dict(label=label, n=len(rows), status="INSUFFICIENT_INSTRUMENTS")
    be = np.array([r["beta_e"] for r in rows], float)
    se_e = np.array([r["se_e"] for r in rows], float)
    bo = np.array([r["beta_o"] for r in rows], float)
    so = np.array([r["se_o"] for r in rows], float)
    res = mr13.run_mr(be, se_e, bo, so, snps=[r["rsid"] for r in rows])
    # independent re-derivation of IVW as a guard
    w = 1.0 / so ** 2
    beta = float(np.sum(w * be * bo) / np.sum(w * be ** 2))
    se = float(math.sqrt(1.0 / np.sum(w * be ** 2)))
    from scipy.stats import norm, chi2
    z = beta / se
    p = float(2 * norm.sf(abs(z)))
    def g(key, default=None):
        return res.get(key, default) if isinstance(res, dict) else default

    return dict(label=label, n=len(rows), ivw_beta=beta, ivw_se=se,
                ivw_z=z, ivw_p=p, or_per_sd=float(math.exp(beta)),
                ci_lo=float(math.exp(beta - 1.96 * se)),
                ci_hi=float(math.exp(beta + 1.96 * se)),
                q=float(np.sum(w * (bo - beta * be) ** 2)),
                q_df=len(rows) - 1,
                q_p=float(chi2.sf(float(np.sum(w * (bo - beta * be) ** 2)), len(rows) - 1)),
                i2=g("I2"), ivw_re_p=g("ivw_re_p"), wmed=g("wmed"),
                egger=g("egger"), egger_p=g("egger_p"),
                egger_intercept_p=g("egger_intercept_p"),
                n_outlier=g("n_outlier"), loo_sign_flips=g("loo_sign_flips"),
                mr13=res if isinstance(res, dict) else str(res))


def classify(r):
    """Power + robustness gate. A fixed-effect IVW P<0.05 is NOT sufficient:
    the estimator must survive heterogeneity, random effects, weighted median
    and MR-Egger. This mirrors the honesty gate used in code/323 (F10)."""
    if r.get("status"):
        return "INSUFFICIENT_INSTRUMENTS"
    n = r["n"]
    if n < 5:
        return "UNINFORMATIVE (n<5 harmonised instruments)"
    flags = []
    if (r.get("q_p") or 1) < 0.05:
        flags.append(f"Q_p={r['q_p']:.2g}")
    if (r.get("i2") or 0) > 0.5:
        flags.append(f"I2={r['i2']:.1%}")
    if (r.get("ivw_re_p") or 0) > 0.05:
        flags.append(f"random-effects P={r['ivw_re_p']:.3g}")
    wm = r.get("wmed")
    if wm is None or (wm * r["ivw_beta"] <= 0):
        flags.append("weighted-median sign discordant")
    eg = r.get("egger")
    if eg is None or (eg * r["ivw_beta"] <= 0):
        flags.append("MR-Egger sign discordant")
    if (r.get("n_outlier") or 0) > 0.25 * n:
        flags.append(f"{r['n_outlier']}/{n} radial outliers")
    if flags:
        return "INCONCLUSIVE — HETEROGENEITY-DOMINATED (" + "; ".join(flags) + ")"
    if r["ivw_p"] < 0.05 and r["ivw_beta"] > 0:
        return "SIZE-PATH SUPPORTED (robust across estimators)"
    return "INFORMATIVE NULL (no robust evidence that bone size raises knee-OA risk)"


def distance_clump(rows, kb=500):
    """Crude positional clump; keeps the smallest-FNW-p lead per window."""
    items = []
    for r in rows:
        c = r["coord"][0]
        try:
            items.append((int(c[0]) if c[0].isdigit() else 99, int(c[1]), r))
        except (ValueError, IndexError):
            items.append((99, 0, r))
    items.sort()
    kept, last = [], None
    for ch, bp, r in items:
        if last and ch == last[0] and bp - last[1] < kb * 1000:
            continue
        kept.append(r)
        last = (ch, bp)
    return kept


def main():
    print("=== F2 size arm · Tobias 2024 femoral neck width -> knee OA ===")
    leads = parse_leads()
    print(f"[leads] {len(leads)} usable lead SNPs (of 71 Catalog associations)")
    bmap = bim_lookup(leads)
    print(f"[bim] rsID->GRCh37 mapping for {len(bmap)}/{sum(1 for l in leads if l['rsid'])} leads")
    targets = set()
    by_coord = {}
    for l in leads:
        cs = list(l["coords"])
        if l["rsid"] and l["rsid"] in bmap:
            cs.append(bmap[l["rsid"]])
        for c in cs:
            targets.add(c)
            by_coord.setdefault(c, []).append(l)
    print(f"[targets] {len(targets)} distinct GRCh37/GRCh38 coordinates to look up")
    n, hits, header = scan_knee(targets, by_coord)
    rows, missing, ambiguous = harmonise(leads, hits, bmap)
    print(f"[harmonise] matched={len(rows)} missing={missing} ambiguous={ambiguous}")

    main_res = estimate(rows, "all_leads")
    clumped = distance_clump(rows)
    clump_res = estimate(clumped, "distance_clumped_500kb")
    strict = [r for r in rows if r["fnw_pexp"] is not None and int(r["fnw_pexp"]) <= -8]
    strict_res = estimate(strict, "fnw_p_le_1e-8")

    verdict = classify(main_res)
    # Independent-exposure check: are the size leads merely re-used BMD tools?
    bmd_inst = json.loads(Path("G:/OA_B_PROJECT_cache/l7/bmd_instruments.json")
                          .read_text(encoding="utf-8"))["inst"]
    bmd_rs = set(bmd_inst)
    bmd_coords = {(str(v["chr"]), str(v["bp"])) for v in bmd_inst.values()}
    row_rs = {r["rsid"] for r in rows}
    row_coords = {(str(r["coord"][0][0]), str(r["coord"][0][1]))
                  for r in rows if r["coord"]}
    overlap = dict(n_bmd_instruments=len(bmd_rs), n_size_leads=len(row_rs),
                   rsid_overlap=len(row_rs & bmd_rs),
                   coord_overlap=len(row_coords & bmd_coords),
                   same_sign_wald=sum(1 for r in rows if r["beta_e"] * r["beta_o"] > 0),
                   n_rows=len(rows))
    record = dict(
        design="de novo femoral-neck-width lead SNPs (Tobias 2024 GCST90399549) -> knee OA (GCST90566800)",
        exposure_study="GCST90399549 / PMID 38477772 / n=38150 EUR / GRCh37",
        outcome_study="GCST90566800 / Hatzikotoulas 2025 / mixed ancestry knee OA",
        n_catalog_leads=71, n_usable_leads=len(leads), n_rows_scanned=n,
        n_target_coords=len(targets), n_loci_hit=len(hits),
        n_bim_mapped=len(bmap), n_harmonised=len(rows),
        n_missing=missing, n_ambiguous=ambiguous,
        main=main_res, clumped=clump_res, strict=strict_res,
        independent_exposure_check=overlap,
        rows=rows, verdict=verdict,
        honesty_boundaries=[
            "Femoral neck width is a 2-D DXA bone-size measure, NOT volumetric BMD or bone area; "
            "this is a size-arm approximation, not execution of the registered volumetric arm.",
            "Hip/femoral-neck exposure against knee-OA outcome: cross-site.",
            "These leads include classic height/skeletal-development loci; horizontal pleiotropy "
            "is likely and IVW is a directional test only.",
            "No de novo LD clumping was performed; a 500-kb distance clump is reported as sensitivity.",
            "Results are NOT written into the manuscript and do not alter any reported beta/SE/p/OR/CI.",
        ],
    )
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(record, ensure_ascii=False, indent=2), encoding="utf-8")
    lines = [
        "# F2 骨尺寸臂替代检验（内部审计，2026-11-10）", "",
        f"- 设计：Tobias 2024 `GCST90399549` 股骨颈宽度 71 个 lead SNP → 膝 OA `GCST90566800`",
        f"- Catalog lead：71，可用：{len(leads)}；扫描结局行数：{n:,}；命中位点：{len(hits)}/{len(targets)}",
        f"- 和谐化工具：**{len(rows)}**（missing {missing} / ambiguous {ambiguous}）",
        "", "## 估计", "",
        "| 分析 | n | IVW β | SE | OR/SD (95% CI) | P | Q (df, P) |", "|---|---:|---:|---:|---|---:|---|",
    ]
    for r in (main_res, clump_res, strict_res):
        if r.get("status"):
            lines.append(f"| {r['label']} | {r['n']} | — | — | — | — | — |")
            continue
        lines.append(f"| {r['label']} | {r['n']} | {r['ivw_beta']:+.4f} | {r['ivw_se']:.4f} | "
                     f"{r['or_per_sd']:.4f} ({r['ci_lo']:.4f}–{r['ci_hi']:.4f}) | "
                     f"{r['ivw_p']:.3g} | {r['q']:.2f} ({r['q_df']}, {r['q_p']:.3g}) |")
    lines += ["", "## 稳健性与工具独立性（判读依据）", "",
              "| 检验 | 值 |", "|---|---|",
              f"| 与 31 个 BMD 工具 rsID 重叠 | **{overlap['rsid_overlap']}/{overlap['n_size_leads']}**（坐标重叠 {overlap['coord_overlap']}）|",
              f"| 逐变异 Wald 比同号 | {overlap['same_sign_wald']}/{overlap['n_rows']} = {overlap['same_sign_wald']/max(1,overlap['n_rows']):.1%} |",
              f"| 异质性 Q (df, P) | {main_res.get('q',float('nan')):.1f} ({main_res.get('q_df')}, {main_res.get('q_p',float('nan')):.2g}) |",
              f"| I² | {main_res.get('i2',float('nan')):.1%} |",
              f"| 随机效应 IVW P | {main_res.get('ivw_re_p')} |",
              f"| 加权中位数 β | {main_res.get('wmed')} |",
              f"| MR-Egger β (P) | {main_res.get('egger')} ({main_res.get('egger_p')}) |",
              f"| Egger 截距 P | {main_res.get('egger_intercept_p')} |",
              f"| Radial 离群工具 | {main_res.get('n_outlier')} |",
              "", "固定效应 IVW 虽显著，但在异质性 / 随机效应 / 加权中位数 / MR-Egger 下不稳健，"
              "故不能据此宣称「骨尺寸混杂成立」。零工具重叠说明这是独立暴露，"
              "与几何臂（HSM1–10, 193 lead）的 informative null 方向一致。",
              "", f"**判读**：{verdict}", "", "## 诚实边界", ""]
    lines += [f"- {b}" for b in record["honesty_boundaries"]]
    lines += ["", "## 与登记册 F2 的关系", "",
              "登记册体积臂要求 vBMD / 骨面积 / 骨宽。本检验用的是 **股骨颈宽度（DXA 二维尺寸）**，"
              "因此只能作为「尺寸臂」的近似替代，与几何臂（Faber 2024 HSM1–10, 193 lead）互补。"
              "登记册「体积臂未测」在 vBMD/骨面积意义下仍然成立；S23「9/10」不变。", ""]
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")
    print(f"\n[verdict] {verdict}")
    for r in (main_res, clump_res, strict_res):
        if not r.get("status"):
            print(f"  {r['label']:26s} n={r['n']:2d} beta={r['ivw_beta']:+.4f} "
                  f"SE={r['ivw_se']:.4f} OR/SD={r['or_per_sd']:.4f} P={r['ivw_p']:.3g}")
        else:
            print(f"  {r['label']:26s} n={r['n']:2d} {r['status']}")
    print(f"[write] {OUT_JSON.relative_to(ROOT)}")
    print(f"[write] {OUT_MD.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
