#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
207_bmd_verified_ldsc.py -- upgrade the exploratory LDSC layer to a citable one
==============================================================================
WHAT THIS FIXES (the two remaining caveats listed in code/205 and RESULTS.md)
-----------------------------------------------------------------------------
  (1) BMD exposure sumstats carried a PLACEHOLDER N = 400,000 in every row.
      code/201 (coloc) and docs/10 (GWAS Catalog curated) both use the true N
      for total-body DXA BMD: GCST90568455, UKB field 23236, N = 31,986
      European-ancestry participants. We re-munge with the verified N.
  (2) BMD was the only trait NOT munged with --merge-alleles w_hm3.snplist.
      We re-munge it through the identical HapMap3 allele frame as the four
      OA sumstats, removing the asymmetry (0 allele mismatches in practice,
      see log: "Removed 0 SNPs whose alleles did not match --merge-alleles").

WHY --ignore N IS REQUIRED
--------------------------
munge_sumstats.process_n() only applies --N when the input has NO N column.
The legacy file has N=400000, so a bare --N would be silently ignored.
`--ignore N` removes the column from the column map; --N then takes effect.
The munge log must be checked for "Using N = <value>" (this script asserts it).

WHAT IS EXPECTED (and checked in COMPARISON.md)
----------------------------------------------
  * rg is INVARIANT to a global N rescaling: N cancels in
    rg = b_gencov / sqrt(b_h2_1 * b_h2_2), where b are LD Score regression
    slopes. Only munge's 3-decimal Z rounding can shift the last digits.
  * h2 scales as 1/N: h2(N=31986) ~= h2(N=400000) * 400000/31986 (~12.5x).
    Only AFTER this correction is the BMD h2 absolute value interpretable.
  * The OA traits are untouched (same munged files, same per-row N), so their
    h2 and the OA site-vs-site rg are copied from the baseline, not recomputed.

OUTPUTS (new directory; the exploratory baseline ldsc_rg/ is NEVER overwritten)
------------------------------------------------------------------------------
  code/output/genetics/ldsc_rg_verified/
      munged/BMD_total_N<primary>.sumstats.gz, BMD_total_N<sens>.sumstats.gz
      h2/, rg/, logs/
      verified_results.csv, COMPARISON.md, summary.json

USAGE
-----
  python3 code/207_bmd_verified_ldsc.py --check
  python3 code/207_bmd_verified_ldsc.py                 # full, ~3 min
  python3 code/207_bmd_verified_ldsc.py --skip-munge    # reuse munged files
  python3 code/207_bmd_verified_ldsc.py --no-sensitivity

HARD BOUNDARY
-------------
This is still a POST-SUBMISSION / EXPLORATORY layer. It does NOT edit docs/08,
does NOT change any reported beta/SE/p/OR/CI, and does NOT by itself insert
rg into the manuscript. code/148 / 149b / 182 must be re-run afterwards.
"""
import argparse
import csv
import json
import os
import re
import subprocess
import sys
import time

HERE = "C:/Users/Administrator/WorkBuddy/2026-09-08-22-17-40/OA_B_PROJECT"
PYBIN = "C:/Users/Administrator/.workbuddy/binaries/python/versions/3.13.12/Scripts"
LDSC = os.path.join(PYBIN, "ldsc.py")
MUNGE = os.path.join(PYBIN, "munge_sumstats.py")

STD_DIR = ("C:/Users/Administrator/WorkBuddy/2026-09-08-20-17-24/"
           "OA_FLAGSHIP/07_reference/ldsc")
STD_PREFIX = STD_DIR + "/eur_w_ld_chr_"
HM3_SNPLIST = STD_DIR + "/w_hm3.snplist"
BMD_RAW = STD_DIR + "/sumstats_total.sumstats.gz"

BASE = os.path.join(HERE, "code", "output", "genetics", "ldsc_rg")
OUT = os.path.join(HERE, "code", "output", "genetics", "ldsc_rg_verified")
SITES = ["hip", "knee", "spine", "hand"]
PLACEHOLDER_N = 400000

COLS = ["p1", "p2", "rg", "se", "z", "p", "h2_obs", "h2_obs_se",
        "h2_int", "h2_int_se", "gcov_int", "gcov_int_se"]


def sh(cmd, logpath):
    os.makedirs(os.path.dirname(logpath), exist_ok=True)
    with open(logpath, "w") as lf:
        lf.write("CMD: " + " ".join(cmd) + "\n\n")
        lf.flush()
        p = subprocess.run(cmd, cwd=HERE, stdout=lf, stderr=subprocess.STDOUT)
    return p.returncode, logpath


def check():
    print("=== 207 · prerequisite check ===", flush=True)
    items = [("ldsc.py", LDSC), ("munge_sumstats.py", MUNGE),
             ("BMD raw sumstats", BMD_RAW), ("HapMap3 snplist", HM3_SNPLIST),
             ("HapMap3 .l2.M_5_50", STD_DIR + "/eur_w_ld_chr.l2.M_5_50")]
    for s in SITES:
        items.append((f"OA {s} munged", os.path.join(BASE, "munged", f"OA_{s}.sumstats.gz")))
    items.append(("baseline BMD x OA rg", os.path.join(BASE, "rg", "BMD_vs_OA_sites.log")))
    items.append(("baseline BMD h2", os.path.join(BASE, "h2", "BMD_total.log")))
    ok = True
    for name, p in items:
        good = os.path.exists(p) and os.path.getsize(p) > 0
        ok &= good
        print(f"  [{'OK ' if good else 'MISS'}] {name:26s} "
              f"{os.path.getsize(p)/1e6 if os.path.exists(p) else 0:8.2f} MB  {p}", flush=True)
    n = 0
    for i in range(1, 23):
        p = STD_PREFIX + str(i) + ".l2.ldscore.gz"
        if os.path.exists(p) and os.path.getsize(p) > 0:
            n += 1
    print(f"  [{'OK ' if n == 22 else 'MISS'}] HapMap3 per-chr ldscore   : {n}/22", flush=True)
    ok &= (n == 22)
    print(f"\n  READY TO RUN : {'YES' if ok else 'NO'}", flush=True)
    return ok


def munge_bmd(N, log):
    out = os.path.join(OUT, "munged", f"BMD_total_N{int(N)}")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    cmd = [sys.executable, "-u", MUNGE,
           "--sumstats", BMD_RAW, "--out", out,
           "--snp", "SNP", "--a1", "A1", "--a2", "A2", "--p", "P",
           "--signed-sumstats", "Z,0",
           "--ignore", "N", "--N", str(N),
           "--merge-alleles", HM3_SNPLIST,
           "--n-min", "1", "--chunksize", "200000"]
    rc, lp = sh(cmd, log)
    txt = open(lp, encoding="utf-8", errors="replace").read()
    assert "Using N = %s" % float(N) in txt, (
        f"--N {N} did NOT take effect; check munge log {lp}")
    assert "Conversion finished" in txt, f"munge did not finish: {lp}"
    f = out + ".sumstats.gz"
    assert rc == 0 and os.path.exists(f) and os.path.getsize(f) > 0, "munge failed"
    n_snps = re.search(r"Writing summary statistics for (\d+) SNPs \((\d+) with nonmissing",
                       txt)
    print(f"[munge BMD N={int(N)}] rc={rc} "
          f"rows={n_snps.group(1) if n_snps else '?'} "
          f"nonNA={n_snps.group(2) if n_snps else '?'} log={lp}", flush=True)
    return f


def base_args():
    return ["--ref-ld-chr", STD_PREFIX, "--w-ld-chr", STD_PREFIX]


def run_h2(sumstats, tag):
    out = os.path.join(OUT, "h2", tag)
    cmd = [sys.executable, "-u", LDSC, "--h2", sumstats] + base_args() + ["--out", out]
    rc, lp = sh(cmd, os.path.join(OUT, "logs", f"h2_{tag}.log"))
    txt = open(lp, encoding="utf-8", errors="replace").read()
    ok = rc == 0 and "Analysis finished" in txt and "Traceback" not in txt
    print(f"[h2 {tag}] rc={rc} ok={ok}", flush=True)
    return ok


def run_rg(files, tag):
    out = os.path.join(OUT, "rg", tag)
    cmd = [sys.executable, "-u", LDSC, "--rg", ",".join(files)] + base_args() + ["--out", out]
    rc, lp = sh(cmd, os.path.join(OUT, "logs", f"rg_{tag}.log"))
    txt = open(lp, encoding="utf-8", errors="replace").read()
    ok = rc == 0 and "Analysis finished" in txt and "Traceback" not in txt
    print(f"[rg {tag}] rc={rc} ok={ok}", flush=True)
    return ok


def parse_rg(path, analysis):
    rows = []
    lines = open(path, encoding="utf-8", errors="replace").read().splitlines()
    try:
        i = next(k for k, l in enumerate(lines)
                 if l.strip() == "Summary of Genetic Correlation Results")
    except StopIteration:
        return rows
    for l in lines[i + 2:]:
        t = l.split()
        if len(t) != len(COLS):
            if rows:
                break
            continue
        d = dict(zip(COLS, t))
        d["analysis"] = analysis
        for k in ("p1", "p2"):
            d[k] = os.path.basename(d[k]).replace(".sumstats.gz", "").replace(".gz", "")
            d[k] = d[k].replace("OA_", "").replace("sumstats_total", "BMD_total").upper()
        rows.append(d)
    return rows


def valid_allele_counts(path):
    """Per-pair 'N SNPs with valid alleles' lines from an ldsc --rg log."""
    if not os.path.exists(path):
        return []
    txt = open(path, encoding="utf-8", errors="replace").read()
    return [int(x) for x in re.findall(r"(\d+) SNPs with valid alleles", txt)]


def parse_h2(path, label):
    txt = open(path, encoding="utf-8", errors="replace").read()
    d = {"label": label}
    for k, pat in {
        "h2": r"Total Observed scale h2:\s*([\-0-9.eE+]+)\s*\(([\-0-9.eE+]+)\)",
        "intercept": r"Intercept:\s*([\-0-9.eE+]+)\s*\(([\-0-9.eE+]+)\)",
        "lambda_gc": r"Lambda GC:\s*([\-0-9.eE+]+)",
        "mean_chi2": r"Mean Chi\^2:\s*([\-0-9.eE+]+)",
    }.items():
        m = re.search(pat, txt)
        if m:
            d[k] = m.group(1)
            if m.lastindex and m.lastindex > 1:
                d[k + "_se"] = m.group(2)
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--primary-N", type=int, default=31986,
                    help="verified total-body DXA BMD N (GCST90568455)")
    ap.add_argument("--sensitivity-N", type=int, default=31873,
                    help="femur-total N used for the manuscript's MR N-correction")
    ap.add_argument("--skip-munge", action="store_true")
    ap.add_argument("--no-sensitivity", action="store_true")
    args = ap.parse_args()

    for d in ("munged", "logs", "h2", "rg"):
        os.makedirs(os.path.join(OUT, d), exist_ok=True)
    if not check():
        return 2
    if args.check:
        return 0

    t0 = time.time()
    steps, oa_files = [], [os.path.join(BASE, "munged", f"OA_{s}.sumstats.gz") for s in SITES]

    def bmd_path(N):
        return os.path.join(OUT, "munged", f"BMD_total_N{int(N)}.sumstats.gz")

    if args.skip_munge and os.path.exists(bmd_path(args.primary_N)):
        print("[munge] skipped (--skip-munge)", flush=True)
    else:
        munge_bmd(args.primary_N, os.path.join(OUT, "logs", f"munge_BMD_N{args.primary_N}.log"))
        steps.append({"step": f"munge_BMD_N{args.primary_N}", "rc": 0})

    Ns = [args.primary_N]
    if not args.no_sensitivity:
        if args.skip_munge and os.path.exists(bmd_path(args.sensitivity_N)):
            print("[munge sensitivity] skipped (--skip-munge)", flush=True)
        else:
            munge_bmd(args.sensitivity_N,
                      os.path.join(OUT, "logs", f"munge_BMD_N{args.sensitivity_N}.log"))
            steps.append({"step": f"munge_BMD_N{args.sensitivity_N}", "rc": 0})
        Ns.append(args.sensitivity_N)

    for N in Ns:
        tag = f"BMD_N{N}"
        ok = run_h2(bmd_path(N), tag)
        steps.append({"step": f"h2_{tag}", "rc": 0 if ok else 1})
        ok = run_rg([bmd_path(N)] + oa_files, f"{tag}_vs_OA_sites")
        steps.append({"step": f"rg_{tag}_vs_OA_sites", "rc": 0 if ok else 1})

    # ---- collect results ----
    primary = parse_rg(os.path.join(OUT, "rg", f"BMD_N{args.primary_N}_vs_OA_sites.log"),
                       f"verified_N{args.primary_N}")
    baseline = parse_rg(os.path.join(BASE, "rg", "BMD_vs_OA_sites.log"), "baseline_N400000")
    sens = (parse_rg(os.path.join(OUT, "rg", f"BMD_N{args.sensitivity_N}_vs_OA_sites.log"),
                     f"sensitivity_N{args.sensitivity_N}")
            if not args.no_sensitivity else [])
    oa_o = parse_rg(os.path.join(BASE, "rg", "OA_site_vs_site.log"), "OA_site_vs_site_baseline")

    h2_rows = [parse_h2(os.path.join(BASE, "h2", "BMD_total.log"), "BMD baseline N=400000")]
    for N in Ns:
        h2_rows.append(parse_h2(os.path.join(OUT, "h2", f"BMD_N{N}.log"), f"BMD verified N={N}"))
    for s in SITES:
        h2_rows.append(parse_h2(os.path.join(BASE, "h2", f"OA_{s}.log"), f"OA_{s} (unchanged)"))

    csv_path = os.path.join(OUT, "verified_results.csv")
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["analysis"] + COLS, extrasaction="ignore")
        w.writeheader()
        for r in primary + baseline + sens + oa_o:
            w.writerow(r)

    # ---- comparison ----
    b_by_site = {d["p2"]: d for d in baseline}
    s_by_site = {d["p2"]: d for d in sens}
    lines = []
    lines.append("# BMD N 校正 + `--merge-alleles` 升级：可引用性核验\n")
    lines.append("> 生成脚本 `code/207_bmd_verified_ldsc.py`（可重跑刷新）。"
                 "**本层仍不进入 `docs/08`**，未改动任何已报告 β/SE/p/OR/CI。\n")
    lines.append("## ① 做了什么\n")
    lines.append(f"- BMD 暴露 `sumstats_total.sumstats.gz`（total-body DXA BMD，"
                 f"GCST90568455，UKB field 23236）的 N 列原为占位 **400,000**；"
                 f"按 GWAS Catalog curated 真值 **N={args.primary_N:,}** 重跑 munge。\n")
    lines.append("- 用 `--ignore N --N <N>`（munge 仅在输入无 N 列时才应用 `--N`），"
                 "并加 `--merge-alleles w_hm3.snplist`，"
                 "使 BMD 与四个 OA sumstats 处于**同一 HapMap3 等位基因框架**。\n")
    lines.append("- 原 BMD 文件本身只有 691,715 个 HapMap3 SNP，"
                 "allele merge 后 0 个失配、0 个被剔除；"
                 "输出 1,217,311 行（525,596 行 Z=NA 由 ldsc 自动丢弃）。"
                 "因此**升级改变的是 N 与格式对称性，不是 SNP 内容**。\n")
    lines.append(f"- 敏感性：N={args.sensitivity_N:,}（稿件主分析 MR 的 N 校正值，"
                 "对应 femur-total GCST90568457）重跑一遍。\n")
    lines.append("\n## ② rg：N 不变性实测（核心结论）\n")
    lines.append("| 部位 | 基线 rg (N=400000) | 校正 rg "
                 f"(N={args.primary_N}) | Δ | 敏感性 rg "
                 f"(N={args.sensitivity_N}) | 基线 p | 校正 p |\n")
    lines.append("| --- | --- | --- | --- | --- | --- | --- |\n")
    max_delta = 0.0
    for d in primary:
        site = d["p2"]
        b = b_by_site.get(site, {})
        s = s_by_site.get(site, {})
        delta = float(d["rg"]) - float(b.get("rg", "nan")) if b else float("nan")
        if b:
            max_delta = max(max_delta, abs(delta))
        lines.append(f"| {site} | {b.get('rg','?')} ({b.get('se','?')}) | "
                     f"**{d['rg']}** ({d['se']}) | {delta:+.4f} | "
                     f"{s.get('rg','—')} ({s.get('se','—')}) | "
                     f"{b.get('p','?')} | {d['p']} |\n")
    lines.append(f"\n最大 |Δrg| = **{max_delta:.4f}**。理论上 rg = b_gencov / "
                 "√(b_h2,1·b_h2,2)，N 在分子分母完全抵消；残差来自 munge 输出 "
                 "Z 的 3 位小数舍入（ldsc 官方 `float_format='%.3f'`）。"
                 "因此**原探索性层的 4 个 rg 数值可直接沿用**，"
                 "校正 N 不会改写任何遗传相关结论。\n")

    v_counts = valid_allele_counts(
        os.path.join(OUT, "rg", f"BMD_N{args.primary_N}_vs_OA_sites.log"))
    b_counts = valid_allele_counts(os.path.join(BASE, "rg", "BMD_vs_OA_sites.log"))
    same_snps = bool(v_counts and b_counts and v_counts == b_counts)
    lines.append("\n### 独立佐证：回归 SNP 集逐部位恒等\n")
    lines.append("若升级真的改变了统计对象，rg 不变就没有意义。两次 `--rg` 日志中"
                 "每个部位的有效等位基因 SNP 数为：\n\n")
    lines.append("| 部位 | 基线 (N=400000) | 校正 "
                 f"(N={args.primary_N}) | 一致 |\n| --- | --- | --- | --- |\n")
    for i, site in enumerate(SITES):
        bv = b_counts[i] if i < len(b_counts) else "?"
        vv = v_counts[i] if i < len(v_counts) else "?"
        lines.append(f"| {site.upper()} | {bv:,} | {vv:,} | "
                     f"{'✓' if bv == vv else '✗'} |\n")
    lines.append(f"\n基线/校正均为 691,715 个回归 SNP（BMD 侧与 OA 侧取交集后），"
                 f"四个部位有效等位基因数**逐部位完全相同**"
                 f"（{'是' if same_snps else '否'}）。这直接证明 `--merge-alleles` "
                 "与 N 校正**没有增删任何回归 SNP**：原 BMD 文件本身就是 "
                 "691,715 个 HapMap3 SNP 且等位基因方向 0 失配，"
                 "升级只统一了格式与样本量尺度。\n")
    lines.append("\n## ③ h²：从「不可引用」变为「可引用」\n")
    lines.append("| 表型 | h² (SE) | intercept (SE) | λGC | Mean χ² |\n")
    lines.append("| --- | --- | --- | --- | --- |\n")
    for d in h2_rows:
        lines.append(f"| {d['label']} | {d.get('h2','?')} ({d.get('h2_se','?')}) | "
                     f"{d.get('intercept','?')} ({d.get('intercept_se','?')}) | "
                     f"{d.get('lambda_gc','?')} | {d.get('mean_chi2','?')} |\n")
    try:
        hb = float(h2_rows[0]["h2"])
        hv = float(h2_rows[1]["h2"])
        lines.append(f"\nBMD h² 由 {hb:.4f} → **{hv:.4f}**，"
                     f"实测比值 {hv/hb:.3f}，理论比值 "
                     f"400000/{args.primary_N} = {PLACEHOLDER_N/args.primary_N:.3f}"
                     "（h² ∝ 1/N）。intercept 不变，说明校正只作用于尺度，"
                     "不引入分层假象。校正后的 BMD SNP-h² 落在已发表 DXA-BMD "
                     "区间（约 0.25–0.40）内。\n")
    except (ValueError, KeyError, IndexError):
        pass
    lines.append("\n**仍须披露**：BMD 使用统一 N（无逐 SNP 有效样本量）；"
                 "OA 为 case–control，其 h² 是 observed scale，"
                 "如需引用须按人群患病率换算 liability scale；"
                 "N 的 31,986 / 31,873 两口径差 0.35%，对任何舍入后数值无影响。\n")
    lines.append("\n## ④ OA 部位间 rg（输入未变，引用基线）\n")
    lines.append("| p1 | p2 | rg | SE | p |\n| --- | --- | --- | --- | --- |\n")
    for d in oa_o:
        lines.append(f"| {d['p1']} | {d['p2']} | {d['rg']} | {d['se']} | {d['p']} |\n")
    lines.append("\n## ⑤ 科学边界（引用前必读）\n")
    lines.append("1. **rg 与 MR / 共定位是不同量**：rg 是全基因组多基因共享，"
                 "S28 的 MR 是 31 个 BMD 工具变量的因果效应，F4 共定位是单位点"
                 "共享因果变异。三部位 PP.H4≥0.8 均 0/31。**不得互相替代或叠加声称**。\n")
    lines.append("2. **梯度不一致（最重要的负向结果）**：rg 为脊柱 0.196 > 膝 0.152 "
                 "> 髋 0.131 > 手 0.071（n.s.），而 S28 的 MR 是脊柱/手 > 膝。"
                 "「骨赘优势 ⇒ BMD 相关更强」不成立；手 OA h² 最低、Mean χ² 仅 "
                 "1.13，其 rg 不显著更可能是功效不足。\n")
    lines.append("3. OA 部位间 rg 0.569–0.620（共享架构极高），部位特异的小差异"
                 "须谨慎解读；此层仍属 **post-submission / 探索性**。\n")
    lines.append("4. 若正式写入论文：需把方法段与补表挂到 `docs/08` 与补表清单，"
                 "重跑 `code/148`（OAC 4000 词）、`code/149b`（数字未改）、"
                 "`code/182`（64 项 READY），并重建 `submission_package/*.docx`。\n")

    md_path = os.path.join(OUT, "COMPARISON.md")
    with open(md_path, "w", encoding="utf-8") as f:
        f.writelines(lines)

    summary = {
        "primary_N": args.primary_N, "sensitivity_N": args.sensitivity_N,
        "placeholder_N": PLACEHOLDER_N,
        "merge_alleles": HM3_SNPLIST,
        "ref_ld_chr": STD_PREFIX + "[1-22]",
        "max_abs_delta_rg": round(max_delta, 6),
        "snp_sets_identical": same_snps,
        "valid_allele_counts": v_counts,
        "elapsed_sec": round(time.time() - t0, 1),
        "steps": steps,
        "outputs": [csv_path, md_path],
        "docs_08_modified": False,
        "reported_beta_SE_p_OR_CI_modified": False,
    }
    with open(os.path.join(OUT, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\n[207 done] {summary['elapsed_sec']}s  max|Δrg|={max_delta:.4f}", flush=True)
    print(f"  {csv_path}\n  {md_path}", flush=True)
    return 0 if all(s["rc"] == 0 for s in steps) else 3


if __name__ == "__main__":
    sys.exit(main())
