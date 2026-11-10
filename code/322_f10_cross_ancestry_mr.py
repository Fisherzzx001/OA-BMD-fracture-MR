# -*- coding: utf-8 -*-
"""
F10 跨祖先 MR · 真实执行（EUR 工具 → 非欧 OA 结局）

设计（诚实边界见底部）
----------------------
登记册 F10 点名 "Hatzikotoulas 2025 Nature 多祖先" 的 EAS 分层膝 OA。复核结论：
  * Hatzikotoulas 2025 只发布了合并多祖先文件（GCST90566800 含 132k EAS 但 1.12M 欧裔主导，
    非干净 EAS 检验）；EAS 分层条目 GCST90134281 仅 n=1704，跑不动 31 工具。
  * Liu TY 2025 干净 EAS OA（GCST90652199, n=221585）**未发布 sumstats**。
因此本脚本执行的是「可防御的跨祖先敏感性」而非「登记册点名的 de novo EAS 复制」：
    暴露 = 稿件锁定的 31 个欧洲 total-body BMD 工具（bmd_instruments.json）
    结局 = 非欧 OA 汇总统计（默认 Genes & Health 南亚 OA excl. spine, GCST90727231, n=44026）
方法 = 按 rsID 匹配（build 无关）→ 等位对齐 → IVW / MR-Egger / LOO（复用 code/13）。
判读 = 若 BMD→OA 方向在非欧人群仍为正且显著 → SOS 方向具跨祖先稳健性（支持，非证明）；
       若反转 → 祖先特异性生物学可能推翻 SOS。
局限 = 工具为欧裔定型（EUR LD 定簇），南亚 LD 较接近欧裔故较有效；结局为 all-site OA excl.
       spine 而非膝 OA 特异；N=44026 中等，功效有限。这些均在报告中明示。

输入：本地 harmonised 文件（默认 G:/.../GCST90727231.h.tsv.gz）。先由 aux_download_gandh.py 下载。
"""
import os
import sys
import json
import gzip
import math
import importlib.util
import numpy as np
from scipy.stats import norm, chi2

ROOT = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(ROOT)
CACHE = "G:/OA_B_PROJECT_cache/l7"
INST = os.path.join(CACHE, "bmd_instruments.json")
DEFAULT_OUTFILE = os.path.join(ROOT, "output", "genetics",
                               "L7_F10_cross_ancestry_mr.json")

# 复用 code/13 的 MR 引擎（与正文/code/314 完全一致）
_spec = importlib.util.spec_from_file_location(
    "mr13", os.path.join(ROOT, "13_mr_analysis.py"))
mr13 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mr13)


def load_instruments():
    d = json.load(open(INST, encoding="utf-8"))["inst"]
    return d  # keyed by rsid


def stream_match_by_rsid(path, rsids):
    """Stream a (possibly large) harmonised .tsv.gz, return rows keyed by rsid."""
    wanted = set(rsids)
    rows = {}
    header = None
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if header is None:
                header = line.rstrip("\n").split("\t")
                idx = {c: i for i, c in enumerate(header)}
                continue
            parts = line.rstrip("\n").split("\t")
            rec = dict(zip(header, parts))
            rs = rec.get("rsid")
            if rs in wanted:
                rows[rs] = rec
                if len(rows) == len(wanted):
                    break
    return rows, header


def harmonize_and_estimate(inst, matched, outcome_label, outcome_n):
    rows = []
    for rs, s in inst.items():
        rec = matched.get(rs)
        if rec is None:
            continue
        exp = str(s["exp_allele"]).lower()
        try:
            ea = str(rec["effect_allele"]).lower()
            oa = str(rec["other_allele"]).lower()
            beta_o = float(rec["beta"])
            se_o = float(rec["standard_error"])
        except (KeyError, TypeError, ValueError):
            continue
        if ea == exp:
            bo = beta_o
        elif oa == exp:
            bo = -beta_o
        else:
            continue  # 等位/链不确定 → 排除该 SNP
        be = float(s["beta_e"])
        se_e = float(s["se_e"])
        if se_e <= 0 or se_o <= 0:
            continue
        rows.append(dict(rsid=rs, beta_e=be, se_e=se_e, beta_o=bo, se_o=se_o,
                         eaf=float(rec.get("effect_allele_frequency", "nan"))
                         if rec.get("effect_allele_frequency") not in (None, "nan", "")
                         else float("nan")))
    if len(rows) < 2:
        return None
    d = {k: np.array([r[k] for r in rows]) for k in
         ("rsid", "beta_e", "se_e", "beta_o", "se_o", "eaf")}
    # IVW (regression, no intercept), weights = 1/se_o^2
    w = 1.0 / d["se_o"] ** 2
    ivw_beta = float(np.sum(w * d["beta_e"] * d["beta_o"]) / np.sum(w * d["beta_e"] ** 2))
    ivw_se = float(np.sqrt(1.0 / np.sum(w * d["beta_e"] ** 2)))
    z = ivw_beta / ivw_se
    p = 2 * norm.sf(abs(z))
    or_sd = math.exp(ivw_beta)
    lo = math.exp(ivw_beta - 1.96 * ivw_se)
    hi = math.exp(ivw_beta + 1.96 * ivw_se)
    # Cochran Q
    pred = ivw_beta * d["beta_e"]
    q = float(np.sum(w * (d["beta_o"] - pred) ** 2))
    k = len(rows)
    q_p = 1 - chi2.cdf(q, max(k - 1, 1))
    # MR-Egger
    eg = mr13.run_mr(d["beta_e"], d["se_e"], d["beta_o"], d["se_o"],
                     snps=d["rsid"].tolist()) if hasattr(mr13, "run_mr") else None
    res = dict(
        outcome_label=outcome_label, outcome_n=outcome_n, n_snps=k,
        ivw_beta=ivw_beta, ivw_se=ivw_se, z=z, p=p,
        or_sd=or_sd, ci_lo=lo, ci_hi=hi, q=q, q_df=k - 1, q_p=q_p,
        egger=eg,  # full run_mr result dict (or None)
        snps=[dict(rsid=r["rsid"], beta_e=r["beta_e"], beta_o=r["beta_o"],
                   se_o=r["se_o"], eaf=r["eaf"]) for r in rows],
    )
    return res


def main():
    # ---- 目标结局（默认 Genes & Health 南亚 OA；可改 GCST90134287 非洲等）----
    target = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        CACHE, "GCST90727231.h.tsv.gz")
    outcome_label = (sys.argv[2] if len(sys.argv) > 2 else
                    "Genes & Health British South Asian OA excl. spine "
                    "(GCST90727231, n=44026)")
    outcome_n = int(sys.argv[3]) if len(sys.argv) > 3 else 44026
    if not os.path.exists(target):
        print(f"[ABORT] target not found: {target}\n"
              f"  download via code/aux_download_gandh.py first")
        return
    inst = load_instruments()
    print(f"[load] {len(inst)} BMD instruments (by rsid)")
    matched, header = stream_match_by_rsid(target, list(inst.keys()))
    print(f"[match] {len(matched)}/{len(inst)} rsids found in outcome")
    out = harmonize_and_estimate(inst, matched, outcome_label, outcome_n)
    os.makedirs(os.path.dirname(DEFAULT_OUTFILE), exist_ok=True)
    if out is None:
        print(f"[FAIL] only {len(matched)}/{len(inst)} rsids matched with "
              f"unambiguous alleles — insufficient for MR")
        json.dump(dict(
            design="F10 cross-ancestry SENSITIVITY feasibility check "
                   "(European BMD instruments -> non-European OA outcome)",
            outcome=outcome_label, outcome_n=outcome_n,
            n_instruments=len(inst), n_rsid_matched=len(matched),
            n_harmonised=0,
            verdict="UNINFORMATIVE: exome-wide sequencing does not capture "
                    "the common regulatory variants used as European BMD "
                    "instruments (1/31 matched). No MR estimate produced; "
                    "does NOT close F10 and must NOT be reported as a result."),
            open(DEFAULT_OUTFILE, "w", encoding="utf-8"),
            ensure_ascii=False, indent=2)
        print(f"[write] {DEFAULT_OUTFILE} (negative feasibility record)")
        return
    print(f"\n=== F10 cross-ancestry MR · {out['outcome_label']} ===")
    print(f"matched SNPs: {out['n_snps']}/{len(inst)}")
    print(f"IVW beta(SE)= {out['ivw_beta']:.5f} ({out['ivw_se']:.5f})  "
          f"OR/SD= {out['or_sd']:.4f} ({out['ci_lo']:.4f}–{out['ci_hi']:.4f})  "
          f"P= {out['p']:.2e}")
    print(f"Cochran Q= {out['q']:.2f} on {out['q_df']} df (P= {out['q_p']:.3f})")
    if out["egger"]:
        eg = out["egger"]
        print(f"Egger slope= {eg.get('egger'):.5f} (SE {eg.get('egger_se'):.5f})")
        print(f"Egger intercept= {eg.get('egger_intercept'):.5f} "
              f"(P= {eg.get('egger_intercept_p'):.3f})")
    json.dump(out, open(DEFAULT_OUTFILE, "w", encoding="utf-8"),
              ensure_ascii=False, indent=2)
    print(f"[write] {DEFAULT_OUTFILE}")


if __name__ == "__main__":
    main()
