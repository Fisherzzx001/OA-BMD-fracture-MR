#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
206_summarise_ldsc_rg.py -- parse code/205 full-run ldsc logs into a table.
==========================================================================
Reads rg/*.log (Summary of Genetic Correlation Results) and h2/*.log, writes
  code/output/genetics/ldsc_rg/rg_results.csv
  code/output/genetics/ldsc_rg/RESULTS.md

Pure stdlib. Re-run after any 205 rerun to refresh the tables.
NOTE: these numbers are EXPLORATORY / POST-SUBMISSION (see RESULTS.md header).
"""
import csv
import os
import re

HERE = "C:/Users/Administrator/WorkBuddy/2026-09-08-22-17-40/OA_B_PROJECT"
BASE = os.path.join(HERE, "code", "output", "genetics", "ldsc_rg")
RG_LOGS = [
    ("BMD_vs_OA_sites", os.path.join(BASE, "rg", "BMD_vs_OA_sites.log")),
    ("OA_site_vs_site", os.path.join(BASE, "rg", "OA_site_vs_site.log")),
]
H2_NAMES = ["BMD_total", "OA_hip", "OA_knee", "OA_spine", "OA_hand"]
COLS = ["p1", "p2", "rg", "se", "z", "p", "h2_obs", "h2_obs_se",
        "h2_int", "h2_int_se", "gcov_int", "gcov_int_se"]


def short(path):
    b = os.path.basename(path.replace("\\", "/"))
    b = b.replace(".sumstats.gz", "").replace(".gz", "")
    if b.startswith("OA_"):
        return b[3:].upper()
    if b.startswith("sumstats_"):
        return "BMD_total"
    return b


def parse_rg(path, analysis):
    rows = []
    with open(path, encoding="utf-8", errors="replace") as f:
        lines = f.read().splitlines()
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
        d["p1"] = short(d["p1"])
        d["p2"] = short(d["p2"])
        rows.append(d)
    return rows


def parse_h2(path, name):
    d = {"phenotype": name}
    txt = open(path, encoding="utf-8", errors="replace").read()
    patterns = {
        "h2": r"Total Observed scale h2:\s*([\-0-9.eE+]+)\s*\(([\-0-9.eE+]+)\)",
        "intercept": r"Intercept:\s*([\-0-9.eE+]+)\s*\(([\-0-9.eE+]+)\)",
        "lambda_gc": r"Lambda GC:\s*([\-0-9.eE+]+)",
        "mean_chi2": r"Mean Chi\^2:\s*([\-0-9.eE+]+)",
    }
    for k, pat in patterns.items():
        m = re.search(pat, txt)
        if m:
            d[k] = m.group(1)
            if m.lastindex and m.lastindex > 1:
                d[k + "_se"] = m.group(2)
    return d


def fmt_p(p):
    try:
        v = float(p)
    except ValueError:
        return p
    return f"{v:.2e}" if v < 1e-3 else f"{v:.4f}"


def main():
    rg = []
    for tag, path in RG_LOGS:
        if os.path.exists(path):
            rg += parse_rg(path, tag)
    h2 = [parse_h2(os.path.join(BASE, "h2", n + ".log"), n)
          for n in H2_NAMES if os.path.exists(os.path.join(BASE, "h2", n + ".log"))]

    csv_path = os.path.join(BASE, "rg_results.csv")
    fieldnames = ["analysis"] + COLS
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(rg)

    out = []
    out.append("# 跨性状 LDSC rg 全量结果（探索性 / 投稿后层）\n")
    out.append("> **本文件不进入稿件**。`docs/08` 是唯一文本真源，本层未改动任何已报告 "
               "β/SE/p/OR/CI。生成脚本：`code/206_summarise_ldsc_rg.py`，"
               "可重跑刷新。\n")
    out.append("## 配置（发表级 HapMap3）\n")
    out.append("- `--ref-ld-chr` / `--w-ld-chr` = Broad "
               "`eur_w_ld_chr_1..22`（1,290,028 SNP，M=1,290,050 / "
               "M_5_50=1,173,569，HLA 已移除）\n")
    out.append("- munge 加 `--merge-alleles w_hm3.snplist`（1,217,311 SNP，"
               "SNP/A1/A2），四个部位各保留约 1,176,000 个有效 HapMap3 SNP\n")
    out.append("- OA GWAS = Hatzikotoulas 2025 部位分层（髋/膝/脊柱/手），"
               "GRCh37 → rsID 映射；N 逐行 = NCASES+NCONTROLS，`--n-min 1`\n")
    out.append("- 全量运行 1198.6 s，13 个步骤 rc 全 0；BMD×OA 回归实际交集 "
               "691,715 SNP\n")
    out.append("\n## 表 1 · BMD × 部位分层 OA 的遗传相关\n")
    out.append("| 分析 | p1 | p2 | rg | SE | z | p | h2_obs(p2) | gencov intercept |\n")
    out.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- |\n")
    for d in rg:
        out.append(f"| {d['analysis']} | {d['p1']} | {d['p2']} | {d['rg']} | "
                   f"{d['se']} | {d['z']} | {fmt_p(d['p'])} | {d['h2_obs']} "
                   f"({d['h2_obs_se']}) | {d['gcov_int']} ({d['gcov_int_se']}) |\n")
    out.append("\n## 表 2 · 单性状 h² 与 intercept（**h² 不可引用**）\n")
    out.append("| 表型 | h² (SE) | LDSC intercept (SE) | λGC | Mean χ² |\n")
    out.append("| --- | --- | --- | --- | --- |\n")
    for d in h2:
        out.append(f"| {d['phenotype']} | {d.get('h2','?')} ({d.get('h2_se','?')}) | "
                   f"{d.get('intercept','?')} ({d.get('intercept_se','?')}) | "
                   f"{d.get('lambda_gc','?')} | {d.get('mean_chi2','?')} |\n")
    out.append("\n## 解读（探索性，不得与 MR 结果互相替代）\n")
    out.append("- **BMD × OA 的 rg 在 4 个部位中 3 个显著为正**：脊柱 0.196 "
               "(p=2.9×10⁻⁵) > 膝 0.152 (p=3.5×10⁻⁶) > 髋 0.131 "
               "(p=5.3×10⁻⁴) > 手 0.071 (p=0.204，不显著)。方向一致支持"
               "「骨形成程序与 OA 存在共享的多基因基础」，但**量级都很小**。\n")
    out.append("- **梯度与 MR 不一致，这是本轮最重要的负向信息**：S28 的部位分层 MR 里"
               "骨赘优势部位（脊柱 OR 1.129、手 OR 1.124）数值上高于软骨优势的膝（1.081）；"
               "但 rg 里脊柱最高、手最低且不显著、髋居中。因此「骨赘优势 ⇒ BMD 相关更强」"
               "**并不成立**。两者测的是不同量：rg 是全基因组多基因共享，MR 是 31 个 BMD "
               "工具变量的因果效应。这一不一致限制了用部位分层替真 F5/T4 的解释力。\n")
    out.append("- **OA 部位间 rg 极高（0.569–0.620，均 p<10⁻⁵⁰）**：四个部位的 OA "
               "遗传架构高度共享，与 F4 共定位三部位 PP.H4≥0.8 均 0/31 所提示的"
               "「弥散多基因」一致；共享背景会稀释部位特异信号，小的部位间差异需谨慎解读。\n")
    out.append("- 手 OA 的 h² 最低（诊断值 0.0087）、Mean χ² 仅 1.13，其 rg 不显著"
               "更可能是**功效不足**而非真实的生物学例外。\n")
    out.append("\n## 诚实边界（引用任何数值前必读）\n")
    out.append("1. **BMD 暴露 sumstats 的 N 列是占位值 400000**（全文件唯一取值）；"
               "`docs/07` 推算真值约 3.2–3.5 万（共定位用 31,986）。"
               "**rg 对 N 全局缩放不变**（分子分母抵消），但 **h² 绝对值与 "
               "intercept 的量级不可引用**。表 2 仅供诊断。\n")
    out.append("2. BMD sumstats 是旧资产，**未经过 `--merge-alleles`**；"
               "`--w-ld-chr` 已把回归限制到 HapMap3，跨性状 A1/A2 由 ldsc 自动对齐，"
               "但这一不对称须在方法段披露。\n")
    out.append("3. 此层是 **post-submission / 探索性**：与 F4 共定位（三部位 PP.H4≥0.8 均 0/31）"
               "及 S28 的 MR 是不同量（全基因组多基因相关 vs 工具变量因果效应），"
               "不得互相替代或叠加声称。\n")
    out.append("4. **2026-11-04 更新：第 1、2 条已解决**。`code/207_bmd_verified_ldsc.py` "
               "已核实 BMD 真实 N = 31,986（GCST90568455, total-body DXA）、"
               "对 BMD 重跑 `--merge-alleles w_hm3.snplist`，h² 校正为 0.3479 (0.0355) "
               "并变为可引用，rg 实测不变（max |Δ| = 0.0002；N=31,873 敏感性一致）。"
               "可提交方法段与结果表见 `docs/Supplementary/Table_S29_ldsc_site_rg.md`；"
               "完整前后对比见 `ldsc_rg_verified/COMPARISON.md`。"
               "本目录**保留为探索性基线**作前后对照：表 2 的 h² 仍不可引用。"
               "是否纳入正文需先做等词压缩（OAC 3999/4000，仅剩 1 词），"
               "再复跑 148 / 149b / 182 三连门禁。\n")

    md_path = os.path.join(BASE, "RESULTS.md")
    with open(md_path, "w", encoding="utf-8") as f:
        f.writelines(out)
    print(f"[206] rg rows={len(rg)} h2 rows={len(h2)}")
    print(f"      {csv_path}")
    print(f"      {md_path}")


if __name__ == "__main__":
    main()
