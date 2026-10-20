# -*- coding: utf-8 -*-
"""
132 · 第十二轮整改 P0-2：显示条目 15 → 8
------------------------------------------------------------------
OAC 官方限制："Full length original articles should amount to no more than 4000 words,
8 figures and tables, and 50 references." 当前主文 10 图 + 5 表 = 15。

处置：
  主文保留 5 图 + 3 表 = 8
    图 1 (r_g 热图)      ← 原 Figure 1
    图 2 (MR 森林图)     ← 原 Figure 3
    图 3 (概念架构)      ← 原 Figure 5
    图 4 (散点图)        ← 原 Figure 6
    图 5 (MRlap 校正)    ← 原 Figure 10
    表 1 (h²/截距)、表 3 (MR)、表 4 (功效)
  移入补充（重编号，避免与既有 S8–S13 冲突）
    原 Figure 2 → Supplementary Figure S13（r_g 森林图）
    原 Figure 4 → S14（h² 条形图）
    原 Figure 7 → S15（PRESSO 漏斗图）
    原 Figure 8 → S16（功效）
    原 Figure 9 → S17（ρ_GE 重叠诊断）
    原 Table 1b → Supplementary Table S14（祖源构成）
    原 Table 2  → Supplementary Table S15（分块 r_g）
同时：正文交叉引用、图注段、补充材料清单、图文件（磁盘）一并同步。
"""
import io
import os
import re
import shutil

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MD = os.path.join(ROOT, "docs", "08_投稿稿_OA_BMD_fracture.md")
FIGS = os.path.join(ROOT, "Figures")
SUPL = os.path.join(ROOT, "Figures_v18_Suppl")
MR_MAIN = os.path.join(ROOT, "Manuscript_Ready", "Figures", "Main")
MR_SUPL = os.path.join(ROOT, "Manuscript_Ready", "Figures", "Supplementary")

t = io.open(MD, encoding="utf-8").read()
L = t.split("\n")
n = 0


def say(s=""):
    print(s)


# ---------------------------------------------------------------- 1. 正文交叉引用
ref_edits = [
    ("(as does the SNP-effect correlation ρ_GE, Fig. 9)",
     "(as does the SNP-effect correlation ρ_GE, Fig. S17)"),
    ("(Figure 1, Table 2)", "(Figure 1, Table S15)"),
]
for old, new in ref_edits:
    c = t.count(old)
    assert c == 1, f"交叉引用 {old!r} 命中 {c} 次"
    t = t.replace(old, new)
    n += 1
    say(f"  OK  交叉引用 {old} → {new}")
t = t.replace("## Table 1b. Population-ancestry composition of the 13 GWAS",
              "## Supplementary Table S14. Population-ancestry composition of the 13 GWAS")
c = t.count("Table 1b")
assert c == 4, f"'Table 1b' 正文引用命中 {c} 次（应为 4）"
t = t.replace("Table 1b", "Table S14")
n += 1
say("  OK  表头重命名 + 4 处正文引用 Table 1b → Table S14")

# ---------------------------------------------------------------- 2. 图注段
def legend_of(num):
    """抽取 '**Figure <num>.** ...' 整行（可能跨行，按以 '**Figure' 或 '**Supplementary' 或 '---' 结束）。"""
    i = next(k for k, ln in enumerate(L) if ln.startswith(f"**Figure {num}.**"))
    j = i + 1
    while j < len(L) and not L[j].startswith("**Figure") and not L[j].startswith(
            "**Supplementary") and L[j].strip() != "---":
        j += 1
    return "\n".join(L[i:j])


leg = {k: legend_of(k) for k in (1, 2, 3, 4, 5, 6, 7, 8, 9, 10)}


def renumber_legend(txt, newnum):
    return re.sub(r"^\*\*Figure \d+\.\*\*", f"**Figure {newnum}.**", txt, count=1)


def renumber_supp(txt, newnum):
    return re.sub(r"^\*\*Figure \d+\.\*\*", f"**Supplementary Figure S{newnum}.**", txt, count=1)


main_legends = "\n\n".join([
    renumber_legend(leg[1], 1),
    renumber_legend(leg[3], 2),
    renumber_legend(leg[5], 3),
    renumber_legend(leg[6], 4),
    renumber_legend(leg[10], 5),
])

# 既有补充图行（S1–S12 + Graphical Abstract）
supp_line = next(ln for ln in L if ln.startswith("**Supplementary Figure S1.**"))
assert "**Graphical Abstract.**" in supp_line
moved = " ".join([
    renumber_supp(leg[2], 13),
    renumber_supp(leg[4], 14),
    renumber_supp(leg[7], 15),
    renumber_supp(leg[8], 16),
    renumber_supp(leg[9], 17),
])
supp_line_new = supp_line.replace("**Graphical Abstract.**", moved + " **Graphical Abstract.**")

i_start = next(k for k, ln in enumerate(L) if ln.strip() == "## Figure legends")
i_end = next(k for k, ln in enumerate(L) if ln.startswith("## Table 1."))
new_block = ("## Figure legends\n\n" + main_legends + "\n\n" + supp_line_new + "\n\n")
L = L[:i_start] + new_block.split("\n") + L[i_end:]
t = "\n".join(L)
n += 1
say("  OK  图注段重写：主文 5 图 + 补充 S1–S12 与新增 S13–S17")

# ---------------------------------------------------------------- 3. 表 1b / 表 2 移入补充
L = t.split("\n")
i_t1b = next(k for k, ln in enumerate(L) if ln.startswith("## Table 1b."))
i_t2 = next(k for k, ln in enumerate(L) if ln.startswith("## Table 2."))
i_t3 = next(k for k, ln in enumerate(L) if ln.startswith("## Table 3."))
blk_t1b = "\n".join(L[i_t1b:i_t2])
blk_t2 = "\n".join(L[i_t2:i_t3])
blk_t1b_new = blk_t1b.replace("## Table 1b. Population-ancestry composition of the 13 GWAS",
                              "## Supplementary Table S14. Population-ancestry composition of "
                              "the 13 GWAS")
blk_t2_new = blk_t2.replace("## Table 2. Genetic correlations by phenotype block",
                            "## Supplementary Table S15. Genetic correlations by phenotype block")
L = L[:i_t1b] + L[i_t3:]
t = "\n".join(L)
n += 2
say("  OK  Table 1b → Supplementary Table S14；Table 2 → Supplementary Table S15（移出主文）")

# 插入到补充表区：Table S13 块之后、"---" 分隔线之前
L = t.split("\n")
i_s9 = next(k for k, ln in enumerate(L) if ln.startswith("## Supplementary Table S9."))
i_sep = next(k for k in range(i_s9, len(L)) if L[k].strip() == "---")
L = L[:i_sep] + ["", blk_t1b_new, "", blk_t2_new, ""] + L[i_sep:]
t = "\n".join(L)
n += 1
say("  OK  S14/S15 插入补充表区")

# ---------------------------------------------------------------- 4. 补充材料清单
L = t.split("\n")
i_s12b = next(k for k, ln in enumerate(L) if ln.startswith("- **Figure S12.**"))
add_figs = [
    "- **Figure S13.** Forest plot of the selected genetic correlations (any-fracture, knee OA "
    "and self-reported OA versus the ten BMD sites; the two OA definitions with each other and "
    "with fracture).",
    "- **Figure S14.** SNP-based heritability (observed scale) of the 13 skeletal phenotypes, "
    "with liability-scale h² for the three binary traits.",
    "- **Figure S15.** MR-PRESSO funnel plots for the four headline pairs.",
    "- **Figure S16.** Statistical power of the two-sample MR design.",
    "- **Figure S17.** Sample-overlap diagnostics (ρ_GE) for the ten MR pairs.",
]
i_s13b = next(k for k, ln in enumerate(L) if ln.startswith("- **Table S13.**"))
add_tabs = [
    "- **Table S14.** Population-ancestry composition of the 13 GWAS and the LDSC reference "
    "used for each.",
    "- **Table S15.** Genetic correlations by phenotype block (n pairs, r_g range, r_g mean, "
    "proportion significant at q<0.05).",
]
L = L[:i_s12b + 1] + add_figs + L[i_s12b + 1:]
i_s13b = next(k for k, ln in enumerate(L) if ln.startswith("- **Table S13.**"))
L = L[:i_s13b + 1] + add_tabs + L[i_s13b + 1:]
t = "\n".join(L)
n += 2
say("  OK  Supplementary Material 清单补 S13–S17 与 Table S14/S15")

io.open(MD, "w", encoding="utf-8").write(t)

# ---------------------------------------------------------------- 5. 磁盘文件
def move_pair(src_stem, dst_stem, dst_dir, src_dir=FIGS):
    for ext in (".png", ".pdf"):
        s = os.path.join(src_dir, src_stem + ext)
        if not os.path.exists(s):
            continue
        d = os.path.join(dst_dir, dst_stem + ext)
        os.makedirs(dst_dir, exist_ok=True)
        shutil.copy2(s, d)
        os.remove(s)
        say(f"  OK  {src_stem}{ext} → {os.path.relpath(d, ROOT)}")


say()
say("---- 磁盘图文件 ----")
# 先搬走（避免与主图重编号冲突）
for src, dst in [("Figure2_forest_key_results", "SupplFig13_rg_forest"),
                 ("Figure4_h2_bar", "SupplFig14_h2_bar"),
                 ("Figure7_MR_funnel", "SupplFig15_presso_funnel"),
                 ("Figure8_MR_power", "SupplFig16_mr_power"),
                 ("Figure9_overlap_pleiotropy", "SupplFig17_overlap_rhoGE")]:
    move_pair(src, dst, SUPL)
    move_pair(src, dst, MR_SUPL, src_dir=MR_MAIN)
# 再重编号主图
for src, dst in [("Figure3_MR_forest", "Figure2_MR_forest"),
                 ("Figure5_causal_network", "Figure3_causal_network"),
                 ("Figure6_MR_scatter", "Figure4_MR_scatter"),
                 ("Figure10_MRlap_overlap_correction", "Figure5_MRlap_overlap_correction")]:
    move_pair(src, dst, FIGS)
    move_pair(src, dst, MR_MAIN, src_dir=MR_MAIN)

say()
say(f"完成：{n} 处文本修改 + 磁盘图文件重编号")
say("主文显示条目 = 5 图 + 3 表 = 8（OAC 上限 8）")
