# -*- coding: utf-8 -*-
"""预审委员会 · 事实核查（供 R2/R3/R4 取证，不修改任何文件）

核查四件事：
  A. 正文实际引用的文献编号 vs 参考文献列表（找"挂名未引"条目）
  B. 关键软件/估计器在 Methods 被描述但是否给了引用编号
  C. 10 个 BMD 位点中进入 MR 的 5 个 —— 选择依据在文件中是否有记录
  D. OA_selfrep 是否进入过任何 MR 分析
"""
import io, os, re, glob, csv

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MS = os.path.join(PROJ, "Manuscript_Ready", "01_Manuscript.md")
md = io.open(MS, encoding="utf-8").read()

# 正文 = 参考文献列表之前
ref_head = md.index("## References")
body = md[:ref_head]

print("=" * 78)
print("A. 引用覆盖：[n] 形式在正文中的出现")
print("=" * 78)
cited = set()
for m in re.finditer(r"\[(\d+)\]", body):
    cited.add(int(m.group(1)))
# 表格/脚注里的 "Refs 1–3,10" 之类
for m in re.finditer(r"Refs?\s+([\d,\u2013\-and ]+)", body):
    for n in re.findall(r"\d+", m.group(1)):
        cited.add(int(n))
listed = set(int(m.group(1)) for m in re.finditer(r"^(\d+)\.\s", md[ref_head:], re.M))
uncited = sorted(listed - cited)
ghost = sorted(cited - listed)
print(f"参考文献列表条目: {len(listed)}   正文中出现过编号的: {len(cited & listed)}")
print(f">> 列表有、正文从未以 [n] 引用的: {uncited}")
print(f">> 正文引用了但列表中没有的(幽灵引用): {ghost}")
print()
for n in uncited:
    m = re.search(rf"^{n}\.\s+(.+?)(?=^\d+\.\s|\Z)", md[ref_head:], re.M | re.S)
    if m:
        title = " ".join(m.group(1).split())[:150]
        print(f"  [{n}] {title}")

print()
print("=" * 78)
print("B. Methods 中被描述但未给引用编号的关键方法")
print("=" * 78)
methods = md[md.index("## 2. Methods"):md.index("## 3. Results")]
for kw, expect in [("weighted median", 7), ("MR-PRESSO", 8), ("coloc", None),
                   ("Cochran Q", None), ("Burgess approximation", 14),
                   ("MRlap", 15), ("Sanderson", 20)]:
    hit = kw in methods
    print(f"  Methods 提到 {kw:<24}: {hit}")

print()
print("=" * 78)
print("C. 进入 MR 的 BMD 位点：文件层面的选择依据")
print("=" * 78)
GEN = os.path.join(PROJ, "code", "output", "genetics")
mr = os.path.join(GEN, "AUDIT_mr_recomputed.csv")
if os.path.exists(mr):
    rows = list(csv.DictReader(io.open(mr, encoding="utf-8")))
    exps = []
    for r in rows:
        e = r.get("exposure") or r.get("exposure_id")
        if e and e not in exps:
            exps.append(e)
    print(f"  AUDIT_mr_recomputed.csv 中出现的暴露: {exps}")
    print(f"  共 {len(rows)} 行")
else:
    print("  [缺失] AUDIT_mr_recomputed.csv")

HARM = os.path.join(GEN, "harmonised")
if os.path.isdir(HARM):
    files = sorted(os.listdir(HARM))
    bmd_as_exp = sorted(set(f.split("__")[0] for f in files if "__" in f))
    print(f"\n  harmonised/ 目录中作为暴露出现过的性状({len(bmd_as_exp)}): {bmd_as_exp}")
    allf = [f for f in files]
    print(f"  harmonised 文件总数: {len(allf)}")

# 十个 BMD 位点的 tophits 是否都存在（即是否有工具）
print("\n  各位点 tophits 文件行数（能否选出工具）:")
for f in sorted(glob.glob(os.path.join(GEN, "tophits_*.csv"))):
    n = sum(1 for _ in io.open(f, encoding="utf-8")) - 1
    print(f"    {os.path.basename(f):<32} {n}")

print()
print("=" * 78)
print("D. OA_selfrep 是否进入过 MR")
print("=" * 78)
hits = [f for f in os.listdir(HARM) if "selfrep" in f or "OA_self" in f] if os.path.isdir(HARM) else []
print(f"  harmonised/ 中含 selfrep 的文件: {hits if hits else '无（OA_selfrep 从未作为 MR 暴露或结局被 harmonise）'}")
if os.path.exists(mr):
    se = [r for r in rows if "selfrep" in (r.get("exposure", "") + r.get("outcome", ""))]
    print(f"  AUDIT_mr_recomputed.csv 中涉及 selfrep 的行: {len(se)}")
