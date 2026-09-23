# -*- coding: utf-8 -*-
import re, os
PROJ = os.path.dirname(os.path.abspath(__file__))
md = open(os.path.join(PROJ, "docs", "08_投稿稿_OA_BMD_fracture.md"), encoding="utf-8").read()

print("=" * 90)
print("稿件引用的数据文件是否存在")
print("=" * 90)
files = set(re.findall(r"`([A-Za-z0-9_./\-]+\.(?:csv|json|py|R|md))`", md))
for f in sorted(files):
    cands = [f, os.path.join(PROJ, f),
             os.path.join(PROJ, "code", "output", "genetics", f),
             os.path.join(PROJ, "code", "output", f),
             os.path.join(PROJ, "code", f),
             os.path.join(PROJ, "docs", f),
             os.path.join(PROJ, "submission_package_v16", f)]
    ok = any(os.path.exists(c) for c in cands)
    print(("OK   " if ok else "MISS "), f)

print()
print("=" * 90)
print("图/表 交叉引用一致性")
print("=" * 90)
fig_refs = sorted(set(re.findall(r"Figure (\d+|S\d+)", md)), key=lambda s: (len(s), s))
fig_defs = sorted(set(re.findall(r"\*\*(?:Supplementary )?Figure (\w+)\.\*\*", md)))
tbl_defs = sorted(set(re.findall(r"## Table (\d+)\.", md)))
tbl_refs = sorted(set(re.findall(r"Table (\d+|S-[A-Za-z]+)", md)))
print("图注定义 :", fig_defs)
print("正文引用 :", fig_refs)
print("表定义   :", tbl_defs)
print("表引用   :", tbl_refs)

missing = [f for f in fig_refs if f not in fig_defs]
print("\n引用了但无图注的图:", missing if missing else "无")
unused = [f for f in fig_defs if f not in fig_refs]
print("有图注但正文未引用:", unused if unused else "无")
mt = [t for t in tbl_refs if t not in tbl_defs and not t.startswith("S-")]
print("引用了但无定义的表:", mt if mt else "无")

print()
print("=" * 90)
print("残留的确定性负向断言 / over-claim 扫描")
print("=" * 90)
PAT = [
    (r"does not raise", "确定性: does not raise"),
    (r"do not raise", "确定性: do not raise"),
    (r"were null", "非显著当零效应: were null"),
    (r"are null", "非显著当零效应: are null"),
    (r"is null\b", "非显著当零效应: is null"),
    (r"no causal", "确定性: no causal"),
    (r"rules? out", "确定性: rules out"),
    (r"does not causally", "确定性: does not causally"),
    (r"not causally", "确定性: not causally"),
    (r"excludes? (?:a|any|the) ", "排除性断言"),
    (r"is not mediated", "确定性: not mediated"),
    (r"no evidence of a causal", "确定性: no evidence of a causal"),
]
for pat, label in PAT:
    for m in re.finditer(pat, md):
        line_no = md[:m.start()].count("\n") + 1
        s = max(0, m.start() - 110)
        e = min(len(md), m.end() + 110)
        ctx = md[s:e].replace("\n", " ")
        print(f"L{line_no:4d} [{label}] ...{ctx}...")

print()
print("=" * 90)
print("关键数值在稿件中出现次数（用于检测跨章节不一致）")
print("=" * 90)
for key in ["0.337", "0.351", "0.366", "0.106", "0.142", "0.092", "0.0354", "0.0258",
            "1.40", "1.42", "1.55", "31,873", "264,973", "1,208,397", "0.303", "0.84"]:
    n = md.count(key)
    print(f"  {key:>12s}: {n} 次")
