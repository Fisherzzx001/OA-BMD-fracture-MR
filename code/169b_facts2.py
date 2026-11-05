# -*- coding: utf-8 -*-
"""取证2: (1) 仅作者年引用而未进编号表的方法/软件文献; (2) 重复表述占词; (3) 既可 housed 的冗余."""
import re

P = "docs/08_投稿稿_OA_BMD_fracture.md"
t = open(P, encoding="utf-8").read()
main, rest = t.split("## References", 1)
refs_block = rest.split("\n## ", 1)[0]

print("=" * 78)
print("E. 正文中的作者年式引用 (Vancouver 编号体例下不应出现)")
print("=" * 78)
pat = re.compile(r"([A-Z][A-Za-z\u00C0-\u017F\-]+(?:\s+et\s+al\.?)?)\s*,?\s*((?:19|20)\d{2})[a-z]?")
found = []
for m in pat.finditer(main):
    s = max(0, m.start() - 55)
    found.append((m.group(1).strip(), m.group(2), main[s:m.end() + 30].replace("\n", " ")))
for a, y, ctx in found:
    print("  %-26s %s   | ...%s" % (a, y, ctx[-95:]))

print()
print("=" * 78)
print("F. 参考文献表中是否登记了上述来源")
print("=" * 78)
for probe in ["Skrivankova", "Giambartolomei", "Hemani", "TwoSampleMR", "PLINK",
              "purge", "1000 Genomes", "1000G Phase 3"]:
    print("  %-16s in reference list: %s" % (probe, probe in refs_block))

print()
print("=" * 78)
print("G. 重复表述占词统计 (用于购买改动所需的词额度)")
print("=" * 78)
probes = {
    "ancestry/population-stratification 警示": r"(?:ancestry-mismatch|population-stratification|ancestry mismatch|European-ancestry)",
    "height does not capture width/geometry": r"(?:bone width|geometry)",
    "not resolved into mass versus geometry": r"(?:not (?:been )?resolved into|not resolved into|not established as bone mass)",
    "multi-estimator concordance": r"(?:estimator concordance|concordance of pleiotropy-robust)",
    "estimates do not support changing prescribing": r"(?:do not support changing|not support changing prescribing)",
}
for label, pat in probes.items():
    hits = list(re.finditer(pat, main, flags=re.I))
    words = sum(len(main[max(0, m.start() - 180):m.end() + 180].split()) for m in hits)
    print("  %-46s 命中 %d 处" % (label, len(hits)))
    for m in hits:
        print("        ...%s..." % main[max(0, m.start()-75):m.end()+75].replace("\n", " ")[-135:])

print()
print("=" * 78)
print("H. 结语/主张强度与证据匹配")
print("=" * 78)
for kw in ["do not support changing prescribing practice", "blood is neither bone nor cartilage",
           "reflects compute availability", "a hypothesis consistent with these data",
           "expectation rather than established positioning"]:
    print("  [%s] -> %s" % (kw, "found" if kw in main else "ABSENT"))
