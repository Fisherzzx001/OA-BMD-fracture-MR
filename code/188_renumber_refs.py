# -*- coding: utf-8 -*-
"""
188 · 参考文献重编号（OAC 顺序编码制）
Cause: 参考文献 [44] (Cui 2022) 与 [45] (Hardcastle 2015) 首次出现在讨论，
       早于原 [40]–[43]，违反 OAC "Guide for Authors" 明文要求的
       "按正文首次出现顺序编号"。
Fix:   把两条新文献插到 40/41，原 40–43 顺延为 42–45（循环置换 -> 两阶段替换）。
Verify: code/182 §C "OAC 顺序编码制（按首次出现编号）"。
"""
import re

F = "docs/08_投稿稿_OA_BMD_fracture.md"
t = open(F, encoding="utf-8").read()

# old -> new
MAP = {40: 42, 41: 43, 42: 44, 43: 45, 44: 40, 45: 41}
SENT = {n: f"\x00OLD{n}\x00" for n in MAP}

# ---- phase 1: in-text citations -> sentinel
def cite_sub(m):
    inner = m.group(1)
    out = []
    for part in inner.split(","):
        part_s = part.strip()
        mm = re.match(r"^(\d+)$", part_s)
        if mm:
            n = int(mm.group(1))
            if n in SENT:
                part_s = SENT[n]
        out.append(part_s)
    return "[" + ", ".join(out) + "]"


t = re.sub(r"\[([0-9,\s]+)\]", cite_sub, t)

# ---- phase 2: sentinels -> new numbers
for n, placeholder in SENT.items():
    t = t.replace(placeholder, str(MAP[n]))

# ---- reference list: renumber list entries then restore ascending order
lines = t.split("\n")
ri = next(i for i, l in enumerate(lines) if l.startswith("## References"))
re_entries = {}
for i in range(ri + 1, len(lines)):
    m = re.match(r"^(\d+)\.\s+(.*)$", lines[i])
    if m:
        n = int(m.group(1))
        body_txt = m.group(2)
        if n in MAP:
            n = MAP[n]
        re_entries[n] = body_txt
        lines[i] = None  # mark for removal and rebuild
    elif re_entries and lines[i].strip() and not lines[i].startswith("#"):
        break

# rebuild block
tail_at = None
for i in range(ri + 1, len(lines)):
    if lines[i] is None:
        continue
    if lines[i].strip() == "":
        continue
    tail_at = i
    break

new_block = [f"{n}. {re_entries[n]}" for n in sorted(re_entries)]
lines = [l for l in lines if l is not None]
insert_at = next(i for i, l in enumerate(lines)
                 if l.startswith("## References")) + 1
lines[insert_at:insert_at] = new_block + [""]

out = "\n".join(lines)
open(F, "w", encoding="utf-8").write(out)
print("renumbered OK")
print("new order tail:")
for n in sorted(re_entries)[-8:]:
    print(f"  {n}. {re_entries[n][:80]}...")
