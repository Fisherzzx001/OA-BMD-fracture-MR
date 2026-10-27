# -*- coding: utf-8 -*-
"""
预审委员会第二轮修复：OSF/体例缺陷 + 参考文献重编号（OAC 顺序编码制）

依据：OAC (Elsevier) Guide for Authors 明确规定
  "Number references in the order they appear in your article."
  "References should be listed in numerical sequence as they are cited in the text."
原稿件编号按字母/主题分组，非首次出现顺序，属会被技术编辑退回的格式缺陷。

修复项：
 A. Highlights 第 3 条缺空格 ":endpoint-specific" -> ": endpoint-specific"
 B. 利益冲突声明重复（两条措辞不同）-> 保留 Elsevier 标准措辞，删除冗余条
 C. 补充材料清单 "F1–F7" -> "F1–F10"（与正文 F1–F10 一致）
 D. STROBE-MR 报告规范仅有作者-年形式且不在文献表 -> 新增 [41] Skrivankova JAMA 2021，
    并替换两处作者-年写法（原写 BMJ 2021 属 elaboration，statement 应在 JAMA）
 E. 参考文献按首次出现顺序重编号（含 "[27–30]" 区间展开与 "(Refs 1–3,10,19,38)" 特例）

安全：每处替换要求出现次数恰好为 1；硬约束不改任何 β/SE/p/OR/CI（本脚本不触碰数值）。
"""
import io
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DOC = ROOT / "docs" / "08_投稿稿_OA_BMD_fracture.md"

NEW_REF = ("41. Skrivankova VW, Richmond RC, Woolf BAR, Yarmolinsky J, Davies NM, Swanson SA, "
           "VanderWeele TJ, Higgins JPT, Timpson NJ, Dimou N, Langenberg C, Golub RM, Loder EW, "
           "Gallo V, Tybjaerg-Hansen A, Davey Smith G, Egger M, Richards JB. Strengthening the "
           "Reporting of Observational Studies in Epidemiology Using Mendelian Randomization: "
           "the STROBE-MR statement. *JAMA*. 2021;326(16):1614–1621. "
           "doi:10.1001/jama.2021.18236. PMID 34698778.")

CITE = re.compile(r"\[(\d+(?:[–\-]\d+)?(?:,\s*\d+)*)\]")


def die(msg):
    print(f"  ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def sub_once(text, old, new, label):
    n = text.count(old)
    if n != 1:
        die(f"{label}: 目标串出现 {n} 次（要求 1 次）-> {old[:60]!r}")
    print(f"  [OK] {label}")
    return text.replace(old, new)


def parse_nums(spec):
    """'27–30' / '16,17' / '1–3,10,19,38' -> [1,2,...]"""
    out = []
    for part in re.split(r",\s*", spec.strip()):
        if "–" in part or "-" in part:
            a, b = re.split(r"[–-]", part)
            out += list(range(int(a), int(b) + 1))
        else:
            out.append(int(part))
    return out


def compress(nums):
    """[5,6,7,10] -> [(5,7),(10,10)]"""
    nums = sorted(set(nums))
    runs, s, p = [], nums[0], nums[0]
    for n in nums[1:]:
        if n == p + 1:
            p = n
        else:
            runs.append((s, p))
            s = p = n
    runs.append((s, p))
    return runs


def fmt(runs):
    return ", ".join(str(a) if a == b else f"{a}–{b}" for a, b in runs)


def main():
    text = io.open(DOC, encoding="utf-8").read()
    print("=== A–D 内容修复 ===")

    # A. Highlights 标点
    text = sub_once(text, "(0/31):endpoint-specific", "(0/31): endpoint-specific",
                    "A. Highlights 补全空格")

    # B. 删除重复的利益冲突声明（保留 Elsevier 标准措辞的那条）
    dup_line = "- **Competing interests.** The authors declare no competing interests.\n"
    text = sub_once(text, dup_line, "", "B. 删除重复的 Competing interests 声明")

    # C. 补充材料 F1–F7 -> F1–F10
    text = sub_once(text, "predictions (P1–P6) and falsifying tests (F1–F7)",
                    "predictions (P1–P6) and falsifying tests (F1–F10)",
                    "C. 补充材料清单 F1–F7 -> F1–F10")

    # D. STROBE-MR 引用
    text = sub_once(text,
                    "STROBE-MR reporting guideline [Skrivankova et al., BMJ 2021]",
                    "STROBE-MR reporting guideline [41]",
                    "D1. 方法段 STROBE-MR 改为编号引用")
    text = sub_once(text,
                    "22-item reporting-compliance table (Skrivankova et al., BMJ 2021)",
                    "22-item reporting-compliance table (Skrivankova et al. 2021 [41])",
                    "D2. 补充材料清单 STROBE-MR 改为编号引用")

    print("\n=== E. 参考文献按首次出现顺序重编号 ===")

    # 切出参考文献区
    m = re.search(r"\n## References\n(.*?)\n## Acknowledgements", text, re.S)
    if not m:
        die("未能定位 References 区块")
    refblock = m.group(1)
    head, tail = text[:m.start(1)], text[m.end(1):]

    # 追加新文献条目
    refblock = refblock.rstrip("\n") + "\n\n" + NEW_REF + "\n"
    print(f"  [OK] 新增 [41] Skrivankova 2021 (JAMA STROBE-MR statement)")

    # 解析条目
    entries = {}
    for line in refblock.splitlines():
        mm = re.match(r"^(\d+)\.\s+(.*)$", line.strip())
        if mm:
            entries[int(mm.group(1))] = mm.group(2).strip()
    if sorted(entries) != list(range(1, len(entries) + 1)):
        die(f"参考文献编号不连续: {sorted(entries)}")
    print(f"  [OK] 解析 {len(entries)} 条参考文献")

    body = head + tail  # 除文献区外的全文

    # 计算首次出现顺序（跳过文献区）
    first = {}
    for cm in CITE.finditer(body):
        for n in parse_nums(cm.group(1)):
            if n not in first:
                first[n] = cm.start()
    missing = sorted(set(entries) - set(first))
    if missing:
        die(f"以下文献从未被引用: {missing}")
    print(f"  [OK] 全部 {len(entries)} 条均被引用")

    old_order = sorted(first, key=lambda n: first[n])
    mapping = {old: new for new, old in enumerate(old_order, 1)}
    print("\n  旧 -> 新 编号映射：")
    for old in old_order:
        print(f"    [{old}] -> [{mapping[old]}]  {entries[old][:58]}")

    if old_order != sorted(first):
        print("  → 需要重排（当前非首次出现顺序）")
    else:
        print("  → 已符合顺序编码制")
        return

    # 替换正文引用
    def repl(m):
        nums = [mapping[n] for n in parse_nums(m.group(1))]
        return "[" + fmt(compress(nums)) + "]"

    new_body, n_cit = CITE.subn(repl, body)
    print(f"\n  [OK] 正文替换 {n_cit} 处引用")

    # "(Refs 1–3,10,19,38)" 特例
    def repl_refs(m):
        nums = [mapping[n] for n in parse_nums(m.group(1))]
        return "Refs " + fmt(compress(nums))

    new_body, n_refs = re.subn(r"Refs (\d+(?:[–\-]\d+)?(?:,\s*\d+)*)", repl_refs, new_body)
    print(f"  [OK] 'Refs ' 特例替换 {n_refs} 处")

    # 重排文献表
    lines = [f"{new}. {entries[old]}" for new, old in enumerate(old_order, 1)]
    new_refblock = "\n".join(lines) + "\n"

    # 重组文件：拼接点必须在【替换之后】重新定位——替换会改变正文长度，
    # 用替换前的偏移量切断会把 "## References" 标题劈坏（第一轮出错即此因）。
    marker = "\n## References\n"
    cut2 = new_body.index(marker) + len(marker)
    assert new_body[cut2:].lstrip().startswith("## Acknowledgements"), "拼接点定位失败"
    final = new_body[:cut2] + "\n" + new_refblock + new_body[cut2:]
    io.open(DOC, "w", encoding="utf-8", newline="\n").write(final)

    # ---- 复核
    print("\n=== 复核 ===")
    t2 = io.open(DOC, encoding="utf-8").read()
    b2 = t2[:t2.index("\n## References")] + t2[t2.index("\n## Acknowledgements"):]
    f2 = {}
    for cm in CITE.finditer(b2):
        for n in parse_nums(cm.group(1)):
            if n not in f2:
                f2[n] = cm.start()
    order2 = sorted(f2, key=lambda n: f2[n])
    ok = order2 == sorted(f2) == list(range(1, len(f2) + 1))
    print(f"  首次出现顺序 = 自然升序且连续 1..{len(f2)}: {ok}")
    print(f"  顺序: {order2}")
    if not ok:
        die("重编号未通过复核")


if __name__ == "__main__":
    main()
