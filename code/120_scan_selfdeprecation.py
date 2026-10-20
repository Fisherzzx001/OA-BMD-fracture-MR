# -*- coding: utf-8 -*-
"""
120_scan_selfdeprecation.py
==========================
从零扫描稿件中的「局限性 / hedging / 自我否定 / 自我贬损」表述。

设计原则（与 100/101/104 系列一致）：
  * 只依赖 docs/08 原文文本，不依赖任何既往审计结论；
  * 输出 = 每一行号 + 原文片段 + 命中类别，可直接对照稿件核；
  * 额外做「重复度」统计：同一语义在全文被复述几次 —— 这是判断
    "过多 limitation / 过度自我否定" 的定量依据。

类别：
  LIM   局限性声明（必要的，但可能被重复）
  HEDGE 模糊限定（may / might / appears / seems / likely / possibly）
  NEG   自我否定（cannot / could not / not established / unclear / unknown / fail）
  DIS   自我贬损（provisional / cautiously / absence of evidence / weakens / underpowered）
  ABS   摘要/亮点区的负面陈述（卖点区放负面 = 结构性问题）
"""

import io
import os
import re
import json
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOC = os.path.join(ROOT, "docs", "08_投稿稿_OA_BMD_fracture.md")

# ---------------------------------------------------------------- 词表
PATTERNS = {
    # --- 自我否定：能力/证据层面的 "做不到" ---
    "NEG": [
        r"\bcannot\b", r"\bcan not\b", r"\bcould not\b", r"\bcouldn't\b",
        r"\bnot able to\b", r"\bunable to\b",
        r"\bdoes not establish\b", r"\bis not established\b", r"\bnot established\b",
        r"\bnot directly comparable\b", r"\bnot a directly measured\b",
        r"\bdo not support\b", r"\bdoes not support\b", r"\bdid not support\b",
        r"\bunsettled\b", r"\bunknown\b", r"\bunclear\b", r"\bunresolved\b",
        r"\bremains? open\b", r"\bno .{0,30}evidence\b", r"\bfails?\b", r"\bfailed\b",
        r"\bnot .{0,20}count(ed)?\b", r"\bnon-significant\b", r"\bnon-significantl?y?\b",
        r"\bnull\b", r"\babsent\b", r"\bnot decompose\b", r"\bnot separate\b",
    ],
    # --- 自我贬损：主动降低自己结论的分量 ---
    "DIS": [
        r"\bprovisional\b", r"\bcautious(ly)?\b", r"\bcaution(s|ed)? against\b",
        r"\babsence of evidence\b", r"\bweakly powered\b", r"\bunderpowered\b",
        r"\bunder-powered\b", r"\bweak\b", r"\bweaken(s|ed)?\b",
        r"\bshifted against\b", r"\blimited\b", r"\blimitation(s)?\b",
        r"\bpreliminary\b", r"\bhypothesis-generating\b",
        r"\bshould be interpreted\b", r"\binterpreted with care\b",
        r"\bresidual\b", r"\bcannot be (fully )?excluded\b",
        r"\bleaves? the causal claim\b", r"\brests? on\b",
        r"\brather than on\b", r"\bis an upper bound\b", r"\bnot a point estimate\b",
        r"\bimperfectly\b", r"\bmarginally short\b", r"\bfell .{0,15}short\b",
        r"\bconservative\b", r"\bdepends on\b", r"\bsensitivity to\b",
        r"\battrition\b", r"\bviolat(e|ing|es|ion)\b", r"\bviolated\b",
    ],
    # --- 模糊限定 ---
    "HEDGE": [
        r"\bmay\b", r"\bmight\b", r"\bcould\b", r"\bpossibly\b", r"\bpossible\b",
        r"\bappears?\b", r"\bseems?\b", r"\blikely\b", r"\bplausibl[ey]\b",
        r"\bpotentially\b", r"\bpresumably\b", r"\bsuggest(s|ing|ed)?\b",
        r"\bconsistent with\b", r"\bcompatible with\b", r"\btentative(ly)?\b",
        r"\bspeculative\b", r"\bwe therefore base\b", r"\bbetter described as\b",
    ],
}

# 语义主题：判断"同一条 limitation 被复述了几次"
THEMES = {
    "mass_vs_size": [
        r"mass[- ]versus[- ]size", r"mass versus size", r"size with mass",
        r"mass or by skeletal size", r"cannot separate", r"not separate",
        r"could not separate", r"not decompose", r"cannot decompose",
        r"size confounding", r"skeletal size", r"confounded by skeletal size",
    ],
    "heterogeneity": [
        r"heterogeneity", r"I\u00b2", r"pleiotropy-robust",
        r"random[- ]effects", r"random effects",
    ],
    "underpowered_reverse": [
        r"underpowered", r"weakly powered", r"power 0\.14", r"achieved power",
        r"0\.14\u20130\.26", r"MDE",
    ],
    "sample_overlap": [
        r"sample overlap", r"overlap", r"MRlap", r"Burgess",
    ],
    "scale_areal_bmd": [
        r"areal[- ]BMD", r"areal BMD", r"per-SD scale", r"effect-scale",
        r"DXA", r"placeholder N",
    ],
    "coloc_null": [
        r"PP\.H4", r"colocalis", r"colocaliz",
    ],
}


def split_sentences(text):
    """按句子切分，保留行号。"""
    out = []
    for ln, line in enumerate(text.split("\n"), start=1):
        if not line.strip():
            continue
        # 不在表格行内做句子切分过细；用简单句号切分
        parts = re.split(r"(?<=[.!?])\s+(?=[A-Z(\u201c])", line)
        for p in parts:
            p = p.strip()
            if p:
                out.append((ln, p))
    return out


def main():
    md = io.open(DOC, encoding="utf-8").read()
    lines = md.split("\n")

    # ---------- 1) 逐句扫描 ----------
    hits = []
    for ln, sent in split_sentences(md):
        cats = defaultdict(list)
        for cat, pats in PATTERNS.items():
            for p in pats:
                m = re.search(p, sent, flags=re.I)
                if m:
                    cats[cat].append(m.group(0))
        if cats:
            hits.append({
                "line": ln,
                "text": sent,
                "cats": sorted(cats.keys()),
                "matches": {k: sorted(set(v)) for k, v in cats.items()},
            })

    # ---------- 2) 主题重复度 ----------
    theme_count = {}
    theme_lines = {}
    for theme, pats in THEMES.items():
        c = 0
        lns = []
        for ln, line in enumerate(lines, start=1):
            for p in pats:
                if re.search(p, line, flags=re.I):
                    c += len(re.findall(p, line, flags=re.I))
                    if ln not in lns:
                        lns.append(ln)
                    break
        theme_count[theme] = c
        theme_lines[theme] = lns

    # ---------- 3) 分区负面密度 ----------
    def section_of(ln):
        if ln <= 17:
            return "Highlights"
        if ln <= 24:
            return "Abstract"
        if ln <= 35:
            return "Introduction"
        if ln <= 53:
            return "Methods"
        if ln <= 71:
            return "Results"
        if ln <= 87:
            return "Discussion"
        if ln <= 113:
            return "Figure legends"
        if ln <= 200:
            return "Tables"
        return "Back matter"

    sec_stat = defaultdict(lambda: {"sent": 0, "hit": 0, "neg": 0, "dis": 0})
    for h in hits:
        s = section_of(h["line"])
        sec_stat[s]["hit"] += 1
        if "NEG" in h["cats"]:
            sec_stat[s]["neg"] += 1
        if "DIS" in h["cats"]:
            sec_stat[s]["dis"] += 1
    for ln, line in enumerate(lines, start=1):
        s = section_of(ln)
        sec_stat[s]["sent"] += len([x for x in re.split(r"(?<=[.!?])\s+", line) if x.strip()])

    # ---------- 4) 输出 ----------
    print("=" * 78)
    print("稿件自我否定 / hedging 扫描报告")
    print("=" * 78)
    print(f"\n总句数(估): {sum(v['sent'] for v in sec_stat.values())}")
    print(f"命中句数:   {len(hits)}")

    print("\n" + "-" * 78)
    print("【A】分区负面密度（sent=句数, hit=命中, NEG=自我否定, DIS=自我贬损）")
    print("-" * 78)
    print(f"{'Section':<16}{'sent':>6}{'hit':>6}{'NEG':>6}{'DIS':>6}{'hit%':>8}")
    for s in ["Highlights", "Abstract", "Introduction", "Methods", "Results",
              "Discussion", "Figure legends", "Tables", "Back matter"]:
        if s in sec_stat:
            v = sec_stat[s]
            pct = 100.0 * v["hit"] / v["sent"] if v["sent"] else 0.0
            print(f"{s:<16}{v['sent']:>6}{v['hit']:>6}{v['neg']:>6}{v['dis']:>6}{pct:>7.1f}%")

    print("\n" + "-" * 78)
    print("【B】主题重复度（同一条 limitation 在全文被复述的频次 / 出现行号）")
    print("-" * 78)
    for t, c in sorted(theme_count.items(), key=lambda x: -x[1]):
        print(f"  {t:<22} mentions={c:<4} lines={theme_lines[t]}")

    print("\n" + "-" * 78)
    print("【C】逐条命中清单（行号 · 类别 · 片段）")
    print("-" * 78)
    for h in hits:
        tag = "+".join(h["cats"])
        snippet = h["text"]
        if len(snippet) > 150:
            snippet = snippet[:147] + "..."
        print(f"L{h['line']:>3} [{tag}] {snippet}")

    out = {
        "n_hits": len(hits),
        "section_density": {k: dict(v) for k, v in sec_stat.items()},
        "theme_counts": theme_count,
        "theme_lines": theme_lines,
        "hits": hits,
    }
    op = os.path.join(ROOT, "code", "output", "AUDIT_120_selfdeprecation.json")
    io.open(op, "w", encoding="utf-8").write(
        json.dumps(out, ensure_ascii=False, indent=1))
    print(f"\n-> {op}")


if __name__ == "__main__":
    main()
