# -*- coding: utf-8 -*-
"""
130 · AIGC 痕迹定量诊断（改前基线 / 改后复核通用）

输出六类指标，全部按「行号 + 片段」可溯源：
  M1 加粗密度及其分类（必要 vs 可去）
  M2 破折号 em dash / en dash
  M3 分号
  M4 句长节奏：均值、标准差、变异系数、最长/最短
  M5 连接词集中度 + AI 高发词
  M6 被动语态比例 + 句首小品词
"""

import io
import os
import re
import sys
import statistics as st
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOC = os.path.join(ROOT, "docs", "08_投稿稿_OA_BMD_fracture.md")

EM = "\u2014"
EN = "\u2013"

# 必要加粗：摘要结构化标签、图注编号、表标题
KEEP_BOLD = re.compile(
    r"^\*\*(Objective|Methods|Results|Conclusion|Abstract|Keywords|Highlights|"
    r"Figure\s|S1|S2|Supplementary Figure|Table\s|Introduction|Discussion)\b")

AI_WORDS = ["robust", "robustness", "coherent", "canonical", "uniformly", "crucial",
            "critically", "exquisite", "underscore", "highlight", "highlighting",
            "landscape", "delve", "Moreover", "Furthermore", "Additionally", "Notably",
            "Importantly", "Taken together", "serves as", "systematically", "pivotal",
            "showcase", "foster", "leverage", "utilize", "groundbreaking", "novel",
            "comprehensive", "multifaceted", "intricate", "interplay", "tapestry",
            "realm", "testament", "enduring", "align with", "key", "vital", "valuable"]

CONNECTORS = ["however", "therefore", "moreover", "furthermore", "additionally",
              "consequently", "indeed", "thus", "nevertheless", "whereas", "although",
              "while", "because", "since", "despite", "accordingly", "hence"]

SENTENCE_START_PARTICLES = [
    "Additionally", "Moreover", "Furthermore", "Consequently", "Notably",
    "Importantly", "Interestingly", "Significantly", "Crucially", "Overall",
    "Taken together", "Collectively", "Together", "Indeed", "Thus", "Hence",
    "Accordingly", "Nevertheless", "Nonetheless", "Specifically", "Particularly",
]

PASSIVE = re.compile(r"\b(?:is|are|was|were|be|been|being)\s+(?:\w+ly\s+)?\w+(?:ed|en)\b",
                     re.I)


def sentences(text):
    out = []
    for ln, line in enumerate(text.split("\n"), 1):
        if line.strip().startswith("|") or line.strip().startswith("#"):
            continue
        if not line.strip():
            continue
        for s in re.split(r"(?<=[.!?])\s+(?=[A-Z(\u201c])", line):
            s = s.strip()
            if s and re.search(r"[A-Za-z]{3}", s):
                out.append((ln, s))
    return out


def main():
    md = io.open(DOC, encoding="utf-8").read()
    lines = md.split("\n")
    body_lines = [l for l in lines if not l.strip().startswith("|")]
    body = "\n".join(body_lines)
    words = re.findall(r"[A-Za-z][A-Za-z'-]*", body)
    nw = len(words)
    sents = sentences(body)
    ns = len(sents)

    print("=" * 78)
    print("AIGC 痕迹定量诊断")
    print("=" * 78)
    print(f"正文词数 {nw}　句数 {ns}")

    # ---------- M1 加粗 ----------
    print("\n" + "-" * 78)
    print("【M1】加粗密度与分类")
    print("-" * 78)
    bold_keep, bold_drop = [], []
    for ln, l in enumerate(lines, 1):
        if l.strip().startswith("|"):
            continue
        for m in re.finditer(r"\*\*([^*]+)\*\*", l):
            txt = m.group(1)
            if KEEP_BOLD.match("**" + txt + "**"):
                bold_keep.append((ln, txt))
            else:
                bold_drop.append((ln, txt))
    print(f"  必要加粗（摘要标签/图注/表题）: {len(bold_keep)}")
    print(f"  可去加粗（正文强调）          : {len(bold_drop)}")
    print(f"  加粗总密度: 1/{nw/max(1,len(bold_keep)+len(bold_drop)):.0f} 词")
    for ln, txt in bold_drop:
        print(f"    L{ln:>3} {txt[:76]}")

    # ---------- M2 破折号 ----------
    print("\n" + "-" * 78)
    print("【M2】破折号 (目标: em dash ≤2 处 / 5000 词)")
    print("-" * 78)
    for ln, l in enumerate(lines, 1):
        if l.strip().startswith("|"):
            continue
        if EM in l:
            segs = [s for s in re.split(EM, l) if s.strip()]
            ctx = (" ... ".join(x.strip()[-38:] + " \u2014 " + y.strip()[:38]
                                for x, y in zip(segs, segs[1:])))
            print(f"    L{ln:>3} x{l.count(EM)}  {ctx[:110]}")
    print(f"  em dash 合计 {body.count(EM)}  -> 1/{nw/max(1,body.count(EM)):.0f} 词")

    # ---------- M3 分号 ----------
    print("\n" + "-" * 78)
    print("【M3】分号（非表格区）")
    print("-" * 78)
    semi = [(ln, l) for ln, l in enumerate(lines, 1)
            if ";" in l and not l.strip().startswith("|")]
    for ln, l in semi[:25]:
        print(f"    L{ln:>3} x{l.count(';')}  {l.strip()[:100]}")
    print(f"  分号合计 {sum(l.count(';') for _, l in semi)} 行内 {len(semi)} 行")
    print(f"  （表格区分号 {sum(l.count(';') for l in lines if l.strip().startswith('|'))}，属数据格式，不计）")

    # ---------- M4 句长节奏 ----------
    print("\n" + "-" * 78)
    print("【M4】句长节奏（burstiness）")
    print("-" * 78)
    lens = [len(re.findall(r"[A-Za-z][A-Za-z'-]*", s)) for _, s in sents]
    if lens:
        mu = st.mean(lens)
        sd = st.pstdev(lens)
        print(f"  均值 {mu:.1f} 词　标准差 {sd:.1f}　变异系数 CV={sd/mu:.3f}")
        print(f"  最短 {min(lens)}　最长 {max(lens)}")
        q = sorted(lens)
        print(f"  四分位 P25={q[len(q)//4]}  中位={q[len(q)//2]}  P75={q[3*len(q)//4]}")
        print(f"  >40 词的长句: {sum(1 for x in lens if x>40)} 句")
        print(f"  <12 词的短句: {sum(1 for x in lens if x<12)} 句")
        print("  参考：人类学术写作 CV 通常 0.45-0.70；AI 文本常 <0.40")

    # ---------- M5 连接词与 AI 词 ----------
    print("\n" + "-" * 78)
    print("【M5】句首小品词 / 连接词 / AI 高发词")
    print("-" * 78)
    ssp = Counter()
    for _, s in sents:
        first = s.split()[0].strip(",;:")
        if first in SENTENCE_START_PARTICLES:
            ssp[first] += 1
    print(f"  句首小品词合计 {sum(ssp.values())} / {ns} 句 "
          f"({100*sum(ssp.values())/max(1,ns):.1f}%)")
    for w, c in ssp.most_common():
        print(f"      {w:<16} {c}")
    low = " ".join(w.lower() for w in words)
    print("  AI 高发词：")
    for w in AI_WORDS:
        c = low.count(" " + w.lower() + " ") + low.startswith(w.lower() + " ")
        if c:
            print(f"      {w:<16} {c}")

    # ---------- M6 被动语态 ----------
    print("\n" + "-" * 78)
    print("【M6】被动语态")
    print("-" * 78)
    npass = sum(1 for _, s in sents if PASSIVE.search(s))
    print(f"  含被动结构句 {npass}/{ns} = {100*npass/max(1,ns):.1f}%")
    print("  参考：学术 Methods 段被动 40-60% 属正常；Results/Discussion >20% 偏高")

    # ---------- 汇总 ----------
    print("\n" + "=" * 78)
    print("汇总（对照 sci-manuscript-deai 目标）")
    print("=" * 78)
    tgt = [
        ("em dash", body.count(EM), "<=2"),
        ("正文强调加粗", len(bold_drop), "<=0"),
        ("句长 CV", round(sd / mu, 3) if lens else 0, ">=0.45"),
        ("句首小品词占比", f"{100*sum(ssp.values())/max(1,ns):.1f}%", "<=8%"),
    ]
    for k, v, t in tgt:
        print(f"  {k:<16} {str(v):>8}   目标 {t}")


if __name__ == "__main__":
    main()
