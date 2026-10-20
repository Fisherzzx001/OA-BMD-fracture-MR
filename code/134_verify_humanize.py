# -*- coding: utf-8 -*-
"""
134 · 降噪后的完整性校验 + AIGC 形式化评分（改前 / 改后对照）

两部分：

【一】数字完整性（证明语言润色没有动数据）
  对改前（.bak-before-humanize131）与改后做**数字多重集比对**。
  任何被删除/新增的数字都会出现在差集中。允许的差异必须逐条说明。

【二】AIGC 形式化评分
  说明：商业检测器（GPTZero / Originality.ai / Sapling）需付费 API key，
  本项目未持有；Sapling 端点实测可达但返回 401 Missing API key。
  因此这里给出的是**可复现的形式化代理评分**，输入正是检测器所依赖的
  表层信号（em dash、加粗、burstiness、句首小品词、AI 词汇、分号、
  被动语态、词汇多样性）。评分权重透明列出，便于复核或替换成真值。

  评分越高 = 越像 AI。总分 100。
"""

import io
import os
import re
import statistics as st
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOC = os.path.join(ROOT, "docs", "08_投稿稿_OA_BMD_fracture.md")
BAK = DOC + ".bak-before-humanize131"

# 注意：不含 robust / robustness —— 在 Mendelian randomisation 文献里
# "pleiotropy-robust estimators" "robustness overview" 是标准术语，不是 AI 痕迹。
AI_WORDS = ["coherent", "canonical", "uniformly", "crucial", "critically",
            "exquisite", "underscore", "highlight", "landscape", "delve", "Moreover",
            "Furthermore", "Additionally", "Notably", "Importantly", "Taken together",
            "serves as", "systematically", "pivotal", "showcase", "foster", "leverage",
            "utilize", "groundbreaking", "comprehensive", "multifaceted", "intricate",
            "interplay", "tapestry", "realm", "testament", "enduring", "vital", "valuable"]

# 句间分号判定：分号后紧跟这些词，说明后面是一个完整子句
CLAUSAL_AFTER = {"we", "it", "this", "which", "the", "however", "these", "they",
                 "such", "all", "both", "no", "each", "our", "that", "for", "under",
                 "across", "when", "because", "although", "whereas", "and", "but",
                 "six", "both"}

SSP = ["Additionally", "Moreover", "Furthermore", "Consequently", "Notably",
       "Importantly", "Interestingly", "Significantly", "Crucially", "Overall",
       "Taken together", "Collectively", "Together", "Indeed", "Thus", "Hence",
       "Accordingly", "Nevertheless", "Nonetheless", "Specifically", "Particularly"]

PASSIVE = re.compile(r"\b(?:is|are|was|were|be|been|being)\s+(?:\w+ly\s+)?\w+(?:ed|en)\b", re.I)
KEEP_BOLD = re.compile(r"^\*\*(Objective|Methods|Results|Conclusion|Abstract|Keywords|"
                       r"Highlights|Figure\s|S\d|Supplementary Figure|Table\s|"
                       r"Introduction|Discussion)\b")


def metrics(path):
    t = io.open(path, encoding="utf-8").read()
    lines = t.split("\n")
    body = "\n".join(l for l in lines if not l.strip().startswith("|"))
    words = re.findall(r"[A-Za-z][A-Za-z'-]*", body)
    nw = max(1, len(words))

    sents = []
    for l in lines:
        if l.strip().startswith("|") or l.strip().startswith("#"):
            continue
        for s in re.split(r"(?<=[.!?])\s+(?=[A-Z(\u201c])", l):
            s = s.strip()
            if s and re.search(r"[A-Za-z]{3}", s):
                sents.append(s)
    ns = max(1, len(sents))
    lens = [len(re.findall(r"[A-Za-z][A-Za-z'-]*", s)) for s in sents] or [1]

    em = body.count("\u2014")
    bold_drop = 0
    for l in lines:
        if l.strip().startswith("|"):
            continue
        for m in re.finditer(r"\*\*([^*]+)\*\*", l):
            if not KEEP_BOLD.match("**" + m.group(1) + "**"):
                bold_drop += 1

    # 分号分两类：列举性（分隔并列统计项，学术 MR 论文必需）与句间（连接两个子句）
    semi = 0
    semi_clausal = 0
    for l in lines:
        if l.strip().startswith("|"):
            continue
        semi += l.count(";")
        for m in re.finditer(r";\s*([A-Za-z']+)", l):
            if m.group(1).lower() in CLAUSAL_AFTER:
                semi_clausal += 1
    low = " ".join(w.lower() for w in words)
    aic = sum(len(re.findall(r"\b" + re.escape(w.lower()) + r"\b", low)) for w in AI_WORDS)
    ssp = sum(1 for s in sents if s.split()[0].strip(",;:") in SSP)
    pas = sum(1 for s in sents if PASSIVE.search(s))
    # TTR 随文本长度单调下降，对 5000 词稿件必然偏低，不能直接设阈值。
    # 改用 MATTR（moving-average TTR，窗口 100 词），长度无关。
    wl = [w.lower() for w in words]
    win = 100
    if len(wl) >= win:
        mattr = st.mean([len(set(wl[i:i + win])) / win
                         for i in range(0, len(wl) - win + 1, 25)])
    else:
        mattr = len(set(wl)) / max(1, len(wl))
    ttr = len(set(w.lower() for w in words)) / nw   # 仅作参考输出

    return dict(nw=nw, ns=ns, em=em, bold=bold_drop, semi=semi, semi_cl=semi_clausal, aiw=aic,
                ssp=100.0 * ssp / ns, pas=100.0 * pas / ns, ttr=ttr, mattr=mattr,
                cv=st.pstdev(lens) / st.mean(lens), mu=st.mean(lens))
    # 单位：em / bold / semi / aiw 均为绝对计数；评分时按每千词归一


def score(m):
    """0-100，越高越像 AI。权重透明。"""
    s, detail = 0.0, []

    def add(name, pts, cap, why):
        nonlocal s
        v = min(pts, cap)
        s += v
        detail.append((name, round(v, 1), cap, why))

    per_k = lambda k: 1000.0 * k / m["nw"]

    # 1 em dash（GPT-4 家族最强表层标志）
    add("em dash 密度", 15 * min(1.0, per_k(m["em"]) / 1.0), 15,
        f"{m['em']} 处 = {per_k(m['em']):.2f}/千词")
    # 2 正文强调加粗
    add("加粗密度", 15 * min(1.0, per_k(m["bold"]) / 5.0), 15,
        f"{m['bold']} 处 = {per_k(m['bold']):.2f}/千词")
    # 3 burstiness（CV 越低越像 AI）
    cv = m["cv"]
    add("低 burstiness", 20 * max(0.0, min(1.0, (0.60 - cv) / 0.25)), 20,
        f"CV={cv:.3f}")
    # 4 句首小品词
    add("句首小品词", 10 * min(1.0, m["ssp"] / 12.0), 10, f"{m['ssp']:.1f}%")
    # 5 AI 高发词
    add("AI 高发词", 10 * min(1.0, per_k(m["aiw"]) / 3.0), 10,
        f"{m['aiw']} 次 = {per_k(m['aiw']):.2f}/千词")
    # 6 分号：只罚句间分号；列举性分号（并列统计项）是学术 MR 写作必需，不计
    add("句间分号密度", 10 * min(1.0, per_k(m["semi_cl"]) / 4.0), 10,
        f"句间 {m['semi_cl']} 处 = {per_k(m['semi_cl']):.1f}/千词"
        f"（总 {m['semi']} 处，其余为列举性）")
    # 7 被动语态偏离（10-25% 视为正常区）
    p = m["pas"]
    dev = 0.0 if 10 <= p <= 25 else (min(abs(p - 10), abs(p - 25)) / 15.0)
    add("被动语态偏离", 10 * min(1.0, dev), 10, f"{p:.1f}%")
    # 8 词汇多样性：用 MATTR（长度无关）。学术英文 MATTR 通常 0.70-0.80。
    add("低词汇多样性", 10 * max(0.0, min(1.0, (0.72 - m["mattr"]) / 0.12)), 10,
        f"MATTR={m['mattr']:.3f} (学术英文参考 0.70-0.80)")

    return round(s, 1), detail


def main():
    print("=" * 78)
    print("【一】数字完整性：改前 vs 改后 数字多重集比对")
    print("=" * 78)
    before = io.open(BAK, encoding="utf-8").read()
    after = io.open(DOC, encoding="utf-8").read()

    def nums(t):
        return Counter(re.findall(r"\d+(?:\.\d+)?", t))

    cb, ca = nums(before), nums(after)
    removed = cb - ca
    added = ca - cb
    print(f"  改前数字总数 {sum(cb.values())}　改后 {sum(ca.values())}")
    print(f"  改后被移除的数字种类: {len(removed)}")
    for k, v in sorted(removed.items()):
        print(f"      - {k} x{v}")
    print(f"  改后新增的数字种类: {len(added)}")
    for k, v in sorted(added.items()):
        print(f"      + {k} x{v}")

    if not removed and not added:
        print("\n  => 数字多重集完全一致：润色未改动任何数值。")
        num_ok = True
    else:
        print("\n  => 存在数字差异，需逐条说明。")
        num_ok = False

    # 引用编号完整性
    rb = Counter(re.findall(r"\[(\d+(?:[\u2013\-]\d+)?(?:,\s*\d+)*)\]", before))
    ra = Counter(re.findall(r"\[(\d+(?:[\u2013\-]\d+)?(?:,\s*\d+)*)\]", after))
    print(f"\n  引用标记: 改前 {sum(rb.values())} 改后 {sum(ra.values())} "
          f"{'一致' if rb == ra else '不一致'}")

    print("\n" + "=" * 78)
    print("【二】AIGC 形式化评分（0-100，越高越像 AI）")
    print("=" * 78)
    mb, ma = metrics(BAK), metrics(DOC)
    sb, db = score(mb)
    sa, da = score(ma)

    print(f"{'指标':<20}{'改前':>12}{'改后':>12}")
    for k, lab in [("nw", "正文词数"), ("ns", "句数"), ("mu", "平均句长"),
                   ("cv", "句长 CV"), ("em", "em dash"), ("bold", "正文强调加粗"),
                   ("semi", "分号(总)"), ("semi_cl", "句间分号"), ("aiw", "AI 高发词"), ("ssp", "句首小品词 %"),
                   ("pas", "被动语态 %"), ("ttr", "TTR(参考)"),
                   ("mattr", "MATTR(长度无关)")]:
        fmt = (lambda v: f"{v:.3f}") if k in ("cv", "ttr", "mattr") else (
            lambda v: f"{v:.1f}") if k in ("mu", "ssp", "pas") else (lambda v: f"{v}")
        print(f"{lab:<20}{fmt(mb[k]):>12}{fmt(ma[k]):>12}")

    print(f"\n改前评分 {sb}/100   分项：")
    for n, v, c, w in db:
        print(f"    {n:<16} {v:>5}/{c}   {w}")
    print(f"\n改后评分 {sa}/100   分项：")
    for n, v, c, w in da:
        print(f"    {n:<16} {v:>5}/{c}   {w}")
    print(f"\n  => 降噪幅度 {sb - sa:+.1f} 分  ({sb} -> {sa})")

    print("\n" + "=" * 78)
    print("结论")
    print("=" * 78)
    print(f"  数字完整性: {'PASS' if num_ok else 'NEEDS REVIEW'}")
    print(f"  AIGC 形式化评分: {sb} -> {sa}")

    op = os.path.join(ROOT, "code", "output", "AUDIT_134_humanize_verify.json")
    io.open(op, "w", encoding="utf-8").write(
        ("{" + '"num_removed": %d, "num_added": %d, "score_before": %s, '
         '"score_after": %s, "metrics_before": %s, "metrics_after": %s' + "}")
        % (len(removed), len(added), sb, sa,
           str(mb).replace("'", '"'), str(ma).replace("'", '"')))
    print(f"\n-> {op}")


if __name__ == "__main__":
    main()
