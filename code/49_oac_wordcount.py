#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
OAC word-count auditor (corrected).

Osteoarthritis and Cartilage Guide for Authors:
    Full length original research article: maximum 4000 words.
    The limit EXCLUDES title page, abstract, tables, figure legends,
    acknowledgements, author contributions, and references.

So the OAC-counted body = Introduction + Methods + Results + Discussion prose
only. Tables/figure legends/refs/title-page/author-contrib are reported as
excluded (informational). Also checks: abstract <=250 words; Highlights
3-5 bullets each <=85 chars.

Usage:  python code/49_oac_wordcount.py [markdown_file]
"""
import re
import sys
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DEFAULT = os.path.join(ROOT, "docs", "08_投稿稿_OA_BMD_fracture.md")

TOKEN = re.compile(r"[A-Za-z0-9][A-Za-z0-9\u00b2\u00b3\-\+/\.]*")


def words(text: str) -> int:
    return len([t for t in text.split() if re.search(r"[0-9A-Za-z]", t)])


def strip_tables(line: str) -> str:
    s = line.strip()
    if s.startswith("|"):
        return ""
    return line


def parse(path: str):
    with open(path, "r", encoding="utf-8") as fh:
        lines = fh.read().split("\n")
    secs, cur, buf = [], None, []
    for ln in lines:
        m = re.match(r"^(#{1,6})\s+(.*)$", ln)
        if m:
            if cur is not None:
                secs.append((cur, buf))
            cur, buf = m.group(2).strip(), []
        else:
            if cur is not None:
                buf.append(ln)
    if cur is not None:
        secs.append((cur, buf))
    return secs


IMRAD_PREFIX = ("1.", "2.", "3.", "4.")


def main(path: str):
    secs = parse(path)
    print("=" * 74)
    print("OAC WORD-COUNT AUDIT  (body = IMRaD prose only; tables/figs/refs excluded)")
    print("file:", os.path.relpath(path, ROOT))
    print("=" * 74)

    body_total = 0
    excluded = 0
    abstract_words = None
    highlights = []

    for title, body in secs:
        prose = "\n".join(strip_tables(l) for l in body)
        n = words(prose)
        if title.lower().startswith("abstract"):
            abstract_words = n
            print(f"{n:>7}  [EXCLUDED] {title[:60]}")
            continue
        if title.lower().startswith("highlights"):
            # collect bullet lines: "- text" or "* text" (single asterisk + space)
            for l in body:
                ls = l.strip()
                if re.match(r"^[-*]\s", ls):
                    txt = re.sub(r"^[-*]\s", "", ls).strip()
                    if txt:
                        highlights.append(txt)
            print(f"{'':>7}  [CHECK] Highlights: {len(highlights)} bullets")
            continue
        is_imrad = any(title.startswith(p) for p in IMRAD_PREFIX)
        if is_imrad:
            body_total += n
            print(f"{n:>7}  [COUNTED ] {title[:64]}")
        else:
            excluded += n
            tag = "figs/tables/refs/title/author"
            print(f"{n:>7}  [EXCLUDED] {title[:58]:<58} ({tag})")

    print("-" * 74)
    print(f"{body_total:>7}  OAC-COUNTED BODY (Intro+Methods+Results+Discussion)")
    print(f"{4000:>7}  OAC full-length limit (4000)")
    print(f"{4000 - body_total:>7}  headroom")
    if abstract_words is not None:
        flag = "OK" if abstract_words <= 250 else "OVER"
        print(f"{abstract_words:>7}  Abstract words (limit 250)  [{flag}]")
    # Highlights check
    print("-" * 74)
    h_ok = True
    for i, h in enumerate(highlights, 1):
        c = len(h)
        ok = c <= 85
        h_ok = h_ok and ok
        print(f"  H{i}: {c:>3} chars {'OK ' if ok else 'OVER'}  {h}")
    hn = len(highlights)
    hn_ok = 3 <= hn <= 5
    print(f"  Highlights count = {hn} (must be 3-5): {'OK' if hn_ok else 'BAD'}")
    print("=" * 74)
    overall = (body_total <= 4000) and (abstract_words is None or abstract_words <= 250) \
        and h_ok and hn_ok
    print("OVERALL OAC COMPLIANT:" , "YES" if overall else "NO — see flags above")
    print("=" * 74)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else DEFAULT)
