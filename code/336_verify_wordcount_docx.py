# -*- coding: utf-8 -*-
"""
code/336 - Independent OAC word-count verification against the ACTUAL docx.

Rationale (2026-11-11, six-item terminal audit):
    code/144 and code/148_oac_compliance shared a narrow tokeniser
    `re.findall(r"[A-Za-z][A-Za-z0-9'\\-]*", s)` that requires a token to START
    with an ASCII letter. It silently drops every numeric / bracketed /
    symbol-initial token (0.075, beta=+0.097, [1];, (OR, 1.40,, 38,472, 13x13),
    under-counting this MR manuscript by ~246 tokens. OAC counts the WORD
    document, where a whitespace-delimited token containing a letter or digit
    counts as one word. This gate counts the real submission docx with a
    Word-like tokeniser and FAILS if the IMRaD body exceeds 4000 words.

    It also emits a numeric-token multiset snapshot so that any compression can
    be verified to leave EVERY reported statistic untouched.

Usage:  python3 code/336_verify_wordcount_docx.py
"""
import re, json, os, sys
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MD = os.path.join(ROOT, "docs", "08_投稿稿_OA_BMD_fracture.md")
DOCXS = [
    os.path.join(ROOT, "submission_package", "Manuscript.docx"),
    os.path.join(ROOT, "submission_package", "Submission_MainText.docx"),
]
LIMIT = 4000


def words_wordlike(s):
    return [t for t in s.split() if re.search(r"[^\W_]", t, re.UNICODE)]


def words_narrow(s):
    return re.findall(r"[A-Za-z][A-Za-z0-9'\-]*", s)


def numeric_tokens(s):
    return Counter(t for t in s.split() if re.search(r"\d", t))


def md_body(text):
    ai = text.find("**Objective.**")
    ci = text.find("**Conclusions.**", ai)
    tail = text[ci:]
    em = re.search(r"\n---", tail)
    ae = ci + (em.start() if em else len(tail))
    leg = text.find("## Figure legends")
    body = text[ae:leg]
    return "\n".join(l for l in body.splitlines() if not l.strip().startswith("|"))


def docx_imrad(path):
    from docx import Document
    from docx.oxml.ns import qn
    from docx.text.paragraph import Paragraph

    doc = Document(path)
    blocks = []
    for child in doc.element.body.iterchildren():
        if child.tag == qn("w:p"):
            blocks.append(("p", Paragraph(child, doc).text))
        elif child.tag == qn("w:tbl"):
            blocks.append(("t", None))
    start = end = None
    for i, (k, tx) in enumerate(blocks):
        s = tx.strip() if tx else ""
        if k == "p" and s.startswith("1. Introduction") and start is None:
            start = i
        if k == "p" and s.startswith("Figure legends"):
            end = i
            break
    if start is None or end is None:
        raise RuntimeError(f"cannot locate IMRaD boundaries in {path}")
    parts = [tx for k, tx in blocks[start:end] if k == "p"]  # tables excluded (OAC)
    return "\n".join(parts)


def main():
    mdtext = open(MD, encoding="utf-8").read()
    mdb = md_body(mdtext)
    wl = len(words_wordlike(mdb))
    nw = len(words_narrow(mdb))
    print("=== code/336 OAC WORD-COUNT VERIFICATION (Word-like) ===")
    print(f"[markdown docs/08] Word-like={wl}  narrow(144/148)={nw}  "
          f"delta=+{wl - nw} numeric tokens dropped by narrow regex")
    ok = wl <= LIMIT
    for p in DOCXS:
        if not os.path.exists(p):
            print(f"[{os.path.basename(p)}] MISSING")
            ok = False
            continue
        b = docx_imrad(p)
        w = len(words_wordlike(b))
        n = len(words_narrow(b))
        print(f"[{os.path.basename(p)}] Word-like={w}  narrow={n}  "
              f"delta=+{w - n}")
        if w > LIMIT:
            ok = False
    snap = {
        "source": "docs/08 IMRaD body",
        "word_like": wl,
        "narrow_144_148": nw,
        "delta_dropped_numeric_tokens": wl - nw,
        "limit": LIMIT,
        "numeric_token_multiset": dict(numeric_tokens(mdb)),
    }
    out = os.path.join(ROOT, "code", "output", "wordcount_numeric_snapshot.json")
    json.dump(snap, open(out, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print(f"numeric snapshot -> {os.path.relpath(out, ROOT)} "
          f"({len(snap['numeric_token_multiset'])} distinct numeric tokens)")
    print("OAC WORD-LIKE COMPLIANT:",
          "YES" if ok else f"NO - IMRaD body exceeds {LIMIT} words")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
