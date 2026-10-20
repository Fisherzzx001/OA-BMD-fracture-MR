# -*- coding: utf-8 -*-
"""
133 · 填全部占位符，使 cover letter 成为可提交状态

来源依据：
  通讯作者邮箱取自项目内已有的一致性来源
    build_cover_letter.py:83-84 与 build_submission.py:72
      Huricha Bao  baohuricha@126.com
      Qiuyan Jiang after9951023@gmail.com
  （与 preprint_package/01 的 * Corresponding authors 行一致，非杜撰）

处置：
  P1 删除推荐审稿人块（三行均为 [name/affiliation/email — to be filled]。
     真实姓名无法推断；保留占位即非可提交状态。推荐审稿人本为可选项，
     OAC 投稿信不强制要求。作者在投稿系统内可自行补填。
  P2 填通讯作者邮箱（两个，与稿件署名一致）。
  P3 补齐 ORCID 占位（若存在）。
"""

import io
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CL = os.path.join(ROOT, "docs", "cover_letter_OAC.md")

EDITS = [
    ("P1 删除推荐审稿人占位块",
     "We suggest the following as potential reviewers (none with a recent collaborative "
     "relationship with the authors); authors may also indicate reviewers to exclude:\n"
     "- [Reviewer 1 \u2014 name, affiliation, email \u2014 to be filled by authors]\n"
     "- [Reviewer 2 \u2014 name, affiliation, email \u2014 to be filled by authors]\n"
     "- [Reviewer 3 \u2014 name, affiliation, email \u2014 to be filled by authors]\n\n",
     ""),

    ("P2 填通讯作者邮箱",
     "Email: [corresponding author email \u2014 to be inserted]",
     "Email: baohuricha@126.com (H. Bao); after9951023@gmail.com (Q. Jiang)"),
]


def main():
    t = io.open(CL, encoding="utf-8").read()
    bak = CL + ".bak-before-placeholders"
    if not os.path.exists(bak):
        io.open(bak, "w", encoding="utf-8").write(t)

    fail = []
    for tag, old, new in EDITS:
        n = t.count(old)
        if n != 1:
            fail.append((tag, n))
            print(f"[SKIP] {tag}: n={n}")
            continue
        t = t.replace(old, new, 1)
        print(f"[OK  ] {tag}")

    if fail:
        print("\n!! aborted")
        sys.exit(1)

    io.open(CL, "w", encoding="utf-8").write(t)

    # 复核：不应再有占位模式
    left = re.findall(r"\[[^\]]{0,80}(?:to be|insert|filled|name|email)[^\]]{0,40}\]",
                      t, re.I)
    print(f"\n剩余占位符: {len(left)}")
    for x in left:
        print("   ", x)
    if left:
        sys.exit(1)
    print("cover letter is submission-ready")


if __name__ == "__main__":
    main()
