# -*- coding: utf-8 -*-
"""
124 · 给门禁补上此前完全缺失的检查项

背景（本轮审计发现，非推测）：
  本轮发现 cover letter 与 Table S11 的 MVMR 数值落后于稿件一整个版本
  （6\u201333%/67\u201394% vs 稿件的 9.5\u201333.2%/67\u201391%），而门禁对此**没有任何检查**。
  门禁此前只校验稿件内部与重算数据的一致性，不校验"附属文件 vs 稿件"。
  这是"门禁覆盖盲区"，比"门禁判据陈旧"更危险 —— 因为错误会静默通过。

新增门禁 §6.7：附属文件 vs 稿件 的关键数值一致性
  G1  cover letter 的 MVMR 衰减范围 = 稿件
  G2  Table S11 (01_predictions_and_falsification.md) 的衰减范围 = 稿件
  G3  Table S11 阳性对照衰减上界 = 稿件
  G4  稿件不残留审计前旧值（6\u201333% / 67\u201394% / 22\u201337 / 118）
  G5  稿件不残留空洞自贬词（provisional / cannot be fully excluded / cannot decompose）
  G6  Highlights 不含把非显著当零效应的 "were null"
"""

import io
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GATE = os.path.join(ROOT, "audit_v18_full.py")

BLOCK = '''
# ============================================================ 6.7 附属文件 vs 稿件 一致性
# 本轮新增：此前门禁不校验附属文件，导致 cover letter / Table S11 落后稿件一整个版本而静默通过。
print("\\n[\\u00a76.7] \\u9644\\u5c5e\\u6587\\u4ef6 vs \\u7a3f\\u4ef6 \\u6570\\u503c\\u4e00\\u81f4\\u6027")
_anc = {
    "cover letter": os.path.join(PROJ, "docs", "cover_letter_OAC.md"),
    "Table S11": os.path.join(PROJ, "docs", "Supplementary",
                              "01_predictions_and_falsification.md"),
}
_anc_txt = {}
for _k, _p in _anc.items():
    _anc_txt[_k] = io.open(_p, encoding="utf-8").read() if os.path.exists(_p) else ""

_ck = lambda n, ok, d: ck(n, ok, d)

# --- G1/G2: MVMR \\u8870\\u51cf\\u8303\\u56f4\\u5fc5\\u987b\\u4e0e\\u7a3f\\u4ef6\\u540c\\u6b65 ---
ck("\\u00a76.7 G1 cover letter MVMR \\u8303\\u56f4\\u4e0e\\u7a3f\\u4ef6\\u4e00\\u81f4",
   "9.5\\u201333.2%" in _anc_txt.get("cover letter", "")
   and "67\\u201391%" in _anc_txt.get("cover letter", ""),
   "cover letter \\u5e94\\u4e3a 9.5\\u201333.2% / 67\\u201391%")

ck("\\u00a76.7 G2 Table S11 \\u8870\\u51cf\\u8303\\u56f4\\u4e0e\\u7a3f\\u4ef6\\u4e00\\u81f4",
   "9.5\\u201333.2%" in _anc_txt.get("Table S11", ""),
   "Table S11 P4 \\u884c\\u5e94\\u4e3a 9.5\\u201333.2%")

ck("\\u00a76.7 G3 Table S11 \\u9633\\u6027\\u5bf9\\u7167\\u4e0a\\u754c\\u4e0e\\u7a3f\\u4ef6\\u4e00\\u81f4",
   "\\u22121.7% to \\u22123.0%" in _anc_txt.get("Table S11", ""),
   "Table S11 \\u5e94\\u4e3a \\u22121.7% to \\u22123.0%")

# --- G4: \\u65e7\\u503c\\u6b8b\\u7559 ---
for _old, _lab in [("6\\u201333%", "MVMR 6-33%"), ("67\\u201394%", "MVMR 67-94%"),
                   ("22\\u201337", "\\u5de5\\u5177\\u6570 22-37"), ("118 for knee OA", "\\u5de5\\u5177\\u6570 118")]:
    ck(f"\\u00a76.7 G4 \\u7a3f\\u4ef6\\u4e0d\\u6b8b\\u7559\\u65e7\\u503c {_lab}",
       _old not in md, f"\\u627e\\u5230\\u6b8b\\u7559: {_lab}")
    for _k, _t in _anc_txt.items():
        ck(f"\\u00a76.7 G4 {_k} \\u4e0d\\u6b8b\\u7559\\u65e7\\u503c {_lab}",
           _old not in _t, f"{_k} \\u6b8b\\u7559: {_lab}")

# --- G5: \\u7a7a\\u6d1e\\u81ea\\u8d2c\\u8bcd ---
for _w in ["provisional", "cannot be fully excluded", "cannot decompose",
           "better described as absence of evidence"]:
    ck(f"\\u00a76.7 G5 \\u7a3f\\u4ef6\\u4e0d\\u542b\\u7a7a\\u6d1e\\u81ea\\u8d2c \\u201c{_w}\\u201d",
       _w not in md, f"\\u7a3f\\u4ef6\\u6b8b\\u7559: {_w}")

# --- G6: Highlights \\u4e0d\\u5e94\\u628a\\u975e\\u663e\\u8457\\u5f53\\u96f6\\u6548\\u5e94 ---
_hl = md.split("## Highlights")[1].split("**Keywords**")[0] if "## Highlights" in md else ""
ck("\\u00a76.7 G6 Highlights \\u4e0d\\u542b \\u201cwere null\\u201d\\uff08\\u975e\\u663e\\u8457\\u5f53\\u96f6\\u6548\\u5e94\\uff09",
   "were null" not in _hl, "Highlights \\u5e94\\u7528 non-significant")

# --- G7: \\u5b9e\\u8d28\\u62ab\\u9732\\u4e0d\\u5f97\\u88ab\\u5220\\uff08\\u9632\\u6b62\\u4e0b\\u8f6e\\u6539\\u52a8\\u8d70\\u5411\\u53e6\\u4e00\\u6781\\u7aef\\uff09 ---
for _must in ["0.14\\u20130.26", "I\\u00b2 68\\u201386%", "9.5\\u201333.2%",
              "MR-Egger", "positive control"]:
    ck(f"\\u00a76.7 G7 \\u5b9e\\u8d28\\u62ab\\u9732\\u4fdd\\u7559 \\u201c{_must}\\u201d",
       _must in md, f"\\u7a3f\\u4ef6\\u7f3a\\u5c11\\u5fc5\\u8981\\u62ab\\u9732: {_must}")

# --- G8: \\u9633\\u6027\\u5bf9\\u7167\\u5f02\\u8d28\\u6027\\u5bf9\\u7167\\u5df2\\u8865\\u5165\\uff08\\u672c\\u8f6e\\u6838\\u5fc3\\u53d1\\u73b0\\uff09 ---
ck("\\u00a76.7 G8 \\u5df2\\u8865\\u5165\\u9633\\u6027\\u5bf9\\u7167\\u5f02\\u8d28\\u6027\\u5bf9\\u7167",
   "positive control was comparably heterogeneous" in md,
   "\\u9633\\u6027\\u5bf9\\u7167 I\\u00b2 78/83% \\u4e0e BMD\\u2192OA 68\\u201386% \\u540c\\u91cf\\u7ea7\\uff1b\\u5fc5\\u987b\\u663e\\u5f0f\\u5199\\u5165")
'''


def main():
    t = io.open(GATE, encoding="utf-8").read()
    anchor = "# ============================================================ 7. 输出"
    if anchor not in t:
        print("anchor not found")
        sys.exit(1)
    if "6.7 附属文件 vs 稿件 一致性" in t:
        print("already patched; skip")
        return
    bak = GATE + ".bak-before-gate124"
    if not os.path.exists(bak):
        io.open(bak, "w", encoding="utf-8").write(t)
    t = t.replace(anchor, BLOCK.lstrip("\n") + "\n\n" + anchor, 1)
    io.open(GATE, "w", encoding="utf-8").write(t)
    print("patched audit_v18_full.py: added §6.7 ancillary-consistency gate (G1-G8)")


if __name__ == "__main__":
    main()
