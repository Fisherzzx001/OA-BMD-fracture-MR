# -*- coding: utf-8 -*-
"""
123 · 同步「比稿件更陈旧」的附属文件 + 修复门禁陈旧判据

审计发现（均为 121/122 之外的独立发现）：
  S1  docs/cover_letter_OAC.md L20 仍为 "6\u201333% / 67\u201394%"，
      稿件自 2026-10-15 审计后已改为 9.5\u201333.2% / 67\u201391%。
      -> cover letter 是编辑看到的第一份文件，数值与稿件冲突 = 直接扣分。
  S2  docs/Supplementary/01_predictions_and_falsification.md (Table S11)：
        P4 行 "OA attenuation 6\u201333%" 陈旧
        P4 行 "fracture positive control moved only \u22121.7% to \u22122.8%" 陈旧（现 \u22123.0%）
        F2 行 "far beyond 33%" 陈旧（现 33.2%）
  S3  audit_v18_full.py L557\u2013561：门禁主判据为旧值 "-0.39"，且提示文字
      "§3.3 正文称 \u22120.39" 已经过时（稿件现为 \u22120.364）。
      -> 门禁靠 `or` 的备用分支通过；其 FAIL 消息文本本身是错的
         （会伪造一个"正文与 Table 3 不一致"的诊断）。
"""

import io
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

TARGETS = []


def T(path, tag, old, new, why):
    TARGETS.append((os.path.join(ROOT, path), tag, old, new, why))


T("docs/cover_letter_OAC.md", "S1a cover letter MVMR 范围",
  "Multivariable MR shows only 6\u201333% of the BMD\u2192OA effect is attributable to "
  "skeletal size (height proxy); 67\u201394% is independent of the size pathway",
  "Multivariable MR shows only 9.5\u201333.2% of the BMD\u2192OA effect is attributable "
  "to skeletal size (height proxy); 67\u201391% is independent of the size pathway",
  "cover letter 落后于稿件：6\u201333%/67\u201394% 为审计前旧值。"),

T("docs/Supplementary/01_predictions_and_falsification.md", "S2a P4 衰减范围",
  "OA attenuation 6\u201333% (vs Table 3); fracture positive control moved only "
  "\u22121.7% to \u22122.8% (mass-transparent)",
  "OA attenuation 9.5\u201333.2% (across the two tabulated baselines); fracture "
  "positive control moved only \u22121.7% to \u22123.0% (mass-transparent)",
  "Table S11 落后于稿件与主表 S9 脚注。"),

T("docs/Supplementary/01_predictions_and_falsification.md", "S2b F2 阈值",
  "If adding a true size phenotype attenuates the effect far beyond 33% \u2192 P4 "
  "residual is not \"bone\"",
  "If adding a true size phenotype attenuates the effect far beyond 33.2% \u2192 P4 "
  "residual is not \"bone\"",
  "F2 阈值应跟随稿件的 33.2% 上界。"),

T("audit_v18_full.py", "S3a 门禁主判据同步为 \u22120.364",
  'print(f"\\n[§3.3] total-body\u2192fracture \u52a0\u6743\u4e2d\u4f4d\u6570\uff1a\u4e25\u683c\u96c6 {wm_strict:+.4f} "\n'
  '      f"(\u00a73.3 \u6b63\u6587\u79f0 \u22120.39\uff1bTable 3 \u62a5 \u22120.364)")\n'
  'ck("\u00a73.3 weighted median \u4e0e Table 3 \u4e00\u81f4", format(wm_strict, ".2f") == "-0.39" or\n'
  '   format(wm_strict, ".3f") == "-0.364",\n'
  '   f"\u6b63\u6587 \u22120.39 \u4e0e Table 3 \u22120.364 \u4e0d\u4e00\u81f4\uff08\u4e25\u683c\u96c6\u91cd\u7b97 {wm_strict:+.4f}\uff09")',
  'print(f"\\n[\u00a73.3] total-body\u2192fracture \u52a0\u6743\u4e2d\u4f4d\u6570\uff1a\u4e25\u683c\u96c6 {wm_strict:+.4f} "\n'
  '      f"(\u7a3f\u4ef6 \u00a73.3 \u4e0e Table 3 \u5747\u62a5 \u22120.364)")\n'
  'ck("\u00a73.3 weighted median \u4e0e Table 3 \u4e00\u81f4", format(wm_strict, ".3f") == "-0.364",\n'
  '   f"\u91cd\u7b97 {wm_strict:+.4f} \u4e0e Table 3 \u22120.364 \u4e0d\u4e00\u81f4")',
  "门禁主判据 \u22120.39 与提示文字均为稿件旧值；现稿件报 \u22120.364。"
  "原 or 分支会掩盖真实失配，且 FAIL 消息会伪造\"正文 vs Table 3 不一致\"。"),


def main():
    failed = []
    for path, tag, old, new, why in TARGETS:
        if not os.path.exists(path):
            failed.append((tag, "file missing"))
            print(f"[MISS] {tag}: {path}")
            continue
        md = io.open(path, encoding="utf-8").read()
        n = md.count(old)
        if n != 1:
            failed.append((tag, f"occurrences={n}"))
            print(f"[SKIP] {tag}: occurrences={n} (expect 1)")
            continue
        bak = path + ".bak-before-sync123"
        if not os.path.exists(bak):
            io.open(bak, "w", encoding="utf-8").write(md)
        io.open(path, "w", encoding="utf-8").write(md.replace(old, new, 1))
        print(f"[OK  ] {tag}")

    if failed:
        print(f"\n!! {len(failed)} failed: {failed}")
        sys.exit(1)
    print(f"\napplied {len(TARGETS)}/{len(TARGETS)} ancillary/gate fixes")


if __name__ == "__main__":
    main()
