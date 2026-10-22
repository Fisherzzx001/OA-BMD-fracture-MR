# -*- coding: utf-8 -*-
"""
OSF DOI 到手后的一次性回写（五步）—— 幂等 + 断言保护

用法：
  python code/184_osf_doi_writeback.py 10.17605/OSF.IO/XXXXX
  python code/184_osf_doi_writeback.py https://doi.org/10.17605/OSF.IO/XXXXX

执行内容：
 1. 校验 DOI 格式
 2. docs/08 与 docs/cover_letter_OAC.md："being deposited" -> "available at <DOI>"
 3. 重建 Cover_Letter.docx / Manuscript.docx / Submission_MainText.docx / Manuscript_Ready/
 4. 复跑 code/182_final_qc.py 终检（49 项）
 5. 更新 docs/128 待办勾选

安全设计：
 - 每处替换要求 count == 1，不符即中止且不写盘
 - 不触碰任何 β/SE/p/OR/CI（改动仅限 OSF 措辞与 DOI）
 - 可用 --dry-run 预览
"""
import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
VENV_PY = Path(r"C:/Users/Administrator/.workbuddy/binaries/python/envs/default/Scripts/python.exe")
PLAIN_PY = Path(r"C:/Users/Administrator/.workbuddy/binaries/python/versions/3.13.12/python.exe")

MANUSCRIPT = ROOT / "docs" / "08_投稿稿_OA_BMD_fracture.md"
COVER = ROOT / "docs" / "cover_letter_OAC.md"
CHECKLIST = ROOT / "docs" / "128_投稿前终检总表.md"

DOI_RE = re.compile(r"^10\.17605/OSF\.IO/[A-Z0-9]{5,10}$", re.IGNORECASE)


def die(msg: str) -> "None":
    print(f"  ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def normalize_doi(raw: str) -> str:
    s = raw.strip().rstrip(".,;")
    s = re.sub(r"^https?://(dx\.)?doi\.org/", "", s, flags=re.IGNORECASE)
    s = re.sub(r"^https?://osf\.io/", "", s, flags=re.IGNORECASE)
    if re.fullmatch(r"[A-Z0-9]{5,10}", s, re.IGNORECASE):  # 只给了 GUID
        s = "10.17605/OSF.IO/" + s
    s = s.replace("10.17605/osf.io/", "10.17605/OSF.IO/")
    if not DOI_RE.match(s):
        die(f"DOI 格式不合法：{raw!r}（应为 10.17605/OSF.IO/XXXXX）")
    return s


def guarded_replace(path: Path, old: str, new: str, dry: bool) -> None:
    text = path.read_text(encoding="utf-8")
    n = text.count(old)
    if n != 1:
        die(f"{path.name}: 目标串出现 {n} 次（要求恰好 1 次），已中止未写盘\n     目标串: {old[:70]}...")
    if not dry:
        path.write_text(text.replace(old, new), encoding="utf-8")
    print(f"  [OK] {path.name}: {old[:48]!r} -> 新措辞")


def run(cmd: list, label: str, cwd: Path, dry: bool) -> None:
    if dry:
        print(f"  [SKIP dry-run] {label}")
        return
    r = subprocess.run([str(c) for c in cmd], cwd=str(cwd), capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-1500:], r.stderr[-800:])
        die(f"{label} 失败（退出码 {r.returncode}）")
    print(f"  [OK] {label}")


def main() -> int:
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    dry = "--dry-run" in sys.argv
    if not args:
        die("缺少参数：python code/184_osf_doi_writeback.py 10.17605/OSF.IO/XXXXX")
    doi = normalize_doi(args[0])
    doi_url = f"https://doi.org/{doi}"

    print("=== OSF DOI 回写 ===")
    print(f"  DOI: {doi}")
    print(f"  模式: {'DRY-RUN（不写盘）' if dry else 'WRITE'}")

    old_ms = "An OSF mirror is being deposited."
    new_ms = f"An OSF mirror is available at {doi_url}."
    old_cl = "An OSF mirror is being deposited."
    new_cl = f"An OSF mirror is available at {doi_url}."

    # ---- step 2: 真源措辞回写
    guarded_replace(MANUSCRIPT, old_ms, new_ms, dry)
    guarded_replace(COVER, old_cl, new_cl, dry)

    # 反向守卫：不得再出现旧措辞
    if not dry:
        for p in (MANUSCRIPT, COVER):
            if "being deposited" in p.read_text(encoding="utf-8"):
                die(f"{p.name}: 仍残留 'being deposited'")
    print("  [OK] 无 'being deposited' 残留")

    # ---- step 3: 重建派生产物
    py = VENV_PY if VENV_PY.exists() else PLAIN_PY
    run([py, ROOT / "build_submission.py"], "Manuscript.docx", ROOT, dry)
    run([py, ROOT / "build_maintext.py"], "Submission_MainText.docx", ROOT, dry)
    run([py, ROOT / "code" / "180_cover_letter_md_to_docx.py"], "Cover_Letter.docx", ROOT, dry)
    run([py, ROOT / "code" / "135_build_manuscript_ready.py"], "Manuscript_Ready/", ROOT, dry)

    # ---- step 4: 终检
    run([PLAIN_PY, ROOT / "code" / "182_final_qc.py"], "code/182 终检 49 项", ROOT, dry)

    # ---- step 5: 勾选表更新
    if CHECKLIST.exists() and not dry:
        t = CHECKLIST.read_text(encoding="utf-8")
        n = t.count("- [ ] **OSF 镜像回写**")
        if n == 1:
            t = t.replace("- [ ] **OSF 镜像回写**", "- [x] **OSF 镜像回写**")
            CHECKLIST.write_text(t, encoding="utf-8")
            print("  [OK] docs/128 待办已勾选")
        else:
            print(f"  [WARN] docs/128 未找到唯一待勾选行（n={n}），请手工勾选")
    elif dry:
        print("  [SKIP dry-run] docs/128 勾选")

    print(f"\n=== 回写完成：{doi} ===")
    if not dry:
        print("  提交建议：git add docs/ code/ Manuscript_Ready/ && git commit")
    return 0


if __name__ == "__main__":
    sys.exit(main())
