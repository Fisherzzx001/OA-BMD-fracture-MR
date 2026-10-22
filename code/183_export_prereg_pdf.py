# -*- coding: utf-8 -*-
"""
docs/118 OSF 预注册 -> PDF（上传 OSF 的保底格式）

为什么需要：OSF 在线渲染 Markdown 表格常错位；注册一经提交不可修改，
上传附件时留一份排版稳定的 PDF 作为 canonical copy。

链路：markdown -> HTML(带 A4 排版 CSS) -> Edge headless print-to-pdf
本机无 pandoc/wkhtmltopdf，Edge 是唯一可靠的 HTML->PDF 通道。

用法： python code/183_export_prereg_pdf.py
"""
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "docs" / "118_OSF预注册_OA_BMD_fracture.md"
OUT = ROOT / "submission_package" / "OSF_Preregistration.pdf"

EDGE_CANDIDATES = [
    Path(r"C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe"),
    Path(r"C:/Program Files/Microsoft/Edge/Application/msedge.exe"),
]

CSS = """
@page { size: A4; margin: 20mm 18mm; }
body {
  font-family: Georgia, "Times New Roman", serif;
  font-size: 10.5pt; line-height: 1.5; color: #111;
}
h1 { font-size: 16pt; margin: 0 0 6pt; line-height: 1.3; }
h2 { font-size: 13pt; margin: 16pt 0 6pt; border-bottom: 1px solid #ccc; padding-bottom: 3pt; }
h3 { font-size: 11.5pt; margin: 12pt 0 4pt; }
p  { margin: 0 0 7pt; }
ul, ol { margin: 0 0 8pt; padding-left: 18pt; }
li { margin-bottom: 3pt; }
code {
  font-family: Consolas, "Courier New", monospace;
  font-size: 9pt; background: #f4f4f4; padding: 1pt 3pt;
}
pre {
  background: #f4f4f4; padding: 8pt; overflow-wrap: break-word;
  font-size: 8.5pt; line-height: 1.35; page-break-inside: avoid;
}
table { border-collapse: collapse; width: 100%; margin: 8pt 0; font-size: 9pt; }
th, td { border: 1px solid #999; padding: 4pt 6pt; text-align: left; vertical-align: top; }
th { background: #eee; font-weight: bold; }
blockquote {
  margin: 8pt 0; padding: 4pt 10pt; border-left: 3px solid #aaa;
  color: #333; font-style: italic;
}
hr { border: 0; border-top: 1px solid #ccc; margin: 14pt 0; }
a { color: #1155aa; text-decoration: underline; word-break: break-all; }
strong { font-weight: bold; }
tr, li, p { page-break-inside: avoid; }
h2, h3 { page-break-after: avoid; }
"""


def build_html(src: Path) -> Path:
    import markdown  # noqa: E402
    text = src.read_text(encoding="utf-8")
    body = markdown.markdown(
        text,
        extensions=["extra", "sane_lists"],
        output_format="html5",
    )
    html = (
        "<!DOCTYPE html>\n<html><head><meta charset=\"utf-8\">"
        f"<title>OSF Preregistration</title><style>{CSS}</style></head>"
        f"<body>\n{body}\n</body></html>"
    )
    tmpdir = tempfile.mkdtemp(prefix="prereg_pdf_")
    html_path = Path(tmpdir) / "prereg.html"
    html_path.write_text(html, encoding="utf-8")
    return html_path


def print_to_pdf(html_path: Path, out_pdf: Path) -> bool:
    edge = next((p for p in EDGE_CANDIDATES if p.exists()), None)
    if edge is None:
        print("  ERROR: 未找到 Edge 可执行文件", file=sys.stderr)
        return False
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    if out_pdf.exists():
        out_pdf.unlink()
    url = html_path.as_uri()
    workdir = tempfile.mkdtemp(prefix="edge_profile_")
    for headless_mode in ("--headless=new", "--headless"):
        cmd = [
            str(edge), headless_mode, "--disable-gpu",
            "--no-pdf-header-footer", "--no-first-run", "--no-extensions",
            f"--user-data-dir={workdir}",
            f"--print-to-pdf={out_pdf.resolve()}", url,
        ]
        try:
            subprocess.run(cmd, timeout=120, check=False,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.TimeoutExpired:
            pass
        if out_pdf.exists() and out_pdf.stat().st_size > 1000:
            return True
    return False


def main() -> int:
    print("=== docs/118 -> PDF ===")
    if not SRC.exists():
        print(f"  ERROR: 源文件不存在 {SRC}", file=sys.stderr)
        return 1

    text = SRC.read_text(encoding="utf-8")
    banned = ["[PENDING", "to be stamped", "as per manuscript", "XXXX", "23400861"]
    hits = [b for b in banned if b in text]
    if hits:
        print(f"  ERROR: 源文件仍含占位符/过期 DOI {hits} —— PDF 会冻结这些内容，先行修复",
              file=sys.stderr)
        return 1
    print("  [OK] 源文件无占位符 / 无过期版本 DOI")

    html_path = build_html(SRC)
    n_tables = len(re.findall(r"<table>", html_path.read_text(encoding="utf-8")))
    print(f"  [OK] HTML 生成，含 {n_tables} 个表格")

    if not print_to_pdf(html_path, OUT):
        print("  ERROR: PDF 生成失败", file=sys.stderr)
        return 1
    size_kb = OUT.stat().st_size / 1024
    with OUT.open("rb") as f:
        head = f.read(5)
    if head != b"%PDF-":
        print("  ERROR: 输出不是合法 PDF", file=sys.stderr)
        return 1
    print(f"  [OK] PDF 已生成：{OUT.relative_to(ROOT)}  ({size_kb:.0f} KB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
