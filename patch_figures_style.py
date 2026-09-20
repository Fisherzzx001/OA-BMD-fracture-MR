# -*- coding: utf-8 -*-
"""patch_figures_style.py
对 build_figures.py / build_figures_mr.py / build_figures_v9.py 做一次性风格补丁，
使其符合 Nature 系列 + Elsevier(Bone) 图片规范：
  (1) CB 调色板统一为 Okabe-Ito 色盲安全配色
  (2) 轴标签/刻度字体统一 >=7pt（原 6.0/6.2/6.5 -> 7.0；5.5/5.6/5.8 -> 6.5）
  (3) save() 增补 300 dpi TIFF 导出（Elsevier 生产首选栅格格式；PDF 仍为矢量）
不改动任何数据读取与统计量。
"""
import os, io

TARGETS = ["build_figures.py", "build_figures_mr.py", "build_figures_v9.py"]

# (1) 颜色 -> Okabe-Ito 色盲安全
COLORS = {
    "#4878CF": "#0072B2",   # blue
    "#D65F5F": "#D55E00",   # red / vermillion
    "#6ACC65": "#009E73",   # green
    "#EE854A": "#E69F00",   # orange
    "#956CB4": "#CC79A7",   # purple / red-violet
    "#82C6E2": "#56B4E9",   # teal / sky blue
    "#8C8C8C": "#999999",   # grey
    "#222222": "#000000",   # ink
    "#D8D8D8": "#DDDDDD",   # faint
}

# (2) 字体 -> >=7pt
FONTS = {
    "fontsize=6.0": "fontsize=7.0",
    "fontsize=6.2": "fontsize=7.0",
    "fontsize=6.5": "fontsize=7.0",
    "labelsize=6.5": "labelsize=7.0",
    "fontsize=5.5": "fontsize=6.5",
    "fontsize=5.6": "fontsize=6.5",
    "fontsize=5.8": "fontsize=6.5",
}

# (3) TIFF 导出：在每处 "name + \".pdf\"))" 后追加 tif 行
TIFF_ADD = 'name + ".pdf"))'

def patch_file(path):
    with open(path, encoding="utf-8") as f:
        s = f.read()
    orig = s
    for a, b in COLORS.items():
        s = s.replace(a, b)
    for a, b in FONTS.items():
        s = s.replace(a, b)
    # 仅在 save() 内增补 TIFF（避免重复添加：若已含 .tif 则跳过）
    if '.tif"), dpi=300' not in s:
        s = s.replace(
            TIFF_ADD,
            TIFF_ADD + '\n    fig.savefig(os.path.join(OUT, name + ".tif"), dpi=300)',
        )
    if s != orig:
        with open(path, "w", encoding="utf-8") as f:
            f.write(s)
        print("patched:", path)
    else:
        print("no change:", path)

if __name__ == "__main__":
    for t in TARGETS:
        if os.path.exists(t):
            patch_file(t)
        else:
            print("MISSING", t)
