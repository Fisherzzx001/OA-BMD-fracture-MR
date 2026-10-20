# -*- coding: utf-8 -*-
# Curated decisive check: every REPORTED statistic value must be byte-identical
# between the committed docs/08 and the working tree. Citation/section/table-ref
# numbers and my intentional "3.6-5.8-fold" line-91 trim are excluded from the
# must-match set (the value 3.6-5.8 is still present at line 69).
import re, subprocess, sys

def load(path):
    return subprocess.check_output(["git", "show", path]).decode("utf-8") if path.startswith("HEAD:") \
        else open(path, encoding="utf-8").read()

head = load("HEAD:docs/08_投稿稿_OA_BMD_fracture.md")
work = load("docs/08_投稿稿_OA_BMD_fracture.md")

# (label, regex) -- first match in each version is compared exactly.
patterns = [
    ("tb beta +0.078", r"β=\+0\.078"),
    ("ls beta +0.095", r"β=\+0\.095"),
    ("arms beta +0.075", r"β=\+0\.075"),
    ("legs beta +0.060", r"β=\+0\.060"),
    ("fn beta +0.058", r"β=\+0\.058"),
    ("fracture IVW beta -0.337", r"β=−0\.337"),
    ("fracture fn beta -0.351", r"β=−0\.351"),
    ("fracture OR 1.40", r"OR 1\.40, 95% CI 1\.34–1\.46"),
    ("fracture fn OR 1.42", r"OR 1\.42, 1\.34–1\.51"),
    ("GEFOS OR 1.55", r"OR 1\.55, 1\.48–1\.63"),
    ("BMD->kneeOA OR range", r"OR 1\.060–1\.099"),
    ("p 3.9e-10", r"p=3\.9×10⁻¹⁰"),
    ("p 7.7e-9", r"p=7\.7×10⁻⁹"),
    ("p 1.7e-6", r"p=1\.7×10⁻⁶"),
    ("p 2.3e-5", r"p=2\.3×10⁻⁵"),
    ("p 6.4e-4", r"p=6\.4×10⁻⁴"),
    ("p 3.7e-51", r"p=3\.7×10⁻⁵¹"),
    ("rg knee +0.14 to +0.19", r"\+0\.14 to \+0\.19"),
    ("rg fracture -0.40 to -0.54", r"−0\.40 to −0\.54"),
    ("rg selfrep +0.84", r"r_g=\+0\.84"),
    ("rg knee-fracture +0.119", r"r_g=\+0\.119"),
    ("PP.H4 0.799", r"PP\.H4=0\.799"),
    ("PP.H4 0.999-1.000", r"PP\.H4=0\.999–1\.000"),
    ("PP.H4 WNT16 0.002", r"PP\.H4=0\.002"),
    ("scale factor 3.389", r"1/√\[K\(1−K\)\] = 3\.389"),
    ("K 0.0964", r"K=0\.0964"),
    ("I2 68-86", r"I²=68–86"),
    ("power 0.14-0.26", r"power 0\.14–0\.26"),
    ("coloc 9/31", r"9/31"),
    ("coloc 0/31", r"0/31"),
    ("Hartley hip OR 1.28", r"OR 1\.28, 1\.05–1\.57"),
    ("Hartley knee OR 1.40", r"OR 1\.40, 1\.20–1\.63"),
    ("Funck-Brentano OR 1.14", r"OR 1\.14, 1\.06–1\.22"),
    ("theta_H +0.145", r"θ_H=\+0\.145"),
    ("N per site 31873", r"N=31,873"),
    ("intercept scale F 3.389 alt", r"3\.389"),
]

flags = 0
for label, pat in patterns:
    mh = re.search(pat, head)
    mw = re.search(pat, work)
    h = mh.group(0) if mh else None
    w = mw.group(0) if mw else None
    if h is None and w is None:
        print(f"[ABSENT both] {label}: pattern {pat!r}")
        flags += 1
    elif h != w:
        print(f"[MISMATCH] {label}: HEAD={h!r}  WORK={w!r}")
        flags += 1
    else:
        print(f"[OK] {label}: {w!r}")

# explicit: 3.6-5.8 must still appear >=1 (line 69 retained)
n658 = len(re.findall(r"3\.6–5\.8", work))
print(f"\n'3.6–5.8' occurrences in working: {n658} (expected >=1; line 91 trim is intentional)")

print("\nRESULT:", "ALL REPORTED VALUES IDENTICAL -> SAFE" if flags == 0 else f"{flags} MISMATCH/ABSENT -> INVESTIGATE")
sys.exit(1 if flags else 0)
