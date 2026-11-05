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
    # NOTE (2026-10-21): these four sub-site betas are written WITHOUT the 'β=' prefix
    # in Results ("lumbar-spine +0.095 (p=7.7x10^-9)"), so they are matched with their
    # site label AND trailing p-value -- this also tightens the check.
    ("ls beta +0.095", r"lumbar-spine \+0\.095 \(p=7\.7×10⁻⁹\)"),
    ("arms beta +0.075", r"arms \+0\.075 \(p=1\.7×10⁻⁶\)"),
    ("legs beta +0.060", r"legs \+0\.060 \(p=2\.3×10⁻⁵\)"),
    ("fn beta +0.058", r"femoral-neck \+0\.058 \(p=6\.4×10⁻⁴\)"),
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

# NEW (2026-10-30): externally sourced values added by the R1 theory-deepening round.
# These are absent from HEAD by construction, so they cannot go through the identical-
# between-versions path above (they would read as MISMATCH). They are therefore
# protected a different way: they MUST be present in the working tree, which catches
# accidental deletion or transcription drift. Both were verified against source:
#   Faber 2025 HMG 34(3):207-217, HSM2 -> hip fracture OR(IVW) 1.27 [1.12-1.44], P=1.79e-4
#   Hardcastle 2014 Arthritis Rheumatol 66(9):2429-39, enthesophytes grade>=1 in HBM:
#       adjusted OR 3.00 [1.96-4.58]
#
# NEW (2026-10-30, R2/R3 review round). Same construction: these are either transcribed
# external values or recomputed summary counts that are absent from HEAD by construction.
#   Cui 2022 Adv Rheumatol 62:22 (PMID 35761320), knee replacement in bisphosphonate
#       users vs non-users: RR 0.71 [0.66-0.77]; adjusted HR 0.76 [0.70-0.83].
#   Colocalisation threshold sensitivity, recomputed from coloc_screen_all_loci.csv by
#       code/187: knee-OA counts 0/31 at PP.H4>=0.8 and 1/31 at >=0.7.
external = [
    ("Faber HSM2 hip-fracture OR", r"HSM2 OR 1\.27, 1\.12–1\.44"),
    ("Hardcastle enthesophyte OR", r"OR 3\.00, 95% CI 1\.96–4\.58"),
    ("Cui bisphosphonate knee-replacement RR", r"RR 0\.71, 0\.66–0\.77"),
    ("coloc 0.7-threshold sensitivity", r"1/31 even at a 0\.7 threshold"),
]
for label, pat in external:
    mw = re.search(pat, work)
    if not mw:
        print(f"[ABSENT external] {label}: pattern {pat!r}")
        flags += 1
    else:
        print(f"[OK external] {label}: {mw.group(0)!r}")

# ======================================================================
# NEW (2026-11-02, Finding D): derived-OR internal consistency.
# Gate 149b previously compared beta/SE/p but NOT the derived OR column,
# which let the Finding-A sign flip (OA_knee->fracture OR=0.973=exp(-beta)
# while beta=+0.027) pass undetected. We now verify, within the WORKING tree,
# that every reported OR equals exp(|IVW beta|) and the 95% CI equals
# exp(|beta| +/- 1.96*SE). This is a self-consistency check, independent of
# HEAD, so it catches a flip even if both versions happened to share it.
import math

def _to_float(s):
    return float(s.strip().replace("−", "-").replace("–", "-").replace("—", "-"))

def _row_or(cells):
    # Table-2 layout after split('|'):
    #   [0]'' [1]exp->out [2]n_snp [3]IVW beta [4]SE [5]p [6]p_rand [7]OR(95% CI) ...
    try:
        beta = _to_float(cells[3])
        se   = _to_float(cells[4])
        orc  = cells[7].strip()
    except (IndexError, ValueError):
        return None
    if orc in ("", "—", "-"):
        return None
    m = re.match(r"([\d.]+)\s*\(([\d.]+)[-]([\d.]+)\)", orc.replace("−", "-").replace("–", "-"))
    if not m:
        return None
    return beta, se, float(m.group(1)), float(m.group(2)), float(m.group(3))

or_flags = 0
print()
print("=" * 96)
print("DERIVED-OR CONSISTENCY (Finding D): OR must equal exp(|IVW beta|), "
      "CI = exp(|beta| +/- 1.96 SE)")
print("=" * 96)
for line in work.splitlines():
    if not line.strip().startswith("|") or "→" not in line:
        continue
    res = _row_or(line.split("|"))
    if res is None:
        continue
    beta, se, or_rep, ci_lo, ci_hi = res
    b = abs(beta)
    or_exp = math.exp(b)
    lo_exp = math.exp(b - 1.96 * se)
    hi_exp = math.exp(b + 1.96 * se)
    label = line.split("|")[1].strip()
    ok = (abs(or_exp - or_rep) <= 0.015 and
          abs(lo_exp - ci_lo) <= 0.030 and
          abs(hi_exp - ci_hi) <= 0.030)
    if not ok:
        print(f"[OR MISMATCH] {label}: beta={beta:+.3f} exp(|b|)={or_exp:.4f} "
              f"reported OR={or_rep} (CI {ci_lo}-{ci_hi} vs exp {lo_exp:.3f}-{hi_exp:.3f})")
        or_flags += 1
    else:
        print(f"[OK] {label}: OR={or_rep} ~ exp(|{beta:+.3f}|)={or_exp:.4f}")
flags += or_flags

# explicit: 3.6-5.8 must still appear >=1 (line 69 retained)
n658 = len(re.findall(r"3\.6–5\.8", work))
print(f"\n'3.6–5.8' occurrences in working: {n658} (expected >=1; line 91 trim is intentional)")

print("\nRESULT:", "ALL REPORTED VALUES IDENTICAL -> SAFE" if flags == 0 else f"{flags} MISMATCH/ABSENT -> INVESTIGATE")
sys.exit(1 if flags else 0)
