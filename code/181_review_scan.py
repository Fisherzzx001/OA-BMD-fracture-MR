#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Pre-review committee defect scan for docs/08 (OA-BMD-fracture MR ms).
Focus: citation integrity, author-year/[n] consistency, duplicate sentences,
and Table 2 vs deposited mr_results_corrected.csv numeric agreement.
No statistic is changed by this script; it only REPORTS.
"""
import re, csv, io, sys, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MS = os.path.join(ROOT, "docs", "08_投稿稿_OA_BMD_fracture.md")
CORR = os.path.join(ROOT, "code", "output", "genetics", "mr_results_corrected.csv")

def load():
    with io.open(MS, encoding="utf-8") as f:
        return f.read()

text = load()

# ---- split body (Abstract Conclusions -> Figure legends) ----
body_start = text.index("## 1. Introduction")
body_end = text.index("## Figure legends")
body = text[body_start:body_end]

# ============ 1. inline numeric citations, order, coverage ============
# capture [n], [n,m], [n–m] ranges etc.
cite_pat = re.compile(r"\[(\d+(?:[–\-]\d+)?(?:,\s*\d+)*)\]")
cites = cite_pat.findall(body)
def expand(c):
    nums=[]
    for part in c.split(","):
        part=part.strip()
        if "–" in part or "-" in part:
            a,b=re.split(r"[–-]",part)
            nums+=list(range(int(a),int(b)+1))
        else:
            nums.append(int(part))
    return nums
all_nums=[]
for c in cites:
    all_nums+=expand(c)
# ascending monotonic check (allow equal/forward)
desc=[]
prev=0
for n in all_nums:
    if n<prev:
        desc.append((prev,n))
    prev=max(prev,n)
print("== 1. CITATION ORDER ==")
print("  total numeric citations (expanded):", len(all_nums))
print("  descending steps (order violations):", desc if desc else "NONE")
# coverage: every 1..38 cited?
cited=set(all_nums)
missing_cited=[n for n in range(1,39) if n not in cited]
print("  refs 1..38 NOT cited anywhere in body:", missing_cited if missing_cited else "NONE (all 38 cited)")
# undefined citations
undefined=[n for n in all_nums if n<1 or n>38]
print("  citations outside 1..38:", undefined if undefined else "NONE")

# ============ 2. References list integrity ============
ref_sec = text[text.index("## References"):text.index("## Acknowledgements")]
ref_lines = [l for l in ref_sec.splitlines() if re.match(r"^\d+\.\s", l.strip())]
ref_nums=[int(re.match(r"^(\d+)\.", l.strip()).group(1)) for l in ref_lines]
print("\n== 2. REFERENCES LIST ==")
print("  count:", len(ref_nums), "sequence 1..N exact:", ref_nums==list(range(1,len(ref_nums)+1)))
# first author of each ref
ref_auth={}
for l in ref_lines:
    m=re.match(r"^(\d+)\.\s+([A-Z][a-z]+)", l.strip())
    if m: ref_auth[int(m.group(1))]=m.group(2)
print("  refs parsed with author:", len(ref_auth), "/", len(ref_lines))

# ============ 3. author-year + [n] consistency ============
print("\n== 3. AUTHOR-YEAR + [n] CONSISTENCY ==")
ay_pat = re.compile(r"([A-Z][a-z]+(?:-[A-Z][a-z]+)?)\s+et al\.\s*\[(\d+)\]")
mism=[]
for m in ay_pat.finditer(body):
    name=m.group(1); num=int(m.group(2))
    if num in ref_auth and ref_auth[num].lower()!=name.lower():
        mism.append((name,num,ref_auth.get(num)))
print("  author-year name vs [n] author mismatches:", mism if mism else "NONE")
# author-year cited WITHOUT bracket
ay_nob = re.compile(r"([A-Z][a-z]+(?:-[A-Z][a-z]+)?)\s+et al\.\s*(?!\[)")
floating=[]
for m in ay_nob.finditer(body):
    # exclude the STROBE guideline inline author-year (Skrivankova)
    nm=m.group(1)
    if nm in ("Skrivankova",): continue
    floating.append(nm)
print("  'X et al.' WITHOUT [n] (floating, excluding Skrivankova):", floating if floating else "NONE")

# ============ 4. duplicate sentences ============
print("\n== 4. DUPLICATE / NEAR-DUP SENTENCES ==")
# strip tables & figure legends & refs for sentence scan in prose sections
prose_sections = [
 ("Intro", text[text.index("## 1. Introduction"):text.index("## 2. Methods")]),
 ("Methods", text[text.index("## 2. Methods"):text.index("## 3. Results")]),
 ("Results", text[text.index("## 3. Results"):text.index("## 4. Discussion")]),
 ("Discussion", text[text.index("## 4. Discussion"):text.index("## Figure legends")]),
]
def norm(s):
    s=re.sub(r"\[[^\]]*\]","",s)        # drop citations
    s=re.sub(r"[^a-z0-9 ]","",s.lower())
    return s
sent_cache={}
for secname,sec in prose_sections:
    sents=re.split(r"(?<=[.)])\s+(?=[A-Z(])", sec)
    for s in sents:
        s=s.strip()
        if len(s.split())<8: continue
        n=norm(s)
        if not n: continue
        if n in sent_cache:
            print(f"  DUP in [{secname}]: \"{s[:90]}...\"")
        else:
            sent_cache[n]=secname
# also cross-section near-dups via token-set jaccard
import itertools
alln=[(secname,s,norm(s)) for secname,sec in prose_sections
      for s in re.split(r"(?<=[.)])\s+(?=[A-Z(])",sec) if len(s.split())>=12]
def toks(n): return set(n.split())
print("  (cross-section jaccard >=0.6):")
for (a,sa,na),(b,sb,nb) in itertools.combinations(alln,2):
    ta,tb=toks(na),toks(nb)
    if not ta or not tb: continue
    j=len(ta&tb)/len(ta|tb)
    if j>=0.6 and a!=b:
        print(f"    [{sa}] <-> [{sb}] j={j:.2f}: \"{sa[:60]}\"")

# ============ 5. Table 2 vs mr_results_corrected.csv ============
print("\n== 5. TABLE 2 vs mr_results_corrected.csv ==")
# parse Table 2 block
t2 = text[text.index("| Exposure → Outcome"):text.index("Univariate instrument strength")]
rows=[l for l in t2.splitlines() if l.strip().startswith("|") and "---" not in l and "Exposure" not in l and "n_snp" not in l]
def parse_row(l):
    c=[x.strip() for x in l.strip().strip("|").split("|")]
    return c
t2rows={}
for r in rows:
    c=parse_row(r)
    # c[0] exposure->outcome (strip **)
    key=re.sub(r"\*\*","",c[0]).strip()
    t2rows[key]=c
# map names
def keymap(exp,out):
    return f"{exp},{out}"
csv_rows={}
with io.open(CORR, encoding="utf-8") as f:
    rd=csv.DictReader(f)
    for row in rd:
        csv_rows[f"{row['exposure']},{row['outcome']}"]=row
# exposure/outcome naming translation
trans = {
 "BMD total → fracture (positive control)":("BMD_total","fracture"),
 "BMD femoral neck → fracture":("BMD_femur_neck","fracture"),
 "BMD total → OA knee":("BMD_total","OA_knee"),
 "Lumbar spine → OA knee":("BMD_lumbar","OA_knee"),
 "Arms → OA knee":("BMD_arms","OA_knee"),
 "Legs → OA knee":("BMD_leg","OA_knee"),
 "Femoral neck → OA knee":("BMD_femur_neck","OA_knee"),
 "OA knee → fracture":("OA_knee","fracture"),
 "OA knee → BMD total":("OA_knee","BMD_total"),
 "OA knee → BMD femoral neck":("OA_knee","BMD_femur_neck"),
}
mism_t2=[]
for label,(e,o) in trans.items():
    t2c=t2rows[label]
    csvr=csv_rows[f"{e},{o}"]
    # table columns: 1 n_snp,2 ivw beta,3 se,4 p,5 re p,6 OR(CI),7 egger beta(se),8 egger p,9 egger int p,10 wmed(ci),11 I2,12 loo
    tb=float(t2c[2].replace("−","-").replace("–","-")); tse=float(t2c[3].replace("−","-").replace("–","-")); tp=t2c[4]; tn=int(t2c[1])
    cb=float(csvr["ivw"]); cse=float(csvr["ivw_se"]); cp=csvr["ivw_p"]; cn=int(csvr["n_snp"])
    def close(a,b,tol=1e-3): return abs(a-b)<=tol
    diffs=[]
    if tn!=cn: diffs.append(f"n_snp {tn}!={cn}")
    if not close(tb,cb,5e-3): diffs.append(f"ivw {tb}!={cb}")
    if not close(tse,cse,5e-3): diffs.append(f"se {tse}!={cse}")
    # p compare in sci
    def tosci(p):
        try: return float(p.replace("−","-").replace("–","-"))
        except: return None
    pta=tosci(tp); pca=tosci(cp)
    if pta and pca and abs(pta-pca)/max(abs(pca),1e-100)>0.02: diffs.append(f"p {tp}!={cp}")
    if diffs: mism_t2.append((label,diffs))
print("  Table2 rows checked:",len(trans))
print("  DISCREPANCIES (Table2 vs corrected csv):", mism_t2 if mism_t2 else "NONE — exact agreement")
