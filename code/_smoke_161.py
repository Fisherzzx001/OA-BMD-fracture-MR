#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
code/_smoke_161.py — 161 的端到端冒烟测试（合成数据，跑完即弃，绝不进稿件）

目的：真实 T1 文件（744/919 MB）到达前，先用合成数据跑通 161 全流程，
      暴露运行期缺陷（缺列 / 键名 / 内存 / 崩溃），避免自动链空跑 5 小时后才发现。

做法：
  1. 从真实 KneeOA 前缀取真实 header（17 列，Chromosome/Position 在第 16/17 列）
  2. 取 BMD 10 位点的全部工具并集作为"必须命中"的位点，写入合成行
  3. 另加随机填充行模拟千万行规模下的匹配行为
  4. 复制 161 -> code/_161_smoke.py，仅改 DATA/OUT 指向临时目录
  5. 运行并捕获 traceback
"""
import gzip, io, os, random, shutil, subprocess, sys, tempfile

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)
import importlib
import numpy as np
import pandas as pd
m09 = importlib.import_module("09_rg_matrix")
m13 = importlib.import_module("13_mr_analysis")

PY = sys.executable
REAL = r"G:/OA_B_PROJECT_cache/KneeOA_sexspecific.tsv.gz"

# ---- 1. 真实 header ----
with gzip.open(REAL, "rt") as f:
    header = f.readline().rstrip("\n")
cols = header.split("\t")
print(f"[smoke] real header cols = {len(cols)}")
print(f"[smoke] Chromosome idx={cols.index('Chromosome')} Position idx={cols.index('Position')}")

# ---- 2. 工具并集 ----
sites = [k for k in m09.PHENOS if k.startswith("BMD_")]
inst_rows, seen = [], set()
for s in sites:
    inst = m13.select_instruments(m13.load_munged(s))
    for _, r in inst.iterrows():
        k = (int(r.CHR), int(r.BP))
        if k in seen:
            continue
        seen.add(k)
        inst_rows.append((int(r.CHR), int(r.BP), str(r.A1).upper(), str(r.A2).upper()))
print(f"[smoke] union instruments = {len(inst_rows)} across {len(sites)} BMD sites")

tmp = tempfile.mkdtemp(prefix="smoke161_")
data = os.path.join(tmp, "data", "gwas_stratified")
out = os.path.join(tmp, "out")
os.makedirs(data, exist_ok=True); os.makedirs(out, exist_ok=True)

rng = np.random.default_rng(20261018)
N_FILLER = 20000

def synth(path, seed, n_inst_use):
    """写合成性别分层文件：真实 header + 工具行 + 随机填充行。"""
    r = np.random.default_rng(seed)
    rows = []
    use = inst_rows[:n_inst_use]
    for (c, p, a1, a2) in use:            # 工具行：给真实可匹配的等位基因与信号
        m_or = float(np.exp(r.normal(0.010, 0.004)))
        f_or = float(np.exp(r.normal(0.014, 0.004)))
        m_p = float(np.clip(r.uniform(1e-6, 0.4), 1e-300, 0.999999))
        f_p = float(np.clip(r.uniform(1e-6, 0.4), 1e-300, 0.999999))
        row = {c_: "" for c_ in cols}
        row.update({"Marker": f"rsSMOKE{c}_{p}", "Chromosome": c, "Position": p,
                    "EffectAllele": a1, "AlternateAllele": a2,
                    "Male.EffectSize": m_or, "Male.Pvalue": m_p,
                    "Female.EffectSize": f_or, "Female.Pvalue": f_p})
        rows.append(row)
    for i in range(N_FILLER):             # 填充行
        c = int(r.integers(1, 23)); p = int(r.integers(1, 250_000_000))
        row = {c_: "" for c_ in cols}
        row.update({"Marker": f"rsFILL{i}", "Chromosome": c, "Position": p,
                    "EffectAllele": r.choice(["A", "C", "G", "T"]),
                    "AlternateAllele": r.choice(["A", "C", "G", "T"]),
                    "Male.EffectSize": float(np.exp(r.normal(0, 0.01))),
                    "Male.Pvalue": float(r.uniform(1e-4, 1)),
                    "Female.EffectSize": float(np.exp(r.normal(0, 0.01))),
                    "Female.Pvalue": float(r.uniform(1e-4, 1))})
        rows.append(row)
    df = pd.DataFrame(rows, columns=cols)
    for c_ in ("Chromosome", "Position"):
        df[c_] = df[c_].astype(int)
    df.to_csv(path, sep="\t", index=False, compression="gzip")
    return len(df)

for label, fn, seed, nuse in (("KneeOA", "KneeOA_sexspecific.tsv.gz", 11, len(inst_rows)),
                              ("AllOA",  "AllOA_sexspecific.tsv.gz",  22, len(inst_rows))):
    n = synth(os.path.join(data, fn), seed, nuse)
    print(f"[smoke] wrote {fn}: {n} rows")

# ---- 3. 复制 161，仅改 DATA/OUT ----
src = open(os.path.join(ROOT, "161_tier1_sex_mr.py"), encoding="utf-8").read()
smoke_path = os.path.join(ROOT, "_161_smoke.py")
sd = src.replace('DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "gwas_stratified")',
                 f'DATA = r"{data}"')
sd = sd.replace('OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "genetics")',
                f'OUT = r"{out}"')
if "DATA = r\"" not in sd or "OUT = r\"" not in sd:
    print("[smoke] !! patch failed — check constant lines"); sys.exit(2)
open(smoke_path, "w", encoding="utf-8").write(sd)

# ---- 4. 运行 ----
print("\n===== RUN SMOKE 161 =====")
pr = subprocess.run([PY, smoke_path], capture_output=True, text=True, cwd=os.path.join(ROOT, ".."))
print(pr.stdout[-4000:])
print("----- STDERR -----")
print(pr.stderr[-4000:])
print(f"===== exit={pr.returncode} =====")

# ---- 5. 产物检查 ----
sj = os.path.join(out, "tier1_sex_mr_summary.json")
if os.path.exists(sj):
    import json
    d = json.load(open(sj))
    print("\n[smoke] summary keys:", list(d.keys()))
    for oc in d:
        print(f"  {oc}: meta keys={list(d[oc].get('meta', {}).keys())}")
    print("SMOKE_RESULT=PASS")
else:
    print("SMOKE_RESULT=FAIL (no summary json)")

# ---- 6. 清理（临时目录 + 冒烟副本）----
shutil.rmtree(tmp, ignore_errors=True)
try:
    os.remove(smoke_path)
except OSError:
    pass
print("[smoke] temp artifacts removed")
