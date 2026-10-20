# -*- coding: utf-8 -*-
"""从重建的 rg_results_long_13.csv 实算 Table S15 全部 5 行 + L65 段落统计"""
import csv
from statistics import mean

with open('code/output/genetics/rg_results_long_13.csv', newline='', encoding='utf-8') as f:
    rows = list(csv.DictReader(f))

print("columns:", list(rows[0].keys()))

def block_of(t1, t2):
    s = {t1, t2}
    if s == {'OA_knee', 'OA_selfrep'}:
        return 'OA x OA'
    if 'OA_knee' in s or 'OA_selfrep' in s:
        if 'fracture' in s:
            return 'OA x Fracture'
        return 'OA x BMD'
    if 'fracture' in s:
        return 'Fracture x BMD'
    return 'BMD x BMD'

blocks = {}
for r in rows:
    t1, t2 = r['trait1'], r['trait2']
    if t1 == t2:
        continue
    b = block_of(t1, t2)
    blocks.setdefault(b, []).append(r)

for b in ['OA x OA', 'OA x BMD', 'Fracture x BMD', 'BMD x BMD', 'OA x Fracture']:
    v = blocks[b]
    rgs = [float(x['rg']) for x in v]
    sig = sum(1 for x in v if float(x['q_fdr']) < 0.05)
    print(f"{b:16s} n={len(v):3d}  range [{min(rgs):+.3f}, {max(rgs):+.3f}]  mean={mean(rgs):+.3f}  q<0.05: {sig}/{len(v)}")

# L65 统计
allp = [r for r in rows if r['trait1'] != r['trait2']]
n_bh = sum(1 for r in allp if float(r['q_fdr']) < 0.05)
n_bonf = sum(1 for r in allp if float(r['p_bonf']) < 0.05)
print(f"\n78 pairs: BH-FDR sig {n_bh}, Bonferroni sig {n_bonf}")
nonsig = [(r['trait1'], r['trait2'], float(r['rg']), float(r['p'])) for r in allp if float(r['q_fdr']) >= 0.05]
print("non-BH pairs:", nonsig)

# OA x Fracture 明细
print("\nOA x Fracture detail:")
for r in blocks['OA x Fracture']:
    print(f"  {r['trait1']} x {r['trait2']}: rg={float(r['rg']):+.4f} se={float(r['se']):.4f} p={float(r['p']):.4g} q={float(r['q_fdr']):.4g}")

# Fracture x BMD 明细（前5+后5）
print("\nFracture x BMD detail (sorted by rg):")
for r in sorted(blocks['Fracture x BMD'], key=lambda x: float(x['rg'])):
    print(f"  {r['trait1']:14s} x {r['trait2']:14s}: rg={float(r['rg']):+.4f} q={float(r['q_fdr']):.3g}")
