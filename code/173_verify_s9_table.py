# -*- coding: utf-8 -*-
"""核验稿件 S9 表全部列 vs mvmr_height_main.csv"""
import csv

with open('code/output/genetics/mvmr_height_main.csv', newline='', encoding='utf-8') as f:
    for r in csv.DictReader(f):
        ref = float(r['ref_uni_beta']); mv = float(r['mvmr_beta'])
        a = 100.0 * (ref - mv) / ref
        print("%-18s -> %-9s uni(T2)=%+.4f mvmr=%+.4f atten=%+.1f%% egger_p=%.3f I2=%.2f "
              "FcB=%.1f FcH=%.1f nB=%s nH=%s mvmr_h=%+.4f mvmr_h_p=%.1e"
              % (r['exposure'], r['outcome'], ref, mv, a,
                 float(r['egger_int_p']), float(r['I2']),
                 float(r['Fc_bmd']), float(r['Fc_h']), r['n_bmd'], r['n_h'],
                 float(r['mvmr_h']), float(r['mvmr_h_p'])))
