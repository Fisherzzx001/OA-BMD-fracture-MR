#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
331_verify_plink_bed_decode.py — REGRESSION GUARD for PLINK1 .bed 2-bit decoding
================================================================================
Background (2026-11-10): a wrong decode mapped het/homA2 dosages backwards
(g==2 -> 2.0, g==3 -> 0.0), silently corrupting every LD matrix / r^2 clump
that used the 1000G EUR reference (code/325 coloc-SuSiE, code/328 Chai clump).
The error was caught only because rs12568930's "A2 frequency" came out 0.276
(= heterozygote rate 2p(1-p)) instead of matching its GWAS eaf 0.8197.

CORRECT PLINK1 bed coding (2 bits/sample, low bits = first sample):
    00 = hom A1  -> A2 dosage 0
    01 = missing -> NaN
    10 = het     -> A2 dosage 1
    11 = hom A2  -> A2 dosage 2

This guard validates BOTH decoding idioms used in the repo:
  (a) numpy unpackbits(bitorder="little"): code = bit_high*2 + bit_low
  (b) pure-python shift:                   code = (byte >> (2*k)) & 3
and asserts they agree SNP-by-SNP, and that each instrument's 1000G A2 allele
frequency equals its GWAS eaf OR its complement (GWAS effect allele may be
1000G A1 or A2), to tolerance. A het/homA2 swap makes A2freq ~= 2p(1-p) and
breaks this for most SNPs, so the guard fails loudly instead of silently.

Usage: python3 code/331_verify_plink_bed_decode.py
Exit 0 = PASS, exit 1 = FAIL (CI gate friendly).
"""
import json, sys
import numpy as np

BED = "G:/REF/1000G_EUR/EUR.bed"
BIM = "G:/REF/1000G_EUR/EUR.bim"
INST = "G:/OA_B_PROJECT_cache/l7/bmd_instruments.json"
N = 503                       # 1000G EUR samples
BY = (N + 3) // 4             # bytes per SNP row
TOL = 0.05                    # 1000G(503) vs GWAS sampling tolerance


def load_bim():
    rs2 = {}
    for i, line in enumerate(open(BIM)):
        p = line.rstrip("\n").split("\t")
        if len(p) < 6:
            continue
        rs2[p[1]] = (i, p[4].upper(), p[5].upper())
    return rs2


def raw_row(fh, row):
    fh.seek(3 + row * BY)
    b = fh.read(BY)
    if len(b) < BY:
        raise IOError(f"short read at row {row}")
    return b


def code_python(raw):
    g = np.zeros(N, dtype=np.uint8)
    for bi in range(BY):
        byte = int(raw[bi])
        for k in range(4):
            j = bi * 4 + k
            if j >= N:
                break
            g[j] = (byte >> (2 * k)) & 3
    return g


def code_numpy(raw):
    bits = np.unpackbits(np.frombuffer(raw, dtype=np.uint8), bitorder="little")
    c = bits[: N * 2].reshape(N, 2)
    return c[:, 1].astype(np.uint8) * 2 + c[:, 0].astype(np.uint8)


def a2_dosage(code):
    d = np.where(code == 3, 2.0, np.where(code == 2, 1.0, 0.0)).astype(float)
    d[code == 1] = np.nan
    return d


def main():
    rs2 = load_bim()
    inst = json.load(open(INST, encoding="utf-8"))["inst"]
    checked = agree = eaf_ok = 0
    fails = []
    with open(BED, "rb") as fh:
        for rs, v in inst.items():
            if rs not in rs2:
                continue
            row, a1, a2 = rs2[rs]
            raw = raw_row(fh, row)
            cp = code_python(raw)
            cn = code_numpy(raw)
            checked += 1
            if not np.array_equal(cp, cn):
                fails.append(f"{rs}: decoders disagree ({np.mean(cp==cn):.3f})")
                continue
            agree += 1
            d = a2_dosage(cn)
            ok = ~np.isnan(d)
            if ok.sum() < 0.9 * N:
                fails.append(f"{rs}: too many missing ({N-ok.sum()})")
                continue
            a2f = float(np.nanmean(d)) / 2.0
            eu = v.get("eaf")
            if eu is None:
                continue
            err = min(abs(a2f - eu), abs(a2f - (1.0 - eu)))
            if err < TOL:
                eaf_ok += 1
            else:
                fails.append(f"{rs}: A2freq={a2f:.3f} vs eaf={eu:.3f}/1-eaf={1-eu:.3f} (min err {err:.3f})")

    print(f"[331] checked {checked} instruments in 1000G EUR ({N} samples)")
    print(f"[331] python vs numpy decoder identical: {agree}/{checked}")
    print(f"[331] A2freq matches eaf or complement (<{TOL}): {eaf_ok}/{checked}")
    if fails:
        print("[331] FAILURES:")
        for x in fails[:20]:
            print("   ", x)
        print("VERDICT: FAIL — bed decoding is WRONG; do NOT trust any LD matrix/clump from it")
        sys.exit(1)
    # Require near-unanimous agreement: a swapped het/homA2 decode fails most SNPs.
    if agree < checked or eaf_ok < 0.95 * checked:
        print("VERDICT: FAIL — decoder agreement or eaf concordance too low")
        sys.exit(1)
    print("VERDICT: PASS — PLINK bed 2-bit decoding verified against GWAS eaf")
    sys.exit(0)


if __name__ == "__main__":
    main()
