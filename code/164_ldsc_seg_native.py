#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
164_ldsc_seg_native.py — F3 组织级机制旁证：LDSC-SEG 细胞类型富集（原生 Python，无需 ldsc.py）
=============================================================================================
目标（docs/118 F3）：BMD 全 GWAS 汇总统计对 53 组织 LD-score 做 S-LDSC 回归，检验 BMD
遗传信号在成骨谱系（MSC/成骨/骨相关组织）是否比软骨谱系/其他组织更富集。
若富集系统性在软骨 → "骨驱动"机制被推翻。

原生实现（不依赖 ldsc.py / python2）—— 复用 docs/08 同族 IRWLS 思想：
    S-LDSC 模型（Finucane 2018 / Viana 2020）：
        χ²_k = a + b·ℓ_base,k + Σ_c τ_c·(N/M)·ℓ_{c,k}
    IRWLS 权重（官方 SEG.weights）：
        het_w = 1/[2·(a + b·(N/M)·ℓ_base + Σ_c τ_c·(N/M)·ℓ_c)²] ,  oc_w = 1/ℓ_base
    τ_c = 组织 c 的富集系数；富集 z_c = τ_c / se(τ_c)（SE 来自 (XᵀWX)⁻¹ 对角，权重即逆方差）。
F3 判别：成骨谱系组织平均 |z| / 平均 τ 是否 > 软骨/其他谱系。

数据（真实、公开，需下载 —— 本会话沙箱不可达，见下）
------------------------------------------------------
· LDSC-SEG 53-tissue 权重（Viana J, et al. Nat Genet 2020;52:986-1000）：
      https://alkesgroup.broadinstitute.org/LDSC_SEG/LDSC_SEG_20180515.tar.gz   (~2.4 GB)
    （data.broadinstitute.org/alkesgroup/LDSC_SEG/ 会 301 重定向到上述 alkesgroup 主机；
      解压后把每个组织的 *.ldscore.gz 放进 ref/ldsc_seg_53tissue/。）
· 基线 LD 计分：FLAGSHIP_LDSC 的 eur_w_ld_chr_*.l2.ldscore.gz（已缓存）。
· BMD 暴露：sumstats_total.sumstats.gz（munge 格式 CHR BP SNP A1 A2 MAF N Z，已缓存）。

诚实边界
--------
本会话实测：alkesgroup.broadinstitute.org 出口不可达（http 000），且沙箱出口带宽
~11 KB/s，2.4 GB 需 ~60 h → F3 **无法在本会话执行**。本脚本提供：
  · 经合成数据自测的原生 SEG 回归（trust-but-verify，证明 τ 恢复与富集 z 正确）；
  · 权重就位后一键运行真实 F3（python 164_ldsc_seg_native.py real）；
  · 权重下载与放置的确切命令（见上），供带宽充足环境复现。
**不虚构任何 F3 结果。**

用法：
  python 164_ldsc_seg_native.py selftest     # 合成数据验证原生 SEG 回归
  python 164_ldsc_seg_native.py real         # 权重就位时跑真实 F3
"""
import os, sys, glob, math
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "output", "genetics")
FLAGSHIP_LDSC = r"C:/Users/Administrator/WorkBuddy/2026-09-08-20-17-24/OA_FLAGSHIP/07_reference/ldsc"
SEG_DIR = os.path.join(ROOT, "ref", "ldsc_seg_53tissue")
BMD_SUMSTATS = os.path.join(FLAGSHIP_LDSC, "sumstats_total.sumstats.gz")

# 成骨/骨谱系组织关键词（SEG 组织名大小写不敏感匹配）；软骨/非骨作为对照
OSTEO_KEYWORDS = ["osteoblast", "bone", "msc", "mesench", "stromal", "chondro"]
CHONDRO_KEYWORDS = ["chondro", "cartilage"]


# ---------------------------------------------------------------- 原生 S-LDSC SEG 回归
def seg_regression(z, L2_base, L2_ctx, N, M, n_iter=12):
    """
    z: (n,) BMD z 分数；L2_base: (n,) 基线 LD 计分；L2_ctx: (n,C) 各组织 LD 计分。
    返回 (a, b, tau[C], se_tau[C], z_tau[C])。
    """
    y = z ** 2
    n, C = L2_ctx.shape
    ones = np.ones(n)
    X = np.column_stack([ones, L2_base, L2_ctx])        # 截距 + 基线 + C 组织
    w = 1.0 / np.maximum(L2_base, 1.0)                  # 初值 oc_w
    beta = None
    for _ in range(n_iter):
        XtWX = X.T @ (X * w[:, None])
        XtWy = X.T @ (w * y)
        beta = np.linalg.solve(XtWX, XtWy)
        a, b = beta[0], beta[1]
        tau = beta[2:]
        pred = a + (b * N / M) * L2_base + (N / M) * (L2_ctx @ tau)
        het_w = 1.0 / (2.0 * np.maximum(pred, 1e-6) ** 2)
        w = het_w / np.maximum(L2_base, 1.0)
    cov = np.linalg.inv(XtWX)
    se = np.sqrt(np.diag(cov))
    # 将回归系数换算为 S-LDSC 富集单位 τ_c（系数 × M/N）；z = τ/se 缩放不变
    tau = beta[2:] * (M / N)
    se_tau = se[2:] * (M / N)
    z_tau = tau / se_tau
    return float(beta[0]), float(beta[1]), tau, se_tau, z_tau


# ---------------------------------------------------------------- 自测（合成数据）
def run_selftest():
    rng = np.random.default_rng(20261010)
    n = 20000
    M = 1_000_000.0
    N = 500_000.0
    L2_base = rng.gamma(1.5, 1.0, n)
    C = 53
    L2_ctx = rng.gamma(1.0, 0.3, size=(n, C))           # 各组织 LD 计分（独立于基线）
    # 在"成骨"组织（索引 0,1）植入真实富集 τ=1.5；其余 τ=0
    tau_true = np.zeros(C)
    tau_true[0] = 1.5
    tau_true[1] = 1.2
    a, b = 1.0, 0.5
    pred = a + (b * N / M) * L2_base + (N / M) * (L2_ctx @ tau_true)
    y = pred + np.sqrt(np.maximum(pred, 1e-6)) * rng.standard_normal(n)   # χ² 方差模型
    z = np.sign(rng.standard_normal(n)) * np.sqrt(np.maximum(y, 0))
    a_h, b_h, tau, se_tau, z_tau = seg_regression(z, L2_base, L2_ctx, N, M)
    ok0 = abs(tau[0] - tau_true[0]) < 0.3 and z_tau[0] > 3
    ok1 = abs(tau[1] - tau_true[1]) < 0.3 and z_tau[1] > 3
    ok_bg = abs(tau[5]) < 3 * se_tau[5] + 0.2           # 背景组织不显著
    print("=" * 78)
    print("F3 NATIVE SEG SELFTEST (synthetic planted osteoblast enrichment)")
    print("=" * 78)
    print("osteoblast#0 planted tau=1.50 recovered=%.3f z=%.1f (ok=%s)" % (tau[0], z_tau[0], ok0))
    print("osteoblast#1 planted tau=1.20 recovered=%.3f z=%.1f (ok=%s)" % (tau[1], z_tau[1], ok1))
    print("background#5  recovered tau=%.3f z=%.1f (ok=%s)" % (tau[5], z_tau[5], ok_bg))
    assert ok0 and ok1 and ok_bg, "native SEG selftest failed"
    print("SELFTEST PASSED — native S-LDSC recovers planted enrichment & correct z")
    return True


# ---------------------------------------------------------------- 真实 F3
def _load_baseline_l2():
    ld = {}
    for c in range(1, 23):
        p = os.path.join(FLAGSHIP_LDSC, f"eur_w_ld_chr_{c}.l2.ldscore.gz")
        if not os.path.exists(p):
            raise FileNotFoundError(f"missing baseline LD: {p}")
        d = pd.read_csv(p, sep=r"\s+", usecols=["SNP", "L2"])
        for s, l in zip(d.SNP, d.L2):
            if l > 0:
                ld[s] = l
    return ld


def _load_munge_bmd():
    d = {}
    import gzip
    with gzip.open(BMD_SUMSTATS, "rt") as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 8:
                continue
            try:
                d[p[2]] = float(p[7])   # Z
            except ValueError:
                continue
    return d


def run_real():
    os.makedirs(OUT, exist_ok=True)
    files = sorted(glob.glob(os.path.join(SEG_DIR, "*.ldscore.gz")) +
                   glob.glob(os.path.join(SEG_DIR, "*.gp.gz")))
    if not files:
        print("LDSC-SEG weights NOT present under %s" % SEG_DIR)
        print("Download (bandwidth-adequate host):")
        print("  curl -O https://alkesgroup.broadinstitute.org/LDSC_SEG/LDSC_SEG_20180515.tar.gz")
        print("  tar xzf LDSC_SEG_20180515.tar.gz ; mv the per-tissue *.ldscore.gz -> %s" % SEG_DIR)
        return
    print("loading baseline LD + BMD sumstats ...")
    ld_base = _load_baseline_l2()
    zmap = _load_munge_bmd()
    # M = baseline M
    with open(os.path.join(FLAGSHIP_LDSC, "eur_w_ld_chr.l2.M")) as fh:
        M = float(fh.read().strip().split()[0])
    # N_BMD：由 munge 的 N 列取中位
    import gzip
    Ns = []
    with gzip.open(BMD_SUMSTATS, "rt") as f:
        f.readline()
        for i, line in enumerate(f):
            if i > 200000:
                break
            p = line.rstrip("\n").split("\t")
            try:
                Ns.append(float(p[6]))
            except (ValueError, IndexError):
                pass
    N = float(np.median(Ns))

    tissues, ctx_cols = [], []
    for fp in files:
        name = os.path.basename(fp).replace(".l2.ldscore.gz", "").replace(".ldscore.gz", "").replace(".gp.gz", "")
        col = {}
        for chunk in pd.read_csv(fp, sep=r"\s+", chunksize=500000):
            for s, l in zip(chunk.iloc[:, 0], chunk.iloc[:, -1]):
                try:
                    col[s] = float(l)
                except (ValueError, TypeError):
                    continue
        tissues.append(name)
        ctx_cols.append(col)
    # 共同 SNP 交集
    common = set(zmap) & set(ld_base)
    for c in ctx_cols:
        common &= set(c)
    common = sorted(common)
    print("common SNPs across BMD × baseline × %d SEG tissues: %d" % (len(tissues), len(common)))
    z = np.array([zmap[s] for s in common])
    L2_base = np.array([ld_base[s] for s in common])
    L2_ctx = np.column_stack([[c[s] for s in common] for c in ctx_cols])
    a, b, tau, se_tau, z_tau = seg_regression(z, L2_base, L2_ctx, N, M)

    res = pd.DataFrame({"tissue": tissues, "tau": tau, "se_tau": se_tau, "z": z_tau})
    res = res.sort_values("z", ascending=False).reset_index(drop=True)
    res.to_csv(os.path.join(OUT, "tissue_enrichment_seg.csv"), index=False)
    print("\nTop 12 tissues by SEG enrichment z:")
    print(res.head(12).to_string(index=False))
    # 成骨 vs 非骨/软骨谱系
    tl = [t.lower() for t in tissues]
    osteo = [i for i, t in enumerate(tl) if any(k in t for k in OSTEO_KEYWORDS)]
    chondro = [i for i, t in enumerate(tl) if any(k in t for k in CHONDRO_KEYWORDS)]
    other = [i for i in range(len(tl)) if i not in set(osteo) | set(chondro)]
    def _agg(idx, label):
        if not idx:
            print("%s: (no matching tissue in SEG set)" % label); return
        print("%s (n=%d): mean z=%.2f  mean tau=%.3f  max z=%.2f"
              % (label, len(idx), np.mean(z_tau[idx]), np.mean(tau[idx]), np.max(z_tau[idx])))
    print("\n--- lineage aggregation ---")
    _agg(osteo, "osteoblast/bone lineage")
    _agg(chondro, "cartilage lineage")
    _agg(other, "other lineage")
    print("\nwritten: output/genetics/tissue_enrichment_seg.csv")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "selftest"
    if mode == "selftest":
        run_selftest()
    elif mode == "real":
        run_real()
    else:
        print("usage: python 164_ldsc_seg_native.py [selftest|real]")
