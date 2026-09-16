#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
17_strict_clump_mrpresso.py
=========================================================================
投稿前整改 #3：
  (A) 严格 PLINK clumping（1000G EUR, r^2<0.001, 10 Mb）命令生成 + 输入文件准备
  (B) MR-PRESSO 全局（distortion）检验 —— 纯 Python 实现，无需 PLINK 即可运行

环境说明：
  · 本机未安装 PLINK，也未挂载 1000G EUR bfile，故 (A) 仅「生成可执行命令 +
    写出 clump 输入文件」，待用户在装有 PLINK 的环境运行；
  · (B) 直接基于现有 1-Mb clumped 仪器集运行，给出 MR-PRESSO 全局检验
    （逐 SNP 离群检测 + distortion test），并与 13_ 的 radial |z|>3 对比；
  · 用严格 clump 后得到的 <exposure>.clumped 覆盖仪器集后，可在本脚本开关控制下
    复跑 MR + MR-PRESSO（见 main 末尾的 STRICT_CLUMP_DIR 选项）。
=========================================================================
"""
import os
import sys
import importlib
import numpy as np
import pandas as pd
from scipy.stats import chi2

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
m13 = importlib.import_module("13_mr_analysis")

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "output", "genetics")

# ---- 可配置：1000G EUR bfile 路径（用户按本地环境修改）------------------
BFILE = "/path/to/1000G_phase3_EUR/1000G_phase3_EUR"

P1 = 1e-5       # --clump-p1 / --clump-p2（与 select_instruments p<5e-8 衔接时可取 5e-8）
R2 = 0.001      # --clump-r2  严格 LD 阈值
KB = 10000      # --clump-kb  10 Mb 窗口

PAIRS = m13.PAIRS

# 严格 clump 结果目录：若此处存在 <exposure>.clumped 则用之覆盖仪器集
STRICT_CLUMP_DIR = None


# ----------------------------------------------------------------- MR-PRESSO
def _ivw(be, bo, so):
    w = 1.0 / (so ** 2)
    return float(np.sum(w * be * bo) / np.sum(w * be ** 2))


def _ivw_se(be, so):
    return float(1.0 / np.sqrt(np.sum(be ** 2 / so ** 2)))


def mr_presso_global(be, se_e, bo, so, bonf=True):
    """
    MR-PRESSO 全局（distortion）检验（Verbanck et al. 2018, Nat Genet）。
      · Q_full     = Cochran Q(IVW)
      · 逐 SNP 离群检验：Q_i = Q_full - Q_{(-i)} ~ χ²(1)，Bonferroni 校正
      · 剔除离群后重估 IVW -> Q_outlier_free
      · χ²_distortion = Q_full - Q_outlier_free, df = m(离群数) -> p = chi2.sf
    返回 dict。
    """
    be = np.asarray(be, float)
    se_e = np.asarray(se_e, float)
    bo = np.asarray(bo, float)
    so = np.asarray(so, float)
    K = len(be)
    b_full = _ivw(be, bo, so)
    se_ivw = _ivw_se(be, so)
    r = bo - b_full * be
    sig = np.sqrt(so ** 2 + (b_full * se_e) ** 2)          # σ_i² = se_Y² + β_IVW²·se_X²
    Q_full = float(np.sum((r / sig) ** 2))
    # 逐 SNP 离群检验
    outliers = []
    thr = (0.05 / K) if bonf else 0.05
    for i in range(K):
        m = np.ones(K, bool)
        m[i] = False
        b_wo = _ivw(be[m], bo[m], so[m])
        r_wo = bo[m] - b_wo * be[m]
        sig_wo = np.sqrt(so[m] ** 2 + (b_wo * se_e[m]) ** 2)
        Q_wo = float(np.sum((r_wo / sig_wo) ** 2))
        Q_i = Q_full - Q_wo
        p_i = float(chi2.sf(Q_i, 1))
        if p_i < thr:
            outliers.append(i)
    m_out = len(outliers)
    keep = np.ones(K, bool)
    keep[outliers] = False
    if int(keep.sum()) >= 2:
        b_out = _ivw(be[keep], bo[keep], so[keep])
        r_out = bo[keep] - b_out * be[keep]
        sig_out = np.sqrt(so[keep] ** 2 + (b_out * se_e[keep]) ** 2)
        Q_out = float(np.sum((r_out / sig_out) ** 2))
    else:
        b_out, Q_out = b_full, Q_full
    chi2_dist = Q_full - Q_out
    df_dist = max(m_out, 1)
    p_dist = float(chi2.sf(chi2_dist, df_dist)) if m_out > 0 else 1.0
    return dict(K=K, Q_full=Q_full, n_outlier=m_out,
                b_ivw_full=b_full, b_ivw_out=b_out,
                Q_out=Q_out, chi2_distortion=chi2_dist, p_distortion=p_dist)


# ----------------------------------------------------------------- PLINK 部分
def prepare_clump_inputs(dfs):
    print("▶ 准备 PLINK clump 输入文件（SNP + P，来自各暴露 GWAS 全基因组成因）...")
    n_written = 0
    for nm, d in dfs.items():
        if "SNP" in d.columns:                       # munge 表型有 SNP 列
            sid = d["SNP"].astype(str)
        else:                                        # fracture 等原始格式无 SNP -> 用 CHR:BP
            sid = d["CHR"].astype(str) + ":" + d["BP"].astype(str)
        sub = pd.DataFrame({"SNP": sid, "P": d["P"]}).dropna()
        sub.to_csv(os.path.join(OUT, f"{nm}_clump_input.txt"), sep="\t", index=False)
        n_written += 1
    print(f"  已写出 {n_written} 个 clump 输入文件 -> {OUT}")


def emit_plink_commands():
    print("=" * 72)
    print("A. 严格 PLINK clumping 命令（需在装有 PLINK + 1000G EUR bfile 的环境运行）")
    print("=" * 72)
    print(f"# bfile 指向 1000G Phase3 EUR：{BFILE}")
    print("# 依赖：plink --bfile 计算 r^2；--clump 输入为上一步写出的 *_clump_input.txt")
    for exp, out, note in PAIRS:
        inp = os.path.join(OUT, f"{exp}_clump_input.txt")
        clumped = os.path.join(OUT, f"{exp}_clumped")
        cmd = (f"plink --bfile {BFILE} --clump {inp} "
               f"--clump-p1 {P1:g} --clump-p2 {P1:g} "
               f"--clump-r2 {R2} --clump-kb {KB} --out {clumped}")
        print(f"\n# [{exp} -> {out}] {note}")
        print(cmd)
        print(f"#   保留 SNP 列表：{clumped}.clumped")


def load_strict_instruments(exp, dfs):
    """若 STRICT_CLUMP_DIR 存在 <exp>.clumped，则据其 SNP 列表返回严格仪器集。"""
    if not STRICT_CLUMP_DIR:
        return None
    fp = os.path.join(STRICT_CLUMP_DIR, f"{exp}_clumped.clumped")
    if not os.path.exists(fp):
        return None
    clumped = pd.read_csv(fp, sep=r"\s+")
    snps = set(clumped["SNP"].astype(str))
    d = dfs[exp]
    return d[d["SNP"].astype(str).isin(snps)]


# ----------------------------------------------------------------- 主流程
def main():
    os.makedirs(OUT, exist_ok=True)
    dfs = m13.load_all()
    prepare_clump_inputs(dfs)
    emit_plink_commands()

    print("\n" + "=" * 72)
    print("B. MR-PRESSO 全局（distortion）检验")
    print("=" * 72)
    if STRICT_CLUMP_DIR and os.path.isdir(STRICT_CLUMP_DIR):
        print(f"★ 使用严格 clump 仪器集：{STRICT_CLUMP_DIR}")
    else:
        print("⚠  PLINK 未运行 / bfile 未挂载：以下基于近似 1-Mb clump 仪器集。")
        print("   严格 clump 后请将 <exposure>.clumped 放入 STRICT_CLUMP_DIR 并复跑。\n")

    inst = {nm: m13.select_instruments(d) for nm, d in dfs.items()}
    rows = []
    for exp, out, note in PAIRS:
        strict = load_strict_instruments(exp, dfs)
        if strict is not None:
            iv = m13.select_instruments(strict)      # 已在显著集内，仅去重窗口
        else:
            iv = inst[exp]
        h = m13.harmonise(iv, dfs[out])
        if len(h) < 3:
            print(f"【{exp} -> {out}】 harmonise 后 IV 不足 ({len(h)})，跳过")
            continue
        mp = mr_presso_global(h.beta_e.values, h.se_e.values,
                              h.beta_o.values, h.se_o.values)
        print(f"【{exp} -> {out}】 {note}")
        print(f"   K={mp['K']}  Q_full={mp['Q_full']:.1f}  离群SNP={mp['n_outlier']}  "
              f"β_IVW(full)={mp['b_ivw_full']:+.4f}  β_IVW(outlier-free)={mp['b_ivw_out']:+.4f}")
        print(f"   χ²_distortion={mp['chi2_distortion']:.1f} (df={max(mp['n_outlier'], 1)})  "
              f"MR-PRESSO 全局P={mp['p_distortion']:.3g}")
        rows.append(dict(exposure=exp, outcome=out, note=note, **mp))

    pd.DataFrame(rows).to_csv(os.path.join(OUT, "mr_presso_results.csv"), index=False)
    print(f"\n✅ MR-PRESSO 结果 -> {OUT}/mr_presso_results.csv")


if __name__ == "__main__":
    main()
