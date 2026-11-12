# -*- coding: utf-8 -*-
"""
code/343 — BMD 效应尺度链独立验证（docs/148 P0-1 铁证 + P0-3 复现入口）
=================================================================================
问题：稿件 L57/L174 写 "no rescaling (is needed / was applied)"，但
  raw  mr_results.csv        BMD_total→fracture IVW = −1.1943
  corr mr_results_corrected  BMD_total→fracture IVW = −0.3371
  比值 = 3.5426 = sqrt(400000 / 31873)
即 deposited 汇总统计是 placeholder N=400,000 尺度，Table 2 是真实 N=31,873 尺度。
本脚本从 deposited sumstats + 交付的 harmonised 逐 SNP 表独立重演这条链，
并确认 harmonised 已是真实尺度（故 Table 2 数值正确，错的是正文披露）。

同时输出交付文件的 SHA256 清单，供 P0-2（唯一权威源声明）与 P0-3（逐 SNP 表纳交付）引用。

复现 Table 2 全部单元格（Egger/Q/I2/WM/LOO）：python3 code/334_data_level_trust_but_verify.py

运行：python3 code/343_scale_chain_verify.py
"""
import gzip
import hashlib
import os
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GEN = os.path.join(ROOT, "code", "output", "genetics")
SS_TOTAL = r"C:/Users/Administrator/WorkBuddy/2026-09-08-20-17-24/OA_FLAGSHIP/07_reference/ldsc/sumstats_total.sumstats.gz"
HARM = os.path.join(GEN, "harmonised", "BMD_total__fracture.csv")
RAW = os.path.join(GEN, "mr_results.csv")
CORR = os.path.join(GEN, "mr_results_corrected.csv")

N_PLACEHOLDER = 400000.0
N_TRUE_BMD = 31873.0
F_BMD = float(np.sqrt(N_PLACEHOLDER / N_TRUE_BMD))


def sha12(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:12]


def load_ss_subset(path, wanted):
    out = {}
    with gzip.open(path, "rt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        iz, im, ip, isnp = hdr.index("Z"), hdr.index("MAF"), hdr.index("P"), hdr.index("SNP")
        for line in f:
            c = line.rstrip("\n").split("\t")
            if c[isnp] in wanted:
                out[c[isnp]] = dict(Z=float(c[iz]), MAF=float(c[im]), P=float(c[ip]))
    return out


def ivw(df, exposure, outcome):
    r = df[(df.exposure == exposure) & (df.outcome == outcome)]
    return float(r.ivw.iloc[0]) if len(r) else None


def main():
    print("=" * 78)
    print("BMD 效应尺度链独立验证 — docs/148 P0-1")
    print("=" * 78)
    print(f"placeholder N={N_PLACEHOLDER:,.0f}; verified N={N_TRUE_BMD:,}; "
          f"sqrt(N_ph/N_true) = {F_BMD:.4f}")

    harm = pd.read_csv(HARM)
    ss = load_ss_subset(SS_TOTAL, set(harm.SNP))
    rec = []
    for _, r in harm.iterrows():
        d = ss.get(r.SNP)
        if not d:
            continue
        f = d["MAF"] * (1 - d["MAF"])
        se_ph = 1 / np.sqrt(2 * N_PLACEHOLDER * f)
        se_tr = 1 / np.sqrt(2 * N_TRUE_BMD * f)
        rec.append(dict(SNP=r.SNP, Zss=abs(d["Z"]), z_harm=abs(r.beta_e / r.se_e),
                        beta_e_harm=abs(r.beta_e), beta_ph=abs(d["Z"]) * se_ph,
                        beta_tr=abs(d["Z"]) * se_tr, se_e_harm=r.se_e,
                        se_ph=se_ph, se_tr=se_tr))
    df = pd.DataFrame(rec)
    print(f"\n[1] 匹配 SNP: {len(df)}/{len(harm)}")
    z_ratio = float((df.z_harm / df.Zss.replace(0, np.nan)).median())
    r_tr = float((df.beta_e_harm / df.beta_tr.replace(0, np.nan)).median())
    r_ph = float((df.beta_e_harm / df.beta_ph.replace(0, np.nan)).median())
    se_ratio = float((df.se_e_harm / df.se_tr).median())
    print(f"    harmonised z / sumstats|Z| 中位 = {z_ratio:.4f}  (≈1 → 同一信号，z 尺度无关)")
    print(f"    harmonised |beta_e| / (|Z|·SE_true[N=31873]) 中位 = {r_tr:.4f}  (≈1 → 已是真实尺度)")
    print(f"    harmonised |beta_e| / (|Z|·SE_ph[N=400000])   中位 = {r_ph:.4f}  (≈{F_BMD:.3f} → placeholder 尺度)")
    print(f"    harmonised SE / SE_true 中位 = {se_ratio:.4f}")

    print("\n[2] raw vs corrected 主结果（IVW β）")
    raw = pd.read_csv(RAW)
    corr = pd.read_csv(CORR)
    ratios = {}
    for e, o in [("BMD_total", "fracture"), ("BMD_total", "OA_knee"),
                 ("BMD_femur_neck", "fracture")]:
        a, b = ivw(raw, e, o), ivw(corr, e, o)
        if a is not None and b is not None:
            ratios[e] = a / b
            print(f"    {e}->{o}: raw={a:+.4f}  corrected={b:+.4f}  ratio={a/b:.4f}")

    print("\n[3] harmonised 逐 SNP 表 → Table 2 IVW 重算")
    # Table 2 的标准 TwoSampleMR IVW：ratio SE = se_o/|beta_e|（只含结局 SE），
    # 权重 w=1/(se_o/beta_e)^2 = 1/se_ratio^2；正是 code/334 的 18/18 复现口径。
    beta_o = harm.beta_o.values
    be = harm.beta_e.values
    so = harm.se_o.values
    se_be = harm.se_e.values
    theta = beta_o / be
    w_std = (be / so) ** 2                      # = 1/(so/be)^2
    ivw_std = float((w_std * theta).sum() / w_std.sum())
    se_std = float(np.sqrt(1 / w_std.sum()))
    # 注记：完整 delta 方差（额外含 beta_e 分母误差）是 code/340 测量误差模型用的
    # V_j 口径，给出略不同的加权估计；它不是 Table 2 报告量，仅作方法学对照。
    var_full = (so / be) ** 2 + (beta_o ** 2 / be ** 4) * se_be ** 2
    w_full = 1 / var_full
    ivw_full = float((w_full * theta).sum() / w_full.sum())
    b2 = ivw(corr, "BMD_total", "fracture")
    print(f"    标准 IVW (w=1/(se_o/beta_e)^2, Table 2 口径) = {ivw_std:+.6f}  "
          f"SE={se_std:.6f}  vs corrected={b2:+.6f}")
    print(f"    完整 delta IVW (w=1/Var_full, 含 beta_e 误差; code/340 口径) = {ivw_full:+.6f}")
    print("    —— 两者权重定义不同，均非错误；Table 2 采标准口径，R1 模型用完整口径。")

    print("\n[4] 交付文件 SHA256（前 12 位）— 唯一权威源 / 逐 SNP 复现数据")
    manifest = []
    hdir = os.path.join(GEN, "harmonised")
    for fn in sorted(os.listdir(hdir)):
        if fn.endswith(".csv"):
            p = os.path.join(hdir, fn)
            manifest.append((f"harmonised/{fn}", sha12(p)))
    for fn in ("mr_results_corrected.csv", "mr_results.csv",
               "_backup_pre_G6/mr_results_corrected.csv"):
        p = os.path.join(GEN, fn)
        if os.path.exists(p):
            manifest.append((fn, sha12(p)))
    for name, h in manifest:
        print(f"    {h}  {name}")

    # ---- 断言：尺度链自洽 ----
    assert abs(r_tr - 1.0) < 0.02, f"harmonised 非真实 N 尺度: {r_tr}"
    assert abs(r_ph - F_BMD) < 0.05, f"placeholder 因子不符: {r_ph} vs {F_BMD}"
    assert abs(se_ratio - 1.0) < 0.02, f"SE 尺度不符: {se_ratio}"
    a, b = ivw(raw, "BMD_total", "fracture"), ivw(corr, "BMD_total", "fracture")
    assert abs(a / b - F_BMD) < 0.01, f"raw/corrected 比值 != 尺度因子: {a/b}"
    assert abs(ivw_std - b) < 1e-6, "harmonised 标准 IVW 重算与 corrected 不符"
    assert abs(se_std - 0.022406) < 1e-5, "标准 IVW SE 与 Table 2 不符"

    print("\n" + "=" * 78)
    print("[结论] deposited sumstats 为 placeholder N=400,000 尺度；交付 harmonised 已是")
    print(f"        真实 N=31,873 尺度；raw mr_results.csv 比 Table 2 高 {F_BMD:.3f}×。")
    print("        Table 2 数值正确，但正文 L57/L174 的 'no rescaling' 披露错误，必须写明")
    print("        BMD β/SE ×sqrt(400000/31873)=3.543（逐 SNP z 不变）这一步。")
    print("        Table 2 全单元格复现 → python3 code/334_data_level_trust_but_verify.py")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
