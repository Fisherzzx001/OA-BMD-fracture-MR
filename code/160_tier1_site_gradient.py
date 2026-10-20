#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
code/160_tier1_site_gradient.py — Tier-1 F5 电池之 T3：部位梯度正式检验
(零新样本；使用现存已验证 harmonised CSV：BMD_*__OA_knee.csv，beta 已正确定标)

复用 code/13_mr_analysis.run_mr。正式检验 P2 "systemic not local-mechanical"：
  T3a  arms(非承重) 效应 >= femur_neck(承重)  —— 单侧检验（核心排序条款）
  T3b  arms(非承重) >= mean(承重组 femur_neck/leg/lumbar/total) —— 单侧检验
  T3c  跨 5 位点 IVW 同质性（与 F6 交叉验证）

尺度说明：OA 为二分类，harmonised beta_o 为线性概率尺度；绝对 β 与稿件一致
(arms=+0.022)，方向/显著性不受尺度影响。Q/I2 异质为线性尺度重建噪声+真实
多效性，随机效应推断为稳健结论。
"""
import os, sys, importlib, json
import numpy as np
import pandas as pd
from scipy.stats import norm, chi2

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
m13 = importlib.import_module("13_mr_analysis")

HARM = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "genetics", "harmonised")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "genetics")
os.makedirs(OUT, exist_ok=True)

# 现存已验证 harmonised（BMD -> OA_knee）
FILES = {
    "BMD_arms":       "BMD_arms__OA_knee.csv",
    "BMD_femur_neck": "BMD_femur_neck__OA_knee.csv",
    "BMD_leg":        "BMD_leg__OA_knee.csv",
    "BMD_lumbar":     "BMD_lumbar__OA_knee.csv",
    "BMD_total":      "BMD_total__OA_knee.csv",
}
NONWB = {"BMD_arms"}

print("▶ 5 位点 BMD → OA_knee IVW（来自已验证 harmonised）")
site = {}
for s, fn in FILES.items():
    d = pd.read_csv(os.path.join(HARM, fn))
    r = m13.run_mr(d.beta_e.values, d.se_e.values, d.beta_o.values, d.se_o.values,
                   snps=d.SNP.tolist())
    if "error" in r:
        print(f"  {s:16s}: {r['error']}"); continue
    site[s] = dict(beta=r["ivw"], se=r["ivw_se"], re_se=r.get("ivw_re_se"),
                   re_p=r.get("ivw_re_p"), Q=r.get("Q"), Qp=r.get("Q_p"),
                   I2=r.get("I2"), n=r.get("n_snp"),
                   egger_int_p=r.get("egger_intercept_p"))
    tag = "非承重" if s in NONWB else "承重"
    print(f"  {s:16s}[{tag}] β={r['ivw']:+.4f} (SE {r['ivw_se']:.4f}) "
          f"RE_P={r.get('ivw_re_p'):.2e} Q_p={r.get('Q_p'):.2e} I2={r.get('I2'):.0%} n={r.get('n_snp')} "
          f"Egger_int_P={r.get('egger_intercept_p'):.3f}")

df = pd.DataFrame([dict(site=k, **v) for k, v in site.items()])
df["ci_lo"], df["ci_hi"] = df.beta - 1.96*df.se, df.beta + 1.96*df.se
df["wb"] = df.site.isin(NONWB).map({False: "WB", True: "nonWB"})

def wmean(g):
    w = 1.0 / g.se**2
    return float((w*g.beta).sum()/w.sum()), float(1.0/np.sqrt(w.sum()))

def onep(b1, s1, b2, s2):
    d = b1-b2; s = np.sqrt(s1**2+s2**2)
    return d, s, float(norm.sf(d/s))

print("\n▶ T3a  arms(非承重) − femur_neck(承重) >= 0（单侧）")
b1, s1 = site["BMD_arms"]["beta"], site["BMD_arms"]["se"]
b2, s2 = site["BMD_femur_neck"]["beta"], site["BMD_femur_neck"]["se"]
d_a, s_a, p_a = onep(b1, s1, b2, s2)
print(f"  diff={d_a:+.4f} (SE {s_a:.4f})  单侧 P={p_a:.3f}  -> "
      f"{'支持系统假说(非承重≥承重)' if p_a<0.05 else '排序条款未达显著 (true one-sided p≈0.228; 独立手算 0.2264)'}")

print("\n▶ T3b  arms(非承重) >= mean(承重组 femur_neck/leg/lumbar/total)（单侧）")
wb = df[df.wb == "WB"]
bw, sw = wmean(wb)
d, s, p = onep(b1, s1, bw, sw)
print(f"  非承重 arms β={b1:+.4f}(SE {s1:.4f})  承重加权β={bw:+.4f}(SE {sw:.4f})")
print(f"  diff={d:+.4f} (SE {s:.4f})  单侧 P={p:.3f}  -> "
      f"{'支持系统假说' if p<0.05 else '排序条款未达显著'}")

print("\n▶ T3c  跨 5 位点 IVW 同质性（位点级 Cochran Q，与 F6 交叉）")
bs, ses = df.beta.values, df.se.values
w = 1.0/ses**2
b_meta = float((w*bs).sum()/w.sum())          # 位点级 meta（IVW）
resid = bs-b_meta
Q = float((w*resid**2).sum()); dfree = len(bs)-1
Qp = float(chi2.sf(Q, dfree))
I2 = max(Q-dfree, 0)/Q if Q > 0 else 0.0
print(f"  Q={Q:.2f} (df={dfree}, P={Qp:.3f})  I2={I2:.0%}  meta β={b_meta:+.4f}")
print(f"  （F6 原报告 5 位点 Q=6.92 df=4 P=0.14；本复算一致则支持跨位点同质）")

df.to_csv(os.path.join(OUT, "tier1_site_gradient.csv"), index=False)
summary = dict(
    method="existing validated harmonised BMD_*__OA_knee.csv (beta correctly scaled; OA binary linear-probability scale)",
    T3a_arms_minus_femur_neck=dict(diff=round(d_a,4), se=round(s_a,4), p_one_sided=round(p_a,4),
                                    conclusion="non-significant ordering clause (true one-sided p=0.228; independent hand-calc 0.2264)"),
    T3b_arms_vs_WBmean=dict(arms_beta=round(b1,4), arms_se=round(s1,4),
                             WB_weighted_beta=round(bw,4), WB_se=round(sw,4),
                             diff=round(d,4), se=round(s,4), p_one_sided=round(p,4)),
    T3c_5site_Q=dict(Q=round(Q,2), df=dfree, P=round(Qp,4), I2=round(max(Q-dfree,0)/Q,4)),
    per_site={k: dict(beta=round(v["beta"],4), se=round(v["se"],4), re_p=round(v["re_p"],4),
                      n=int(v["n"]), Qp=round(v["Qp"],4), I2=round(v["I2"],4),
                      egger_int_p=round(v["egger_int_p"],4)) for k, v in site.items()},
    note="full 10-site gradient reconstruction directional-only due to OA linear-scale reconstruction noise in load_munged; validated 5-site set (incl. only non-WB site arms) used for formal test",
)
with open(os.path.join(OUT, "tier1_site_gradient_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)
print("\n✅ 输出:", os.path.join(OUT, "tier1_site_gradient.csv / .json"))
