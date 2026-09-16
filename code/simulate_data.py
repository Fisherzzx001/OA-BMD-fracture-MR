"""
simulate_data.py — 生成结构与 CHARLS / UKB 一致的**模拟**数据
================================================================================
用途：在拿到真实数据之前，端到端验证分析管线（匹配→描述→网络→推断→出图）能否跑通，
并确认管线确实能检出"中国代谢驱动 / 英国骨-关节驱动"的预设结构。

⚠️ 严禁将本模块产生的任何数字写入论文。运行产物一律标注 SYNTHETIC。

关键设计：用二分法标定每个疾病的截距，使**边际患病率**落在真实的流行病学区间。
这一点很重要——O/E 有数学上限 1/max(p_i, p_j)，若患病率不真实（如高血压 67%），
超额共现会被结构性压低，网络退化，检验全部失效。
================================================================================
"""

import numpy as np
import pandas as pd

from config import (DKEYS, DISEASES, DMODULE, SIM, SEED)

# 目标边际患病率（45+ 人群，CHARLS/UKB 量级）
TARGET_PREV = {
    "hypertension": 0.34, "diabetes": 0.11, "dyslipidaemia": 0.22,
    "heart_disease": 0.10, "stroke": 0.045, "asthma": 0.035, "copd": 0.060,
    "liver_disease": 0.030, "kidney_disease": 0.035, "digestive": 0.10,
    "cancer": 0.025, "depression": 0.11, "dementia": 0.020,
    "back_pain": 0.30, "neck_pain": 0.20, "multisite_pain": 0.15,
}


def _expit(x):
    return 1.0 / (1.0 + np.exp(-x))


def _calibrate(u, target, lo=-14.0, hi=10.0, iters=40):
    """二分法求截距 c，使 mean(expit(u + c)) ≈ target"""
    for _ in range(iters):
        mid = (lo + hi) / 2
        if _expit(u + mid).mean() < target:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def _make_arm(arm, rng, n_mc=20000):
    cfg = SIM[arm]
    n = cfg["n"]

    age = np.clip(rng.normal(cfg["age_mean"], cfg["age_sd"], n), 45, 95).round()
    sex = rng.binomial(1, 0.52, n)                      # 1 = 女
    bmi = np.clip(rng.normal(cfg["bmi_mean"], cfg["bmi_sd"], n), 15, 55)
    stratum = rng.integers(0, 10, n)                   # 波次 / 招募中心
    weight = rng.uniform(0.5, 2.5, n)
    edu = rng.choice([0, 1, 2], n, p=[0.45, 0.40, 0.15])
    smoke = rng.choice([0, 1, 2], n, p=[0.55, 0.25, 0.20])
    alcohol = rng.choice([0, 1, 2], n, p=[0.45, 0.15, 0.40])

    z_age = (age - 60) / 10
    z_bmi = (bmi - 25) / 5

    p_oa = _expit(cfg["oa_prob_base"] + cfg["oa_b_age"] * z_age
                  + cfg["oa_b_bmi"] * z_bmi + cfg["oa_b_sex"] * sex
                  + 0.06 * z_age * z_bmi)
    oa = rng.binomial(1, p_oa)

    mods = ["metabolic", "cardiovascular", "respiratory",
            "organ_other", "mental", "musculoskeletal"]
    F = {m: rng.normal(0, 1, n) for m in mods}
    # 全局"多病倾向"因子：模拟就医密度/健康行为导致的跨模块共现，
    # 这是真实多病网络的主体部分，缺了它网络会稀疏到无法分析。
    F_mm = rng.normal(0, 1, n)

    # --- 用于标定的独立蒙特卡洛样本（与正式样本同分布）---
    z_age_mc = rng.normal(0, 1, n_mc)
    z_bmi_mc = rng.normal(0, 1, n_mc)
    sex_mc = rng.binomial(1, 0.52, n_mc)
    F_mc = {m: rng.normal(0, 1, n_mc) for m in F}
    F_mm_mc = rng.normal(0, 1, n_mc)

    cols = {}
    for key, label, mod in DISEASES:
        lam = cfg["factor_load"][mod]
        target = TARGET_PREV[key] * cfg["prev_scale"][mod]

        u_mc = (0.35 * z_age_mc + 0.12 * z_bmi_mc + 0.10 * sex_mc
                + lam * F_mc[mod] + cfg["global_load"] * F_mm_mc
                + 0.20 * rng.normal(0, 1, n_mc))
        c = _calibrate(u_mc, min(max(target, 0.005), 0.85))

        u = (c + 0.35 * z_age + 0.12 * z_bmi + 0.10 * sex
             + lam * F[mod] + cfg["global_load"] * F_mm * (1 + 0.30 * oa)
             + 0.20 * rng.normal(0, 1, n)
             + cfg["oa_extra_by_module"][mod] * oa)
        cols[key] = rng.binomial(1, _expit(u))

    df = pd.DataFrame(cols)
    df.insert(0, "pid", [f"{arm.upper()[:3]}{i:07d}" for i in range(n)])
    df["age"] = age
    df["sex"] = sex
    df["bmi"] = np.round(bmi, 1)
    df["stratum"] = stratum
    df["weight"] = weight
    df["education"] = edu
    df["smoking"] = smoke
    df["alcohol"] = alcohol
    df["oa"] = oa
    return df


def generate(seed=SEED):
    rng = np.random.default_rng(seed)
    return {"ukb": _make_arm("ukb", rng), "charls": _make_arm("charls", rng)}


if __name__ == "__main__":
    d = generate()
    for k, v in d.items():
        print(f"{k}: n={len(v)}, OA n={int(v['oa'].sum())} "
              f"({100 * v['oa'].mean():.1f}%)")
        print("   患病率:", {c: round(float(v[c].mean()), 3) for c in DKEYS[:6]})
