"""
01_build_cohort.py — 队列构建与 1:1 匹配
================================================================================
输入：data/charls.csv, data/ukb.csv（或 --simulate 生成模拟数据）
输出：output/cohort_{arm}.csv（匹配后长表）、output/table1_baseline.csv、
      output/table_smd.csv（匹配质量）

匹配变量（预注册）：年龄 ±1 岁、性别（精确）、波次/招募中心（精确）
对应 R 实现：MatchIt::matchit(method="exact"+nearest, caliper)
================================================================================
"""

import os
import sys
import numpy as np
import pandas as pd
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import DKEYS, PARAMS, OUT_DIR, DATA_DIR


# ----------------------------------------------------------------------------
def load_arm(arm, use_real=False):
    if use_real:
        path = os.path.join(DATA_DIR, f"{arm}.csv")
        if not os.path.exists(path):
            raise FileNotFoundError(f"缺少真实数据文件：{path}")
        df = pd.read_csv(path)
        print(f"[{arm}] 读取真实数据 {df.shape}")
        return df
    from simulate_data import generate
    print(f"[{arm}] ⚠️ 使用【模拟数据】——仅用于管线验证")
    return generate()["charls" if arm == "charls" else "ukb"]


def apply_criteria(df):
    """纳排标准（预注册 §4）。返回 (df, 流失计数)"""
    log = {}
    n0 = len(df)
    df = df[df["age"] >= PARAMS["min_age"]].copy()
    log["age>=45"] = n0 - len(df)
    n1 = len(df)
    df = df.dropna(subset=["age", "sex", "oa"] + DKEYS)
    log["关键字段完整"] = n1 - len(df)
    return df, log


def match_1to1(df, age_window=1, seed=20260908):
    """
    按 性别 × 波次/中心 分层，年龄 ±age_window 岁的 1:1 无放回最近邻匹配。
    返回 matched DataFrame（长表，含 pair_id、is_case）
    """
    rng = np.random.default_rng(seed)
    cases = df[df["oa"] == 1].copy()
    ctrls = df[df["oa"] == 0].copy()

    pool = defaultdict(lambda: defaultdict(list))
    for idx, a, s, st in zip(ctrls.index, ctrls["age"], ctrls["sex"], ctrls["stratum"]):
        pool[(s, st)][int(a)].append(idx)
    for k in pool:
        for a in pool[k]:
            rng.shuffle(pool[k][a])

    used = set()
    pairs = []
    order = rng.permutation(cases.index)
    for ci in order:
        ca, cs, cst = int(cases.at[ci, "age"]), cases.at[ci, "sex"], cases.at[ci, "stratum"]
        picked = None
        for d in sorted(range(-age_window, age_window + 1), key=abs):
            lst = pool[(cs, cst)].get(ca + d, [])
            while lst:
                cand = lst.pop()
                if cand not in used:
                    picked = cand
                    break
            if picked is not None:
                break
        if picked is None:
            continue
        used.add(picked)
        pairs.append((ci, picked))

    m = pd.concat([
        df.loc[[c for c, _ in pairs]].assign(is_case=1),
        df.loc[[p for _, p in pairs]].assign(is_case=0),
    ])
    m["pair_id"] = np.concatenate([np.arange(len(pairs))] * 2)
    return m.sort_values(["pair_id", "is_case"], ascending=[True, False]), len(pairs)


# ----------------------------------------------------------------------------
def smd(x1, x0):
    """标准化均数差（二分类）"""
    p1, p0 = np.mean(x1), np.mean(x0)
    denom = np.sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2)
    return 0.0 if denom == 0 else (p1 - p0) / denom


def smd_cont(x1, x0):
    m1, m0 = np.mean(x1), np.mean(x0)
    s1, s0 = np.std(x1, ddof=1), np.std(x0, ddof=1)
    denom = np.sqrt((s1 ** 2 + s0 ** 2) / 2)
    return 0.0 if denom == 0 else (m1 - m0) / denom


def build_table1(m, arm):
    cases = m[m.is_case == 1]
    ctrls = m[m.is_case == 0]
    w = m.get("weight", pd.Series(np.ones(len(m)), index=m.index))
    wc, w0 = w[cases.index], w[ctrls.index]

    def wmean(x, wt):
        return float(np.average(x, weights=wt))

    def wprop(x, wt):
        return float(np.average(x, weights=wt)) * 100

    rows = []
    rows.append(dict(Variable="N", OA=f"{len(cases):,}", Control=f"{len(ctrls):,}", SMD="—"))
    for var, name, fn in [
        ("age", "Age, mean (SD)", "cont"),
        ("bmi", "BMI, mean (SD)", "cont"),
        ("sex", "Female, %", "bin"),
        ("smoking", "Current smoking, %", "bin"),
        ("alcohol", "Current drinking, %", "bin"),
    ]:
        if var not in m.columns:
            continue
        xc, x0 = cases[var].values, ctrls[var].values
        if fn == "cont":
            rows.append(dict(Variable=name,
                             OA=f"{wmean(xc, wc):.1f} ({np.std(xc, ddof=1):.1f})",
                             Control=f"{wmean(x0, w0):.1f} ({np.std(x0, ddof=1):.1f})",
                             SMD=f"{smd_cont(xc, x0):.3f}"))
        else:
            rows.append(dict(Variable=name,
                             OA=f"{wprop(xc, wc):.1f}",
                             Control=f"{wprop(x0, w0):.1f}",
                             SMD=f"{smd(xc, x0):.3f}"))
    for d in DKEYS:
        if d not in m.columns:
            continue
        xc, x0 = cases[d].values, ctrls[d].values
        rows.append(dict(Variable=d, OA=f"{wprop(xc, wc):.1f}",
                         Control=f"{wprop(x0, w0):.1f}", SMD=f"{smd(xc, x0):.3f}"))
    t = pd.DataFrame(rows)
    t.insert(0, "Cohort", arm.upper())
    return t


# ----------------------------------------------------------------------------
def main(use_real=False):
    os.makedirs(OUT_DIR, exist_ok=True)
    all_t1, all_log = [], []
    for arm in ["charls", "ukb"]:
        df = load_arm(arm, use_real)
        df, log = apply_criteria(df)
        m, npairs = match_1to1(df, PARAMS["match_age_window"])
        print(f"[{arm}] 匹配对 = {npairs:,}；纳排：{log}")
        m.to_csv(os.path.join(OUT_DIR, f"cohort_{arm}.csv"), index=False)
        t1 = build_table1(m, arm)
        all_t1.append(t1)
        all_log.append(dict(cohort=arm, n_pairs=npairs, **log))

    pd.concat(all_t1).to_csv(os.path.join(OUT_DIR, "table1_baseline.csv"), index=False)
    pd.DataFrame(all_log).to_csv(os.path.join(OUT_DIR, "cohort_flow.csv"), index=False)
    print("✅ 队列构建完成 →", OUT_DIR)


if __name__ == "__main__":
    main(use_real="--real" in sys.argv)
