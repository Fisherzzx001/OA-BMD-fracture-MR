"""
03_network_analysis.py — 描述层 → 疾病对 → 共病网络 → 统计推断
================================================================================
输出（均在 output/）：
  table2_prevalence.csv      单病加权患病率 + McNemar（匹配数据的正确检验）
  table2_pairs.csv           Top-10 二元 / 三元组合与 O/E
  table3_network_metrics.csv 五类网络指标 + 置换检验 P 值
  table4_hubs.csv            度中心性 Top-5 枢纽病（中英对比）
  table5_cross_country.csv   跨库 Δ 平均度 / Δ 模块内平均度 + bootstrap 95%CI

对应 R：survey（加权） / igraph（网络） / NetworkComparisonTest（推断）
================================================================================
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import DKEYS, DLABEL, DMODULE, MODULES, PARAMS, OUT_DIR
from lib_network import (edge_stats, select_edges, network_metrics, module_degree,
                         degree_centrality, permutation_test_pairs,
                         bootstrap_ci_metric, bootstrap_cross_country)

SYNTHETIC = "--real" not in sys.argv


# ---------------------------------------------------------------- 描述层
def weighted_prop(x, w):
    return float(np.average(x, weights=w)) if len(x) else np.nan


def mcnemar_paired(m, dkey):
    """1:1 匹配数据的正确检验：McNemar。返回 (OA%, 对照%, OR, P, b, c)"""
    sub = m[["pair_id", "is_case", dkey]].pivot(index="pair_id",
                                                columns="is_case", values=dkey)
    if sub.shape[1] < 2:
        return np.nan, np.nan, np.nan, np.nan, 0, 0
    a, b0 = sub[1].values.astype(int), sub[0].values.astype(int)
    b = int(((a == 1) & (b0 == 0)).sum())
    c = int(((a == 0) & (b0 == 1)).sum())
    orr = b / c if c > 0 else np.inf
    p = stats.binomtest(b, b + c, 0.5).pvalue if (b + c) > 0 else np.nan
    return 100 * a.mean(), 100 * b0.mean(), orr, p, b, c


def table2_prevalence(cohorts):
    rows = []
    for arm, m in cohorts.items():
        w = m.get("weight", pd.Series(np.ones(len(m)))).values
        case_mask = (m.is_case == 1).values
        for d in DKEYS:
            x = m[d].values
            p_oa = 100 * weighted_prop(x[case_mask], w[case_mask])
            p_ct = 100 * weighted_prop(x[~case_mask], w[~case_mask])
            _, _, orr, p, b, c = mcnemar_paired(m, d)
            rows.append(dict(cohort=arm.upper(), disease=DLABEL[d],
                             prev_oa=round(p_oa, 2), prev_control=round(p_ct, 2),
                             matched_OR=(round(orr, 3) if np.isfinite(orr) else np.inf),
                             p_mcnemar=(f"{p:.3g}" if np.isfinite(p) else "NA"),
                             discordant_b=b, discordant_c=c))
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(OUT_DIR, "table2_prevalence.csv"), index=False)
    return t


def table2_pairs(cohorts, topn=10):
    rows = []
    for arm, m in cohorts.items():
        X = m.loc[m.is_case == 1, DKEYS].values.astype(float)
        n, p_ = X.shape
        st = edge_stats(X)
        iu = np.triu_indices(p_, k=1)
        pairs = []
        for t in range(len(iu[0])):
            i, j = iu[0][t], iu[1][t]
            cnt = st["co"][i, j]
            if cnt < 5:
                continue
            pairs.append((DLABEL[DKEYS[i]] + " + " + DLABEL[DKEYS[j]],
                          100 * cnt / n, st["oe"][i, j], np.exp(st["logor"][i, j])))
        pairs.sort(key=lambda z: -z[1])
        for r in pairs[:topn]:
            rows.append(dict(cohort=arm.upper(), combination=r[0],
                             prevalence_oa=round(r[1], 2), OE=round(r[2], 2),
                             OR=round(r[3], 2), type="dyad"))

        # 三元组合
        tri = []
        for i in range(p_):
            for j in range(i + 1, p_):
                for k in range(j + 1, p_):
                    cnt = int(np.sum((X[:, i] > 0) & (X[:, j] > 0) & (X[:, k] > 0)))
                    if cnt < 5:
                        continue
                    pi, pj, pk = X[:, i].mean(), X[:, j].mean(), X[:, k].mean()
                    oe = (cnt / n) / (pi * pj * pk)
                    tri.append((" + ".join([DLABEL[DKEYS[i]], DLABEL[DKEYS[j]],
                                            DLABEL[DKEYS[k]]]), 100 * cnt / n, oe))
        tri.sort(key=lambda z: -z[1])
        for r in tri[:topn]:
            rows.append(dict(cohort=arm.upper(), combination=r[0],
                             prevalence_oa=round(r[1], 2), OE=round(r[2], 2),
                             OR=np.nan, type="triad"))
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(OUT_DIR, "table2_pairs.csv"), index=False)
    return t


# ---------------------------------------------------------------- 网络层
def build_arm(cohorts, arm, top_k=None):
    m = cohorts[arm]
    X_oa = m.loc[m.is_case == 1, DKEYS].values.astype(float)
    X_ct = m.loc[m.is_case == 0, DKEYS].values.astype(float)
    prev_oa = X_oa.mean(axis=0)
    st_oa = edge_stats(X_oa)
    st_ct = edge_stats(X_ct)
    adj_oa, q_oa, keep_oa = select_edges(st_oa, len(DKEYS),
                                         min_prev=PARAMS["min_prevalence"],
                                         prev=prev_oa, top_k=top_k)
    adj_ct, q_ct, keep_ct = select_edges(st_ct, len(DKEYS),
                                         min_prev=PARAMS["min_prevalence"],
                                         prev=X_ct.mean(axis=0), top_k=top_k)
    return dict(X_oa=X_oa, X_ct=X_ct, adj_oa=adj_oa, adj_ct=adj_ct,
                prev_oa=prev_oa, keep_oa=keep_oa, keep_ct=keep_ct,
                st_oa=st_oa, st_ct=st_ct)


def equal_density_k(built):
    """
    共同 K = 两臂【OA 网络】显著边数的较小值。

    ⚠️ 只按 OA 网络定 K，不把对照网络纳入：
      对照组网络只用于库内比较（H1），而 H1 不做稀疏化；
      若把对照网络也算进来，K 会被样本量较小的 CHARLS 对照组拖到很低，
      白白牺牲跨国比较的分辨率。
    """
    ks = [int(np.triu(built[arm]["adj_oa"], 1).sum()) for arm in ["charls", "ukb"]]
    return max(1, int(np.min(ks)))


def table3_metrics(built, k):
    rows = []
    for arm in ["charls", "ukb"]:
        b = built[arm]
        for tag, grp, kmask in (("adj_oa", "OA", "keep_oa"),
                                ("adj_ct", "Matched control", "keep_ct")):
            mt = network_metrics(b[tag], b[kmask])
            rows.append(dict(cohort=arm.upper(), group=grp,
                             n_edges=int(mt["n_edges"]),
                             average_degree=round(mt["average_degree"], 3),
                             density=round(mt["density"], 4),
                             modularity=round(mt["modularity"], 3)
                             if np.isfinite(mt["modularity"]) else np.nan,
                             centralization=round(mt["centralization"], 3)
                             if np.isfinite(mt["centralization"]) else np.nan,
                             clustering=round(mt["clustering"], 3)
                             if np.isfinite(mt["clustering"]) else np.nan))
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(OUT_DIR, "table3_network_metrics.csv"), index=False)
    return t


def table4_hubs(built):
    rows = []
    for arm in ["charls", "ukb"]:
        dc = degree_centrality(built[arm]["adj_oa"])
        keep = built[arm]["keep_oa"]
        items = [(DLABEL[DKEYS[i]], dc.get(i, 0.0), DMODULE[DKEYS[i]])
                 for i in range(len(DKEYS)) if keep[i]]
        items.sort(key=lambda z: -z[1])
        for rank, (name, c, mod) in enumerate(items[:5], 1):
            rows.append(dict(cohort=arm.upper(), rank=rank, hub_disease=name,
                             degree_centrality=round(c, 3), module=mod))
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(OUT_DIR, "table4_hubs.csv"), index=False)
    return t


# ---------------------------------------------------------------- 推断层
def inference(built, k, built_k):
    """
    ⚠️ 关键：H1 与 H2/H4 使用【不同的】网络构建方式。
      H1（同库内 OA vs 对照）→ **不做**等密度稀疏化。
        原因：等密度会强制两网边数相等，平均度差异恒为 0，检验失去意义。
        同库内两组的节点可得性完全一致，直接比较即为有效比较。
      H2/H4（跨库中国 vs 英国）→ **必须**等密度稀疏化到共同 K。
        原因：两库的患病率水平、样本量、诊断口径不同，不稀疏化则比较的是
              "谁的网络更密"这一由规模决定的伪差异。
    """
    out = {}
    # H1：库内 OA vs 对照（配对标签置换，不稀疏化）
    h1 = []
    for arm in ["charls", "ukb"]:
        b = built[arm]
        delta, p, null = permutation_test_pairs(b["X_oa"], b["X_ct"],
                                                n_perm=PARAMS["n_permutation"],
                                                top_k=None)
        h1.append(dict(cohort=arm.upper(), hypothesis="H1 (OA vs control)",
                       delta_average_degree=round(delta, 3),
                       p_permutation=f"{p:.4g}",
                       n_permutation=PARAMS["n_permutation"]))
        print(f"  [{arm}] H1 Δ平均度 = {delta:+.3f}, 置换 P = {p:.4g}")
    pd.DataFrame(h1).to_csv(os.path.join(OUT_DIR, "table5_H1_permutation.csv"),
                            index=False)
    out["H1"] = h1

    # H2 / H4：跨库比较（两样本 bootstrap，**用等密度稀疏化后的网络**）
    cc = bootstrap_cross_country(built_k["charls"]["X_oa"], built_k["ukb"]["X_oa"],
                                 n_boot=PARAMS["n_bootstrap"], top_k=k)
    rows = []
    for key, v in cc.items():
        label = "Average degree (whole network)" if key == "average_degree" \
            else f"Module: {key}"
        rows.append(dict(metric=label,
                         delta_CN_minus_UK=round(v["delta"], 3),
                         ci_low=round(v["ci"][0], 3), ci_high=round(v["ci"][1], 3),
                         p_bootstrap=f"{v['p']:.4g}",
                         significant="YES" if v["p"] < 0.05 else "no"))
        print(f"  {label:38s} Δ = {v['delta']:+.3f} "
              f"[{v['ci'][0]:+.3f}, {v['ci'][1]:+.3f}] P = {v['p']:.4g}")
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(OUT_DIR, "table6_H2H4_cross_country.csv"), index=False)
    out["H2H4"] = rows
    return out


# ---------------------------------------------------------------- main
def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    cohorts = {arm: pd.read_csv(os.path.join(OUT_DIR, f"cohort_{arm}.csv"))
               for arm in ["charls", "ukb"]}
    print("▶ 描述层")
    table2_prevalence(cohorts)
    table2_pairs(cohorts)

    print("\n▶ 网络层")
    built = {arm: build_arm(cohorts, arm) for arm in ["charls", "ukb"]}
    k = equal_density_k(built)
    print(f"  未稀疏化边数：CHARLS OA/CT = "
          f"{int(np.triu(built['charls']['adj_oa'], 1).sum())}/"
          f"{int(np.triu(built['charls']['adj_ct'], 1).sum())}，"
          f"UKB OA/CT = {int(np.triu(built['ukb']['adj_oa'], 1).sum())}/"
          f"{int(np.triu(built['ukb']['adj_ct'], 1).sum())}")
    print(f"  跨库等密度稀疏化 K = {k}")
    built_k = {arm: build_arm(cohorts, arm, top_k=k) for arm in ["charls", "ukb"]}

    t3 = table3_metrics(built, k)
    table4_hubs(built)
    print("▶ 推断层")
    inference(built, k, built_k)
    print("✅ 分析完成")
    return built, k, built_k


if __name__ == "__main__":
    main()
