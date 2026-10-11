"""
lib_network.py — 共病网络核心算法（全向量化）
================================================================================
对应 R 实现：igraph / qgraph / NetworkComparisonTest
本模块的所有函数均对"个体 × 疾病"的 0/1 矩阵 X 操作，保证 10,000 次置换可承受。
================================================================================
"""

import numpy as np
import networkx as nx
from scipy import stats
from scipy.sparse import csr_matrix

from config import DKEYS, DMODULE, MODULES, PARAMS


# ----------------------------------------------------------------------------
# 1. 超额共现（O/E）与边的选择
# ----------------------------------------------------------------------------
def cooccurrence(X):
    """cooccur[i,j] = 同时患 i 与 j 的人数（含对角线的 n_i）。"""
    XtX = X.T @ X
    return np.asarray(XtX, dtype=float)


def oe_matrix(X, eps=1e-9):
    """
    超额共现比 O/E = (n * O_ij) / (n_i * n_j)
    在"两病独立"的零假设下 O/E 的期望为 1。
    """
    n, _ = X.shape
    co = cooccurrence(X)
    ni = np.diag(co).copy()
    expected = np.outer(ni, ni) / max(n, 1)
    return co / (expected + eps)


def edge_stats(X, eps=0.5):
    """
    为每一个疾病对计算：共现数、O/E、log(OR)、SE、z、双侧 P。
    使用 OR 的 Woolf 对数标准误（向量化），避免 153 次 Fisher 精确检验的开销。
    返回 dict of (p×p) 矩阵。
    """
    n, p = X.shape
    co = cooccurrence(X)
    ni = np.diag(co).copy()
    a = co.astype(float)
    b = ni[:, None] - a
    c = ni[None, :] - a
    d = n - ni[:, None] - ni[None, :] + a

    # Haldane–Anscombe 修正，避免 0 单元格
    aa, bb, cc, dd = a + eps, b + eps, c + eps, d + eps
    expected = np.outer(ni, ni) / n
    oe = a / (expected + 1e-9)
    log_or = np.log((aa * dd) / (bb * cc))
    se = np.sqrt(1.0 / aa + 1.0 / bb + 1.0 / cc + 1.0 / dd)
    z = log_or / se
    pval = 2.0 * stats.norm.sf(np.abs(z))
    return {"co": a, "oe": oe, "logor": log_or, "se": se, "z": z, "p": pval}


def bh_fdr(pvals):
    """Benjamini–Hochberg 校正，返回与输入同形状的 q 值。"""
    p = np.asarray(pvals, dtype=float).ravel()
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out.reshape(np.asarray(pvals).shape)


def select_edges(stats_dict, n_nodes, min_prev=None, prev=None,
                 oe_threshold=1.0, alpha=0.05, top_k=None):
    """
    边的选择规则（与预注册一致）：
      1. 仅保留患病率 > min_prev 的节点
      2. O/E > oe_threshold
      3. BH-FDR < alpha
      4. 若指定 top_k，只保留 O/E 最大的 top_k 条（等密度稀疏化）
    返回 (adjacency 0/1 矩阵, q 值矩阵, 有效节点 mask)
    """
    if min_prev is not None and prev is not None:
        keep_nodes = prev > min_prev
    else:
        keep_nodes = np.ones(n_nodes, dtype=bool)

    # FDR 只在【唯一疾病对】（上三角）上做。
    # 若把整个 p×p 矩阵（含对角线与重复）丢进 BH，等于把多重度放大约 2 倍，过度保守。
    q = np.zeros((n_nodes, n_nodes))
    q[np.triu_indices(n_nodes, k=1)] = bh_fdr(
        stats_dict["p"][np.triu_indices(n_nodes, k=1)])
    q = q + q.T
    np.fill_diagonal(q, 1.0)

    adj = (stats_dict["oe"] > oe_threshold) & (q < alpha)

    iu = np.triu_indices(n_nodes, k=1)
    mask = keep_nodes[:, None] & keep_nodes[None, :]
    adj = adj & mask
    adj = np.triu(adj, 1)
    adj = adj + adj.T

    if top_k is not None:
        oe = stats_dict["oe"].copy()
        oe[~np.triu(np.ones_like(oe, dtype=bool), 1)] = -np.inf
        vals = oe[iu]
        valid = np.isfinite(vals) & (vals > oe_threshold) & (q[iu] < alpha) \
                & (keep_nodes[:, None] & keep_nodes[None, :])[iu]
        idx = np.where(valid)[0]
        if len(idx) > top_k:
            take = idx[np.argsort(-vals[idx])[:top_k]]
            newadj = np.zeros_like(adj, dtype=bool)
            for t in take:
                i, j = iu[0][t], iu[1][t]
                newadj[i, j] = newadj[j, i] = True
            adj = newadj
    return adj.astype(float), q, keep_nodes


# ----------------------------------------------------------------------------
# 2. 网络指标
# ----------------------------------------------------------------------------
def network_metrics(adj, node_mask=None):
    """
    报告 5 项指标（预注册要求：不能只报平均度）。
      average_degree : 平均度（主指标）
      density        : 密度
      modularity     : 模块度 Q（按预先固定的病因模块，非数据驱动）
      centralization : 度中心化指数
      clustering     : 平均聚类系数
    """
    G = nx.from_numpy_array(adj)
    n = G.number_of_nodes()
    if n == 0 or G.number_of_edges() == 0:
        return dict(average_degree=0.0, density=0.0, modularity=np.nan,
                    centralization=0.0, clustering=0.0, n_edges=0)

    deg = np.array([d for _, d in G.degree()], dtype=float)
    avg_deg = deg.mean()
    density = nx.density(G) if n > 1 else 0.0
    try:
        clustering = nx.average_clustering(G)
    except Exception:
        clustering = np.nan
    centralization = (deg.max() - deg).sum() / ((n - 1) * (n - 2)) if n > 2 else np.nan

    communities = None
    if node_mask is None:
        node_mask = np.ones(n, dtype=bool)
    try:
        comms = []
        for i in range(n):
            if node_mask[i]:
                comms.append(DMODULE.get(DKEYS[i], "organ_other"))
            else:
                comms.append("__excluded__")
        communities = nx.algorithms.community.modularity(G, _groups(comms))
    except Exception:
        communities = np.nan

    return dict(average_degree=float(avg_deg), density=float(density),
                modularity=float(communities), centralization=float(centralization),
                clustering=float(clustering), n_edges=int(G.number_of_edges()))


def _groups(comms):
    d = {}
    for i, c in enumerate(comms):
        d.setdefault(c, set()).add(i)
    return list(d.values())


def module_degree(adj, node_mask=None):
    """模块内平均度：每个模块内部子图的平均度。"""
    G = nx.from_numpy_array(adj)
    n = G.number_of_nodes()
    if node_mask is None:
        node_mask = np.ones(n, dtype=bool)
    out = {}
    for m in MODULES:
        idx = [i for i in range(n) if node_mask[i]
               and DMODULE.get(DKEYS[i], "") == m]
        if len(idx) < 2:
            out[m] = np.nan
            continue
        sub = G.subgraph(idx)
        degs = [d for _, d in sub.degree()]
        out[m] = float(np.mean(degs)) if degs else 0.0
    return out


def module_strength(stats_dict, node_mask=None):
    """
    模块内聚集强度 = 该模块内所有疾病对的 O/E 均值。

    ⚠️ 为什么不用"模块内平均度"：
      平均度是计数型指标，会饱和。3 节点模块上限=2，4 节点上限=3。
      实测中英两国的网络都撞到这个上限，导致 Δ 恒为 0，假设 H2 无法检验。
      改用 O/E 的均值后指标连续、不饱和、且直接回答
      "该模块内部的疾病共现超出随机期望多少倍"，正是 H2 要问的问题。
    """
    n = stats_dict["oe"].shape[0]
    if node_mask is None:
        node_mask = np.ones(n, dtype=bool)
    out = {}
    for m in MODULES:
        idx = [i for i in range(n) if node_mask[i] and DMODULE[DKEYS[i]] == m]
        vals = [stats_dict["oe"][i, j]
                for a, i in enumerate(idx) for j in idx[a + 1:]]
        out[m] = float(np.mean(vals)) if vals else np.nan
    return out


def degree_centrality(adj):
    G = nx.from_numpy_array(adj)
    return nx.degree_centrality(G)


# ----------------------------------------------------------------------------
# 3. 推断：置换检验与 bootstrap
# ----------------------------------------------------------------------------
def _rebuild_metric(X_oa, X_ct, prev_min=PARAMS["min_prevalence"], top_k=None):
    """给定两组的疾病矩阵，重建网络并返回 OA 组平均度与模块内平均度。"""
    res = {}
    for tag, X in (("oa", X_oa), ("ct", X_ct)):
        prev = X.mean(axis=0)
        st = edge_stats(X)
        adj, q, keep = select_edges(st, X.shape[1], min_prev=prev_min, prev=prev,
                                    oe_threshold=PARAMS["oe_threshold"],
                                    alpha=PARAMS["alpha"], top_k=top_k)
        res[tag] = (network_metrics(adj, keep)["average_degree"],
                    module_degree(adj, keep))
    return res


def permutation_test_pairs(X_oa, X_ct, n_perm=PARAMS["n_permutation"],
                           top_k=None, seed=20260908):
    """
    【假设 H1】OA vs 匹配对照的平均度差异 —— 配对内标签置换检验。
    零假设：在同一匹配对内，谁是 OA 谁是对照与共病结构无关。
    返回 obs_delta, p_value, null_distribution
    """
    rng = np.random.default_rng(seed)
    prev_min = PARAMS["min_prevalence"]

    obs = _rebuild_metric(X_oa, X_ct, prev_min, top_k)
    obs_delta = obs["oa"][0] - obs["ct"][0]

    n = X_oa.shape[0]
    flip = rng.random((n_perm, n)) < 0.5
    null = np.empty(n_perm)
    for b in range(n_perm):
        f = flip[b][:, None]
        A = np.where(f, X_oa, X_ct)
        C = np.where(f, X_ct, X_oa)
        r = _rebuild_metric(A, C, prev_min, top_k)
        null[b] = r["oa"][0] - r["ct"][0]

    p = (np.sum(np.abs(null) >= abs(obs_delta) - 1e-12) + 1) / (n_perm + 1)
    return obs_delta, float(p), null


def bootstrap_ci_metric(X, n_boot=PARAMS["n_bootstrap"], top_k=None,
                        seed=20260908, stat="average_degree"):
    """
    单库指标（网络整体平均度 / 各模块内平均度）的参与者级 bootstrap 95%CI。
    """
    rng = np.random.default_rng(seed)
    n = X.shape[0]
    prev_min = PARAMS["min_prevalence"]
    vals, mods = [], {m: [] for m in MODULES}
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        Xb = X[idx]
        prev = Xb.mean(axis=0)
        st = edge_stats(Xb)
        adj, q, keep = select_edges(st, X.shape[1], min_prev=prev_min, prev=prev,
                                    oe_threshold=PARAMS["oe_threshold"],
                                    alpha=PARAMS["alpha"], top_k=top_k)
        vals.append(network_metrics(adj, keep)["average_degree"])
        md = module_degree(adj, keep)
        for m in MODULES:
            mods[m].append(md[m])
    out = {stat: (float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5)))}
    for m in MODULES:
        arr = np.array([v for v in mods[m] if np.isfinite(v)])
        out[m] = (float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5))) \
            if len(arr) else (np.nan, np.nan)
    return out


def bootstrap_cross_country(X_cn, X_uk, n_boot=PARAMS["n_bootstrap"],
                            top_k=None, seed=20260908):
    """
    【假设 H2/H4】跨库比较：两样本 bootstrap 得到 Δ（中国 − 英国）的 95%CI。
    分别对整体平均度与各模块内平均度计算。CI 不含 0 判为显著。
    """
    rng = np.random.default_rng(seed)
    prev_min = PARAMS["min_prevalence"]
    dc, mc = [], {m: [] for m in MODULES}
    du, mu = [], {m: [] for m in MODULES}
    for _ in range(n_boot):
        for X, (dstore, mstore) in ((X_cn, (dc, mc)), (X_uk, (du, mu))):
            idx = rng.integers(0, X.shape[0], X.shape[0])
            Xb = X[idx]
            prev = Xb.mean(axis=0)
            st = edge_stats(Xb)
            adj, q, keep = select_edges(st, X.shape[1], min_prev=prev_min, prev=prev,
                                        oe_threshold=PARAMS["oe_threshold"],
                                        alpha=PARAMS["alpha"], top_k=top_k)
            dstore.append(network_metrics(adj, keep)["average_degree"])
            ms = module_strength(st, keep)
            for m in MODULES:
                mstore[m].append(ms[m])
    res = {}
    dc, du = np.array(dc), np.array(du)
    d = dc - du
    res["average_degree"] = _delta_summary(d)
    for m in MODULES:
        a = np.array(mc[m]); b = np.array(mu[m])
        ok = np.isfinite(a) & np.isfinite(b)
        if ok.sum() < 10:
            res[m] = dict(delta=np.nan, ci=(np.nan, np.nan), p=np.nan)
            continue
        res[m] = _delta_summary(a[ok] - b[ok])
    return res


def _delta_summary(d):
    """
    两样本 bootstrap 差异的汇总。p 为双侧 bootstrap P 值（裁剪至 [0,1]）。
    当 bootstrap 分布完全退化为 0（罕见，说明指标不随样本变动）时返回 P=1。
    """
    d = np.asarray(d, dtype=float)
    if np.allclose(d, 0):
        return dict(delta=0.0, ci=(0.0, 0.0), p=1.0)
    p = 2.0 * min((d <= 0).mean(), (d >= 0).mean())
    return dict(delta=float(np.mean(d)),
                ci=(float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))),
                p=float(min(1.0, max(0.0, p))))
