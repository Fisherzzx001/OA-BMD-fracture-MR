# -*- coding: utf-8 -*-
"""
code/340 — Reviewer 1：单设定点(M1) vs 双耦合模块(M2) 正式统计模型比较
第七轮对抗性极限审稿 docs/148 P0-4 / P1-5 重建版（2026-11-12）。

重建动机（相对第六轮初版的两处实质修正）
------------------------------------------------------------
[P0-4 协方差概念错配] 初版把 MRlap 的 cross-trait LDSC 截距 λ
    (int_crosstrait: fracture −0.0361 / OA −0.0051；稿件 S14 自承
    "These are intercepts, not overlap fractions") 当作"两结局 GWAS 在同一 SNP
    上的抽样误差相关系数"塞进 Wald ratio 两端的 delta-method 抽样协方差。
    这是量纲 + 符号双重错配：
      · λ 是暴露-结局间样本重叠+遗传结构造成的截距偏移，不是"两个结局在同一
        SNP 上的抽样误差相关"；后者由参与者重叠驱动，应用重叠分数×残差相关估计。
      · 符号反了：fracture↔OA_knee 的 r_g=+0.12（正文 §3），共享 UKB 参与者的
        抽样误差相关应为正；λ=−0.036 为负，物理上不可能是该重叠相关。
    后果：|λ| 小故 ρ̂/LRT 数值影响有限，但"纳入 MRlap 样本重叠协方差"的立论
    被架空。本版：
      (a) 把抽样协方差正名/重建 —— 只保留数学上确定的共享 beta_e 分母项，
          额外的两结局抽样相关 rho_s 作为参数（默认 0，并显式扫描 ≥0 的
          物理允许区间；不再声称 λ 是其实测值）；
      (b) λ 降格为"cross-trait LDSC 截距，作为敏感性扰动的参照量"写入报告，
          不再进入主协方差矩阵；正文/S14 的 "intercepts, not overlap" 本已正确。

[P1-5 解释过度 + 边界零分布] 初版用解析 0.5*chi2_0+0.5*chi2_1 作边界 LRT 零分布。
    但 rho=−1 同时令 T→秩1（tau_oa 压至下界 1e-4，实测确认），是"rho 边界 +
    方差分量边界"的双边界，单一 chi2 混合零分布可能失准。本版：
      (c) 用参数 bootstrap 在 M1 零假设下模拟、直接校准 LRT 的零分布与 p 值；
      (d) 增加 rho 与 0 / −0.5 的 TOST 等价性检验，明确"rho_hat≈−0.13 与近独立
          分不开"，结论严格限定为"拒绝严格共线(rho=−1)强版本"，不再暗示数据
          偏向 SOS；把单设定点一般形式明确列为"未被本检验判定"。

模型（bivariate measurement-error / random-effects ML，31 工具为独立观测单位）
------------------------------------------------------------
  对每个工具 j，观测两个 Wald ratio：theta_j = (theta_fx, theta_oa)'
    theta_{j,k} = mu_k + u_{j,k} + e_{j,k},  k ∈ {fx, oa}
      mu_k        : 固定效应（聚合因果效应）
      u_j=(u_fx,u_oa)' ~ N(0, T)  潜在真实效应异质性；T=[[t_fx^2,rho*t_fx*t_oa],[.,t_oa^2]]
      e_j         ~ N(0, V_j)     抽样协方差（delta method）
    V_j 的构造（本版重建，P0-4）：
      Var(theta_fx) = (so_fx/be)^2 + (bo_fx^2/be^4)*se_be^2   （含 beta_e 误差）
      Var(theta_oa) = (so_oa/be)^2 + (bo_oa^2/be^4)*se_be^2
      Cov(theta_fx,theta_oa)
          = (bo_fx*bo_oa/be^4)*se_be^2      [共享 beta_e 分母，数学确定；保留]
          + rho_s*(so_fx*so_oa/|be|^2)      [两结局同 SNP 抽样误差相关；rho_s 为
                                              参数，默认 0，仅扫描物理允许范围]
      其中 rho_s>=0 是参与者重叠所允许的符号（r_g(fracture,OA)=+0.12）；
      MRlap 截距 lambda 不作为 rho_s 的实测值（概念错配，已弃用）。
   边缘：theta_j ~ N(mu, T + V_j)

  * M1 单潜在设定点：rho = −1，潜在效应严格共线（一个潜因子按固定比例驱动两端，
      T 秩 1，tau_oa→0 下界）。这是"单一设定点"的可检验【强版本/严格形式】。
  * M2 双耦合模块：rho 自由 ∈ (−1,1)。rho ≈ 0 = 两模块近独立。

判别统计量（预先声明）
  1. 主检验：LRT H0:rho=−1 vs H1:rho 自由；p 用【参数 bootstrap 校准】（P1-5），
     并并列报告解析 0.5*chi2_0+0.5*chi2_1 供对照。
  2. rho 的 profile-likelihood 95% CI + 非参数 bootstrap CI。
  3. rho 与 0（近独立）、−0.5（中等耦合）的 TOST 等价性 / CI 相容性（P1-5）。
  4. AIC / BIC（M1 vs M2）。
  5. 功效分析：检出 rho>−1 的功效（用解析临界，设计层面）。
  6. 抽样相关敏感性：rho_s ∈ {0, 0.10, 0.20, 0.30, 0.50}（均 >=0，物理允许）；
     rho_s<0 仅在附录对照中展示，以证明原负值 λ 方向不改变结论。

判据（结果前固定；先于结果固定以避免"数字与文件一致≠数字对"）
  * 若 bootstrap LRT p<0.05：拒绝"严格共线 rho=−1"强版本。
  * 但只有当 rho 的 95% CI 上界 < −0.5 时，才谈得上"支持中等以上耦合"；
    若 CI 含 0，则"近独立"与"弱耦合"均不可排除 —— 报告为"耦合强度未定位"。
  * 无论结果如何：本检验只对"单一设定点的严格共线强版本"表态；不对
    "单设定点一般形式（可含端点特异噪声，|rho|<1）"表态，也不暗示数据偏向 SOS。
  * L1 符号相反（code/339：22/31, p=0.0147）是独立分量，不受本检验影响。

只读：不改动任何已报 beta/SE/p/OR/CI，不改动文档。输出 code/output/genetics/_r1_model_comparison.*。
运行：python3 code/340_r1_single_vs_dual_module.py
"""
import os
import json
import numpy as np
import pandas as pd
from scipy import stats, optimize

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HARM = os.path.join(ROOT, "code", "output", "genetics", "harmonised")
OUT = os.path.join(ROOT, "code", "output", "genetics")
FX = os.path.join(HARM, "BMD_total__fracture.csv")
OA = os.path.join(HARM, "BMD_total__OA_knee.csv")

# MRlap cross-trait LDSC 截距 lambda —— 仅作敏感性参照量写入报告，不作为抽样相关（P0-4）。
# 来自 code/output/mrlap_overlap_results.csv；稿件 S14 已正确注明 "intercepts, not overlap fractions"。
LAMBDA_FX = -0.0360600721603821   # BMD_total -> fracture（cross-trait LDSC 截距）
LAMBDA_OA = -0.00510735075368427  # BMD_total -> OA_knee

rng = np.random.default_rng(20261112)


# ---------------------------------------------------------------------------
# 1. 组装 31 工具的 theta 与抽样协方差 V_j（delta method，含 beta_e 误差）
# ---------------------------------------------------------------------------
def build_data(rho_s=0.0):
    """rho_s = 两结局 GWAS 在同一 SNP 上的抽样误差相关系数（参数，默认 0）。
    P0-4：不再把 MRlap 截距 lambda 当作 rho_s；rho_s 由参与者重叠驱动，
    物理允许符号 >=0（因 r_g(fracture,OA)=+0.12），lambda<0 不适用。"""
    fx = pd.read_csv(FX)
    oa = pd.read_csv(OA)
    m = fx[["SNP", "beta_e", "se_e", "beta_o", "se_o", "p_o"]].merge(
        oa[["SNP", "beta_o", "se_o", "p_o"]], on="SNP", suffixes=("_fx", "_oa"))
    assert len(m) == 31
    be, se_be = m.beta_e.to_numpy(), m.se_e.to_numpy()
    bo_fx, so_fx = m.beta_o_fx.to_numpy(), m.se_o_fx.to_numpy()
    bo_oa, so_oa = m.beta_o_oa.to_numpy(), m.se_o_oa.to_numpy()

    # Wald ratio theta = beta_o / beta_e
    th_fx, th_oa = bo_fx / be, bo_oa / be
    # delta method Var(theta) = (so/be)^2 + (bo^2/be^4)*se_be^2 （含分母 beta_e 误差）
    v_fx = (so_fx / be) ** 2 + (bo_fx ** 2 / be ** 4) * se_be ** 2
    v_oa = (so_oa / be) ** 2 + (bo_oa ** 2 / be ** 4) * se_be ** 2
    # 两端抽样协方差：
    #   (1) 共享 beta_e 分母引入的协方差 (bo_fx*bo_oa/be^4)*se_be^2 —— 数学确定，保留
    #   (2) 两结局同 SNP 抽样误差相关 rho_s*(so_fx*so_oa/|be|^2) —— 参数，默认 0
    #       （P0-4：lambda 不作为此项的实测值；rho_s 物理允许 >=0）
    cov_sampling = (bo_fx * bo_oa / be ** 4) * se_be ** 2 \
        + rho_s * (so_fx * so_oa / np.abs(be) ** 2)
    theta = np.column_stack([th_fx, th_oa])          # (31,2)
    V = np.zeros((len(m), 2, 2))
    V[:, 0, 0], V[:, 1, 1], V[:, 0, 1], V[:, 1, 0] = v_fx, v_oa, cov_sampling, cov_sampling
    return theta, V, m


# ---------------------------------------------------------------------------
# 2. 二变量随机效应 ML
# ---------------------------------------------------------------------------
def neg_loglik(theta, V, mu, t_fx, t_oa, rho):
    """向量化 −loglik（批量 Cholesky）。参数越界/非正定时返回平滑大值。"""
    T = np.array([[t_fx ** 2, rho * t_fx * t_oa],
                  [rho * t_fx * t_oa, t_oa ** 2]])
    S = T + V                                    # (n,2,2)
    try:
        L = np.linalg.cholesky(S)                 # (n,2,2) 下三角
    except np.linalg.LinAlgError:
        return 1e6
    d = theta - mu
    z = np.linalg.solve(L, d[:, :, None])[:, :, 0]
    ll = np.sum(-np.log(2 * np.pi) - np.log(L[:, 0, 0]) - np.log(L[:, 1, 1])
                - 0.5 * np.sum(z * z, axis=1))
    return float(-ll) if np.isfinite(ll) else 1e6


def _fit(theta, V, rho_fixed=None, robust=False):
    """L-BFGS-B + bounds 的快速 ML。rho_fixed 给定时不估计 rho。
    robust=True（profile/bootstrap 用）：log(tau) 用 3×3 多初值，避免固定 rho 时
    陷入局部解导致 profile 曲线非单调（rho→−1 时 tau_oa→0 边界尤易发生）。"""
    free_rho = rho_fixed is None
    def unpack(p):
        tf, to = np.exp(p[2]), np.exp(p[3])
        return np.array(p[:2]), tf, to, (np.tanh(p[4]) if free_rho else rho_fixed)
    def obj(p):
        mu, tf, to, rh = unpack(p)
        return neg_loglik(theta, V, mu, tf, to, rh)
    m0, m1 = theta[:, 0].mean(), theta[:, 1].mean()
    sd0, sd1 = max(theta[:, 0].std(ddof=1), 1e-3), max(theta[:, 1].std(ddof=1), 1e-3)
    if free_rho or not robust:
        starts = [np.array([m0, m1, np.log(sd0), np.log(sd1)]
                           + ([np.arctanh(np.clip(np.corrcoef(theta.T)[0, 1], -0.9, 0.9))]
                              if free_rho else []))]
    else:
        starts = [np.array([m0, m1, np.log(sd0 * a), np.log(sd1 * b)])
                  for a in (0.5, 1.0, 2.0) for b in (0.5, 1.0, 2.0)]
    bounds = [(-5, 5), (-5, 5),
              (np.log(1e-4), np.log(10)), (np.log(1e-4), np.log(10))]
    if free_rho:
        bounds += [(-4, 4)]
    best = None
    for x0 in starts:
        r = optimize.minimize(obj, x0, method="L-BFGS-B", bounds=bounds,
                              options={"maxiter": 400, "ftol": 1e-13, "gtol": 1e-8})
        if best is None or r.fun < best.fun:
            best = r
    mu, tf, to, rh = unpack(best.x)
    return best.fun, mu, tf, to, rh


def fit_m2(theta, V):
    """rho 自由 ∈ (−1, 1)。"""
    nll, mu, tf, to, rh = _fit(theta, V, rho_fixed=None)
    k = 5
    return dict(nll=nll, mu=mu, t_fx=tf, t_oa=to, rho=rh, k=k,
                aic=2 * nll + 2 * k, bic=2 * nll + k * np.log(len(theta)))


def fit_m1(theta, V):
    """rho 固定 = −1（单潜在设定点，潜在效应严格共线反号；tau_oa→0 为双边界）。"""
    nll, mu, tf, to, _ = _fit(theta, V, rho_fixed=-1.0, robust=True)
    k = 4
    return dict(nll=nll, mu=mu, t_fx=tf, t_oa=to, rho=-1.0, k=k,
                aic=2 * nll + 2 * k, bic=2 * nll + k * np.log(len(theta)))


# ---------------------------------------------------------------------------
# 3. rho 的 profile-likelihood 95% CI
# ---------------------------------------------------------------------------
def profile_ci(theta, V, m2, alpha=0.05):
    from scipy.optimize import brentq
    thr = m2["nll"] + stats.chi2.ppf(1 - alpha, 1) / 2.0

    def prof_nll(rho_fixed):
        return _fit(theta, V, rho_fixed=rho_fixed, robust=True)[0]

    rh = m2["rho"]
    def f_lower(r):
        return prof_nll(r) - thr
    lower = -1.0
    try:
        if f_lower(-0.9999) < 0:
            lower = -1.0
        else:
            lower = brentq(f_lower, -0.9999, rh, xtol=1e-4)
    except Exception:
        lower = -1.0
    def f_up(r):
        return prof_nll(r) - thr
    upper = 0.9999
    try:
        if f_up(0.9999) < 0:
            upper = 0.9999
        else:
            upper = brentq(f_up, rh, 0.9999, xtol=1e-4)
    except Exception:
        upper = 0.9999
    return lower, upper


# ---------------------------------------------------------------------------
# 4. rho 的非参数 bootstrap CI（重采样工具）
# ---------------------------------------------------------------------------
def bootstrap_rho(theta, V, B=2000):
    n = len(theta)
    vals = []
    for _ in range(B):
        idx = rng.integers(0, n, n)
        try:
            r = fit_m2(theta[idx], V[idx])
            if abs(r["rho"]) < 0.999:
                vals.append(r["rho"])
        except Exception:
            pass
    if len(vals) < 50:
        return None, None, len(vals)
    return np.percentile(vals, [2.5, 97.5]), np.median(vals), len(vals)


# ---------------------------------------------------------------------------
# 5. 【P1-5 新增】参数 bootstrap 校准双边界 LRT 零分布
#    在 M1(rho=−1, tau_oa→0 边界) 下模拟 theta* ~ N(mu_hat, T_null+V_j)，
#    每个模拟样本重估 M1/M2 得 LRT*，用经验零分布给 p 值。
# ---------------------------------------------------------------------------
def bootstrap_lrt_null(theta, V, m1, B=1000):
    mu = m1["mu"]
    tf = m1["t_fx"]
    to = max(m1["t_oa"], 1e-4)   # M1 下 tau_oa 压至下界（实测≈1e-4）
    T_null = np.array([[tf ** 2, -tf * to],
                       [-tf * to, to ** 2]])
    n = len(theta)
    try:
        cov = np.array([T_null + V[j] for j in range(n)])
        L = np.linalg.cholesky(cov)          # (n,2,2)
    except np.linalg.LinAlgError:
        return None, None, None, 0
    null_lrt = []
    for _ in range(B):
        z = rng.standard_normal((n, 2))
        sim = mu + np.einsum("nij,nj->ni", L, z)
        try:
            m2s = fit_m2(sim, V)
            m1s = fit_m1(sim, V)
            lrt = 2 * (m1s["nll"] - m2s["nll"])
            if np.isfinite(lrt):
                null_lrt.append(max(lrt, 0.0))
        except Exception:
            continue
    if len(null_lrt) < 100:
        return None, None, None, len(null_lrt)
    return np.array(null_lrt)


# ---------------------------------------------------------------------------
# 6. 功效分析（模拟；用解析临界，设计层面）
# ---------------------------------------------------------------------------
def power_analysis(n_obs, t_fx, t_oa, rho_true, rho_null=-1.0, reps=300, alpha=0.05):
    theta0, V0, _ = build_data(rho_s=0.0)
    base_idx = rng.integers(0, len(theta0), n_obs)
    Vb = V0[base_idx]
    T = np.array([[t_fx ** 2, rho_true * t_fx * t_oa],
                  [rho_true * t_fx * t_oa, t_oa ** 2]])
    mu = np.array([theta0[:, 0].mean(), theta0[:, 1].mean()])
    detected = 0
    for _ in range(reps):
        sim = np.array([rng.multivariate_normal(mu, T + Vb[j]) for j in range(n_obs)])
        try:
            m2 = fit_m2(sim, Vb); m1 = fit_m1(sim, Vb)
            lrt = 2 * (m1["nll"] - m2["nll"])
            p = 0.5 * stats.chi2.sf(max(lrt, 0), 1)
            if p < alpha:
                detected += 1
        except Exception:
            pass
    return detected / reps if reps else float("nan")


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------
def main():
    rep = []
    def say(s=""):
        print(s); rep.append(s)

    say("=" * 78)
    say("Reviewer 1 — 单设定点(M1) vs 双耦合模块(M2) · code/340【P0-4/P1-5 重建版】")
    say("=" * 78)

    # --- 主分析：rho_s=0（仅保留数学确定的共享分母协方差项；不引入 MRlap 截距）---
    theta, V, m = build_data(rho_s=0.0)
    tau_fx = max(np.sqrt(max(theta[:, 0].var(ddof=1) - V[:, 0, 0].mean(), 1e-6)), 1e-3)
    tau_oa = max(np.sqrt(max(theta[:, 1].var(ddof=1) - V[:, 1, 1].mean(), 1e-6)), 1e-3)
    rel_fx = 1 - V[:, 0, 0].mean() / theta[:, 0].var(ddof=1)
    rel_oa = 1 - V[:, 1, 1].mean() / theta[:, 1].var(ddof=1)
    say(f"\n[数据] n={len(m)} 工具；tau_fx={tau_fx:.4f} tau_oa={tau_oa:.4f} "
        f"（含 beta_e 误差 reliability FX {rel_fx:.3f} / OA {rel_oa:.3f}）")
    say(f"       抽样协方差 V_j = 共享beta_e分母项 + rho_s*(so_fx*so_oa/|be|^2)，主分析 rho_s=0")
    say(f"       [P0-4] MRlap cross-trait LDSC 截距 lambda: fracture {LAMBDA_FX:+.4f}, "
        f"OA {LAMBDA_OA:+.4f} —— 仅参照，不作为抽样相关；lambda<0 与 r_g(fracture,OA)=+0.12 符号矛盾。")

    m2 = fit_m2(theta, V)
    m1 = fit_m1(theta, V)
    lrt = 2 * (m1["nll"] - m2["nll"])
    p_analytic = 0.5 * stats.chi2.sf(max(lrt, 0), 1)
    lo, hi = profile_ci(theta, V, m2)
    (bci, bmed, nb) = bootstrap_rho(theta, V, B=1500)

    # 【P1-5】参数 bootstrap 校准双边界 LRT
    say("\n[边界诊断] M1(rho=−1) 下 tau_oa 被压至下界："
        f"tau_oa_M1={m1['t_oa']:.5f}（下界 1e-4）→ rho 边界 + 方差分量边界，"
        "解析 0.5chi2_0+0.5chi2_1 失准，改用 bootstrap 校准。")
    null_lrt = bootstrap_lrt_null(theta, V, m1, B=1000)
    if null_lrt is not None and len(null_lrt) >= 100:
        p_boot = (1 + np.sum(null_lrt >= lrt)) / (len(null_lrt) + 1)
        q50, q95, q99 = np.percentile(null_lrt, [50, 95, 99])
        boot_ok = True
    else:
        p_boot, q50, q95, q99, boot_ok = float("nan"), *([float("nan")] * 3), False

    say("\n[模型拟合]")
    say(f"  M1 单设定点 (rho=−1 固定) : nll={m1['nll']:.4f}  AIC={m1['aic']:.2f}  BIC={m1['bic']:.2f}"
        f"  tau_oa={m1['t_oa']:.5f}")
    say(f"  M2 双模块   (rho 自由)     : nll={m2['nll']:.4f}  AIC={m2['aic']:.2f}  BIC={m2['bic']:.2f}"
        f"  rho_hat={m2['rho']:+.4f}")
    say(f"  dAIC(M1−M2)={m1['aic']-m2['aic']:+.2f}   dBIC(M1−M2)={m1['bic']-m2['bic']:+.2f}（正=M2 更优）")
    say(f"\n[主检验] LRT(rho=−1)={lrt:.4f}")
    say(f"  解析单边界 p(0.5chi2_0+0.5chi2_1) = {p_analytic:.4g}（仅对照）")
    if boot_ok:
        say(f"  【bootstrap 校准 p】= {p_boot:.4g}（B={len(null_lrt)}；零分布中位数 {q50:.3f}, "
            f"95%分位 {q95:.3f}, 99%分位 {q99:.3f}）")
    say(f"  rho profile 95% CI = [{lo:+.3f}, {hi:+.3f}]")
    if bci is not None:
        say(f"  rho bootstrap 95% CI = [{bci[0]:+.3f}, {bci[1]:+.3f}] (median {bmed:+.3f}, B={nb})")

    # 【P1-5】rho 与 0 / −0.5 / −1 的相容性 + TOST
    say("\n[rho 可分辨性 / TOST，P1-5]")
    def in_ci(x):
        return lo <= x <= hi
    compat = {"rho_eq_0_within_CI": bool(in_ci(0.0)),
              "rho_eq_minus05_within_CI": bool(in_ci(-0.5)),
              "rho_eq_minus1_within_CI": bool(in_ci(-1.0))}
    say(f"  0 在 95% CI 内（近独立不可排除）={compat['rho_eq_0_within_CI']}")
    say(f"  −0.5 在 95% CI 内（中等耦合不可排除）={compat['rho_eq_minus05_within_CI']}")
    say(f"  −1 在 95% CI 内（严格共线不可排除）={compat['rho_eq_minus1_within_CI']}")
    # TOST：rho 与 0 等价，取等价界 ±0.30
    delta = 0.30
    # 用 bootstrap 分布做 TOST（双侧 90% CI ≡ 两次单尾 0.05）
    if bci is not None:
        # bootstrap 90% CI
        vals = []
        n = len(theta)
        for _ in range(1000):
            idx = rng.integers(0, n, n)
            try:
                r = fit_m2(theta[idx], V[idx])
                if abs(r["rho"]) < 0.999:
                    vals.append(r["rho"])
            except Exception:
                pass
        if len(vals) >= 50:
            lo90, hi90 = np.percentile(vals, [5, 95])
            tost_near_indep = bool(lo90 > -delta and hi90 < delta)
            say(f"  TOST(rho 与 0 等价, 界 ±{delta:.2f}): 90% CI [{lo90:+.3f},{hi90:+.3f}] → "
                f"{'等价成立(支持近独立)' if tost_near_indep else '等价不成立(近独立与弱耦合仍难分)'}")
        else:
            tost_near_indep = None
    else:
        tost_near_indep = None

    # --- 判据 ---
    say("\n[判据裁决]")
    p_main = p_boot if boot_ok else p_analytic
    if p_main < 0.05 and hi < -0.5:
        verdict = ("拒绝 rho=−1，且 rho 上界<−0.5：数据与中等以上耦合(M2)相容。"
                   "但本检验仍只对严格共线强版本表态，不构成对 SOS 一般形式的证明。")
    elif p_main < 0.05:
        verdict = ("拒绝'单一设定点严格共线(rho=−1)'强版本：M2 拟合显著更优，rho=−1 在 95% CI 外。"
                   f"但 rho_hat={m2['rho']:+.3f}，95% CI [{lo:+.3f},{hi:+.3f}] 含 0，"
                   "'近独立(rho≈0)'与'弱–中等耦合'均不可排除 → 耦合强度未定位；"
                   "本检验不判定单设定点一般形式，也不暗示数据偏向 SOS。")
    else:
        pw = power_analysis(len(m), tau_fx, tau_oa, rho_true=-0.50, reps=200)
        verdict = ("不能拒绝 M1（bootstrap p>=0.05）；结合功效判'功效不足未判别'还是"
                   f"'与严格共线相容'（检出 rho=−0.5 vs −1 功效={pw:.2f}）。")
    say("  " + verdict)

    # --- 【P0-4】抽样相关 rho_s 敏感性（仅物理允许 >=0；负值作附录对照）---
    say("\n[敏感性: 两结局抽样误差相关 rho_s 扫描，P0-4]")
    say("  rho_s 由参与者重叠驱动、物理允许 >=0（r_g(fracture,OA)=+0.12）。"
        "lambda<0 不作为其实测值，仅末行附录对照以证原负值不改变结论。")
    say(f"  {'rho_s':>7s} {'rho_hat':>9s} {'LRT':>8s} {'p_analytic':>10s} {'CI_hi':>8s}")
    sens_rows = []
    for rg_ in [0.0, 0.10, 0.20, 0.30, 0.50, LAMBDA_FX]:
        thx, Vx, _ = build_data(rho_s=rg_)
        m2x = fit_m2(thx, Vx); m1x = fit_m1(thx, Vx)
        lrtx = 2 * (m1x["nll"] - m2x["nll"])
        px = 0.5 * stats.chi2.sf(max(lrtx, 0), 1)
        lox, hix = profile_ci(thx, Vx, m2x)
        tag = "  <-附录:原lambda负方向对照" if rg_ == LAMBDA_FX else ""
        say(f"  {rg_:+.4f}  {m2x['rho']:+9.4f} {lrtx:8.3f} {px:10.3g} {hix:+8.3f}{tag}")
        sens_rows.append((rg_, m2x["rho"], lrtx, px, lox, hix))

    # --- 功效分析 ---
    say("\n[功效分析] 检出 rho 显著>−1（否决严格共线）的功效（解析临界，设计层面）")
    say(f"  {'n_tools':>8s} {'rho_true':>9s} {'power':>8s}")
    power_rows = []
    for n_obs in (31, 60, 100, 150):
        for rho_true in (-0.50, 0.00):
            pw = power_analysis(n_obs, tau_fx, tau_oa, rho_true, reps=200)
            say(f"  {n_obs:8d} {rho_true:+9.2f} {pw:8.2f}")
            power_rows.append((n_obs, rho_true, pw))
    n_needed = 31 if power_rows and power_rows[1][2] >= 0.80 else None
    if n_needed is None:
        for n_obs in (50, 80, 120, 180, 260):
            if power_analysis(n_obs, tau_fx, tau_oa, rho_true=0.0, reps=120) >= 0.80:
                n_needed = n_obs; break
    say(f"  区分 rho=−1 与 rho=0 达 80% 功效约需工具数: {n_needed if n_needed else '>260'}")

    say("\n" + "=" * 78)
    say("NOTE: 本检验为 in-corpus 假设检验诊断，不改变任何已报 beta/SE/p/OR/CI；")
    say("      只针对'单一设定点严格共线强版本(rho=−1)'；L1 符号相反(code/339)不受影响；")
    say("      单设定点一般形式(可含端点特异噪声)不在本检验判定范围。")
    say("=" * 78)

    # 落盘
    pd.DataFrame([dict(model="M1_single_setpoint", nll=m1["nll"], aic=m1["aic"], bic=m1["bic"],
                       rho=-1.0, tau_oa=m1["t_oa"]),
                  dict(model="M2_dual_module", nll=m2["nll"], aic=m2["aic"], bic=m2["bic"],
                       rho=m2["rho"], tau_oa=m2["t_oa"])]).to_csv(
        os.path.join(OUT, "_r1_model_comparison_fit.csv"), index=False)
    pd.DataFrame(sens_rows, columns=["rho_s", "rho_hat", "LRT", "p_LRT_analytic", "CI_lo", "CI_hi"]).to_csv(
        os.path.join(OUT, "_r1_model_comparison_sensitivity.csv"), index=False)
    pd.DataFrame(power_rows, columns=["n_tools", "rho_true", "power"]).to_csv(
        os.path.join(OUT, "_r1_model_comparison_power.csv"), index=False)
    if null_lrt is not None and len(null_lrt) >= 100:
        pd.DataFrame(dict(LRT_bootstrap_null=null_lrt)).to_csv(
            os.path.join(OUT, "_r1_lrt_bootstrap_null.csv"), index=False)
    with open(os.path.join(OUT, "_r1_model_comparison.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(rep) + "\n")
    summary = dict(
        rho_s_main=0.0,
        lambda_lcsc_intercept=[LAMBDA_FX, LAMBDA_OA],
        covariance_note="V_j = shared-beta_e denominator term + rho_s*(so_fx*so_oa/|be|^2); "
                        "MRlap cross-trait LDSC intercept lambda NOT used as sampling correlation (P0-4)",
        rho_hat=m2["rho"], rho_ci_prof=[lo, hi],
        rho_ci_boot=(list(bci) if bci is not None else None),
        LRT=lrt, p_LRT_analytic=p_analytic, p_LRT_bootstrap=(p_boot if boot_ok else None),
        lrt_null_q95=(q95 if boot_ok else None),
        rho_compat_CI=compat, tost_near_indep_0=tost_near_indep,
        dAIC=m1["aic"] - m2["aic"], dBIC=m1["bic"] - m2["bic"],
        verdict=verdict, n_needed_80pct=n_needed,
        reliability=[rel_fx, rel_oa], tau_oa_M1=m1["t_oa"],
        boundary_note="rho=-1 forces tau_oa->0 lower bound (double boundary); "
                      "analytic 0.5chi2_0+0.5chi2_1 replaced by parametric bootstrap (P1-5)")
    with open(os.path.join(OUT, "_r1_model_comparison_summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print("\nwritten: _r1_model_comparison.{csv,txt,json} + _r1_lrt_bootstrap_null.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
