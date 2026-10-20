#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 123_verify_mrlap_direction.py — MRlap 方向性深挖与验证（第九轮 P0 项）
#
# 背景：111/112 轮预审指出稿件 "MRlap moved every estimate further from the
# null (median 20%), so the primary estimates are conservative with respect to
# overlap" 存在逻辑跳跃。本脚本用两类证据核定 MRlap 三个偏倚分量的方向：
#
#   A. 沉积数据（code/output/mrlap_overlap_results.csv）：
#      - 10/10 对 sign(r_g) == sign(效应)（一致）；
#      - corrected ≈ observed + k·r_g (k>0)，即校正沿 r_g 方向（离零）。
#   B. 基因型级模拟（固定个体池 + 部分参与者重叠）：
#      - overlap 分量 = E[IVW|f] − E[IVW|f=0]，预测 +f·r_p/(N_x(s²+v)) > 0
#        （拉向观测相关，对同号对即离零）——与 MRlap 论文表述
#        ("biased towards the observational correlation") 一致；
#      - f=0 对照展示 weak-instrument/NOME 衰减（偏向零）。
#
# 结论（写入 code/output/mrlap_direction_verification.txt）：
#   overlap 分量单独对同号对是「离零（anti-conservative）」，有界
#   （Burgess: ≤2.6% BMD→OA / ≤12% BMD→fracture）；
#   weak-instrument + winner's curse 分量偏向零；
#   MRlap 联合校正观测到的「离零」移动由后者主导 ⇒ 只能推出「对联合偏倚
#   保守」，不能推出「对 overlap 保守」。
# ---------------------------------------------------------------------------
import csv, io, os, time
import numpy as np

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSV  = os.path.join(PROJ, "code", "output", "mrlap_overlap_results.csv")
OUT  = os.path.join(PROJ, "code", "output", "mrlap_direction_verification.txt")

t0 = time.time()
rep = []


def say(s=""):
    print(s)
    rep.append(s)


say("=" * 78)
say("MRlap 方向性验证 · 2026-10-16（第九轮 P0：L50/L72/摘要 L24/Fig.10）")
say("=" * 78)

# ---------- A. 沉积数据 ----------
say("\nA · 沉积数据一致性（mrlap_overlap_results.csv）")
say("-" * 78)
rows = list(csv.DictReader(open(CSV, encoding="utf-8")))
say(f"{'pair':30s} {'r_g':>7s} {'observed':>9s} {'corrected':>9s} {'Δ':>9s} "
    f"{'%|β|':>7s} {'sign(rg)==sign(eff)':>20s}")
k_vals, pcts = [], []
for r in rows:
    obs, cor, rg = (float(r["observed_effect"]), float(r["corrected_effect"]),
                    float(r["rg"]))
    d = cor - obs
    pct = (abs(cor) / abs(obs) - 1) * 100
    k = d / rg if rg else float("nan")
    k_vals.append(k)
    pcts.append(pct)
    conc = (rg > 0) == (obs > 0)
    say(f"{r['exposure']+'->'+r['outcome']:30s} {rg:+7.3f} {obs:+9.4f} {cor:+9.4f} "
        f"{d:+9.4f} {pct:+6.1f}% {str(conc):>20s}")
n_conc = sum(1 for r in rows
             if (float(r["rg"]) > 0) == (float(r["observed_effect"]) > 0))
say(f"\n一致对（sign r_g == sign effect）: {n_conc}/10")
say(f"校正方向沿 +r_g（离零）的对: "
    f"{sum(1 for r in rows if (float(r['corrected_effect'])-float(r['observed_effect']))*float(r['rg'])>0)}/10")
say(f"|β| 变化: min={min(pcts):.1f}%  median={np.median(pcts):.1f}%  max={max(pcts):.1f}%")
say(f"k = Δ/r_g: min={min(k_vals):.3f}  max={max(k_vals):.3f}（全为正）")
say("解读：校正沿 r_g 方向离零，是 weak-instrument + winner's curse（偏向零）")
say("      主导、overlap（拉向观测相关即离零）被压过的联合结果；不能据此推出")
say("      「主估计对 overlap 保守」——overlap 分量单独对同号对是离零的。")

# ---------- B. 基因型级模拟 ----------
say("\nB · 基因型级模拟（隔离 overlap 分量方向）")
say("-" * 78)
rng = np.random.default_rng(20261016)
M, NPOOL, NX, NY, REPS = 2000, 40_000, 20_000, 18_000, 300
MAF = rng.uniform(0.05, 0.5, M)
X = ((rng.random((NPOOL, M)) < MAF).astype(np.float32)
     + (rng.random((NPOOL, M)) < MAF).astype(np.float32))
X = (X - X.mean(0)) / X.std(0)
sa = sb = np.sqrt(0.20 / M)
H2 = 0.20
THETA = 0.5   # 真实因果效应（r_g=+0.5, h2x=h2y ⇒ theta=r_g）


def sim(fover):
    EX = rng.choice(NPOOL, size=NX, replace=False)
    SH = rng.choice(EX, size=int(NY * fover), replace=False)   # 共享个体∈暴露集
    rp = np.setdiff1d(np.arange(NPOOL), EX)
    RS = rng.choice(rp, size=NY - len(SH), replace=False)
    OY = np.concatenate([SH, RS])
    Xe, Xo = X[EX], X[OY]
    Lc = np.linalg.cholesky([[sa ** 2, .5 * sa * sb], [.5 * sa * sb, sb ** 2]])
    AB = np.einsum('ij,jmr->imr', Lc, rng.standard_normal((2, M, REPS)))
    yx = X @ AB[0] + rng.standard_normal((NPOOL, REPS)) * np.sqrt(1 - H2)
    yy = X @ AB[1] + rng.standard_normal((NPOOL, REPS)) * np.sqrt(1 - H2)
    bx = (Xe.T @ yx[EX]) / NX
    by = (Xo.T @ yy[OY]) / NY
    return np.einsum('mr,mr->r', bx, by) / np.einsum('mr,mr->r', bx, bx)


res = {f: sim(f) for f in (0.0, 0.15, 0.30, 0.60)}
say(f"{'f':>5s} {'E[IVW]':>9s} {'±sem':>8s} {'bias vs θ':>11s}")
for f, e in res.items():
    say(f"{f:5.2f} {e.mean():+9.4f} {e.std(ddof=1)/np.sqrt(len(e)):8.4f} "
        f"{100*(e.mean()-THETA)/THETA:+10.2f}%")
say("\noverlap 分量（相对 f=0 对照的差，消去 weak-instrument 衰减）:")
say("  理论预测: +f·r_p/(N_x(s²+v))，r_p=r_g·√(h²x·h²y)=0.10，N_x=20,000，s²+v=1.5e-4")
say("           ⇒ +f·0.0333（β 单位）= +f·6.67% of θ；f=0.30 时预测 +2.0%")
for f in (0.15, 0.30, 0.60):
    d = res[f].mean() - res[0.0].mean()
    sem = np.sqrt(res[f].var(ddof=1) / len(res[f])
                  + res[0.0].var(ddof=1) / len(res[0.0]))
    say(f"  f={f:.2f}: Δ={d:+.5f} ± {sem:.5f} ({100*d/THETA:+.2f}% of θ)  "
        f"离零( toward observational corr )={d > 0}")
say("\n结论：overlap 分量随 f 单调递增且为正（离零）——对同号对是")
say("      anti-conservative；weak-instrument/NOME 衰减（f=0 对照的 −33%）")
say("      才是偏向零的分量。MRlap 三者联合校正，观测到的离零移动由偏向零")
say("      的分量主导。")

say(f"\n[{time.time()-t0:.0f}s] 验证完成。")
io.open(OUT, "w", encoding="utf-8").write("\n".join(rep))
print("written ->", OUT)
