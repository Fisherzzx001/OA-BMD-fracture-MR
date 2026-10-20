# 106 · 数据级 trust-but-verify 全量终审报告

**审计对象**：`docs/08_投稿稿_OA_BMD_fracture.md`（OAC 版，284 行）+ 全部底层产出
**审计日期**：2026-10-16
**审计原则**：不复用项目既有结论；不信任任何"我检查过了"；每条声称必须回答
**「对应哪个文件的第几行 / 用什么独立手段验证过 / 验证结论是什么」**。
**独立脚本**：`code/100_audit_independent_recompute.py`、`code/101_audit_variants.py`
（IVW、Cochran Q、I²、MR-Egger、加权中位数、LOO、power 全部从零实现，未 import 任何项目代码）

---

## 0. 结论速览

| 等级 | 数量 | 摘要 |
|---|---|---|
| 🔴 **致命/必须在投稿前修** | 2 | 尺度校正叙事与真实数据矛盾；表 S10 衰减列与表 S9 冲突且违反自身列定义 |
| 🟠 **重要/强烈建议修** | 4 | 门禁自身 bug 污染 4 条 WARN 中的 3 条；随机效应约定未声明且结论依赖它；加权中位数估计器未声明；工具数叙事三处不一致 |
| 🟡 **次要/说明性** | 4 | GEFOS 差距用比值尺度弱化；h² 的 N 精度；MRlap/PRESSO 依赖原 flags；n_h=0 退化基线未解释 |

**核心结论**：**稿件里报告的全部 β / SE / p (固定效应) / MR-Egger 斜率 / I² / LOO 翻转 / 功效 / OR 与 CI，全部被独立重算精确复现，数字是对的。**
**真正的问题不在数字本身，而在三处"数字对、文件错 / 文本错 / 门禁错"：**

1. 正文宣称的 `×0.2823` 校正**从未作用于报告值**（报告值 = 未经乘算的原始 IVW），而且**如果按正文照做，阳性对照会崩**；
2. 补充表 **S10 的衰减列算错/换错基线**，与 S9 在同一配置上直接冲突；
3. **门禁脚本本身有 bug**，其 4 条 WARN 里有 3 条是 bug 产物，反而掩盖了一个真实问题。

---

## 1. 逐声称溯源表（claim → 证据位置 → 独立验证手段 → 结论）

### 1.1 Table 3 · MR 主表（全部数值）

**数据源**：`code/output/genetics/harmonised/<exposure>__<outcome>.csv`（SNP 级，`beta_e/se_e/beta_o/se_o`）
**验证手段**：重算 Wald 比值 β_i=β_o/β_e、SE_i=|SE_o/β_e|，再做逆方差加权；Eg​ger 用加权回归 + 离散度校正；I² 用 (Q−df)/Q；功效用 Φ(|z|−1.96)。

| Table 3 断言 | 独立重算值 | 结论 |
|---|---|---|
| 全部 10 对 **n_snp** (31/25/31/21/16/24/25/69/69/69) | 与 harmonised 文件行数逐一相等 | ✅ 一致 |
| 全部 10 对 **IVW β** | 误差全部 ≤0.0004（纯属印刷取整） | ✅ 精确复现 |
| 全部 10 对 **固定效应 p** | 全部复现（如 3.86e-10 vs 3.9e-10） | ✅ 精确复现 |
| **MR-Egger 斜率** 10 个 | 全部复现（如 −0.34258 vs −0.343） | ✅ 精确复现 |
| **Egger 截距 p** 10 个 | 接近但不完全相等（0.5626 vs 0.558 等，差 0.001–0.01） | 🟡 离散度处理细节差异，**不改变任何结论**（0.018 与 0.027 都显著） |
| **I²** 10 个 | 全部复现（78.3/83.3/80.9/82.8/85.7/84.7/68.4/48.6/48.5/30.0） | ✅ 精确复现 |
| **LOO 符号翻转** 全为 0 | 全为 0 | ✅ 精确复现 |
| **OR / CI**：1.40(1.34–1.46)、1.42(1.34–1.51) | 1.4009(1.341–1.464)、1.4200(1.338–1.507) | ✅ 精确复现 |

**数据完整性旁证**：文件自带的 `ratio` 列与我重算的 β_o/β_e 偏差 ≤3e-14 → 文件内部自洽，未被篡改。

### 1.2 Table 4 · 功效表

| 断言 | 独立重算 | 结论 |
|---|---|---|
| \|z\| = 15.05 / 11.55 / 6.26 / 5.77 / 4.78 / 4.23 / 3.42 / 0.85 / 0.98 / 1.31 | 15.046 / 11.546 / 6.260 / 5.775 / 4.785 / 4.234 / 3.416 / 0.851 / 0.981 / 1.312 | ✅ 全部复现 |
| MDE = 2.8×SE | 全部复现 | ✅ |
| 达成功效 = Φ(\|z\|−1.96) = 0.14 / 0.17 / 0.26 | 0.134 / 0.164 / 0.258 | ✅ 复现（公式已反推确认） |

### 1.3  liabilities scale h²（Table 1 脚注）

用 Lee et al. 2011 公式 **h²_liab = h²_obs × K(1−K)/z²**（z=密度在阈值处）独立验算：

| trait | K | 换算系数 K(1−K)/z² | h²_obs | 算得 h²_liab | 稿件 | 结论 |
|---|---|---|---|---|---|---|
| OA_knee | 0.12 | 2.6373 | 0.0354 | 0.0934 | 0.094 | ✅ |
| OA_selfrep | 0.083 | 3.2367 | 0.0223 | 0.0722 | 0.073 | ✅ |
| fracture | 0.143 | 2.4009 | 0.0258 | 0.0619 | 0.062 | ✅ |

---

## 2. 🔴 致命项

### 🔴 F1 · `×0.2823` 尺度校正文案与真实数据矛盾（会把审者引向"你算错了 3.54 倍"）

**稿件位置**：Methods「Effect-scale reconstruction」（第 49 行）、Limitations（第 85 行）

> 第 49 行："Munged BMD statistics carried a placeholder N=400,000 … shrinking reconstructed β 3.54-fold; forward MR β **were corrected by ×0.2823** (÷3.5426)"
> 第 85 行："absolute MR magnitudes are **provisional** because exposure β were reconstructed from a placeholder N=400,000 … and scaled down 3.54-fold (×0.2823)"

**我的验证**（`code/101_audit_variants.py` 的 A 段与 F 段）：

**(a) 报告值到底乘没乘 0.2823？** —— 没有。
全部 10 对：|报告值 − 原始 IVW| ≤ 0.0004；|报告值 − 原始 IVW×0.2823| ≥ 0.012。
→ **Table 3 报告的就是未经任何后乘的原始 IVW。**

**(b) 那原始 IVW 本身对不对？** —— 对，而且已经用真 N 标准化。
对每个 SNP 用闭式解核对：标准化效应的标准误应满足 **SE = 1/√(2·N·f·(1−f))**（f=MAF）。

```
BMD_total__OA_knee.csv  rs12568930  se_e=0.010303   N=31,873 → 理论 0.010303  ratio=1.0000 ✅
                                                    N=400,000 → 理论 0.002908  ratio=3.5426
(rs2566755, rs4846607 等全部同样：对 N=31,873 精确等于 1.0000)
```

→ **文件里的暴露 β 是用真 N=31,873 标准化的，丝毫没有"被占位 N 缩小 3.54 倍"。**
（若真如正文所说用了 N=400,000，观测 SE 应当是 1/3.5426 ≈ 0.2823 倍，实测否。）

**(c) 如果真的照正文乘 0.2823 会怎样？** —— 阳性对照当场崩掉。

| | 不乘（=报告值） | 乘 0.2823 |
|---|---|---|
| BMD total → fracture β | −0.337 | −0.0952 |
| 对应 OR（每低 1 SD） | **1.40** | **1.10** |
| 对照 GEFOS 1.55 | 比值尺度低 5.3% | 比值尺度低 **29.2%** |

即：**正文描述的修正一旦落到数字上，就会把论文赖以自证的阳性对照从"接近文献（1.40 vs 1.55）"变成"荒谬地弱（1.10）"。**

**真实状态**：这份 munge 早就在上游按真 N 重标化了，因此×0.2823 不需要再乘（乘了才是错）。**数字是对的，文案是错的。**

**危害**：审者只要按第 49 行照做一次乘法，就会判定"全部效应量高估 3.54 倍"——这正是审稿人最容易动手验算的一步。同时 Limitations（第 85 行）还主动告诉对方"这些绝对值只是暂定值"，把一个根本不存在的不确定性写成了自我否定。

**修复方案**（改文案，不动任何数字）：
- 第 49 行改为陈述已核查的事实而非一个待施加的修正，例如：
  "BMD betas were re-standardised to per-SD units using the verified site N (31,873). We confirmed the reconstruction empirically: per-SNP standard errors satisfy SE = 1/√(2·N·f(1−f)) with N=31,873 to within 0.1% at every instrument, and are inconsistent with the placeholder N=400,000 (which would predict SE 3.54-fold smaller). **No post-hoc rescaling was therefore applied to the reported estimates.**"
- 第 85 行删除"magnitudes are provisional because … ×0.2823"，改为：
  "absolute magnitudes depend on the areal-BMD standardisation, which we verified against the external GEFOS positive control (see §Effect-scale validation)."

---

### 🔴 F2 · 表 S9 与表 S10 在同一配置上互相矛盾，且 S10 违反自己的列定义

**文件**：`docs/08` 的 Table S9（第 185–198 行）vs `docs/Supplementary/02_mvmr_height_grid.md`（Table S10）

S10 自述定义：`Atten% vs Table 3 = 100×(β_uni(Table3) − β_MVMR)/β_uni(Table3)`

同一配置（EUR, n_h=100, primary）：

| 表 | β_uni(Table3) | β_MVMR | 声称衰减 | **按 S10 自身公式重算** |
|---|---|---|---|---|
| **S9**（第 191 行）BMD total→OA knee | +0.0231 | +0.0187 | **19.1** | 19.07 ✅ |
| **S10**（第 8 行）BMD total→OA knee | 0.02312 | 0.01871 | **9.8** | **19.07 ❌** |
| **S9**（第 196 行）BMD total→fracture | −0.3371 | −0.3427 | **−1.7** | −1.66 ✅ |
| **S10**（第 48 行）BMD total→fracture | −0.33713 | −0.34274 | **−2.8** | **−1.66 ❌** |

S10 里写的 9.8 其实是**拿 β_uni(MVMR set)=0.02074 当分母**算出来的（(0.02074−0.01871)/0.02074 = 9.79）；
S10 里写的 −2.8 同样是拿 β_uni(MVMR set)=−0.33339 当分母算出来的（= −2.80）。
→ **S10 的 `Atten% vs Table 3` 列实际填的是"vs MVMR-set 基线"的值，与其表头定义不符。**

**顺带暴露的选择性呈现问题**：S10 网格里存在**负衰减**（BMD total→OA knee，EUR n_h=0 时 **−9.0**；EUR_noUKB n_h=0 时 **−9.2**），即加height变量后效应反而变大 9%。而正文 headline「6–33%」只取正值区间，未披露全网格范围。

**修复**：重算 S10 两列并按 .rows 对齐 S9；在该表脚下注明"negative attenuation indicates the estimate moved *away* from the null after conditioning"，并在正文说明 headline 6–33% 的取值规则。
（**好消息**：正文 headline 的 13.5–33.2% 取自 S9，而 S9 是算对的，所以 headline 本身不受影响。）

---

## 3. 🟠 重要项

### 🟠 F3 · 门禁脚本自身有 bug —— 正好是你要求查的"门禁比稿件更错"

**门禁**：`audit_v18_full.py`（自称 CHECKS=430 / FAILS=0 / WARNS=4）
**运行现状**：我重跑了一次，仍是 430/0/4 —— 但这 4 条 WARN 里 **3 条是 bug 的产物**。

门禁 §3.6 在比较"移除 PRESSO 离群后 β"时，**对同源数据错误地把"移除后"的值再乘了 0.2823，而"全样本"值保持原尺度**，于是凭空造出 3 条大幅度负向变化：

| 反向 MR 对 | 移除数 | 全样本 β | **正确**移除后 β / 变化 | **门禁报的** β / 变化 |
|---|---|---|---|---|
| OA_knee → fracture | 2 | +0.0922 | **+0.2001 / +117%** | +0.2001 / +117.2% ✅ 碰巧对 |
| OA_knee → BMD_total | 1 | +0.1061 | **+0.1714 / +61.5%** | +0.0484 / **−54.4%** ❌ |
| OA_knee → BMD_femur_neck | **0** | +0.1419 | **+0.1419 / 0.0%**（没移除任何东西，不可能变） | +0.0401 / **−71.8%** ❌ |

复算脚本直接证明：第三对一个离群都没移，门禁却报出 −71.8% 的变化。**这是门禁的尺度串档 bug。**

**最讽刺的后果**：bug 把"OA_knee→BMD_total 移除 1 个离群后从 +0.106 跳到 +0.171（+61.5%）"这个**真实存在的不稳定**按反方向压了下去，导致历轮都没人正视它。
→ **建议**：修正门禁；并在稿件补一句（落在 §3.4 反向 MR 段）：移除 1–2 个 PRESSO 离群可使反向 MR 点估计移动 61%–117%（仍不显著），即反向无效不仅有功效问题，还有对工具变量选择敏感的问题。这与你已披露的"三对无效 MR 功效仅 0.14–0.26 且对工具变量选择敏感"是一致的，应把具体幅度写进去。

### 🟠 F4 · 随机效应约定未声明，且"只有 2/5 存活"依赖该约定

我试了三种标准随机效应实现，稿件用的是 **multiplicative（SE 乘 √φ，φ=max(1,Q/df)，点估计不变）**，十对全部精确复现（如 0.0062/0.0166/0.0699/0.0976/0.0547）。

但换成更常见的 **DerSimonian-Laird 加性 τ²** 后：

| 对 | 稿件(multiplicative) | DL | 是否改变结论 |
|---|---|---|---|
| **Arms → OA knee** | **0.0699（不显著）** | **0.0333（显著）** | **❗直接改变 2/5 → 3/5** |
| Legs → OA knee | 0.0976 | 0.0735 | 否 |
| Total → OA knee | 0.0062 | 0.0048 | 否 |

→ 摘要/结论反复强调的"随机效应下 5 对中仅 2 对仍显著"，**换一种同样标准的实现就变成 3 对**。这不是错误（multiplicative 合法），但必须声明，否则审稿人用 TwoSampleMR 默认设置复核会与你不吻合并被指 selective convention。

**修复**：Methods 补一句 "random-effects IVW was implemented with multiplicative heterogeneity inflation of the standard error (SE_random = SE_fixed × √max(1, Q/df)); we note that a DerSimonian–Laird additive τ² implementation additionally retains significance for the arms BMD exposure (p=0.033)".

### 🟠 F5 · 加权中位数估计器未声明 → 该列不可由论文复现

我用三种标准定义（step 加权中位 / 核平滑 / 简单中位）都无法一致复现 Table 3 的该列：

| 对 | 稿件 | step w=1/se² | 核平滑 | 简单中位 |
|---|---|---|---|---|
| BMD total → fracture | −0.364 | −0.3926 | −0.3398 | −0.3559 |
| Femur-neck → fracture | −0.365 | **−0.3646** ✅ | −0.3417 | −0.3646 |
| **Arms → OA knee** | **0.031** | **0.0076** | 0.0245 | 0.0377 |
| **Legs → OA knee** | **0.025** | **0.0077** | 0.0179 | 0.0162 |
| Lumbar → OA knee | 0.037 | 0.0338 | 0.0313 | 0.0338 |

arm/leg 两对在标准定义下只有 0.008 左右，稿件报 0.031/0.025（反而高于 IVW）。不同实现之间可差 4 倍。**STROBE-MR 要求指定估计器**，Methods 目前只写"weighted median"。
→ **修复**：在 Methods 写明具体估计器与权重定义（例如 "Bowden weighted median, kernel-smoothed, as implemented in <pkg/version>"），或在脚注给出重现脚本路径。

### 🟠 F6 · 工具数量叙事三处不一致

| 出处 | 声称 | 实测（`*_plink.clumped` 行数） |
|---|---|---|
| Methods 第 45 行 "retained 20–35 instruments per BMD exposure" | 20–35 | **22 / 26 / 28 / 30 / 37**（total=**37 > 35**） |
| Methods 第 45 行 "116 for knee OA" | 116 | **118** |
| Table 3 实际分析用 n_snp | 16/21/24/25/31；knee OA **69** | — |

**额外发现**：strict clump 到进入分析，**knee OA 从 118 掉到 69（流失 42%）**，这一损耗在 Methods 完全未提。STROBE-MR 要求报告工具筛选各阶段的流失。
→ **修复**：改为实际数字并补一句"harmonisation and allele-matching against the outcome reduced the knee-OA set from 118 to 69 instruments"。

---

## 4. 🟡 次要 / 说明性

- **G1 · GEFOS 差距只给比值尺度（建议同时给 log 尺度）**。"8.4% / 9.6% below GEFOS"用的是 OR 比值尺度，**这两个数是算对的**（重算 OR 1.4009/1.4200：(1.55−1.42)/1.55=8.39%、(1.55−1.4009)/1.55=9.62%，与稿件一致）。但在 log-odds 尺度（效应比较的自然尺度）上差距是 **20.0% / 23.1%**。建议两者都给，避免被指挑弱化差距的尺度。
  - ⚠️ 自查更正：审计脚本初版曾把 GEFOS 参考值误取为 1.48（下置信限）而非 1.55，导致本节首版 log 尺度数字偏低；已更正。**稿件原 8.4%/9.6% 无误**，此项属可选增强而非错误。
  - ✅ **已修**（2026-10-16）：稿件 §3.7 现已同时给出比值尺度与 log-odds 尺度差距。
- **G2 · Table 1 声称 N=31,873，实际按 N=32,000 计算 h²**（门禁自己 WARN 过），差 0.40%，属于透明度小瑕疵。
- **G3 · 离群剔除与 MRlap 依赖证件自带的 flags**：`presso_outlier` 列与 `mrlap_overlap_results.csv` 是直接采用的项目产出，我只做了"接受 flags 前提下的复算"，未独立重跑 MR-PRESSO/MRlap 选择过程。
- **G4 · S10 中 n_h=0 是退化行**：没有任何 height 工具却仍有 β_MVMR ≠ β_uni（且给出 −9.0% 衰减），概念上可疑，应在该表脚注解释，否则看起来像 bug。

---

## 5. 框架级复核：研究完成度与补强空间

**已完成且经独立验证（可采信）**
LDSC 13×13 矩阵（r_g 区块、h²、截距、观察/易患性尺度换算）· 10 对两样本 MR 全部估计量 · 阳性对照 · 反向 MR · MR-PRESSO · 留一法 · 功效计算 · MVMR 主配置（S9 正确）· coloc（含 31 位点全表与 BMD–fracture 管道对照）· height MVMR 网格。

**真正的补强空间（按性价比排序）**

1. **先修上面 2 个致命项** —— 这是"数字对但文本/补充材料错"的典型，审稿人一查就中。
2. **把"随机效应 2/5 或 3/5"的约定敏感性写进 limitation**，比被动挨问好得多。
3. **补反向 MR 的离群敏感性数值**（+61.5% / +117%），让"反向无效"的定性限定更有说服力。
4. **未验部分建议至少说清**：`Figure 4 h²`、`Figure 8 power`、`表 S8 coloc` 我这次只验了数字层，未逐图重绘比对（`Figures_v17_Nature/plotted_data.json` 存在，可用 `audit_figures_v17.py` 起 gate 补一轮）。
5. **工具变量流失 118→69** 按 STROBE-MR 补齐，这是审稿人常问的 attrition。

---

## 6. 本次审计使用的方法与证据位置（便于复核）

| 步骤 | 脚本 / 命令 | 产出 |
|---|---|---|
| SNP 级数据完整性 + 全估计量重算 | `code/100_audit_independent_recompute.py` | `code/output/genetics/AUDIT_100_independent_recompute.json` |
| 尺度判定 / RE 三实现比对 / WM 三定义比对 / OR-CI / 功效 / SE 闭式解核 N | `code/101_audit_variants.py` | 控制台 A–F 六段 |
| 工件计数核对 | `wc -l code/output/genetics/harmonised/*.csv`、`*_plink.clumped` | 见 §1.1、F6 |
| 门禁复跑 | `python code/49_oac_wordcount.py`、`python audit_v18_full.py` | 字数门禁 PASS；主审计 430/0/4，其中 3 条 WARN 系 bug |
| 反向 MR 离群敏感性独立重算 | 见 §3 F3 表 | 证明门禁 bug |
