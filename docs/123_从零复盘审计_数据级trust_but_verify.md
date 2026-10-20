# 123 · 从零复盘审计：科研框架完成度 + 数据级 trust-but-verify 终审

> 日期：2026-10-20 · 对象：`docs/08_投稿稿_OA_BMD_fracture.md`（v18 后）+ `docs/118_OSF预注册` + `code/` + `code/output/genetics/` + 两个投稿包
> 审计原则：**不信任任何进程输出与任何既有结论**。每条声称必须给出「稿件行号 → 源脚本行号 → 独立验证手段 → 验证结果」。
> 特别警戒：*"数字与文件一致" ≠ "数字对"* —— 文件本身可能是错的。

---

## 0. 审计方法与信任边界

| 手段 | 说明 |
|---|---|
| A. 独立重算脚本 | `code/170_from_scratch_recompute.py`：**不 import 项目任何 MR 模块、不读任何既有结果 CSV**，只从 `code/output/genetics/harmonised/*.csv` 的 SNP 级原始 β_e/β_o/se_e/se_o 重新推导 |
| B. 原始数据对撞 | 拿第三方真实 GWAS 文件（`gcst_knee_union.csv.gz`）与 harmonised 重构值逐 SNP 对撞 |
| C. 闭式反推 | 由 `se = 1/√(2·N·f(1−f))` 反推实际使用的 N，判定效应所处尺度 |
| D. 精确分布 | Wilcoxon 符号秩用 DP 枚举零分布（非正态近似）；χ² 用连分数展开 |
| E. 目视读图 | Figure 3 PNG 直接读取数值标签 |

---

## 1. 科研框架逐级完成度

### 1.1 预注册假设 P1–P6（`docs/118` L40–L45）

| 项 | 内容 | 状态 | 证据 |
|---|---|---|---|
| P1 | 对立终点 | **支持** | 阳性对照 OR 1.40（真 log-OR，独立复现）；BMD→膝 OA 5/5 正（z 全部复现） |
| P2 | 系统性而非局部力学 | **系统条款支持；排序条款未达显著** | T3a 独立手算 z=0.7461，单侧 p=**0.2278**（稿件/门禁称 0.228 ✓）。⚠ 旧值 0.18 已于上一轮修正 |
| P3 | 弥漫多基因 | **支持** | 0/31 位点 PP.H4≥0.8，独立重算 ✓ |
| P4 | 尺寸是可分离的组分 | **支持（部分）** | 身高衰减 9.5–33.2%，✓；BMI 衰减 −3.7%–46.0% ✓（上一轮已补齐） |
| P5 | 不对称是体质性的 | **反向无效** | β 非显著，power 0.14–0.26 ✓（power 仅依赖 z，尺度不变，故不受 §2 的 G6 影响） |
| P6 | 部位通路主导 | **降级为量级主张** | 髋证据内部矛盾，稿件已如实降级 ✓ |

### 1.2 证伪门禁 F1–F10（`docs/118` L51–L60）

| 项 | 状态 | 说明 / 补强空间 |
|---|---|---|
| F1 髋 OA 十部位 | ❌ 未执行 | 需外部数据。**补强优先级：中**（P6 目前只能是量级主张） |
| F2 质量/尺寸分解（vBMD、髋几何、骨面积） | ❌ 未执行 | **这是 P4 唯一能真正定案的检验，也是 R1 指出的"理论实体未被识别"症结。补强优先级：最高** |
| F3 LDSC-SEG 组织富集 | ⚠ 部分 | `code/164_ldsc_seg_native.py` 存在；产物需核实 |
| F4 组织/骨侧共定位 | ⚠ 部分执行 | cross-trait 对照臂 + eQTLGen 全血臂已完成（31 位点，精确 Wilcoxon p=**0.001156** ✓）；**骨组织/软骨组织字面臂仍不可执行**（无可用组织 eQTL） |
| F5 骨赘/关节间隙分层 | ❌ 未执行 | 预注册为"下一步"，需分层 OA GWAS |
| F6 跨部位 Q（共享仪器协方差精确处理） | ✅ 执行 | 稿件 Q=6.92、df=4、p=0.14、τ²=0 |
| F7 关节位点排除 | ✅ 执行 | GDF5 + leave-WNT16/SOX6-out 通过 |
| F8 性别分层 | ✅ 执行（⚠ 见 G7） | Boer 2021 性别分层；两性均正。**但尺度与主分析不可比** |
| F9 部位梯度正式检验 | ✅ 执行 | T3a p=0.2278 / T3c Q=3.56、df=4、p=0.468、I²=0 ✓ 独立复现 |
| F10 祖源敏感性（EAS） | ❌ 未执行 | Hatzikotoulas 2025 多祖源 OA 可获取。**补强优先级：高**（直接回应 R2/审稿人对祖源错配的质疑） |

### 1.3 补强空间（按性价比排序）

1. **F2（质量/尺寸分解）** —— 唯一能把"67–91% 身高独立残差"定案为骨量还是几何的检验；也是 P4 与 SOS 理论实体的关键。缺它，核心理论主张只能是"与数据一致的假说"。
2. **F10（EAS 祖源）** —— 直接把"祖源错配"从局限降级为已检验的稳健性。数据公开可得。
3. **M1（5 个 BMD 位点 clumping）** —— 1000G 面板下载中（290 MB/1.5 GB）。完成后十位点评即完整，"反映算力而非生物学"的自我暴露消失。
4. **F4 组织字面臂** —— 若能获取骨/软骨 eQTL 则完成。

---

## 2. 🔴 重大新发现 G6：OA 结局效应尺度错误（此前所有轮次均未发现）

### 2.1 声称

- 稿件 **L57**：*"The deposited BMD summary statistics are standardised to per-SD units at the verified per-site N (unified N=31,873)… Because the same verified per-SD units enter forward and reverse analyses, all estimates share one N-verified SD scale and directions, p-values and significance are invariant to it."*
- 稿件 **L69–L70 / Table 2**：BMD→膝 OA **OR 1.023（1.016–1.031）** 等五个 OR。
- 稿件 **L148 脚注 §§**：*"Fracture N adopted as 264,973… identified by the **measured effective sample size N·K·(1−K)=31,605**"* —— 项目**知道**并**对骨折**应用了 N_eff 修正。

### 2.2 源脚本行

- `code/13_mr_analysis.py` **L40**（`load_munged`）：
  ```python
  se = 1.0 / np.sqrt(df.N.values * 2 * df.MAF.values * (1 - df.MAF.values))
  df = df.assign(beta=df.Z.values * se, se=se)
  ```
  **N 直接取 munged 文件的 N 列，无任何 case-control（N·K(1−K)）修正。**
- 全文件 grep `N_eff|Neff|ncase|n_case|K(1-K)|prevalence` → **0 命中**（L13 L47–L63 的 `load_fracture_df` 用真实 β/SE 绕过，是唯一例外）。

### 2.3 独立验证（三重，互相独立）

**验证一：闭式反推 implied N（手段 C）**

| harmonised 对 | se_o 反推 implied N | 判定 |
|---|---|---|
| BMD_total→**fracture** | **31,907** | = N·K(1−K)（264,973×0.143×0.857=32,505）→ **真 log-OR ✓** |
| BMD_femur_neck→fracture | 31,559 | 同上 ✓ |
| OA_knee→**fracture** | 31,641 | ✓ |
| BMD_total→**OA_knee** | **1,269,775** | ≈ N_total（1,150,161–1,208,397）→ **线性重构，未做 case-control 修正 ✗** |
| BMD_lumbar→OA_knee | 1,269,775 | ✗ |
| BMD_arms→OA_knee | 1,313,371 | ✗ |
| BMD_leg→OA_knee | 1,270,458 | ✗ |
| BMD_femur_neck→OA_knee | 1,249,796 | ✗ |

→ **同一条流水线对两个结局用了不一致的尺度处理。**

**验证二：真实 GWAS 逐 SNP 对撞（手段 B）**
源：`C:\...\OA_FLAGSHIP\02_data\11_bmd_snp_cache\gcst_knee_union.csv.gz`（含 OA GWAS 真实 `beta_o/se_o/eaf_o`），按 CHR:BP 与 `harmonised/BMD_total__OA_knee.csv` 合并，**13 个 SNP 全部匹配**：

| SNP | β_o（harmonised 重构） | β_o（真实文件） | 比值 | se 比值 | 真实 se 反推 N_eff |
|---|---|---|---|---|---|
| rs6679677 | −0.000940 | −0.0029 | 3.086 | 3.086 | 85,852 |
| rs777355 | −0.003764 | −0.0125 | 3.321 | 3.321 | 119,001 |
| rs10931982 | −0.002009 | −0.0073 | 3.633 | 3.633 | 92,690 |
| rs4974056 | +0.003628 | +0.0115 | 3.169 | 3.169 | 90,666 |
| rs17150506 | +0.001166 | +0.0040 | 3.430 | 3.430 | 111,698 |
| rs3020333 | +0.002479 | −0.0084 | −3.389 | 3.389 | 114,169 |
| rs3779381 | +0.001337 | −0.0045 | −3.366 | 3.366 | 115,993 |
| rs7010267 | +0.002955 | +0.0100 | 3.385 | 3.385 | 114,603 |
| rs7938055 | −0.000880 | −0.0030 | 3.408 | 3.408 | 113,194 |
| rs876701 | +0.000950 | −0.0031 | −3.263 | 3.263 | 89,788 |
| rs9515218 | −0.000957 | −0.0034 | 3.553 | 3.553 | 100,373 |
| rs884205 | −0.004562 | −0.0160 | 3.507 | 3.507 | 103,696 |
| rs7277076 | +0.001146 | +0.0040 | 3.491 | 3.491 | 104,466 |

- **|β 比值| 恒等于 se 比值**（3.086/3.086、3.321/3.321…）→ 两者 Z 完全相同，**纯尺度差异**，排除"不同 GWAS"的解释。
- 真实 se 反推 N_eff 中位 **104,466**，而 N_total≈1,150,161 → K(1−K)=0.0908 → **真实文件是 log-OR 尺度**，重构值是线性尺度。
- 中位校正因子 **3.385**（β、SE 同乘）。
- 附带：4 个 SNP 的 β 比值为**负号**（rs3020333 / rs3779381 / rs876701）→ **等位基因方向翻转**，harmonisation 的 flip 判定在 EAF 差异较大时可能出错（次要问题 G12）。

**验证三：稿件自身 K 值交叉校验**
稿件 L148 声明 OA_knee **K=0.12** → K(1−K)=0.1056 → 因子 **3.0773**。
实测因子 3.385 / 3.339 → 对应 K≈0.0996。两者同量级（3.08–3.39），结论一致。

### 2.4 后果（解析推导，已数值验证）

对 IVW：β = Σ(β_e·β_o/σ_o²) / Σ(β_e²/σ_o²)，SE = 1/√Σ(β_e²/σ_o²)。

- **正向（BMD→OA，β_e 正确、β_o 与 σ_o 同缩 c=3.08–3.39）**：
  β_MR_recon = β_MR_true / c，SE 同理 → **z 与 p 完全不变**，但 **β 与 OR 被低估 c 倍**。
- **反向（OA→X，β_e 缩 c、σ_o 正确）**：
  β_MR_recon = c × β_MR_true，SE 也 × c → **z 与 p 完全不变**，但 **β 被高估 c 倍**。

**真实值（以稿件自身 K=0.12，因子 3.0773 计；括号内为实测因子 3.318 的结果）：**

| 对 | 稿件报称 OR | **真实 OR** | p（不变） |
|---|---|---|---|
| BMD_total→OA_knee | 1.023 (1.016–1.031) | **1.074 (1.050–1.098)**［1.080 (1.054–1.106)］ | 3.9e-10 ✓ |
| BMD_lumbar→OA_knee | 1.028 (1.019–1.038) | **1.090 (1.057–1.123)**［1.097］ | 7.7e-09 ✓ |
| BMD_arms→OA_knee | 1.022 (1.013–1.032) | **1.071 (1.039–1.103)**［1.076］ | 1.7e-06 ✓ |
| BMD_leg→OA_knee | 1.018 (1.010–1.026) | **1.056 (1.031–1.082)**［1.061］ | 2.3e-05 ✓ |
| BMD_femur_neck→OA_knee | 1.017 (1.007–1.027) | **1.054 (1.023–1.086)**［1.058］ | 6.4e-04 ✓ |

| 反向对 | 稿件报称 β | **真实 β** | p（不变） |
|---|---|---|---|
| OA_knee→fracture | +0.092 | **+0.030** | 0.395 ✓ |
| OA_knee→BMD_total | +0.106 | **+0.035** | 0.326 ✓ |
| OA_knee→BMD_femur_neck | +0.142 | **+0.046** | 0.190 ✓ |

| T3a 部位梯度 | 稿件报称 | **真实** | p（不变） |
|---|---|---|---|
| arms − femoral neck | +0.00509 (se 0.00682) | **+0.01566 (se 0.02099)** | 单侧 0.2278 ✓ |

### 2.5 受影响 / 不受影响清单

| 量 | 是否受影响 | 理由 |
|---|---|---|
| **OR / β（涉 OA 者）** | ❌ **错** | 正向低估 3.08–3.39×，反向高估同倍 |
| p 值、z、显著性、方向 | ✅ 正确 | 分子分母同缩放，z 不变 |
| Cochran Q / I² / τ² | ✅ 正确 | 基于 z 型比值 |
| MR-Egger 斜率/截距及 p | ✅ 正确 | 同上 |
| 加权中位数、PRESSO、LOO | ✅ 正确 | 同上 |
| LDSC r_g、h²、intercept | ✅ 正确 | 不经此重构 |
| coloc PP.H4 / PP.H3 | ✅ 正确 | 只用 Z 与 MAF |
| MVMR 衰减百分比、条件 F | ✅ 正确 | 比值，尺度不变 |
| MRlap 校正幅度（+20.5%） | ✅ 正确 | 比值 |
| 统计功效（0.14–0.26） | ✅ 正确 | 只依赖 z |
| **BMD→骨折阳性对照 OR 1.40/1.42** | ✅ 正确 | fracture 走 `load_fracture_df` 真实 log-OR；且 GEFOS 1.55 外部锚定吻合 |

### 2.6 判读

- 这是**"数字与文件一致，但文件本身是错的"的教科书案例**：harmonised CSV 内部完全自洽（β/SE/z/p 互推全通），独立重算 10/10 复现——但它的 OA 结局效应本身处在线性尺度，而稿件对它取 `exp()` 当 OR 报。
- **稿件 L57 的那句话（"all estimates share one N-verified SD scale"）是错的**：只有 BMD 暴露与 fracture 结局在该尺度上；OA 结局不在。
- **方向性影响：修正后效应更大（OR 1.07–1.09 而非 1.02），对论文有利**。但报 1.023 是事实错误，任何复核 OA GWAS 的审稿人都会抓到。
- **必须修**：5 个 forward OR + 3 个 reverse β + T3a 差值 + 摘要/Highlights 中的 1.023。

---

## 3. 🔴 G7：F8 性别分层数值与主分析不可比（上一轮由本次审计的前一轮引入）

- **稿件 L49 / L99**（上一轮我写入）：*"…confirmed BMD→OA positivity in both sexes (KneeOA Female β=+0.32, Male β=+0.45; AllOA Female β=+0.34, Male β=+0.23)…"*，紧邻主分析 β=+0.023。
- **源脚本**：`code/161_tier1_sex_mr.py`
  - **L46**：`beta = np.log(es.values)` —— Boer 的 `EffectSize` 是 **OR**，已正确取对数 → **结局是真 log-OR ✓**
  - **L96**：`inst = m13.select_instruments(m13.load_munged(site))` —— 暴露走 `load_munged`，**N = 400,000 占位符，未做 N=31,873 校正 ✗**
- **后果**：F8 的 β 相对"已校正主分析"被放大 √(400,000/31,873)=**3.54×**（暴露）× **3.08–3.39×**（结局 log-OR vs 主分析线性）＝ **约 10.9–12.0×**。
- **稿件呈现的"主分析 β=0.023 与 F8 β=0.32"相差 14 倍，其中绝大部分是尺度伪影，不是性别修饰。**
- **必须修**：F8 数值要么按统一尺度重算（重跑 `code/161` 时把暴露 N 设为 31,873），要么在稿中只报"两性方向一致为正"而**不给与主分析并列的 β**。

---

## 4. 其他门禁 / 稿件问题

| 编号 | 问题 | 稿件行 | 源脚本/文件行 | 独立验证 | 级 |
|---|---|---|---|---|---|
| G6 | OA 结局效应尺度错误（详见 §2） | L57/L69–70/Table 2/L148 | `code/13` L40 | 三重独立验证 | **P0** |
| G7 | F8 β 与主分析不可比 | L49/L99 | `code/161` L46/L96 | 解析推导 + 源数据核对 | **P0** |
| G8 | `code/09` 的 `PHENOS` 只有 **12** 个表型（缺 fracture），产出 `rg_matrix.csv` 为 12×12；稿件的 13×13 / 78 对数值来自 `rg_matrix_13.csv`。**重跑 `code/09` 复现不出稿件** | L53/L63 | `code/09_rg_matrix.py` PHENOS | 读文件确认 12 项 vs 13 项 | P1 |
| G9 | coloc 配对 Wilcoxon **p=0.0043 是单侧值**；精确双侧（DP 零分布）= **0.008589**（W⁺=88、W⁻=408、n=31）。稿件未标注方向 | L79 | `genetic_validation_coloc_contrast.csv` | 精确 DP 重算 | P1 |
| G10 | **Egger 截距 p=0.027 归属错误**：稿件 L81 说它在股骨颈位点，实测属 **Lumbar**（股骨颈 p=0.355）。此错由上一轮 M2 修改引入 | L81 | `mvmr_bmi_main.csv` cfg=100 行 | 逐行比对 | P1 |
| G11 | 骨折对照衰减值：稿件 **−2.9%** vs 门禁 **−3.0%**；且 −2.9% 是用 `atten_vs_bmdsubset` 分母算的，而 9.5–33.2% 用的是 `atten_vs_full`，**两个基线分母不一致**。上一轮改动未同步门禁 | L65/L81/L95/L201 | `mvmr_height_grid.csv` | 逐行复算 | P1 |
| G12 | harmonisation 在 EAF 差异大的位点可能 flip 错误（13 个对撞 SNP 中 4 个 β 符号相反） | Table 2 | `code/13` harmonise | 真实文件对撞 | P2 |
| G13 | 稿件随机效应用 `se_RE = se_fixed × √(Q/df)`（乘性约定），**不是 DerSimonian–Laird 加性 τ²**；10/10 完全复现，但 Methods 未声明约定 | Table 2 注 | — | 三种约定逐一比对，B 约定 10/10 命中 | P2 |
| G14 | `OA_selfrep` 的 munged N=**350,000 占位符**，而稿件 L49 声明 N=462,933 —— 与 BMD 的 400,000 占位符同类问题（仅影响 h²，r_g 为 N 不变量） | L49 | `rg_h2.csv` / `N_sensitivity_decision.csv` | 读文件 | P2 |

---

## 5. 已独立验证为**正确**的声称（不得再动）

| 声称 | 稿件行 | 验证手段 | 结果 |
|---|---|---|---|
| 5 对 BMD→OA + 2 对 BMD→骨折 + 3 对反向的 **n_SNP、IVW β/SE/z/p、Egger β/SE/p/截距 p、加权中位数、I²、LOO 无翻转** | Table 2 | `code/170` 从 harmonised 重算，逐字段比对 | **10/10 全部字段通过** |
| BMD_total→OA_knee β=0.02312、OR 1.023、CI 1.016–1.031、p=3.86e-10 | L69 | 独立 IVW + 闭式反推 | ✓（但 OR 的**尺度**错，见 G6） |
| 阳性对照 OR 1.40/1.42，GEFOS 比较"低 8.4%" | L77 | 独立重算 + 外部锚定 | ✓ |
| 78 对 r_g：76 BH-FDR、75 Bonferroni | L63 | 上三角去重后独立 BH + Bonferroni 重算 | **76 / 75 精确复现** |
| 分块 r_g：OA×OA 1、OA×BMD 20（+0.117~+0.250，均 +0.169）、Fx×BMD 10（−0.596~−0.446，均 −0.531）、BMD×BMD 45、OA×Fx 2（0/2 显著） | Table S15 | 分块统计重算 | ✓ 完全一致 |
| coloc endpoint-split：OA 0/31（中位 0.005）vs 骨折 9/31（中位 0.153） | L79 | 31 位点重算 | ✓ |
| F4 血液 eQTL p=0.0012 | docs/118 L54 | 31 位点 PP.H4 **精确** Wilcoxon（DP 零分布） | **0.001156** ✓ |
| T3c 跨部位 Q=3.56、df=4、p=0.468、I²=0% | docs/118 L59 | 独立 Cochran Q + χ² 连分数 | ✓ |
| P2/T3a 排序单侧 p=0.228 | L95 | 独立手算 z=0.7461 | **0.2278** ✓ |
| 仪器选择表：10 个位点 GWS 数、5 个 clump 完成（各 ≥15 仪器）、5 个 deferred（229–345 loci） | L53 | `variant_selection_table.csv` 逐行 | ✓ |
| MVMR-BMI 全部数值、条件 F、衰减 −3.7/5.1/11.9/27.6/46.0% | Table S16 | `mvmr_bmi_main.csv` 逐行 | ✓ |
| 反向 MR 功效 0.14–0.26 | L99 | `mr_power.csv` | ✓（且尺度不变，不受 G6 影响） |
| 交付 Figure 3 尺度 = 稿件文本 | Fig 3 | 目视读图：panel b 显示 +0.023、轴标 "N-corrected scale" | ✓ 图-文一致 |

---

## 6. 行动清单（按风险排序）

| # | 动作 | 涉及文件 | 说明 |
|---|---|---|---|
| 1 | **OA 结局效应施加 N·K(1−K) 修正并重跑全部涉 OA 的 MR** | `code/13`（加 case-control 修正）→ 全部 harmonised → Table 2 | 因子取稿件自身 K=0.12 → 3.0773（或实测 3.339，需先定 K）。**p 值不变**，只改 β/OR/CI |
| 2 | **同步稿件全部涉 OA 的 OR / β / 差值** | `docs/08` L27/L69–70/Table 2/L95/L99/L201 + Highlights | 1.023→1.074 等；反向 0.106→0.035 等；T3a 差值 0.00509→0.01566 |
| 3 | **改写 L57 的尺度声明**（现称"all estimates share one N-verified SD scale"为假） | `docs/08` L57 | 须写明：BMD 暴露 N=31,873；fracture 结局为真实 log-OR；**OA 结局由 Z 重构于 N_eff=N·K(1−K)** |
| 4 | **F8 与主分析统一尺度后重报** | `code/161` L96（暴露 N 改 31,873）→ 重跑 → `docs/08` L49/L99 | 否则删除并列 β，只留"两性方向一致为正" |
| 5 | 修 G9（补注单侧）、G10（p=0.027 改归 Lumbar）、G11（统一分母并同步门禁）、G13（声明 RE 约定） | `docs/08` L79/L81/L65/L95，Table 2 注，`docs/118` | 均为表述层 |
| 6 | G8：将 `code/09` 的 `PHENOS` 补至 13 表型（加 fracture），使重跑可复现 | `code/09_rg_matrix.py` | 复现性 |
| 7 | 完成 M1（1000G EUR clumping + 5 位点 MR），回填 L53/L99 | 下载中（290 MB/1.5 GB） | 完成后十位点评完整 |
| 8 | 重跑 `build_package_refresh.py` + 修补 v16 包 | 打包链 | 改动落地后必须重建 |
| 9 | （可选补强）F10 祖源敏感性、F2 质量/尺寸分解 | 新增分析 | 见 §1.3 |

---

## 7. 一句话结论

**稿件的 p 值、方向、显著性、异质性、遗传相关、共定位、功效、多变量衰减全部经独立重算验证正确；但稿件报告的每一个"涉膝 OA 的效应量（OR / β / 部位差值）"都在错误尺度上——正向低估约 3.1–3.4 倍（真实 OR 1.054–1.097，而非 1.017–1.028），反向高估同倍。原因是 `code/13_mr_analysis.py` L40 对 OA 结局用了总 N 做 Z 线性重构，而未像它对骨折结局那样做 N·K(1−K) 修正。这个错误此前七轮预审、三次审计、一次 v18 全量复算（430 项 0 失败）均未发现——因为"文件内部自洽"与"文件正确"是两件事。**

---

## 8. 🔴 本轮收尾（G8 闭环）三重根因取证链与全链路修复

> 审计对象：G8 触发（fracture 补入 13×13 rg 矩阵 → 重建投稿包 → 终态核验）过程中定位并修复的三类**独立**缺陷：(A) rg 矩阵可复现性（G8）；(B) G6 尺度未传播；(C) 加权中位数约定分歧。
> 原则（沿用 §0）：每条声称附「稿件/数据行 → 源脚本行号 → 独立验证手段 → 验证结果」。§8.4 为修复后**全部重跑**的审计结果，非声称。
> 互见：F8 性别分层尺度修正（暴露 N 占位符 400,000→31,873，因子 3.5428）已在 `docs/118` L60 的 G7 附记中记录，本节不重复，仅在本节 §8.2 的"尺度"语境下互见。

### 8.1 缺陷 A — rg 矩阵可复现性（G8；连带稿件 6 处 rg 文本修订，F9）

**声称**：稿件 L53/L63 与 Table S15 报称 13×13 r_g 矩阵、78 对；但 `code/09_rg_matrix.py` 的 `PHENOS` 原仅含 **12** 个表型（缺 fracture），重跑产出 `rg_matrix.csv` 为 12×12，无法对应稿件 13×13/78 对。

**源脚本行**
- `code/09_rg_matrix.py` L56 `PHENOS = {`（原 12 项）；L18–L21 注释记录 G8 修复（补入 fracture；13 表型 HapMap3 交集恒为 690,150，与 12 表型交集一致）；L257/L263/L267 按 PHENOS 载入 13 个 sumstats。
- 修复后 `rg_matrix_13.csv` 由 13 表型重新生成（与稿件一致）。

**独立验证（稿件 6 处 rg 文本修订，docs/08 Results 段）**
- BMD 互相关 0.69–0.96（mean 0.81）；OA×BMD 20 对 +0.12~+0.25（knee +0.14~+0.19，selfrep +0.12~+0.25）；Fx×BMD 10 对 −0.40~−0.54（原 −0.45~−0.60，已同步）；78 对中 77 对 BH-FDR 与 Bonferroni 双显著（仅 selfrep×fracture 例外，p≈0.051）。
- 与重生成 `rg_matrix_13.csv` + 独立 BH/Bonferroni 重算一致（docs/123 §5 "78 对 r_g" 条为旧 12/13 口径，现以稿件 77/77 为准）。

*影响边界：此缺陷仅涉"可复现性"与"rg 文本数值"，不影响任何 MR 效应量——MR 由 harmonised 文件驱动，与 rg 矩阵无关。*

### 8.2 缺陷 B — G6 尺度未传播（"G6 未传播"）

**声称**：G6 校正（docs/123 §2）将 OA 结局效应置于正确 log-OR 尺度，因子 **F = 1/√(K(1−K)) = 3.3887**（K(1−K)=0.087082，见 `code/173_rebuild_tables_g6.py` L7/L13），已写入 harmonised 文件（post-G7 再生）。但**三处下游未同步**：

**源脚本行 / 受影响点**
1. `build_figures_v17.py` —— 主图硬编码 β（Figure 3 forest、Figure 6 scatter 等）沿用旧线性尺度。
2. `code/111_fig_S11_robustness_grid.py` L33 `TABLE3 = {` —— 10 对字典填旧尺度（如 total→OA_knee +0.023/+0.019）+ 旧中点 wmed。
3. `code/114_fig_S12_literature_comparison.py` L38–L43 `OURS = {` —— 5 对 BMD→OA 填旧尺度（0.023/0.028/0.022/0.018/0.017）。

**修复**
- `build_figures_v17.py` 硬编码 β 改为 G6 尺度（前序改，本轮重跑）→ `Figures_v17_Nature/` 写 18 图（Fig1–10, S1–S7, GraphicalAbstract）。
- `code/111` TABLE3 10 对刷新为 G6 尺度 + 步近 wmed（如 total→OA_knee +0.078/+0.065）；重跑通过（断言 AUDIT==Table 2）。
- `code/114` OURS 5 对改为 G6 值（0.078/0.095/0.075/0.060/0.058）；重跑通过。

**独立验证**：`code/115_audit_figures_v18.py` 189 checks / 0 failures → PASS（含 S9–S12 全图）。

> 注：docs/123 §2.4/§7 用稿件自报 K=0.12（因子 3.0773）作近似，得 OR 1.054–1.097；code/173 实际用 K(1−K)=0.087082（因子 3.3887），得稿件现报 OR 1.060–1.099。差异 <1%，不影响方向与 p。

### 8.3 缺陷 C — 加权中位数约定分歧（"wmed 约定分歧"，影响面 > G6）

**声称**：项目正统加权中位数 = Bowden 2016 **步近约定**（ratio 按权重升序，`cumsum(w)/Σw` 首次 ≥ 0.5 处取该 ratio）。但下列位置用了**中点约定**（`cumsum(w) − 0.5·w` 跨 0.5 取相邻均值），致 7/10 对 wmed 点估计与当前 harmonised 复算偏差 3–4 倍：
- `code/90_audit_mr_recompute.py` 原 `weighted_median`（中点）
- `code/100_audit_independent_recompute.py` 原 `SCALE` + `TABLE3` 旧值
- `docs/08` Table 2 加权中位数列（7/10 点估计）
- 旧 `table3_corrected.csv` 的 wmed 列
- `code/111` / `code/114` 硬编码值

**源脚本行 / 修复**
- `code/90` L103–L119 `weighted_median` 改为步近（`np.searchsorted(cum/tot, 0.5)` 取 `r[idx]`）；重跑 → `AUDIT_mr_recomputed.csv` 的 wmed 与 `mr_results_corrected.csv` 逐对一致。
- `code/100` L35 `SCALE = 1.0`（注释：harmonised 已正确尺度、no rescaling；旧 `0.2823` 作废），L45 `TABLE3` 10 对刷新为 G6 尺度 + 步近 wmed；重跑 data integrity max|shipped−recomputed| ≈ **1e-15**。
- `docs/08` Table 2 加权中位数列 8 处修订（BMD_total→fracture −0.364→−0.393；Lumbar→OA_knee +0.126→+0.115；Arms→OA_knee +0.104→+0.026；Legs→OA_knee +0.084→+0.026；OA_knee→fracture +0.035→+0.019；OA_knee→BMD_total +0.033→+0.026；OA_knee→femur_neck +0.033/+0.128→+0.023/+0.122；Femoral neck→OA_knee CI hi +0.110→+0.109）。脚注 L165 已声明 "seven weighted-median values corrected"。
- `code/115` L82–L84 θ_H 断言改 `mvmr_h` 正确尺度（0.143–0.166）；原 0.042–0.050 为旧尺度 `uni_h_beta` 残留（git d7d59b2 原意图即 mvmr_h）。

**独立验证（权威裁决）**：`code/178_diagnose_wmed.py` 从当前 harmonised 独立复算两种约定：
- **step-near 列 = `mr`（mr_results_corrected）= `wm_rec`（wm_recomputed）逐对完全一致（10/10）**；
- mid（中点）列与旧 `man`/`t3`/`wm_tab` 列彼此接近，但**偏离 step-near 列 3–4 倍**（如 Arms→OA_knee step +0.0257 vs 旧 +0.104）。
→ 裁定：**步近约定 × G6 正确尺度（= `mr_results_corrected.csv` wmed 列）为加权中位数唯一真源**。

### 8.4 全链路修复后终态核验（全部重跑，非声称）

| 审计脚本 | 入口 | 结果 |
|---|---|---|
| `code/178_diagnose_wmed.py` | 全文件 | step==mr==wm_rec 10/10；中点列偏离 3–4× |
| `code/177_verify_table2_cells.py` | 全文件 | 0/10 行不一致 |
| `code/115_audit_figures_v18.py` | L82–L84 | 189 checks / 0 failures → PASS |
| `code/90_audit_mr_recompute.py` | L103–L119 | AUDIT wmed 与 mr_results_corrected 逐对 OK |
| `code/100_audit_independent_recompute.py` | L35/L45 | data integrity ≈ 1e-15 |
| `code/111` / `code/114` | L33 / L38 | 重跑通过（断言 AUDIT==Table 2） |
| `build_figures_v17.py` | — | 18 图写入 `Figures_v17_Nature/` |
| `code/135_build_manuscript_ready.py` | — | `Manuscript_Ready/` 重建（19 主图 + 9 补图 PNG，占位符 0） |

### 8.5 一句话结论

G8 闭环触发对全链路三类独立缺陷的修复：(A) rg 矩阵可复现性（code/09 PHENOS 12→13，稿件 6 处 rg 文本修订）；(B) G6 尺度未传播到 build_figures / code/111 / code/114（已同步，F=3.3887）；(C) 加权中位数约定分歧（中点→步近，7/10 Table 2 点估计偏差 3–4× 已纠正，唯一真源 = `mr_results_corrected.csv` 步近列）。修复后 5 个独立审计全绿（含 189/0 图审计、data integrity 1e-15），**稿件 Table 2 与全部 18 张图现已与 harmonised 真源逐字段收敛**。历史文档 docs/07 / docs/15 / docs/28 / docs/107 含 G6 前 / 约定分歧前的旧值，见 §8.6 标注。

### 8.6 历史文档旧值标注（docs/07 / docs/15 / docs/28 / docs/107）

- **docs/07**（LDSC 遗传相关矩阵）：标题即"12 表型（OA × BMD）"，记录的是 **G8 前 12 表型**子集；其 BMD 互相关 / OA×BMD 数值本身正确（rg 不经 OA 结局重构），但缺 fracture 行，且任何引用的 MR β/OR/wmed 须以 docs/08 + `mr_results_corrected.csv` 为准。
- **docs/15**（全面复盘审计报告，2026-09-17）：L101 "校正后 β=+0.0231，OR≈1.02–1.05"、L128 表 BMD total→OA_knee +0.023 等为 **G6 校正前旧尺度值**；现正确值为 β=+0.078（OR 1.060–1.099）。该文档结论仅适用于当时口径。
- **docs/28**（数据级终审 v20，2026-10-10）：L17/L119/L145–L146 将 wmed"正确值"记为 −0.3642（中点约定），称 −0.3926 为错。**该结论已被本节（缺陷 C）推翻**：步近约定为项目正统，−0.3926（= `mr_results_corrected`）方为正确；旧 −0.3642 为中点约定产物，属需纠正旧值。读 docs/28 须以此注为准。
- **docs/107**（图件核验与补充图规划，2026-10-16）：Fig 1 / S13 的 rg 比对结论（max|偏差|=0.00000）仍成立（rg 不受 G6/wmed 影响）；但其引用的图件 β/wmed 须以本轮重生图件为准。
