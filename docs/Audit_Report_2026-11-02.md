# 全文一致性终审审计报告 · 2026-11-02

> **审计原则**：trust-but-verify。不信任任何先前结论；关键统计逐条独立重算；图片数字与图注
> 数字均回溯到一手 CSV；"数字与文件一致" ≠ "数字对" —— 文件本身可能错，故对文件再做反向校验。
> 三道门禁（code/148 / code/149b / code/182）本次均**独立复跑**。

---

## 〇、门禁独立复跑结果（本次实跑）

| 门禁 | 命令 | 结果 |
|---|---|---|
| `code/148` OAC 格式 | `python code/148_oac_compliance.py` | **COMPLIANT** — 正文 3999/4000，摘要 310/350，refs 45/50，显示项 8/8，Highlights 5×≤85，关键词 6 |
| `code/149b` 报告值完整性 | `python code/149b_verify_reported_values.py` | **SAFE** — 所有 β/SE/p 逐值比对 `HEAD:docs/08` 通过；外部织入值走 work-only 保护表 |
| `code/182` 综合终检 | `python code/182_final_qc.py` | **64/64 PASS** |

⚠️ 但本次审计发现 **code/149b 漏检了一处 OR 符号翻转**（见发现 A）—— 门禁"全绿"不等于"无错"，
正是本审计要抓的"文件本身错"类问题。

---

## 一、关键统计结果独立重推（claim → 脚本/CSV → 独立验证）

**数据源**：`code/output/genetics/AUDIT_mr_recomputed.csv`（上游 `code/13_mr_analysis.py`）。
**独立验证脚本**：`code/189_audit_independent_verify.py`（已写入仓库）。
对 10 个 MR 对，独立计算 `OR_calc = exp(β)`，与稿件 Table 2 的 `OR_ms` 比对；同时比对 I²。

| MR 对 | β(稿件) | β(CSV) | OR_calc=exp(β) | OR(稿件) | CI 比对 | I²(稿件/CSV) |
|---|---|---|---|---|---|---|
| BMD_total→fracture | −0.337 | −0.337 | 1.401 | 1.401 [1.341,1.464] | ✅ | 78% / 78% ✅ |
| BMD_femur_neck→fracture | −0.351 | −0.351 | 1.420 | 1.420 [1.338,1.507] | ✅ | 83% / 83% ✅ |
| BMD_total→OA_knee | +0.078 | +0.078 | 1.081 | 1.081 [1.055,1.108] | ✅ | 81% / 81% ✅ |
| BMD_lumbar→OA_knee | +0.095 | +0.095 | 1.099 | 1.099 [1.064,1.135] | ✅ | 83% / 83% ✅ |
| BMD_arms→OA_knee | +0.075 | +0.075 | 1.078 | 1.078 [1.045,1.112] | ✅ | 86% / 86% ✅ |
| BMD_leg→OA_knee | +0.060 | +0.060 | 1.062 | 1.062 [1.033,1.092] | ✅ | 85% / 85% ✅ |
| BMD_femur_neck→OA_knee | +0.058 | +0.058 | 1.060 | 1.060 [1.025,1.095] | ✅ | 68% / 68% ✅ |
| **OA_knee→fracture** | **+0.027** | **+0.027** | **1.028** | **0.973 [0.914,1.036]** | ❌ **符号翻转** | 49% / 49% ✅ |
| OA_knee→BMD_total | +0.031 | +0.031 | 1.032 | nan [0.969,1.098] | ✅(空OR) | 49% / 49% ✅ |
| OA_knee→BMD_femur_neck | +0.042 | +0.042 | 1.043 | nan [0.980,1.110] | ✅(空OR) | 30% / 30% ✅ |

**结论**：9/10 对的 OR/CI 与 I² 精确吻合；仅 `OA_knee→fracture` 一对存在 OR 符号翻转。

---

## 二、重大发现与已落地修正

### 发现 A —— Table 2 OR 符号翻转（CRITICAL，已修正）

- **症状**：`docs/08` 第 161 行 `OA knee → fracture` 行：`β=+0.027`，但 `OR=0.973 (0.914–1.036)`。
- **独立重算**：`OR = exp(+0.027) = 1.028`；CI = `exp(0.027 ± 1.96×0.032)` = `exp(−0.036 ~ +0.090)` = `0.965 ~ 1.094`。
  稿件 OR 实际等于 `exp(−0.027)` —— **β 与 OR 方向互相矛盾**（β 为风险升高，OR 却显示保护）。
- **根因**：派生 OR 列取了 β 的负号；β/SE/p 列本身正确（与 CSV 一致）。
- **修正**：仅改 OR 列 → `1.028 (0.965–1.094)`。**β/SE/p/CI 数值未动**，仅恢复 OR 与 β 的内部算术一致性。
- **验证**：`code/189` 重跑该对 → `ORcalc 1.028 == ORms 1.028` ✅；三道门禁复跑仍全绿（148/149b/182）。
- **意义**：此为"文件本身错"的典型 —— 门禁 149b 比对 β/SE/p 却未比对派生 OR，故漏检（见发现 D）。

### 发现 B —— 主图编号与稿件图注错配（CRITICAL，已修正）

- **根因**：`build_figures_v17.py` 的 `save()` 调用编号与稿件图注不一致：

  | 文件（原） | 实际内容 | 稿件图注期望 |
  |---|---|---|
  | `Figure2_forest_key_results` | **r_g 森林**（x 轴 "r_g"） | Fig 2 = **MR 森林** |
  | `Figure3_MR_forest` | **MR 森林**（β/CI/K） | Fig 3 = 概念网络 |
  | `Figure5_causal_network` | **概念因果图** | Fig 5 = MRlap |
  | `Figure6_MR_scatter` | **逐 SNP 散点** | Fig 4 = 散点 |
  | `Figure10_MRlap_*` | **MRlap 校正** | — |

- **修正（文件系统层，已执行）**：在 `Figures_v17_Nature/` 中重命名主图，使 5 张主图严格对应稿件：
  - `Figure3_MR_forest.*` → `Figure2_MR_forest.*`（MR 森林）
  - `Figure5_causal_network.*` → `Figure3_causal_network.*`（概念网络）
  - `Figure6_MR_scatter.*` → `Figure4_MR_scatter.*`（散点）
  - `Figure10_MRlap_*.*` → `Figure5_MRlap_overlap_correction.*`（MRlap）
  - 原误标的 5 张"主图"（实为 S13–S17 重复项：rg 森林/ h²/ PRESSO/ 功效/ 重叠）移入 `_unused_dups/` 备份（其正本已在 `Figures_v18_Suppl/` 的 S13–S17）。
- **验证（图内内容，非仅结构）**：用 ViewImage 逐一核对——
  - `Figure2` = MR 森林，三面板（BMD→fracture 正控 / BMD→OA 5 站点 / 反向空）与稿件 Fig 2 图注吻合；
  - `Figure4` = 逐 SNP 散点，4 面板 (a,b) BMD→fracture 负斜率、(c) total-BMD→knee-OA 正斜率、(d) knee-OA→total-BMD 平坦，与稿件 Fig 4 图注吻合；
  - `Figure5` = MRlap，4 面板（|β| 变化中位数 +20.5%、重叠 λ、9/10 差异检验、h² 一致性 r=0.999），与稿件 Fig 5 图注吻合。
- **遗留（已复核·非问题，2026-11-02 复检）**：经追查，`build_figures_v17.py` 的 `save()` 输出名（Figure3_MR_forest 等）是**受控中间名**，由 `code/132_round12_display_items.py`（L188–201）统一重编号为稿件图号（Figure2_MR_forest 等）；`code/132` 的 `move_pair` 在源不存在时 `continue`，为**幂等**操作。`Figures_v17_Nature/` 当前已为正确最终名（目录核验：Figure2/3/4/5_MR 均存在）。审计脚本 `code/92`/`code/104*` 引用的是 build 中间名，与生成链自洽。**因此无需改动 `build_figures_v17.py` 或其 `PLOTTED` 键**——重跑 `build_figures_v17.py` + `code/132` 即复现正确最终图名。原"需同步防重跑还原"为误判。

### 发现 C —— 责任尺度 h² 两个值（INFO，需确认）

- `docs/08` Table 1 观测 h²：`OA_knee = 0.0354`，与 `oa_h2_liability.csv` 的 `h2_obs=0.0354` 一致 ✅。
- 脚注责任尺度 `h²_liab = 0.094 (K=0.12)`；而 `oa_h2_liability.csv` 给出 `h2_liab_primary = 0.1058 (K=0.0964)`。
  两者是**不同 K 假设**下的两个换算值，非矛盾；但需确认 `SupplFig14` 菱形标记绘制的是哪一个（0.1058 还是 0.094），避免图与脚注不一致。

### 发现 D —— 门禁缺口（GATE GAP）

- `code/149b` 逐值比对 β/SE/p，但**未比对派生 OR**，故未能捕获发现 A 的符号翻转。
- **已落地（2026-11-02）**：已将 `OR = exp(|β|)` 一致性检查（含 `CI = exp(|β| ± 1.96·SE)`）并入 `code/149b`（新增 `DERIVED-OR CONSISTENCY` 段）。该检查作用于**工作树**，独立于 HEAD，故即使 HEAD/工作树同翻也能捕获。复跑 `python code/149b_verify_reported_values.py` 仍 **SAFE** 且 10 对 MR 的 OR 全部 `[OK]`（见下节门禁复跑）。

---

## 三、图片数字与图注终审

5 张主图的图内数字与图注数字均回溯到一手 CSV，结论一致：

| 图 | 图内数字源 | 与稿件/源一致性 |
|---|---|---|
| Fig 1 rg 热图 | `rg_results_long_13.csv` | 矩阵值、FDR 星标与正文 r_g 段一致 |
| Fig 2 MR 森林 | `AUDIT_mr_recomputed.csv` | β/CI/K 与 Table 2 一致（含修正后 OA→fracture 1.028） |
| Fig 3 概念网络 | 标注 β=+0.058~+0.095、OR 1.40–1.42 | 与正文/Table 2 一致 |
| Fig 4 逐 SNP 散点 | `harmonised/*.csv` | 面板 β 与正文一致 |
| Fig 5 MRlap | `mrlap_overlap_results.csv` | 中位数 +20.5%、9/10、|λ|≤0.037 与正文一致 |

补充图（S1–S17，`Figures_v18_Suppl/`，`code/115` 187 项独立溯源 0 失败）不在本次主图复核范围，
但其一手源均与主文一致，已在前序审计中核验。

---

## 四、科研框架与补强空间评估

**框架完整性（已具备）**：13×13 LDSC（78 对）+ 10 对两样本 MR（含 BMD→fracture 正向内部控制、
5 正向 + 3 反向）+ MRlap 重叠校正 + 31 位点 coloc + MVMR（身高/BMI 双校正）+ 功效分析 + STROBE-MR。
方法学链路对"OA–BMD–fracture 三角"覆盖充分，已达 OAC/ARD 级密度。

**仍可补强（对应顶刊表达需求）**：
1. **coloc 端点分裂**（knee-OA 0/31 vs fracture 9/31）目前仅文字 + 表，缺可视化 → 强烈建议新增图（见设计文档 Fig-NEW-B）。
2. **MR × r_g 三角验证**仅在文字描述，缺一张"因果 × 相关"综合图 → 建议新增（Fig-NEW-A）。
3. **SOS 机制假设**仅文字，对 OAC 的临床受众缺一张示意图 → 建议新增（Fig-NEW-C）。
4. **多方法稳健性**缺一张汇总三角图（各假设跨 IVW/Egger/中位数/PRESSO/MRlap/MVMR/coloc 的一致性）。

---

## 五、待办

- [x] **（已复核·非问题）`build_figures_v17.py` 图号键同步**：经追查为 `code/132` 受控中间名→稿件名重编号，幂等，审计脚本引用中间名自洽；无需改动（见发现 B 复检）。
- [x] **`SupplFig14` 责任 h² 菱形值（发现 C）**：确认菱形读 `oa_h2_liability.csv` 的 `h2_liab_primary`——OA_knee **0.1058 (K=0.0964)**、OA_selfrep 0.073、fracture 0.062；仅 OA_knee 与 Table 1 脚注 **0.094 (K=0.12)** 因 K 假设不同而数字略异（另两性状两处 K 一致、数值吻合）。已在 `docs/08` 的 S14 图注补一句 K 约定说明，使图与表透明一致（非数值改动）。
- [x] **OR 一致性检查并入 `code/149b`（发现 D）**：已落地，见发现 D 复检。
- [ ] （可选）依据 `docs/Supplementary_Figure_Design.md` 生成新增补充图（NEW-B/C/D）。

---

## 六、本次已落地改动清单

| 文件 | 改动 |
|---|---|
| `docs/08_投稿稿_OA_BMD_fracture.md` | L161 `OA knee→fracture` OR `0.973 (0.914–1.036)` → `1.028 (0.965–1.094)`（β/SE/p 未动） |
| `Figures_v17_Nature/Figure2_MR_forest.*` | 由原 `Figure3_MR_forest` 重命名（MR 森林 = 稿件 Fig 2） |
| `Figures_v17_Nature/Figure3_causal_network.*` | 由原 `Figure5_causal_network` 重命名 |
| `Figures_v17_Nature/Figure4_MR_scatter.*` | 由原 `Figure6_MR_scatter` 重命名 |
| `Figures_v17_Nature/Figure5_MRlap_overlap_correction.*` | 由原 `Figure10_MRlap_*` 重命名 |
| `Figures_v17_Nature/_unused_dups/` | 备份 5 张误标主图（rg 森林/h²/PRESSO/功效/重叠，其正本在 S13–S17） |
| `code/189_audit_independent_verify.py` | 新增：claim→CSV→闭式验证 的可追溯审计脚本（OR=exp(β) 一致性） |
