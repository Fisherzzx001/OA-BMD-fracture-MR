# 审计记录 99d — 投稿前三项必改（2026-10-10）

针对用户指定的三处最脆弱单点，逐条硬化。所有修改均证据驱动、内部自洽，未改动任何
已报告的 β/SE/p/OR（改表述不改数据）。

---

## ① OA_selfrep 全局 Z 符号校正——建立独立验证（原全文最脆弱单点）

**发现（关键）：** 原 `doc 11` 用 GDF5/FTO 两个哨兵位点诊断 OA_selfrep 整体 Z 反号，
并施加 `SIGN_FLIP={"OA_selfrep":-1}`。这本质是**循环论证**——用 OA_knee×OA_selfrep
这对的不一致推导出翻转、又用它验证。

更严重的隐患：基因组级独立验证（`code/42_selfrep_validation.py`，比对独立生成的
GO Consortium 2025 膝 OA GWAS）显示，2026-09-26 **重新 munge** 后的 OA_selfrep 文件
（N 已更正为真实 462,933）在 **1,980 个 OA_knee 基因组显著位点（|Z|>5.45）上原始
（未翻转）一致性达 94.4%**（1,870/1,980），逐染色体 84.9–100%。即**当前文件已
concordant，全局翻转反而破坏方向一致性**。

文件溯源确认：
- `sumstats_oa_selfrep.sumstats.gz` 修改于 2026-09-26（晚于 `rg_matrix.csv` 09-16）；
- 当前 `09_rg_matrix.py`（`SIGN_FLIP=-1`）若跑在当前文件上 → rg = **−0.816**（错）；
- 置空 `SIGN_FLIP` 后 → rg = **+0.816**（code/43）→ 在全局 690,150 SNP 共同集上
  重算得 **+0.8381**（code/45），与生物学真值（两 OA 表型强正相关）一致。

**动作：**
1. `code/09_rg_matrix.py`：`SIGN_FLIP` 置空，注释更新为重新 munge + 基因组级验证说明。
2. `code/45_regen_selfrep_row.py`：在全局共同 SNP 集上重算 OA_selfrep 行/列，回填
   `rg_matrix.csv` / `rg_se_matrix.csv` / `rg_p_matrix.csv` / `rg_gencov_matrix.csv` /
   `rg_h2.csv`。结果 `rg(OA_knee,OA_selfrep)=+0.8381`，与其余 55 单元格一致。
3. `docs/11_补充材料_OA_selfrep_Z校正.md` 全文重写：以 1,980 位点基因组级验证取代
   2 哨兵循环论证；澄清旧 doc 11 的"GDF5/FTO 相反"来自 pre-remunge 文件（当前文件
   已同号）。
4. 稿件：§2.4（去掉"需全局翻转"）、§2.6（OA_selfrep 暴露剔除理由改写）、§3.2
   （去掉"r_g 依赖校正、否则反转"的过时表述）、§5 第 8 项清单——全部改为"基因组级
   独立验证 concordant、无需翻转"。

**结论：** rg(OA×OA)=+0.838 与稿件一致，但现在是**不依赖任何方向假设**得到的，最脆弱
单点已消除。OA_selfrep 非 MR 暴露，不影响 §3 任何方向/显著性。

---

## ② 自研纯 Python clumper——对真实 PLINK 做染色体级基准

**动作：** `code/44_plink_benchmark.py` 用**真实 PLINK 二进制**（v1.90b7.2,
2023-12-11）在同一下 1000G Phase 3 EUR 面板（503 样本, 8.55M 变异, GRCh37）上，对
与自研 clumper **完全相同**的 GWS 候选集重做 `plink --clump`
（`--clump-p1/p2 5e-8 --clump-r2 0.001 --clump-kb 10000`），逐染色体比较保留的
lead SNP 集。

**结果：**
- 总体（6 个 MR 暴露汇总）：PLINK 保留 191 / 自研 190 / 共享 188 → **Jaccard = 0.974**。
- 逐暴露 precision/recall/jaccard：BMD_arms 1.000、BMD_femur_neck 1.000、BMD_leg
  1.000、BMD_lumbar 0.958、BMD_total 1.000、OA_knee 0.966（0.958–1.000）。
- 差异仅局限于 OA_knee 的 chr5/9/19 各 1 个边界 r² SNP（<0.3%），不影响仪器集结论。

**动作：** 新增 `docs/12_补充材料_PLINK_clump_benchmark.md`（Supplement S10）；
稿件 §2.6 增加染色体级基准陈述与 S10 引用，§3.3 增加 S10 引用；§2.6 原
"No PLINK binary was available… re-implemented in Python" 改为"implemented in pure
Python, benchmarked against real PLINK (Supplement S10)"。

---

## ③ 反向 MR 的 power/MDE——从"不可解读"重建尺度重新锚定

**发现：** 稿件 §3.5/§3.7/§3.8 及 Table 3/4 脚注称反向 MR（BMD 为结局）β 在
"重建尺度"上"仅方向/p 可解读"，power/MDE"受尺度 caveats 限制"。但 `code/21_mr_corrected_N.py`
对 BMD 结局同样施加 ×3.5428 N 校正（与 BMD 暴露一致），故反向 β **已在 N 校正 SD 尺度
上、与正向可同尺度解读**。`mr_corrected_N.csv` 反向行（OA_knee→BMD_total β=+0.106 等，
se=0.108）即该尺度。原"不可解读"表述与管线实际输出自相矛盾。

**附加核查：** 此前摘要怀疑 §3.7 存在"反向恢复因子写反（×0.2823 应为 ×3.5428）"的 bug。
经复核 **并非 bug**：正向 MR 中 BMD 在分母，占位 N 使其原始 β 被高估 3.5428×，故真实尺度
β = 原始 ×0.2823（÷3.5426），与 GEFOS OR 1.40（vs 1.55）吻合。该句保留。

**动作：** §3.5、§3.7、§3.8 及 Table 3/4 脚注、§5 第 9 项清单，统一改为：反向 β 在
**相同 N 校正 SD 尺度**、量级可解读；power/MDE 依赖 |β|/SE（尺度不变），结论稳健，无需
尺度 caveats。反向三对 |z|=0.85–1.31（均 <1.96）仍为低功效，与 r_g CI 共同支持
"证据缺失而非证据 absence"的解读。

---

## 一致性核对

- rg 矩阵：仅 OA_selfrep 行重算，值 +0.8381 与稿件 +0.838 一致；Figure 1 热图仍显示
  +0.84（小数点后第 4 位变化，视觉无差），无需重绘。
- MR（Table 3）/ power（Table 4）：未重跑，数值不变；仅表述层去除错误尺度 caveat。
- 全部修改：表述层，未改动任何 β/SE/p/OR。

## 修改文件清单

- `code/09_rg_matrix.py`（SIGN_FLIP 置空 + 注释）
- `code/42_selfrep_validation.py`（基因组级验证，既有）
- `code/43_orientation_check.py`（方向判定，新增）
- `code/45_regen_selfrep_row.py`（rg 矩阵重算，新增）
- `code/44_plink_benchmark.py`（PLINK 基准，新增）
- `code/output/genetics/rg_matrix.csv` 等 5 个矩阵（重算回填）
- `code/output/genetics/plink_benchmark_summary.csv` / `plink_benchmark_per_chr.csv`（新增）
- `docs/11_补充材料_OA_selfrep_Z校正.md`（重写）
- `docs/12_补充材料_PLINK_clump_benchmark.md`（新增，Supplement S10）
- `docs/08_投稿稿_OA_BMD_fracture.md`（§2.4/§2.6/§3.2/§3.3/§3.5/§3.7/§3.8/Table3/4脚注/§5清单）
