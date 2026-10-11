# R4 逻辑 7.0 · 主链自洽审计：≥80%「属 areal BMD」的概念偷换

## 问题（用户指出的脆弱点）
主链自洽偏脆：稿件把「BMD→OA 效应中**未被身高解释**的部分」等同为「**由骨量（bone mass）解释**」，
在 §3.10 / §4.2 / §4.3 / 摘要 / S9 表注中落点为「≥80% 属 areal BMD = 骨量效应」。
这是**概念偷换（equivocation）**：aBMD = BMC ÷ 投影面积，是一个「骨量 + 骨尺寸」的复合体；
「height 未解释」是「骨量解释」的**超集**，而非等同。

## 为什么原论证站不住（深度拆解）
1. **可识别性（identifiability）缺失**：MVMR 只把 height 当 size proxy。
   height 捕获的是「长骨身高的线性尺寸」，并不捕获 aBMD 残差里仍混有的其他尺寸分量
   （骨横截面积、皮质厚度、股骨/关节几何）。因此 aBMD|height 残差 = 真实体积骨密度(residual BMC) + height 未捕获的尺寸。
   「未被身高解释」⊃「由骨量解释」—— 骨量只是残差的**子集**，不是全部。
2. **正向对照只证明「兼容」而非「等同」**：BMD→fracture 经 height 校正后衰减 −2.8%，
   说明 height 校正对「已知骨量效应」是 *mass-transparent*（不侵蚀骨量效应）。
   这能证明 BMD→OA 中观测到的 9.5–21.2% 衰减是**真实尺寸信号**（不是骨量被稀释），
   但 mass-transparency ⇒ 残差**兼容**骨量，并不 ⇒ 残差**就是**骨量。兼容 ≠ 同一。
3. **数值推论本身也偏松**：衰减跨位点 9.5–21.2%，残差 = 78.8–90.5%；
   原文「at least 80%」把最小值 78.8% 上舍入到 80%，且把残差直接读作「骨量份额」——
   正确读法：残差的「骨量份额」是**上界（≤≈79–90%）**，其下界由未测尺寸决定，而非由身高决定。
4. **主链下游放大**：§4.3 把 BMD→OA 默认当成「系统性骨量蓄积 → 骨软骨单元」机制，
   但该机制只有在「残差是骨量介导」时才成立；若残差实为未测尺寸（关节力学），机制就变成尺寸→力学，
   二者本设计无法区分。「≥80% 骨量」一旦被推翻，§4.3 的机制叙事就悬空 → 主链自洽破裂。

## 数据溯源（确认数值本身无误，只是解读越界）
- §3.10 的 9.8–21.2% / θ_H=+0.0426 / 正向对照 −2.8% 精确对应 `mvmr_height_grid.csv` 的
  `EUR|n_h=100` 行（`atten_vs_full` 列：total 9.79 / femur 14.58 / lumbar 9.54 / arms 10.96 / leg 21.23）。
  即正文与其引用的 S9 表**数值自洽**；`mvmr_height_main.csv` 是另一 config（stale），非正文来源，无关。
- 残差精确范围 78.8–90.5% → 正文中统一改为「≈79–90%，独立于 height–size 通路」。

## 已落地的修订（仅改表述，未动任何 β/SE/p/OR）
- **§3.10（line 205）**：残差改为「≈79–90% 独立于 height–size 通路（未被身高解释）」；
  显式声明「未被身高解释 ≠ 由骨量解释」，骨量份额为**未解决的上界（≤≈79–90%）**，非点估计。
- **§4.2（line 219）**：删除「at least 80% attributable to the areal-BMD phenotype … not to skeletal size alone」
  及「predominantly … rather than of stature」的骨量化落点；改为「≈79–90% 独立于 height–size 通路」，
  并展开 mass-transparent = 兼容而非证明；明确骨量份额 = **上界**，下界由未测尺寸设定；
  分离质量需 vBMD/髋几何 GWAS 作第二暴露做 conditional MR。
- **§4.3（line 225）**：骨软骨单元机制改为**条件句**——「若残差确为骨量介导（mass-transparent 兼容但未证），
  则系统性骨量蓄积经骨软骨轴升高 OA 风险」；并列替代假设「残差由 height 未捕获尺寸驱动、经关节力学而非骨形成」，
  本设计无法区分。
- **§4.7（line 241）**：结论「BMD→OA 效应属于临床实测 aBMD 表型，其质量-尺寸拆分未能分离（§4.2）」。
- **S9 表注（line 417）**：「at least 80% attributable to the areal-BMD phenotype itself」→
  「≈79–90% 独立于 height–size 通路；残差究竟为骨量还是未测尺寸未解决」。
- **摘要（line 16 结论 / line 24 引言）**：临床启示改为「升高 aBMD 的治疗（经骨量或经尺寸）应记录 Incident OA」，
  并注明本设计无法分离 aBMD 哪一分量驱动该效应。

## 一致性终审
- 全文 grep：「at least 80% / ≥80% / attributable to the areal-BMD phenotype itself / not to skeletal size alone」→ 0 处。
- 所有「bone mass」用法的语境现均为：正向对照（fracture=骨量，正确）/  hedged 假设 / 条件机制 / 上界声明。
- 主链 §2.10 → §3.10 → §4.2 → §4.3 → §4.4(已 hedged) → §4.7 → 摘要 现自洽：
  「未被身高解释」与「由骨量解释」被严格区分；骨量份额为未解决的上界。
- 中心结论「genetically higher BMD raises knee OA risk」不受影响（univariable + MVMR 均正向、正向对照干净）。
