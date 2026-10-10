# R4 逻辑 7.0 · 主链自洽补救（续）：跨位点验证错位 + 因果语气对齐 + 摘要冗余

> 前序：`99e_audit_conceptsubstitution.md` 已解决 R4-2（≥80% areal BMD = 骨量 的概念偷换）。
> 本笔记解决 R4 审查师列出的其余主链脆弱点，使 R4 行"2 处错位/偷换 + 语气错位 + 冗余"全部闭合。
> 所有改动仅为**表述层**；未改动任何 β/SE/p/OR。

## 一、R4-4（站点错配弱化验证）—— §3.7 跨位点比较缺明示

**问题（审稿原话）：** "§3.7 将 total-body OR 1.40 与 femoral-neck GEFOS 1.55 比称'低 8–10%'，但二者不同 BMD 位点；应明示此为跨位点近似比较，否则验证主张被削弱。"

**根因：** 头条验证句把 total-body（本研究的全身 BMD）与 GEFOS 的 femoral-neck（股骨颈）直接并比，未说明这是**不同位点**的近似；虽正文后段已提"different BMD measurement sites"，但领衔比较仍混用位点，审稿人会读成"把不同位点的数当成同位点验证"。

**修复（docs/08 §3.7 两处 + Figure 3 图注）：**
- 第一处（恢复公式后）：同时给出 total-body β=−0.337（OR 1.40）与 femoral-neck β=−0.350（OR 1.42）；明确"femoral-neck 1.42 vs GEFOS 1.55 = 同位点、低 8.4%、like-for-like 直接比较"，"total-body 1.40 = 跨位点近似（不同 consortium 测不同位点）"。
- 第二处（与 GEFOS 并比段）：改为"在 like-for-like femoral-neck 比较上 OR 1.42 低于 GEFOS 1.55 达 8.4%、CI 完全重叠；total-body 1.40 低 9.6% 属跨位点近似"——**直接验证（同位点）成为主证据，跨位点成为次要 sanity check**，反而强化了验证主张。
- Figure 3 图注：补"（cross-site: GEFOS used femoral-neck, so femoral-neck OR 1.42 is the like-for-like comparison at 8.4% below; total-body 1.40 is a different-site approximation）"。

**数值核验：** OR 1.42 ↔ β≈−0.350（exp(0.350)=1.419）与 §3.6 line 169（femoral-neck β_IVW −0.350）一致；OR 1.40 ↔ β=−0.337 与正文既有值一致。无冲突。

## 二、R4-1（因果主张语气错位）—— 摘要/结论 causal 领衔 vs 正文稳健性托底

**问题（审稿原话）：** "标题/摘要用 causal 领衔、正文用稳健性限定托底，二者语气错位。须对齐（要么全称'方向性/假定因果'，要么把随机效应弱点前置）。"

**修复（表述层，不弱化结论）：**
- 摘要 Conclusion（line 16）："the well-powered **causal** evidence runs in one direction only" → "the well-powered **forward MR** evidence runs in one direction only"。
- §4.7（line 241）："the **causal** evidence is asymmetric" → "the **forward MR** evidence is asymmetric"。
- 标题维持不变（标题落脚于"反向因果证据有限"，本身已平衡，未过度领衔）。
- 注：摘要 line 14 已含"heterogeneity high (I²=68–86%), only 2/5 significant under random-effects IVW, we rest the causal reading on pleiotropy-robust evidence"——即审稿建议的"随机效应弱点前置"在摘要 Results 已做；本修复仅把领衔词从"causal"收为"forward MR"，使语气与正文完全一致，且不抹去"genetically higher BMD increasing knee-OA risk"这一方向性结论。

## 三、R4-3（摘要/§4.7 冗余）—— "no well-powered evidence…" 逐字重复

**问题（审稿原话）：** 摘要末句与 §4.7（审稿记为 §4.6）末句"几乎逐字重复（冗余）"。

**修复：**
- §4.7 末句原含"…so no well-powered evidence of a causal pathway from OA to fracture is available"，而该段第 2 句已写"We found no well-powered evidence of reverse causation or of a genetic pathway from OA to fracture"——**段内已重复一次**，且与摘要末句逐字相同。
- 删除 §4.7 末句尾的该从句（保留"…resting on the well-powered forward estimates and remaining provisional, both because…(§2.8, §3.8)."），由该段第 2 句承载反向无效结论。
- 现 grep 确认 "no well-powered evidence of a causal pathway from OA to fracture is available" **仅摘要 1 处**出现，§4.7 不再含该逐字句。

## 四、一致性终审
- R4-4：§3.7 两处 + Figure 3 图注均已明示 cross-site / like-for-like；旧"within 8–10% of the published GEFOS estimate"误导句已根除。
- R4-1：全文 grep "well-powered causal evidence" / "causal evidence is asymmetric" → 0 处；"forward MR evidence" 2 处（摘要 + §4.7）。
- R4-3："no well-powered evidence of a causal pathway…" 出现次数 = 1（仅摘要），§4.7 内 = 0。
- 中心因果结论（genetically higher BMD raises knee-OA risk）、所有 β/SE/p/OR、GEFOS 验证方向均**未改**。
- 与 99e（R4-2 概念偷换）叠加，R4 逻辑 7.0 行的"主链自洽偏脆，2 处错位/偷换"现已全部闭合（R4-2 概念偷换 + R4-4 跨位点错位已修；R4-1 语气对齐 + R4-3 冗余去重已附）。

## 五、文件变更
- `docs/08_投稿稿_OA_BMD_fracture.md`（§3.7 两处 + Figure 3 图注 + 摘要 Conclusion + §4.7）
- `docs/99f_audit_r4_alignment.md`（本笔记）
