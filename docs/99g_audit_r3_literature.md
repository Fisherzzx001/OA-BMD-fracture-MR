# 审计说明 · R2 三条必改核验 + R3 文献对话补强（99g）

## 背景
审查文档 `docs/99_学术预审委员会_全面审稿.md` 的 R2（方法论，3 必改）与 R3（文献对话，5 缺口）在本轮推进。
结论：**R2 三条必改已于前序提交（a915e4a / 581f74c）落地，本轮核验其真实落位；R3 三个命名缺口（身高→OA、骨文献、首个 OA→fracture 核验）+ 两点增强（R3-④ 骨尺寸锚点、R3-⑤ coloc p12 理由）本轮闭合。**

---

## 一、R2 三条必改——核验真实落位（无编辑，确已存在）

| 必改项 | 稿件落位证据 | 结论 |
|---|---|---|
| ① OA_selfrep 全局 Z 符号校正独立验证 | §2.4/§3.2 已写"genome-wide validation vs GO Consortium 2025 knee-OA（1,980 位点，94.4% 一致，no flip）"；`code/09` 中 `SIGN_FLIP = {}`；rg(OA_knee,OA_selfrep)=+0.8381（code/45） | ✅ 已落（commit a915e4a + 99d） |
| ② PLINK clumper 染色体级基准 | §2.6 已写"benchmarked chromosome-by-chromosome against real PLINK binary (v1.90b7.2)…Jaccard 0.974 overall (per-exposure 0.958–1.000)…Supplement S10" | ✅ 已落（commit a915e4a + 99d + 补充 S10） |
| ③ 反向 MR power/MDE 尺度重锚 | §3.8/Table 4 脚注已写"N-corrected SD scale / scale-invariant"；反向行 mr_corrected_N.csv 已 N 校正（F_BMD=3.5428 双向） | ✅ 已落（commit a915e4a + 581f74c） |

> 说明：R2 三项在"三处必改"段已解决，本轮仅做 grep 核验确认无残留"不可解读/非独立"表述。grep 显示"non-interpretable / 不可解读 / 旧口径 global Z-sign"均为 0 处。

---

## 二、R3-① 身高→OA 流行病学（已引用 + 本轮显式回链）

**核验（web 检索，确认非编造）：**
- [24] Welling et al. 2017, *Int J Rheum Dis* 20:1095–1104, PMID 28296269 — Northern Finland Birth Cohort 1966：最高身高四分位 knee-OA HR 男性 **2.5**(95% CI 1.4–4.5)、女性 **1.8**(1.0–3.1)。✅ 与稿件一致。
- [25] Hunter et al. 2005, *Arthritis Rheum* 52:1418–1423, PMID 15880346 — Beijing Osteoarthritis Study：greater knee (leg) height → knee-OA prevalence via joint torque。✅ 真实。
- [26] Sun et al. 2023, *PLoS One* 18:e0279198, PMID 36716300 — height→OA univariable OR **1.115**(95% CI 1.054–1.180)，MVMR 中调整后趋零（仅腰围留存）。✅ 与稿件一致。

**本轮编辑（§3.10 line 203）：** θ_H=+0.0426 处显式回引"established taller-stature→knee-OA epidemiology [24–26]"，使结果层即与已知流行病学对话，直接闭合 R3-①"与已知流行病学一致却未引用"。

---

## 三、R3-② 生物学 OA-骨文献（已引用，本轮核验真实性）

**核验（web 检索）：**
- [27] Goldring & Goldring 2016, *Nat Rev Rheumatol* 12:632–644, PMID 27652499, doi:10.1038/nrrheum.2016.148 — "Changes in the osteochondral unit during osteoarthritis"（骨-软骨单元权威综述，IF≈32）。✅
- [28] Burr & Gallant 2012, *Nat Rev Rheumatol* 8:665–673 — "Bone remodelling in osteoarthritis"。✅
- [29] Zhen & Cao 2013, *Nat Med* 19:704–712 — "Targeting TGFβ signaling in subchondral bone"（软骨下骨 TGFβ）。✅
- [30] Findlay & Kuliwaba 2016, *Bone Res* 4:16028 — "Bone–cartilage crosstalk"。✅

§4.3 已用上述文献展开 osteochondral unit / subchondral bone / osteophyte / TGFβ–WNT crosstalk 机制对话。R3-② 实质已满足（前轮写入，本轮确认文献真实、权威，非编造）。

---

## 四、R3-③ 首个 OA→fracture 直接检验——核验后限定（本轮核心编辑）

**检索发现（关键）：** 2025 年 **Tong et al.**, *Aging Clin Exp Res* 37:305, doi:10.1007/s40520-025-03219-w, PMID 41148443 已做 OA→**足/踝骨折** MR（"knee OA elevates the risk of… foot fractures"）。此为**部位特异性骨折亚型**，与本文 GEFOS **any-fracture** 表型不同。广义检索（"OA as exposure → fracture MR"）未见 OA→任一骨折（非足踝亚型）的先例；Alhassan 2024 / Liu 2025 等 OA-MR 综述仅在共病反向因果中提及 OA 作暴露，未测 OA→fracture。

**原稿风险（R3-③ 警示的 overclaim）：** "to our knowledge the OA→fracture direction has not previously been tested directly" 未区分骨折表型，易被挑刺。

**本轮编辑（§4.4 line 229）：**
> "to our knowledge no prior MR has tested OA as a direct genetic exposure for any-fracture (the GEFOS phenotype examined here). A 2025 MR examined OA→foot/ankle fracture subtypes [31] but used a site-specific fracture endpoint and did not address the systemic BMD–fracture axis investigated here, so the present analysis is the first direct test of OA→any-fracture."

→ 既坐实了"首个 any-fracture 直接检验"的主张，又诚实标注最近邻先例（足踝亚型），消除 overclaim 风险。新增 ref [31]。

---

## 五、R3-④ 骨尺寸混淆锚点增强（本轮）

§4.2 line 219 原仅 [18] Kröger 1995。审查建议补更强 areal-vs-volumetric / bone-size 实证文献。
**新增 [32]** Looker AC, Beck TJ, Orwoll ES. *J Bone Miner Res* 2001;16(7):1291–1299, doi:10.1359/jbmr.2001.16.7.1291（实证 aBMD 受骨尺寸混淆，web 核验真实）。改为 [18,32] 同引。

---

## 六、R3-⑤ coloc p12=5×10⁻⁶ 选择理由（本轮）

§2.9 line 109 原仅给出 p12=5×10⁻⁶ 而未说明理由（审查 flagged）。补充：
> "p12 设略高于常规 1×10⁻⁶ 下界，因 BMD–OA 区域多基因且 LD 密集，过严 p12 易漏掉真实共享变异；更宽松的 p12 只使 H4 更易达成，故我们报告的 null colocalisation 在此选择下反而保守；±250 kb–±1 Mb / p12 1×10⁻⁵–1×10⁻⁷ 敏感性扫描确认 null 非先验人为产物。"

（注：本文 coloc 结论为 null，宽松 p12 使 null 更保守，理由自洽。）

---

## 七、一致性终审

- grep "has not previously been tested directly" → 0 处（旧 overclaim 已清除）
- [31] 出现于 §4.4 + 参考文献；[32] 出现于 §4.2 + 参考文献
- 引用编号 1–32 连续，无断号
- 所有 β/SE/p/OR 未动，仅补文献与限定表述
- R3-①/② 现有引用经 web 核验**全部真实**（PMID 匹配），无编造风险

## 八、提交
- commit （见 git log）：R3 文献对话补强（docs/08 + 99g 审计）
