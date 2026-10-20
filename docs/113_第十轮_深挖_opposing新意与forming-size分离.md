# 第十轮深挖记录 · opposing-endpoints 新意锚定 + bone-forming/size 未分离

> 对象：第九轮预审报告（`docs/112_第九轮预审报告_OAC.md`）R1「理论贡献质询师」致命漏洞 2 与 3
> 日期：2026-10-16 · 落地补丁：`code/126_anchor_novelty_and_sizemass.py`
> 验证：`audit_v18_full.py` 482 项 0 失败（新增 G10/G11 门禁）；`code/115_audit_figures_v18.py` 189 项 0 失败

---

## 一、漏洞 3：「opposing endpoints」新意被高估

### 核验（外部文献）

| 待核验命题 | 结论 | 依据 |
|---|---|---|
| BMD–OA 正相关是已知概念 | **是，40+ 年** | Foss & Byers 1972（文中 [11]）；Dequeker 2003（文中 [12]）综述原话："While predisposing to articular cartilage loss, bone alterations associated with GOA may provide protection against osteoporotic fractures" —— 「同一骨表型既促软骨流失又防骨质疏松骨折」的对立端点框架在文献中已存在四十年 |
| BMD–fracture 负相关是已知概念 | **是** | GEFOS 2018（文中 [10]）及全部骨质疏松遗传学基础文献 |
| 流行病学量级 | **早已确立** | Rotterdam 研究（Bone 2005）：最高 FN-BMD  quartile 膝 ROA 发病 OR 2.8（1.2–6.8）、LS-BMD OR 4.7（2.1–10.7） |
| 是否有先例用**同一套 BMD 遗传工具同时估计 OA 与 fracture** | **未发现** | 系统检索（PubMed/网络）无对应 MR 研究；唯一相关沉积为本人前期 Zenodo five-layer funnel（H1–H5，eBMD 口径），非同行评议、非三角联合设计 |

**结论**：审稿人判断正确——「对立端点」概念本身零新意；论文的真实增量是三件可检验的事：**三角联合估计**（同一工具集、同一校准管线、骨折端点为内部阳性对照）、**位点特异性**（10 位点）、**定量分解**（size/mass 部分划分）。文中 "the OA–BMD–fracture triangle has not been analysed jointly" 的 gap 声明成立且已被检索支持。

### 落地

- **Intro gap 3**（L36）：补 "although each association is individually familiar [11,12,33]" —— 明示概念古老、联合分析才是增量。
- **Discussion L82** 段首改为锚定句："The opposing-endpoint pattern is an old clinical observation [11,12]; what is new is its demonstration for a single genetic instrument set within one calibrated pipeline."（删去原 "The paradox resolves…" 的概念叙事开头及与 L88 重复的 subchondral/osteophyte 句）
- 摘要/H1 的结果式陈述（"protects against fracture yet raises knee-OA risk"）保留——那是结果陈述而非概念声明。

---

## 二、漏洞 2：「bone-forming」与「bone-size/shape」set-point 未分离

### 可行性判定（为什么没有补定量鉴别检验）

审稿人建议方向为严格术语 + Faber [35] 定性支撑；本轮进一步核查了「用现有数据补鉴别性检验」的可能，结论是**不可行，且原因是结构性的**：

1. **数据面**：本地仅有 height GWAS 全量汇总统计（Yengo 2022, 79 MB）；BMD/OA/fracture 均无全基因组 sumstat（只有 clump 输入、instrument 级与 tophits 文件）→ 无法估计 r_g(BMD_site, height)，也无法做 height→OA 的中介分解。
2. **识别面（更根本）**：MVMR 残差是「条件化 height 后的 areal BMD」。areal BMD = BMC/投影面积，**条件化身高后仍编码骨宽与几何**（身高不投影骨横断面形态）。因此该残差在定义上就是 size-or-mass 的混合体，**任何本地数据都无法将其鉴别为 mass**；能做鉴别的只有 F2（以 vBMD/骨面积/髋形状模式为第二暴露条件化）。
3. **伪检验排除**：曾考虑「跨位点衰减率 vs 位点-身高遗传相关」相关检验——已论证其为 MVMR 代数必然（β_uni = θ_BMD + θ_H·γ_BH，两假设下衰减都 ∝ r_g(BMD,height)），非鉴别性，故**拒绝补做**（补假检验正是审稿人批评的叙事包装）。

Faber [35] 已核实（Hum Mol Genet 2025;34(3):207–217, PMID 39574169）：203 个髋形状独立信号；MR 显示**形状对髋 OA 无因果效应**，而 HSM2（颈干角更钝）/HSM4（股骨颈更宽）对骨折有因果效应——同时支撑两点：几何是可分离遗传轴（身高只是其粗糙代理），且几何不导致 OA（与残差非纯几何一致，但不构成证明）。

### 落地（术语 + 边界 + 定性支撑）

| 位置 | 修改 |
|---|---|
| Highlights H4 | "only 9.5–33.2% explained by skeletal size" → "9.5–33.2% … attributable to the **height-proxied** size pathway"（82 字符，合规） |
| Highlights H6 | "explains the opposing endpoint consequences" → "remains a **hypothesis consistent with** these data"（第九轮 P0-3 遗留项同句落地，82 字符） |
| 摘要 | "independent of the height–size pathway" → "height-independent, a residual **not resolved into bone mass versus geometry**" |
| Results L76 | 残差改述 + "height indexes stature but not bone width or geometry [35] … not identified as bone mass; the mass-versus-geometry split is resolved only by test F2" |
| Discussion L82 | 9.5–33.2% 明确为 **size 归因下界**；残差非 mass；阳性对照 height-transparent + Faber 形状→骨折/不→OA 作定性支撑 |
| Discussion L88 | SOS operationalisation 句重写："the name states the hypothesis, not an established mechanism"；身高划分双向不完整 |
| 局限性 L94 | "bounds only the height-proxied share of size … a lower bound, since bone width and geometry form a separable genetic architecture [35]" |
| 结论 L96 | "independent of skeletal size" → "independent of height (not resolved into bone mass versus skeletal geometry)" |
| S10 图（c） | 标题 "Share of the BMD effect explained by skeletal size" → "removed by conditioning on height"；右列表头 → "Height-independent residual"（已重生成，图审通过） |
| Table S9 脚注 | "Size-attributable" → "Height-attributable … a lower bound for total size"；"(mass effect intact)" → "(control intact)" |
| S11（Table S11） | 核心命题补「命名即假说」段；P4 判定改为「 Supported for the height-proxied share; residual not identified as mass」；F2 补几何候选第二暴露（QCT vBMD/骨面积/髋形状模式 [35]）；诚实边界补 height 粗糙代理句 |
| 封面信 finding 3 | 同义严格化（保留 G1 门禁要求的 9.5–33.2% / 67–91% 字样） |

---

## 三、防回退门禁（audit_v18_full.py §6.7 新增）

- **G10 过度断言禁用**（稿件 + cover letter + Table S11 三文件）：`mass-attributable`、`explains the opposing`、`independent of skeletal size`、`explained by skeletal size`、`mass effect intact`、`true mass effect survives` —— 18 项检查
- **G11 必备披露**：`height-independent` 术语、`bone width`/`geometry` 披露、`lower bound` 下界表述、S11 同步 —— 4 项检查
- 基线 430 → 482 项检查，0 失败；`code/115_audit_figures_v18.py` 189 项 0 失败

## 四、已知遗留（既有问题，未越界处理）

1. **正文 4734 词 vs OAC 上限 4000（超 734；本轮 +53）**。建议下一步走 `manuscript-word-limit-reduction`；L86 Hartley 量级段（约 330 词）可按 R3 建议压 ~1/3（约 −110 词）抵偿本轮增量——该项属第九轮 P2，待用户裁决后执行。
2. **摘要 261 词（限 250）**；**Highlights 6 条（应 3–5）且 H3（145 字符）/H5（97 字符）超 85 上限**——H3/H5 的压缩方案已备好，待裁决。
3. `Submission_Package_v16.zip` / `Figures_Package_v17.zip` / `submission_package/Cover_Letter.docx` 仍为旧版（封面信 docx 的 TITLE 与现行稿件不一致），投稿前统一重打。
