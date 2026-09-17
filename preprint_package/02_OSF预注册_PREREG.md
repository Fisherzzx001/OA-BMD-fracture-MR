# OSF 预注册文档（Retrospective Preregistration）

> 用途：上传至 OSF（osf.io）取得时间戳与 DOI，作为研究优先权与透明度的佐证。
> **诚实性声明（必须保留）**：本分析在撰写本预注册时**已经完成**。本文档是对已完成分析计划的**回溯性预注册（retrospective preregistration）**，用于公开完整分析口径、支持结果复核与研究优先权，**不具有验证性（confirmatory）注册的效力**。稿件中所有基于本计划的结果均已在正文中完整报告。
> ★ 标注条目为草拟，请作者确认。

---

## 1. Title

Genetic architecture and causal direction of the osteoarthritis–bone mineral density–fracture triangle: a cross-trait LDSC and two-sample Mendelian randomisation study of 13 skeletal phenotypes

## 2. Contributors ★

Zhenxin Zheng¹, Song Chen¹, Shenbo Huang¹, Chao Qian¹, Hang Yuan¹, Huricha Bao¹\*, Qiuyan Jiang²\*
（¹ Orthopedics, ² Surgery — The Quzhou Affiliated Hospital of Wenzhou Medical University, Quzhou People's Hospital, Quzhou, China；\* corresponding contributors: baohuricha@126.com / after9951023@gmail.com）

## 3. Description（研究问题）

Osteoarthritis (OA), bone mineral density (BMD) and fracture are epidemiologically entangled, but whether their shared genetic architecture reflects causal pathways is unresolved. We ask three questions:

1. Do OA, BMD and fracture share genetic architecture (LDSC genetic correlations)?
2. Does genetically predicted BMD causally influence fracture risk and OA risk (two-sample MR)?
3. Does genetically predicted OA causally influence BMD or fracture (reverse MR; directionality)?

## 4. Hypotheses（预先声明）

- **H1（阳性对照）：** 遗传预测的较高 BMD 降低任意骨折风险（BMD→fracture β < 0；外部基准：GEFOS 每降 1 SD 股骨颈 BMD，OR=1.55，95%CI 1.48–1.63）。
- **H2：** OA 与 BMD 正向遗传相关（r_g > 0）。
- **H3：** 遗传预测的较高 BMD 因果性升高膝 OA 风险（BMD→OA β > 0）。
- **H4（方向性检验）：** 反向 MR（OA 作为暴露）对 BMD 与骨折均无因果效应（β 不显著）。
- **H5（补充）：** OA 与骨折之间无显著遗传相关（r_g ≈ 0，预期零结果）。

## 5. Data sources（公开汇总统计，全部可得）

| 表型 | 来源 | Accession / 标识 | N |
|---|---|---|---|
| Knee OA（影像学） | GO Consortium 2025（*Nature* 641:1217–1224） | PMID 40205036 | 1,208,397（膝 OA 子集） |
| Self-reported OA | IEU OpenGWAS | `ukb-b-14486` | 462,933（38,472/424,461） |
| BMD ×10 部位（DXA） | Qian 2025（*Genomics Proteomics Bioinformatics*） | GWAS Catalog GCST90568448–GCST90568458 | 30,449–32,017（欧洲 UKB，逐部位见稿件 §2.1） |
| Any-fracture | GEFOS 2018（Trajanoska, *BMJ* k3225） | PMID 30158200 | 264,973（37,857/227,116） |

祖先限定：欧洲血统。LD 参考面板：1000 Genomes Phase 3 EUR。

## 6. Analysis plan（分析口径，与稿件完全一致）

1. **QC 与 munge**：统一构建（CHR/BP/A1/A2/MAF/N/Z/P）；SNP 键合采用 chromosome + base-pair position（不使用 rs_id）；MAF 0.01–0.99。
2. **遗传相关**：LDSC（S-BOLT/LDSC 实现为 ldsc.py 流程），78 个表型对（13×13 去自相关），多重校正采用 Bonferroni 与 Benjamini–Hochberg（FDR q<0.05）；截距 1.01–1.07 视为无残余分层。
3. **工具变量**：暴露 GWAS 中 p<5×10⁻⁸；**主分析采用严格 PLINK 等价 clumping（1000G EUR，r²<0.001，10 Mb 窗口）**；1-Mb 近似窗口作为敏感性分析。
4. **两样本 MR**：IVW（固定效应，结局方差加权）为主估计；MR-Egger（截距检验水平多效性）、加权中位数、Cochran's Q/I²、radial 离群剔除（|z|>3）、留一法为稳健性诊断；**MR-PRESSO 全局（distortion）检验**（Bonferroni α=0.05/K 逐 SNP 离群检验）。
5. **效应量标度重构**：munge 数据仅含 Z/N/MAF，按 SE=1/√(N·2p(1−p))、β=Z·SE 重构。**N 标度校验为预注册时已知的必要步骤**（见 §7 偏差说明）。
6. **外部校验**：BMD→fracture 估计与 GEFOS 基准（OR=1.55/SD）比对，一致性带宽预设为"方向一致且 OR 量级在 2 倍以内"；实测落入 8–10% 差异。

## 7. Deviations from the plan（偏差披露，预注册诚实性关键项）

1. **N 标度校正在分析过程中发现并实施**：BMD munge 文件携带占位 N=400,000，真实逐部位 N 为 30,449–32,017（GWAS Catalog 已核验）。重构 β 因此被低估 √(N真/400,000) 倍（≈3.54×），已统一校正（对 headline OR 的影响 <0.3%；lumbar 最大 2.3%）。该校正只改变效应量绝对值，z/p/方向/显著性为 N 不变量。
2. **OA_selfrep N 校正**：munge 文件占位 N=350,000 → 真实 462,933（因子 ×0.8695），同样只影响效应量绝对值。
3. 除上述标度事项外，无其他偏离：clumping 策略、MR 方法组合、多重校正与假设集均按 §6 执行。

## 8. Exploratory（未预设的探索性分析，已在稿件中标注）

- 图 4 的易感性尺度 h² 转换（Lee 2011 Eq 23，K 取文献病例对照比例）；
- 图 5 因果网络概念图；
- 逐 SNP 散点/漏斗/留一法诊断图（图 6、图 7、补充图 S1–S2）。

## 9. Timeline ★

- 计划上传 OSF 并取得 DOI：稿件投稿（*Bone* / *Osteoporosis International*）之前或同日；
- 计划上传 medRxiv 预印本：OSF 注册后 7 日内（取得优先权时间戳）。
