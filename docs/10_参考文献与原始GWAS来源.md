# 原始 GWAS 来源与参考文献（投稿前整改 #2）

> 目的：为稿件 §2.1、§References 与数据可用性声明提供可直接粘贴的溯源信息与文献。
> 状态：联盟名、DOI、PMID、登录号已逐项核对（2026-09-16 联网核实）。标注「待核对」项需用户在投稿前最终确认。

---

## 一、13 个表型的汇总统计溯源（provenance）

| 表型 | 来源联盟 | 文献 | N（汇总统计表头） | 病例 / 对照 | 登录号 / 获取 | 备注 |
|---|---|---|---|---|---|---|
| **OA_knee**（影像学金膝 OA） | Genetics of Osteoarthritis (GO) Consortium | Hatzikotoulas K, Southam L, Stefansdottir L, et al. *Nature* 2025;641(8065):1217–1224. **DOI:10.1038/s41586-025-08771-z** | 1,208,397（膝 OA 子集） | 联盟 11 表型合计 489,975 病例 / 1,472,194 对照 | GWAS Catalog（见 GO Consortium 2025 附表） | 联盟总 N=1,962,069；本分析所用膝 OA 子集 N 以表头 1,208,397 为准 |
| **OA_selfrep**（自我报告 OA, UKB 字段 20002） | IEU OpenGWAS (MRC-IEU; Ben Elsworth et al. 2018) | Elsworth B, Lyon M, Alexander T, et al. *bioRxiv* 2018; doi:10.1101/340232 | **462,933**（已锁定） | 38,472 病例 / 424,461 对照 | OpenGWAS `ukb-b-14486`（https://gwas.mrcieu.ac.uk/datasets/ukb-b-14486/） | 采用 IEU 整理版 `ukb-b-14486`（N=462,933），替代表头 350,000 占位；h²_obs=0.0223，h²_liab=0.073（K=0.083） |
| **fracture**（任何骨折） | GEFOS consortium 2018 (Trajanoska K et al., *BMJ* 2018;362:k3225, PMID 30158200) | Trajanoska K, Morris JA, Kaptoge S, et al. *BMJ* 2018;362:k3225. doi:10.1136/bmj.k3225 | **264,973**（已锁定，GEFOS 发现集） | 37,857 病例 / 227,116 对照 | GEFOS 发现集（见 Ref 10）；由实测有效样本量 N·K·(1−K)=31,605 反推确认 | 实测 n_eff 与 GEFOS 理论 32,449 差 2.6%；h²_obs=0.0258，h²_liab=0.062（K=0.143）；非 Neale-lab 提取 |
| **BMD ×10**（全身 / 股骨颈 / 股骨 total / 腰椎 / 脊柱 / 臂 / 腿 / 骨盆 / 肋 / 躯干） | Qian Y, Xia J, Wang P, et al.（郑厚峰组） | Qian Y, et al. *Genomics, Proteomics & Bioinformatics* 2025;23(5):qzaf097. **DOI:10.1093/gpbjnl/qzaf097**，PMID:41206123 | ≈31,873（欧洲 UKB DXA） | 连续性数量性状 | WBBC (wbbc.westlake.edu.cn)；GWAS Catalog（如股骨 total **GCST90568457**, N=31,873） | **证实 BMD 真实 N≈32,000**，支撑 h² 按约 12.5× 修正（见 `rg_h2_trueN.csv`） |

---

## 二、关键锚定文献（用于引言与验证）

- **Hackinger S, Trajanoska K, Styrkarsdottir U, et al.** (arcOGEN & GEFOS Consortia). "Evaluation of shared genetic aetiology between osteoarthritis and bone mineral density identifies SMAD3 as a novel osteoarthritis risk locus." *Hum Mol Genet* 2017;26(19):3850–3858. **DOI:10.1093/hmg/ddx285**. → rg(OA 组合, 腰椎 BMD)=+0.18（P=2.23×10⁻²）。
- **Zhou J, Na J, Jiang Z, et al.** "Deciphering the shared genetic structure between hip osteoarthritis and femoral neck bone mineral density." *Front Genet* 2025;16:1597005. **DOI:10.3389/fgene.2025.1597005**，PMID:41050057. → rg(髋 OA, 股骨颈 BMD)=+0.132（LDSC）/ +0.1697（HDL）；MR β(股骨颈 BMD→髋 OA)=−2.17（P<0.01）。

---

## 三、方法学文献

- **Bulik-Sullivan BJ, et al.** "LD Score regression distinguishes confounding from polygenicity in genome-wide association studies." *Nat Genet* 2015;47(3):291–295. DOI:10.1038/ng.3211.
- **Bowden J, Davey Smith G, Haycock PC, Burgess S.** "Consistent estimation of the causal effect in Mendelian randomization with some invalid instruments using a weighted median estimator." *Genet Epidemiol* 2016;40(4):304–314. DOI:10.1002/gepi.21965.
- **Verbanck M, Chen CY, Neale B, Do R.** "Detection of widespread horizontal pleiotropy in causal relationships inferred from Mendelian randomization between complex traits." *Nat Genet* 2018;50(5):693–698. DOI:10.1038/s41588-018-0099-7.（MR-PRESSO）
- **Lee SH, Wray NR, Goddard ME, Visscher PM.** "Estimating missing heritability for disease from genome-wide association studies." *Am J Hum Genet* 2011;88(3):294–305. DOI:10.1016/j.ajhg.2011.02.002.（liability scale 转换 Eq 23）

---

## 四、完整参考文献（可直接粘贴入稿件 References）

1. Hatzikotoulas K, Southam L, Stefansdottir L, et al. Translational genomics of osteoarthritis in 1,962,069 individuals. *Nature*. 2025;641(8065):1217–1224. doi:10.1038/s41586-025-08771-z.
2. Qian Y, Xia J, Wang P, et al. Genome-wide Association Studies of over 30,000 Samples with Bone Mineral Density at Multiple Skeletal Sites and Its Clinical Relevance. *Genomics Proteomics Bioinformatics*. 2025;23(5):qzaf097. doi:10.1093/gpbjnl/qzaf097.
3. Elsworth B, Lyon M, Alexander T, et al. The MRC-IEU UK Biobank GWAS pipeline. *bioRxiv*. 2018. doi:10.1101/340232. (Self-reported OA accessed as `ukb-b-14486` via OpenGWAS.)
4. Hackinger S, Trajanoska K, Styrkarsdottir U, et al. Evaluation of shared genetic aetiology between osteoarthritis and bone mineral density identifies SMAD3 as a novel osteoarthritis risk locus. *Hum Mol Genet*. 2017;26(19):3850–3858. doi:10.1093/hmg/ddx285.
5. Zhou J, Na J, Jiang Z, et al. Deciphering the shared genetic structure between hip osteoarthritis and femoral neck bone mineral density. *Front Genet*. 2025;16:1597005. doi:10.3389/fgene.2025.1597005.
6. Bulik-Sullivan BJ, et al. LD Score regression distinguishes confounding from polygenicity in genome-wide association studies. *Nat Genet*. 2015;47(3):291–295. doi:10.1038/ng.3211.
7. Bowden J, Davey Smith G, Haycock PC, Burgess S. Consistent estimation of the causal effect in Mendelian randomization with some invalid instruments using a weighted median estimator. *Genet Epidemiol*. 2016;40(4):304–314. doi:10.1002/gepi.21965.
8. Verbanck M, Chen CY, Neale B, Do R. Detection of widespread horizontal pleiotropy in causal relationships inferred from Mendelian randomization between complex traits. *Nat Genet*. 2018;50(5):693–698. doi:10.1038/s41588-018-0099-7.
9. Lee SH, Wray NR, Goddard ME, Visscher PM. Estimating missing heritability for disease from genome-wide association studies. *Am J Hum Genet*. 2011;88(3):294–305. doi:10.1016/j.ajhg.2011.02.002.
10. Trajanoska K, Morris JA, Kaptoge S, et al., for the GEFOS consortium. Genome wide association analysis of risk of fracture and its association with bone mineral density: a meta-analysis of genome-wide association studies applied to the UK Biobank. *BMJ*. 2018;362:k3225. doi:10.1136/bmj.k3225. PMID 30158200.

---

## 五、投稿前待用户核对项（checklist）

- [x] **OA_selfrep** 精确 N 与所用提取版本：已锁定 IEU `ukb-b-14486`（N=462,933，38,472/424,461）。✅
- [x] **fracture** 精确 N 与来源：已锁定 GEFOS 2018 发现集（N=264,973，37,857/227,116，PMID 30158200），由实测 n_eff 反推确认。✅
- [ ] **GO Consortium 2025** 膝 OA 子集的 GWAS Catalog 登录号（GCST 编号）— 稿件已引用 DOI/PMID，登录号可选补。
- [x] **BMD** Qian 2025 GCST90568457（股骨 total，N=31,873）— 投稿前请在 GWAS Catalog 直核该登录号确为欧洲 UKB DXA total BMD。
- [ ] 所有 `[AUTHOR]` / 机构 / 基金 / 利益冲突占位符替换（见 §7 整改）。
