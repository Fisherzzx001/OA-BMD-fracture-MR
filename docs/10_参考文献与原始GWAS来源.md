# 原始 GWAS 来源与参考文献（投稿前整改 #2）

> 目的：为稿件 §2.1、§References 与数据可用性声明提供可直接粘贴的溯源信息与文献。
> 状态：联盟名、DOI、PMID、登录号已逐项核对（2026-09-16 联网核实）。标注「待核对」项需用户在投稿前最终确认。

---

## 一、13 个表型的汇总统计溯源（provenance）

| 表型 | 来源联盟 | 文献 | N（汇总统计表头） | 病例 / 对照 | 登录号 / 获取 | 备注 |
|---|---|---|---|---|---|---|
| **OA_knee**（影像学金膝 OA） | Genetics of Osteoarthritis (GO) Consortium | Hatzikotoulas K, Southam L, Stefansdottir L, et al. *Nature* 2025;641(8065):1217–1224. **DOI:10.1038/s41586-025-08771-z** | 1,208,397（膝 OA 子集） | 联盟 11 表型合计 489,975 病例 / 1,472,194 对照 | GWAS Catalog（见 GO Consortium 2025 附表） | 联盟总 N=1,962,069；本分析所用膝 OA 子集 N 以表头 1,208,397 为准 |
| **OA_selfrep**（自我报告 OA, UKB 字段 20002） | Neale Lab UK Biobank Round 2 | Neale Lab. "UK Biobank Round 2 GWAS", 2018（发布页 nealelab.is/uk-biobank） | 350,000（表头） | 待核定 | OpenGWAS `ukb-d-*`；IEU 整理版 `ukb-b-14486`（N=462,933，病例 38,472） | 表头 N=350,000 与 IEU 版 462,933 不一致 → **投稿前核对所用提取版本** |
| **fracture**（任何骨折） | Neale Lab UKB Round 2（原始格式） | Neale Lab, 2018 | 400,000（表头，待核定） | 待核定 | OpenGWAS `ukb-d-*`；备选 GEFOS consortium 大样本骨折 GWAS（N≈416,795，病例 53,184） | λ_gc≈1.056 提示真实 N 约 6–18 万（中心 ~9 万），待核对；见 §5 |
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
3. Neale Lab. UK Biobank Round 2 GWAS. 2018. https://nealelab.is/uk-biobank (summary statistics via OpenGWAS, `ukb-d-*`).
4. Hackinger S, Trajanoska K, Styrkarsdottir U, et al. Evaluation of shared genetic aetiology between osteoarthritis and bone mineral density identifies SMAD3 as a novel osteoarthritis risk locus. *Hum Mol Genet*. 2017;26(19):3850–3858. doi:10.1093/hmg/ddx285.
5. Zhou J, Na J, Jiang Z, et al. Deciphering the shared genetic structure between hip osteoarthritis and femoral neck bone mineral density. *Front Genet*. 2025;16:1597005. doi:10.3389/fgene.2025.1597005.
6. Bulik-Sullivan BJ, et al. LD Score regression distinguishes confounding from polygenicity in genome-wide association studies. *Nat Genet*. 2015;47(3):291–295. doi:10.1038/ng.3211.
7. Bowden J, Davey Smith G, Haycock PC, Burgess S. Consistent estimation of the causal effect in Mendelian randomization with some invalid instruments using a weighted median estimator. *Genet Epidemiol*. 2016;40(4):304–314. doi:10.1002/gepi.21965.
8. Verbanck M, Chen CY, Neale B, Do R. Detection of widespread horizontal pleiotropy in causal relationships inferred from Mendelian randomization between complex traits. *Nat Genet*. 2018;50(5):693–698. doi:10.1038/s41588-018-0099-7.
9. Lee SH, Wray NR, Goddard ME, Visscher PM. Estimating missing heritability for disease from genome-wide association studies. *Am J Hum Genet*. 2011;88(3):294–305. doi:10.1016/j.ajhg.2011.02.002.
10. Morris JA, et al. / GEFOS consortium. Genome-wide association study of fracture. *[待补 GEFOS 骨折 GWAS 精确题名与 DOI]* — 如改用 GEFOS 大样本骨折数据重跑，请补此条。

---

## 五、投稿前待用户核对项（checklist）

- [ ] **OA_selfrep** 精确 N 与所用提取版本（Neale R2 `ukb-d` vs IEU `ukb-b-14486`）：影响 h² 绝对值，不影响 r_g / MR 方向。
- [ ] **fracture** 精确 N 与来源（Neale R2 vs GEFOS）：决定骨折 h²（当前 0.017–0.213 区间，见 `fracture_h2_sensitivity.csv`）。
- [ ] **GO Consortium 2025** 膝 OA 子集的 GWAS Catalog 登录号（GCST 编号）。
- [ ] 是否将 fracture 改用 **GEFOS** 大样本（N 更大、更权威）重跑并替换当前 Neale-lab 提取。
- [ ] 所有 `[AUTHOR]` / 机构 / 基金 / 利益冲突占位符替换（见 §7 整改）。
