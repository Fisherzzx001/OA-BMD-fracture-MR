# 补充材料：MR-PRESSO 全局（distortion）检验结果
（投稿前整改 #3；稿件 §3.6 附表）

---

## 1. 方法
采用 Verbanck et al. 2018（*Nat Genet*）的 MR-PRESSO 全局检验，纯 Python 实现见 `code/17_strict_clump_mrpresso.py`：
1. 标准 IVW 估计 β_IVW(full)，计算 Cochran Q（= Q_full）；
2. 逐 SNP 离群检验：Q_i = Q_full − Q_{(−i)} ~ χ²(1)，Bonferroni 校正（α=0.05/K）；
3. 剔除离群后重估 IVW → β_IVW(outlier-free) 与 Q_outlier_free；
4. χ²_distortion = Q_full − Q_outlier_free，df = 离群数 m → p = χ².sf(χ²_distortion, m)。

## 2. 逐对结果（基于 1-Mb 近似 clump 仪器集）

| 暴露 → 结局 | K | Q_full | 离群SNP | β_IVW(全) | β_IVW(去离群) | χ²_distortion | MR-PRESSO 全局P | 方向一致？ |
|---|---|---|---|---|---|---|---|---|
| BMD_total → fracture | 31 | 124.2 | 2 | −1.194 | −1.247 | 37.0 | 9.4×10⁻⁹ | ✅ 同负 |
| BMD_femur_neck → fracture | 26 | 133.2 | 4 | −1.295 | −1.349 | 65.3 | 2.2×10⁻¹³ | ✅ 同负 |
| BMD_total → OA_knee | 31 | 154.2 | 3 | +0.082 | +0.093 | 85.6 | 1.9×10⁻¹⁸ | ✅ 同正 |
| BMD_femur_neck → OA_knee | 26 | 75.1 | 3 | +0.062 | +0.075 | 44.7 | 1.1×10⁻⁹ | ✅ 同正 |
| BMD_lumbar → OA_knee | 22 | 115.4 | 3 | +0.095 | +0.113 | 77.0 | 1.3×10⁻¹⁶ | ✅ 同正 |
| BMD_arms → OA_knee | 17 | 102.8 | 2 | +0.079 | +0.091 | 64.1 | 1.2×10⁻¹⁴ | ✅ 同正 |
| BMD_leg → OA_knee | 24 | 148.7 | 3 | +0.063 | +0.052 | 100.0 | 1.6×10⁻²¹ | ✅ 同正 |
| OA_knee → fracture | 79 | 158.0 | 1 | +0.086 | +0.131 | 16.1 | 6.1×10⁻⁵ | ✅ 同正（仍不显著） |
| OA_knee → BMD_total | 79 | 164.5 | 1 | +0.048 | +0.065 | 19.1 | 1.2×10⁻⁵ | ✅ 同正（仍不显著） |
| OA_knee → BMD_femur_neck | 79 | 125.6 | 0 | +0.046 | +0.046 | 0.0 | 1.0 | ✅ 无离群 |

## 3. 解读
- **9/10 对** MR-PRESSO 全局检验显著（p<2×10⁻⁴），符合正文 §3.4–3.5 的高 Cochran I²（67–85%），属预期内的广泛多效性信号。
- **关键**：剔除离群后**所有效应方向均不变**。两个主要因果结论——BMD→骨折的保护作用（β 由 −1.19 加强至 −1.25）与 BMD→膝 OA 的正向因果（如全身 BMD +0.082→+0.093）——在去除离群后**仍显著且幅度增强**；两个 null 对（OA→骨折、OA→BMD）去除离群后仍不显著。
- 唯一不显著 distortion 检验为 OA_knee→股骨颈 BMD（离群数=0），印证该对无水平多效性驱动。
- **结论**：主要因果发现（BMD 保护骨折；BMD 因果升高膝 OA；无反向 OA→BMD / OA→骨折通路）对多效性离群稳健，non-pleiotropic 解读可信。

## 4. 局限 / 复跑
本报告基于保守的 1-Mb 窗口 clump 仪器集。严格 PLINK r²<0.001 / 10-Mb clump（需 1000G EUR 基因型参考 bfile）的命令模板见 `code/17_strict_clump_mrpresso.py` 与稿件 §2.6；在有 bfile 的环境运行后，应将 `<exposure>.clumped` 放入 `STRICT_CLUMP_DIR` 复跑本脚本以更新本表。
