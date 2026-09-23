# OA–BMD–骨折 跨性状遗传研究 · 项目导航

> 题目：*Bone Mineral Density and Osteoarthritis Share Genetic Architecture, but
> Evidence for an Osteoarthritis-to-Fracture Causal Pathway Remains Limited:
> A Cross-Trait Analysis of 13 Skeletal Phenotypes*
>
> 目标期刊：*Bone* / *Osteoporosis International*（现实档），*Osteoarthritis and Cartilage*（冲刺档）

## 1. 研究设计（一句话）

13 个骨骼表型（2 种 OA 定义 + 10 个区域 BMD + 任意骨折）：
**LDSC** 估 78 对遗传相关 → **双向两样本 MR** 判因果方向（含 BMD→骨折内部阳性对照）
→ **MRlap** 校正样本重叠与赢家诅咒 → 异质性/多效性/功效诊断。

## 2. 目录结构（2026-09-23 清理后）

```
OA_B_PROJECT/
├── README.md                    ← 本文件
├── docs/
│   ├── 08_投稿稿_OA_BMD_fracture.md   ← 唯一文本真源（手稿 + 图注 + 表 + 参考文献）
│   ├── 09–14                          ← 数据下载实操、参考文献与 GWAS 来源、补充材料 S1–S4
│   ├── 18–21                          ← 第三/四轮预审报告与修复确认
│   └── 22–24                          ← 第五/六/七轮审稿与从零复盘审计报告
├── code/                        ← 分析管线（01–28）+ 审计脚本（90–92）
│   └── output/genetics/         ← 全部结果 CSV（harmonised/ 为逐 SNP 数据）
├── Figures_v17_Nature/          ← 18 张图的 PDF/PNG/TIFF（含 Figures_Original_v14/ 原图）
├── submission_package_v16/      ← 投稿包工作目录（docx + 图 + 表 + 补充材料 + MANIFEST）
├── archive/                     ← 已归档的历史版本与临时产物（见 archive/README.md，可删）
├── build_*.py                   ← 构建脚本（submission / maintext / cover_letter / figures / package）
├── audit_*.py                   ← 独立重算审计脚本（v18 全量 430 项，0 失败）
└── Submission_Package_v16.zip   ← 最终投稿包（18 MB，98 个条目）
    Figures_Package_v17.zip      ← 图件单独下载包（14.5 MB，含原图）
    Preprint_Prereg_Package.zip  ← 预印本 + OSF 回溯预注册材料
```

## 3. 一键重建

```bash
PY=<venv python>
OUTDIR=submission_package_v16 $PY build_submission.py      # Manuscript.docx
OUTDIR=submission_package_v16 EMBED_FIG=1 $PY build_submission.py
                                                            # + Manuscript_with_Figures.docx
OUTDIR=submission_package_v16 $PY build_maintext.py        # 盲法主文本
OUTDIR=submission_package_v16 $PY build_cover_letter.py    # 投稿信
$PY build_figures_v17.py                                    # 重制 18 张图
$PY build_package_v16.py                                    # 组装投稿包并打包
$PY audit_v18_full.py                                       # 全量独立复算审计
```

## 4. 审稿状态

- 学术预审委员会第七轮（从零全面复盘）：**430 项独立重算 0 失败**，发现并修复 16 处表述层问题，
  核心统计量（β/SE/p/OR/r_g/h²/功效）全部独立复现，无数据错误。详见 `docs/24_*`。
- 图件：845 项图—数据一致性检查通过（`Figures_v17_Nature/AUDIT_REPORT.txt`）。

## 5. 投稿前待办

- [ ] 补 ref 16（Alhassan 2023）与 ref 18（Kröger 1995）的卷/期/页
- [ ] 选定目标期刊并改 `build_cover_letter.py` 中 `JOURNAL`，重跑投稿信
- [ ] GitHub 仓库公开 + Zenodo DOI 回填 Data availability（PAT / token 步骤见 `code/README.md`）
- [ ] 确认无误后删除 `archive/`（约 340 MB 历史产物）

## 6. 已如实披露的局限（不要弱化）

1. 三个无效 MR 对（OA→BMD×2、OA→骨折）功效仅 0.14–0.26，且对工具变量选择敏感
   （p 在 0.013 / 0.097 / 0.326 之间变动），只能排除大效应，不能排除中等效应。
2. 随机效应 IVW 下五对 BMD→OA 只有 2 对仍 p<0.05，因果主张不依赖随机效应模型。
3. BMD 效应尺度为重建尺度（源 GWAS 用占位 N=400,000，真实 31,873），绝对量级已按
   √(31,873/400,000) 校正，方向与显著性不受影响；与 GEFOS 基准相差 8–10%。
4. UKB 内部样本重叠无法完全剥离，已用 MR-Egger / 加权中位数 / MR-PRESSO / ρ_GE / MRlap 五种方式吸收。
