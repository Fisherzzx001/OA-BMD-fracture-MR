# OA 共病网络双库研究 · B 档项目导航

> 目标：复刻蒸馏论文的"双库互补 + 匹配对照 + 共病网络平均度"三板斧到 OA，
> 补上"可比性四件套"等方法学升级，发表 5–10 分 SCI（Osteoarthritis and Cartilage / ACR / BMC Medicine）。

## 0. 当前进度（截止 2026-09-08）

- [x] 立项 / 预注册文本草稿（`docs/01_*`）
- [x] 变量协调映射表（`docs/02_*`）
- [x] 完整可运行分析管线 + 模拟数据验证（`code/`）
- [x] 端到端跑通：H1 显著（UKB P=0.003，CHARLS P=0.086 边缘），H2 检出代谢模块 Δ O/E +0.066, P=0.011
- [x] 8 项敏感性分析代码实现并跑通（S1–S4/S8 落地，S5–S7 待真数据）
- [ ] **真数据接入**（CHARLS 下载 + UKB 申请）
- [ ] OSF 正式提交
- [ ] 论文初稿

## 1. 目录结构

```
OA_B_PROJECT/
├── README.md                           ← 本文件
├── docs/
│   ├── 01_研究方案与OSF预注册.md       ← 预注册文本（先提交这版）
│   ├── 02_变量协调映射表_CHARLSxUKB.md ← 取数时按此表对号入座
│   ├── 03_结果图表模板与论文写作框架.md
│   ├── 04_执行时间表与风险清单.md
│   └── 05_UKB申请材料.md               ← 可直接粘贴的 UKB 申请文本
└── code/
    ├── config.py            ← 改这一处即可切换分析口径
    ├── fields.py            ← 字段候选清单 + ICD-10 编码集（取数时回填）
    ├── simulate_data.py     ← 患病率标定后的模拟数据生成器
    ├── lib_network.py       ← O/E、FDR、等密度、置换、bootstrap 全向量化
    ├── 01_build_cohort.py   ← 纳排 + 1:1 匹配 + SMD
    ├── 02_harmonise.py      ← 原始数据协调（CHARLS .dta / UKB .tsv）
    ├── 03_network_analysis.py  ← 描述+网络+推断
    ├── 04_make_outputs.py   ← 论文级图件
    ├── 05_sensitivity.py    ← 8 项预注册敏感性分析
    ├── run_all.py           ← 一键跑通
    └── output/              ← 跑通后的表格 / 图件 / 中间 CSV（模拟数据）
```

## 2. 真数据上线（最短路径）

1. **CHARLS**：注册并下载 2011/2013/2015/2018/2020 五波 + 2011/2015 血检
   - 注册：[charls.charlsdata.com](https://charls.charlsdata.com)
   - 字段对照按 `docs/02_变量协调映射表_CHARLSxUKB.md` 中 ⚠️ 项逐一核实
2. **UK Biobank**：以学生身份提交申请（£500 通道），3–6 周回执
   - 申请入口：[ams.ukbiobank.ac.uk](https://ams.ukbiobank.ac.uk)
3. **把协调后的两表 CSV 放到 `code/data/`**：
   - `charls.csv` 与 `ukb.csv`，列名见 `config.py` 中 `REQUIRED_COLUMNS`
   - 或运行 `python 02_harmonise.py --charls path/to/charls.dta --ukb path/to/ukb.tsv`
     自动按 `fields.py` 候选模式匹配并生成字段核对报告
4. **运行真实数据**：
   ```bash
   cd code
   python run_all.py --real
   python 05_sensitivity.py --real
   ```
5. **对照模拟结果**：`code/output/` 下的表格/图件以同样的文件名覆盖。差异处逐条核查。

## 3. 真数据 vs 模拟数据的预期差异（提前预判）

| 指标 | 模拟 | 真数据预期 | 处理 |
|---|---|---|---|
| CHARLS H1 P | 0.086 | 0.001–0.05（n 更大） | 直接替换 |
| UKB H1 P | 0.003 | < 0.001 | 直接替换 |
| CHARLS OA 节点数 | 33 边 | 25–60 边（取决于匹配后 n） | K 自动重算 |
| 跨库代谢 Δ | +0.07, P=0.011 | **方向不变**，量级待定 | 写论文 |
| 肌骨模块 Δ | +0.00, P=0.93 | 可能在英国侧显著，取决于 CHARLS 是否有 osteoporosis | 若 CHARLS 无 osteoporosis，见 §4 |

## 4. CHARLS 缺骨质疏松字段的应急方案

`02_*` 注 2 已说明：若 CHARLS 实际无 osteoporosis 字段，主分析 CNS 降为 16 节点，
骨-关节模块由"慢性背颈痛 + 多部位关节痛"代表。
讨论中必须把"中国老年人骨质疏松诊断率低"作为**解释性发现**（而非"数据缺陷"）写进 Limitations。

## 5. 跑通验证日志（模拟数据）

```
H1 OA vs control:
  CHARLS  Δ = +2.50, permutation P = 0.08596
  UKB     Δ = +2.88, permutation P = 0.002999
等密度 K = 33 edges
跨库模块差异（O/E 强度，CN − UK）:
  metabolic  Δ = +0.066 [+0.015, +0.117] P=0.011  ✓ 检出
  cardiovascular / respiratory / organ_other / mental / musculoskeletal
            Δ 不显著（CI 含 0）
```

## 6. 命令速查

```bash
# 模拟数据跑通（验证管线）
python run_all.py

# 真实数据
python run_all.py --real

# 仅做描述层（草稿阶段反复试）
python 01_build_cohort.py --real
python 03_network_analysis.py

# 提升生产配置（修改 config.py 中 n_permutation=10000, n_bootstrap=5000）
```

## 7. 预注册提交前必做清单

- [ ] 与统计合作者过一遍 SAP（`01_*` 的 §6）
- [ ] 锁定第一版共同节点集（16 或 17 节点）—— **必须**预先确定，OSF 提交后不得改
- [ ] 锁定假设 H1–H4 的精确统计描述
- [ ] 锁定 8 项敏感性分析清单
- [ ] 在 OSF 创建项目 → 上传本目录的 docs/01 + docs/02
- [ ] 取 Git tag → 上传代码至 GitHub
- [ ] **同周内**提交 UKB 申请（OSF 优先权时间戳 = 抢发保险）
