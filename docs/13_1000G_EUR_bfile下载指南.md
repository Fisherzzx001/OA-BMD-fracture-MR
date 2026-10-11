# 13 · 1000G EUR PLINK bfile 下载与严格 clump 实操指南

> 用途：投稿前整改 #3 的**严格 PLINK clumping（r²<0.001 / 10 Mb）**前置条件。
> 本机（agent 环境）无 PLINK、也无 1000G EUR bfile，故此处给出**可在你本地 Windows / 服务器执行的完整路径与步骤**。
> 坐标系：**GRCh37 / hg19** —— 与本项目 GO Consortium 2025、Qian 2025、Neale-lab、UK Biobank 等汇总统计一致，无需 liftover。

---

## 0. 磁盘与工具前置检查

- **可用空间**：EUR 单群体 `.bed` 约 **10–12 GB**（解压峰值含 `.tgz` 约 14–15 GB）。请确认目标盘 ≥ 20 GB 余量。
  - 本 agent 机的 `C:` 盘当前剩余约 29 GB，刚好够；建议下载到**独立目录**（如 `C:/REF/1000G_EUR/`），不要混进仓库。
- **下载工具**：本机 Git Bash **无 `wget`**，但有 **`curl`**（Windows 自带）。以下步骤统一用 `curl`。
- **PLINK**：需自行下载（见 §2）。

---

## 1. 下载 1000G EUR 参考面板（GRCh37）

### 1.1 主源（推荐）：MR-Base / OpenGWAS 同款面板 `1kg.v3.tgz`

这正是 `TwoSampleMR::ld_clump` 在线 API 背后用的 LD 参考，结果可直接与已发表 MR clump 对齐。

```
# 建议目录
mkdir -p "C:/REF/1000G_EUR"
cd "C:/REF/1000G_EUR"

# 下载（约 3–4 GB，支持断点续传 -C -）
curl -L -C - -o 1kg.v3.tgz "http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz"

# 解压（Git Bash 的 tar 可用；含 EUR/AFR/AMR/EAS/SAS 五个群体）
tar -xzvf 1kg.v3.tgz
```

解压后得到（我们只用 EUR）：
```
C:/REF/1000G_EUR/EUR.bed
C:/REF/1000G_EUR/EUR.bim
C:/REF/1000G_EUR/EUR.fam
```
> 若不需要其他群体，可删掉 `AFR.* / AMR.* / EAS.* / SAS.*` 省空间（各约数 GB）。

### 1.2 备用源：Zenodo 6614170（hg19，仅 EUR，含 rsid）

若主源不可达，用此镜像（按 super-population 拆分好的 EUR 单文件，`*.bim` 带 rsid）：

```
cd "C:/REF/1000G_EUR"
curl -L -C - -o 1000G_EUR.bed "https://zenodo.org/records/6614170/files/1000G_EUR.bed"
curl -L -C - -o 1000G_EUR.bim "https://zenodo.org/records/6614170/files/1000G_EUR.bim"
curl -L -C - -o 1000G_EUR.fam "https://zenodo.org/records/6614170/files/1000G_EUR.fam"
```

### 1.3 第三选项（LDSC 全套，含多祖先）：Broad 1000G Phase3

若同时也做 LDSC 权重，可下全祖源版再取 EUR 子集：
`https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_plinkfiles.tgz`
（本任务无需，仅作备注。）

---

## 2. 下载 PLINK 1.9（Windows）

- 官网页面：<https://www.cog-genomics.org/plink/1.9/>
- 直链（plink 1.9 稳定版 win64，若 404 改从官网页面获取最新包）：
  `https://s3.amazonaws.com/plink1-assets/plink_win64_20231025.zip`
- 解压后得到 `plink.exe`，放到例如 `C:/REF/plink_win64/plink.exe`。

> 注：PLINK 2.0 也支持 `--clump`，但输出格式略有差异；严格 clump 建议用 **1.9** 以与 MR-Base 对齐。

---

## 3. 校验 bfile 可加载

```
cd "C:/REF/1000G_EUR"
"C:/REF/plink_win64/plink.exe" --bfile EUR --freq --out eur_check
```

预期：报告样本数 ≈ **503（EUR）**、变异数数千万级，无报错即说明面板正常。
（生成的 `eur_check.*` 可删。）

---

## 4. 接入本项目的严格 clump 流水线

脚本 `code/17_strict_clump_mrpresso.py` 已支持「自动跑 PLINK + 用严格仪器集复跑 MR-PRESSO」。

### 4.1 设置环境变量（Windows cmd / PowerShell）

```bat
set BFILE_1KG_EUR=C:/REF/1000G_EUR/EUR
set PLINK_BIN=C:/REF/plink_win64/plink.exe
```

（Git Bash 下用 `export BFILE_1KG_EUR=...` / `export PLINK_BIN=...`）

### 4.2 一键运行

```bat
cd "C:/Users/Administrator/WorkBuddy/2026-09-08-22-17-40/OA_B_PROJECT/code"
"C:/Users/Administrator/.workbuddy/binaries/python/envs/default/Scripts/python.exe" 17_strict_clump_mrpresso.py
```

脚本会依次：
1. 写出各暴露的 `*_clump_input.txt`（SNP + P，全基因组成因）；
2. 调用 PLINK 对每个暴露做 `--clump --clump-r2 0.001 --clump-kb 10000 --clump-p1 5e-8 --clump-p2 5e-8`，生成 `<exposure>_clumped.clumped`；
3. 自动检测 clump 结果，用之覆盖仪器集，复跑 MR + MR-PRESSO，更新 `mr_presso_results.csv`。

> 等价 PLINK 命令（脚本内部即生成并执行）：
> ```
> plink --bfile C:/REF/1000G_EUR/EUR --clump <exp>_clump_input.txt \
>       --clump-p1 5e-8 --clump-p2 5e-8 --clump-r2 0.001 --clump-kb 10000 \
>       --out <exp>_clumped
> ```

### 4.3 无 PLINK 时的降级行为

若未设置 `PLINK_BIN` 或不可达，脚本仅**打印**上述命令（不报错），并继续用现有近似 1-Mb 仪器集跑 MR-PRESSO——即当前已交付的行为，结果不变。

---

## 5. 预期产出与后续

- 严格 clump 后，暴露仪器数通常会比 1-Mb 近似版**更少**（r²<0.001 比 1-Mb 窗口剔除更彻底），MR 估计应更稳定。
- 对比 `mr_presso_results.csv` 严格版 vs 当前近似版：关注 BMD→骨折、BMD→膝 OA 的 β 方向与显著性是否一致（预期稳健，不翻转）。
- 将严格版结果回填稿件 §2.6 / §3.6 / Table 3，并注明「PLINK r²<0.001 / 10 Mb，1000G Phase3 EUR」。

---

## 6. 一句话总结

下载 `1kg.v3.tgz`（MR-Base 同款 EUR 面板，GRCh37）→ 解压取 `EUR.bed/.bim/.fam` → 下载 PLINK 1.9 → 设 `BFILE_1KG_EUR` / `PLINK_BIN` → 跑 `17_strict_clump_mrpresso.py`，严格 clump + MR-PRESSO 全自动完成。
