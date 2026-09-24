# 获取真实 GitHub 仓库地址与 Zenodo DOI 的详细方法

> 目的：把 `docs/08` 数据可用性声明里的两处占位符换成真实值：
> - `https://github.com/&lt;owner&gt;/OA-BMD-fracture-MR` → 真实仓库 URL
> - `doi:10.5281/zenodo.XXXXXXX` → 真实 Zenodo DOI
>
> **原则：绝不编造。** DOI 只能由 Zenodo 在你真正归档后下发；拿到真实值后告诉我，我来回填并重建投稿包。

---

## 一、GitHub 仓库地址（替换 `&lt;owner&gt;`）

### 步骤
1. 注册 / 登录 <https://github.com>（用你自己的账号，或课题组机构账号）。
2. 右上角 **New repository**：
   - Repository name：`OA-BMD-fracture-MR`（与稿件一致，便于审稿人核对）
   - **Visibility：Public**（数据可用性必须可公开访问，私有仓库 Zenodo 无法抓取）
   - 初始化选项：**不要**勾选 Add a README / .gitignore / License（本地已有仓库结构，避免冲突）
   - 点 **Create repository**
3. 本地关联并推送（在项目根目录 `OA_B_PROJECT/` 执行）：
   ```bash
   git remote add origin https://github.com/<你的用户名>/OA-BMD-fracture-MR.git
   git add code/ docs/10_*.md docs/11_*.md docs/12_*.md docs/14_*.md   # 只推分析与说明，勿推投稿包/大图
   git commit -m "Analysis code for OA-BMD-fracture cross-trait MR"
   git push -u origin main        # 或 master，视默认分支名而定
   ```
4. 推送成功后，仓库地址即：
   `https://github.com/<你的用户名>/OA-BMD-fracture-MR`
   这就是要回填到 `docs/08` 第 366 行的值。

### 注意事项
- **不要**把 `submission_package_v16/`、`Figures_v17_Nature/`（大图）、`archive/` 推上 GitHub。
  本仓库根目录已有 `.gitignore` 会忽略 `submission_package_v16/`，但请确认 `Figures_v17_Nature/` 也被忽略（或手动 `git add` 时不含它们），否则仓库会过大。
- 如确实需要把图也公开，建议用 **Git LFS** 或大图只放 Zenodo，GitHub 仅放代码与关键 CSV。
- 确保 `code/` 下稿件引用的文件真实存在：`code/prep_mrlap_gwas.py`、`code/fix_mrlap_inputs.py`、`code/run_mrlap.R`，以及输出 `code/output/mrlap_overlap_results.csv`。

---

## 二、Zenodo DOI（替换 `XXXXXXX`）

Zenodo DOI 有两种获取路径，**推荐路径 A（GitHub 联动）**，因为稿件写的是"archived at Zenodo"。

### 路径 A：GitHub Release → Zenodo 自动归档（首选）
1. 登录 <https://zenodo.org>（点 **Log in** → 选 **GitHub**，用同一个 GitHub 账号 OAuth 授权）。
2. 顶部菜单 **GitHub** → 找到 `OA-BMD-fracture-MR` 仓库，**打开开关**（Authorize / flip on）。
   Zenodo 从此会对该仓库的每个 **Release** 自动抓取并归档。
3. 回到 GitHub 仓库页面 → **Releases** → **Draft a new release**：
   - Tag version：`v1.0.0`
   - Release title：`OA-BMD-fracture cross-trait MR analysis v1.0.0`
   - 描述：简述分析内容与对应稿件
   - 点 **Publish release**
4. 几分钟后回到 Zenodo → 你的 **Uploads** 里会出现一条新记录（标题同 release）。
   - Zenodo 会为该版本下发一个 **DOI**，形如 `10.5281/zenodo.1234567`。
   - 同时有一个不变的 **concept DOI**（用于引用整个项目）。
5. 复制该版本 DOI，回填 `docs/08` 第 366 行的 `doi:10.5281/zenodo.XXXXXXX` → `doi:10.5281/zenodo.1234567`。
6. 之后每次发新 release，Zenodo 会下发**新版本号 DOI**（如 `.1234568`），旧 DOI 永久保留——适合论文接收后再发一版"accepted manuscript"对应 release。

### 路径 B：Zenodo 直接上传（不走 GitHub）
1. 登录 <https://zenodo.org> → 右上 **Upload** → **New upload**。
2. 把分析代码 + 关键 CSV（或打个 `code.zip`）拖入上传区。
3. 填写元数据：Title、Authors（与稿件作者一致）、Description、Keywords（Mendelian randomization; Osteoarthritis; Bone mineral density）、License（代码用 **MIT**，数据用 **CC-BY-4.0**）。
4. 点 **Publish** → Zenodo 立即下发 DOI。
5. 回填稿件占位符即可。

### 注意事项
- Zenodo 要求**公开**归档才能下发 DOI； embargo（禁阅）状态下无 DOI。
- DOI 一旦下发**不可删除、不可修改**，只能发新版本。所以先确认代码无误再 Publish。
- 稿件引用的 `code/output/mrlap_overlap_results.csv` 等若要公开，确保包含在归档内（路径 A 会随 release 整体抓取；路径 B 需手动加入）。

---

## 三、回填后由我重建投稿包

拿到两个真实值后，把以下信息发给我：
- GitHub URL：`https://github.com/<用户名>/OA-BMD-fracture-MR`
- Zenodo DOI：`10.5281/zenodo.1234567`

我会：
1. 改 `docs/08` 第 366 行两处占位符；
2. 重跑 `build_submission.py`（OUTDIR=submission_package_v16）+ `build_package_v16.py` 重建 `Manuscript.docx` 与 `Submission_Package_v16.zip`；
3. 提交并交付。

---

## 四、当前进度核对

| 占位符 | 状态 |
|---|---|
| `https://github.com/&lt;owner&gt;/OA-BMD-fracture-MR` | 待你建仓后提供真实 URL |
| `doi:10.5281/zenodo.XXXXXXX` | 待你发 GitHub release（或 Zenodo 直传）后提供真实 DOI |
| 参考文献 16/17/18 卷期页 | ✅ 已补（Crossref/PubMed 核实：76(3):366–375 / 57(1):2533428 / 17(2):157–159） |
