# 预印本与预注册包 · Preprint & Preregistration Package

**目的**：在期刊投稿前/投稿同时取得**研究优先权时间戳**，防御"网络结构/跨库复杂度比较"等创新点被抢发（抢发风险评估见项目方案 §抢发风险）。

**顺序建议**：OSF 预注册（当天，~30 分钟）→ medRxiv 预印本（1–3 个工作日过审）→ 期刊投稿（*Bone* / *Osteoporosis International*）。

---

## 包内文件

| 文件 | 用途 |
|---|---|
| `01_medRxiv_标题页与必填声明.md` | medRxiv 提交表单所需全部文本（标题页、资助、COI、伦理、数据/代码可用性、CRediT 贡献、致谢）。★条目需作者确认 |
| `02_OSF预注册_PREREG.md` | OSF 预注册全文（假设 H1–H5、数据 accession、分析口径、**偏差披露**）。标注为回溯性预注册，诚实性条款勿删 |
| 稿件正文 | 使用投稿包 `submission_package_v2/Manuscript.docx`（medRxiv 接受 docx 或 PDF，建议 PDF） |

## 操作步骤

### 第 1 步 · OSF 预注册（osf.io）

1. 注册/登录 osf.io → **Create project**，标题用 `02_OSF预注册_PREREG.md` §1 的标题；
2. Contributors 添加全部作者（OA = view only 或 admin，按意愿）；★贡献者排序请作者确认；
3. 把 `02_OSF预注册_PREREG.md` 内容粘入 wiki（或作为 PDF 上传至 Files）；
4. 建议直接用 **OSF Prereg 模板（"Prereg Challenge"或"Standard Pre-Data Collection Registration"）**时，在开头保留"回溯性预注册"声明；
5. 项目设为 **Public**（或 Embargo 至预印本上线同日），生成 DOI。

### 第 2 步 · medRxiv 预印本（submit.medrxiv.org）

1. 注册账号（通讯作者邮箱 after9951023@gmail.com 或 baohuricha@126.com 之一，需 ORCID——**请提前为两位通讯作者绑定 ORCID**）；
2. 选择 Category：Epidemiology（或 Genetics）；
3. 上传 `Manuscript.docx`（建议另存 PDF）；标题页/声明按 `01_medRxiv_标题页与必填声明.md` 逐项粘贴；
4. Funding 一栏填 "None"（medRxiv 允许无资助声明）；License 建议 CC BY-NC-ND；
5. medRxiv 为**筛审（moderated preprint）**，通常 1–3 个工作日上线，上线即得预印本 DOI；
6. 上线后：把预印本 DOI 写回期刊稿件 §Data availability（"A preprint of this article is available at medRxiv: doi:xxx"）并同步更新 Cover Letter。

### 第 3 步 · 代码公开（GitHub + Zenodo，★可与预印本同步）

1. 新建 GitHub 仓库，推送 `code/`（01–26 全部脚本）与 `build_figures*.py`；README 用 `code/README.md`；
2.Zenodo 连接 GitHub 仓库 → release 打 tag（如 `v1.0-submission`）→ 取得代码 DOI；
3. 把代码 DOI 填入 medRxiv 稿件 Code availability（替换"on request"表述）——**强烈建议**，审稿人与编辑对"代码可得"的权重很高。

## 注意事项

- **回溯性预注册的定位**：`02` 文件已明确声明分析在注册前完成，不删除、不改写该声明；这是诚实做法，且不影响作为优先权与透明度证据的价值。
- medRxiv 不接受与已发表论文完全重复的内容本包不涉及；预印本与期刊稿可同存。
- 预印本上线后如需修订，medRxiv 支持 version 2（保留原 DOI，追加版本号）。
- 若期刊（*Bone*/*Osteoporosis International*/*OAC*）有 preprint 政策疑问：三刊均为 preprint-friendly（Springer Nature 政策允许预印本）。
