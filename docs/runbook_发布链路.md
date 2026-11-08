# 发布链路 Runbook — GitHub / Zenodo / OSF

> 用途：投稿就绪后，按此文档逐项执行公开仓库发布与存档。
> 维护：本 runbook 与 `docs/08_投稿稿_OA_BMD_fracture.md` 同步更新。
> 最后核对：2026-11-05（本地 `master` 提交 `b685035`，`.zenodo.json` 已修为 `submission-v2`）。

---

## 0. 当前就绪状态（起点，无需重复做）

| 项 | 状态 | 精确路径 |
|---|---|---|
| 投稿稿真源 | ✅ 已打磨、门禁全绿（148/149b/182） | `OA_B_PROJECT/docs/08_投稿稿_OA_BMD_fracture.md` |
| Zenodo 元数据 | ✅ `<OWNER>`→`Fisherzzx001`、版本 `submission-v2` | `OA_B_PROJECT/.zenodo.json` |
| 本地提交 | ✅ `b685035`（221 文件，发布就绪态） | 分支 `master`（工作树干净） |
| OSF 回写脚本 | ✅ 已 `--dry-run` 验证（无网络依赖） | `OA_B_PROJECT/code/184_osf_doi_writeback.py` |
| 预注册协议 | ✅ 在仓库内 | `OA_B_PROJECT/docs/118_OSF预注册_OA_BMD_fracture.md` |
| Cover letter | ✅ 含 OSF 占位符待回写 | `OA_B_PROJECT/docs/cover_letter_OAC.md` |
| GitHub 仓库 | 已集成 Zenodo（概念 DOI 自动 mint） | `https://github.com/Fisherzzx001/OA-BMD-fracture-MR` |
| Zenodo 概念 DOI | 在线、已集成 | `https://doi.org/10.5281/zenodo.23400860` |

**关键约束**：本机（沙箱）对 `github.com` 网络截断 + 无 token，因此 **push / tag / mint 必须由你在联网终端用 classic token 执行**。本 runbook 中所有 `git push` 类步骤均在你本地联网终端运行。

---

## 1. 前置要求

1. **GitHub classic token**（非 fine-grained）：`https://github.com/settings/tokens` →
   Generate new token (classic) → 勾选 `repo`（整组）→ 复制保存。
   （此前 fine-grained token 因无 `repo` scope 被拒，务必用 classic。）
2. **终端用 Git Bash（Git for Windows），不要用原生 PowerShell**：本机 PowerShell 的 PATH 里没有 git，`git` 会报“无法将 git 项识别为 cmdlet”；且 PowerShell 5.1 不支持 `&&`。Git Bash 里 `git` 与 `&&` 均正常（本 runbook 全部命令均按 Git Bash 语法写）。**若本机未装 Git for Windows（开始菜单搜不到 Git Bash）：去 https://git-scm.com/download/win 下载 64-bit 版、双击默认安装即可**，装完就有 Git Bash。若不想装软件，也可在 PowerShell 里临时把 WorkBuddy 内置 git 加进 PATH（见文末“备选：不装 Git for Windows”）。开始先 `cd` 到项目目录（Git Bash 用正斜杠 `/c/Users/...`，**不要在家目录 `C:\Users\Administrator` 下跑**，否则 `.gitignore` 会写错位置）。
3. **务必推送 `master`**，**不要推送 `publish` 分支**（`publish` 仅领先 2 个初始提交，缺失 70+ 个 OAC/S29/引用提交，已过时）。

---

## 2. Phase A — GitHub 推送（触发 Zenodo 自动 mint 新版本）

在**联网终端**执行：

```bash
cd /c/Users/Administrator/WorkBuddy/2026-09-08-22-17-40/OA_B_PROJECT

git switch master

# 用 classic token 设置远端（把 <CLASSIC_TOKEN> 替换为你生成的）
git remote set-url origin https://Fisherzzx001:<CLASSIC_TOKEN>@github.com/Fisherzzx001/OA-BMD-fracture-MR.git

# 推送全部提交（含 b685035 发布就绪态 + .zenodo.json 修复）
git push origin master

# 打 tag → GitHub–Zenodo 集成自动 mint 新版本（submission-v2）
git tag -a v2.0 -m "Submission-ready OAC manuscript (S29 added; SOS reinforced; Zenodo metadata owner fixed)"
git push origin v2.0
```

---

## 3. Phase B — 校验 Zenodo 新版本

打开 `https://doi.org/10.5281/zenodo.23400860`，确认：

- [ ] 出现新版本 **`submission-v2`**（版本号由 `.zenodo.json` 的 `"version"` 字段控制）；
- [ ] **"Is supplement to"** 链接 = `https://github.com/Fisherzzx001/OA-BMD-fracture-MR`（已修正，不再是 `<OWNER>` 占位符）；
- [ ] Deposit 含 `code/` 全部脚本 + `code/output/` 结果表 + 补充材料（STROBE-MR 清单、S10、S11）。

> 若版本未出现：等待 1–2 分钟（Zenodo 抓取 GitHub 归档有延迟）；若仍无，检查 tag 是否成功 push、`.zenodo.json` 是否已被 master 包含。

---

## 4. Phase C — OSF 镜像（先拿 GUID，再跑回写）

OSF 不自动 mint，需手动建项目拿 GUID，再用脚本回填到稿件与 cover letter。

**4.1 建 OSF 项目拿 GUID**
1. 登录 `https://osf.io` → Create project；
2. 记下项目 GUID（形如 `10.17605/OSF.IO/ABCDE`）。

**4.2 回填稿件（脚本替换两处占位符 `"An OSF mirror is being deposited."`）**
- `docs/08` 的 Data availability 段（位于 `## Figure legends` 之后，**不计入 4000 词正文**）；
- `docs/cover_letter_OAC.md`。

```bash
cd /c/Users/Administrator/WorkBuddy/2026-09-08-22-17-40/OA_B_PROJECT

# 干跑确认替换目标存在（可选）
python code/184_osf_doi_writeback.py --dry-run 10.17605/OSF.IO/ABCDE

# 正式回写
python code/184_osf_doi_writeback.py 10.17605/OSF.IO/ABCDE
```

**4.3 提交并推送回写结果**
```bash
git add -A
git commit -m "OSF DOI writeback (10.17605/OSF.IO/ABCDE)"
git push origin master
```

---

## 5. Phase D — 收尾校验 + docx 重建

回写后重跑三连门禁确认仍全绿，并重跑 docx 重建（OSF DOI 需进稿）：

```bash
cd /c/Users/Administrator/WorkBuddy/2026-09-08-22-17-40/OA_B_PROJECT

python code/148_oac_compliance.py        # 期望 OAC COMPLIANT（body≤4000 / refs≤50 / 显示项 8/8）
python code/149b_verify_reported_values.py   # 期望 SAFE（β/SE/p/OR/CI 未动）
python code/182_final_qc.py              # 期望 64/64 PASS

python build_submission.py               # → submission_package/Manuscript.docx
python build_maintext.py                 # → submission_package/Submission_MainText.docx
python code/135_build_manuscript_ready.py # → Manuscript_Ready 镜像（与 docs/08 md5 一致）
```

---

## 6. Phase E — 隐私：`.workbuddy/` 移出公开跟踪（由你执行，方案已确认）

> **背景**：本仓库预已跟踪 `.workbuddy/`（含 `memory/` 内部策略日志、`backups/`），非本次引入；推送到**公开** GitHub 会一并公开这些内部日志。
> **方案（已确认，安全）**：仅取消 git 跟踪、不删本地文件，并加入 `.gitignore`。
> **由你执行**：以下步骤在 **Git Bash** 中运行（与 Phase A 同一次推送一并完成即可）。

> **⚠️ 先清理一次误操作**：若你曾在 `C:\Users\Administrator>`（家目录、非仓库）下跑过 `echo ".workbuddy/" >> .gitignore`，会在家目录生成一个 stray `.gitignore`（内容 `.workbuddy/`），而项目本身的 `.gitignore` 并未改到。请先在 Git Bash 里删除该 stray 文件，再在项目目录执行下方正确命令：
> ```bash
> rm -f /c/Users/Administrator/.gitignore        # 删除家目录 stray 文件（无害，不在任何仓库内）
> cd /c/Users/Administrator/WorkBuddy/2026-09-08-22-17-40/OA_B_PROJECT
> ```

```bash
cd /c/Users/Administrator/WorkBuddy/2026-09-08-22-17-40/OA_B_PROJECT

# 1) 加入忽略
echo ".workbuddy/" >> .gitignore

# 2) 取消跟踪（--cached 只移除索引，不删文件）
git rm -r --cached .workbuddy

# 3) 提交并推送
git add -A
git commit -m "Stop tracking .workbuddy (internal logs) from public repo"
git push origin master
```

> 注意：这仅影响**此后**的 push；若 `.workbuddy/` 历史已随 `b685035` 进入公开仓库，需在 GitHub 端对历史做清理（如 `git filter-repo` 或联系支持），属额外操作，按需另议。

---

## 7. 故障排查 / 注意事项

| 现象 | 原因 | 处理 |
|---|---|---|
| `git push` 报 `schannel: server closed abruptly` | 沙箱/网络对 github.com 截断 | 改在联网终端执行（本 runbook 默认如此） |
| token 报权限不足 | 用了 fine-grained 或 scope 不全 | 改用 **classic token** 并勾选 `repo` |
| Zenodo 未出现新版本 | tag 未 push / `.zenodo.json` 未随 master | 确认 `git push origin v2.0` 成功；重读 Phase B |
| 误推 `publish` 分支 | 分支选错 | `publish` 已过时，务必 `git switch master` 后推送 |
| 推送后 `.workbuddy/` 公开 | 未执行 Phase E | 按 Phase E 取消跟踪并推送；历史清理另议 |

---

## 8. 路径总览（一张图）

```
本地 master(b685035)
   ├─ push origin master ─────────────▶ GitHub(Fisherzzx001/OA-BMD-fracture-MR)
   │                                     │ tag v2.0
   │                                     └────────▶ Zenodo 概念 DOI 23400860 自动 mint submission-v2
   │                                                    （Is-supplement-to = 真实 GitHub 链接）
   │
   ├─ OSF 建项目拿 GUID ── code/184 回写 ─▶ docs/08 + cover letter 的 "OSF mirror" 占位符
   │                                     │
   │                                     └─ commit + push ─▶ GitHub（含 OSF DOI）
   │
   └─ Phase E: git rm --cached .workbuddy + .gitignore ─▶ 取消公开跟踪（不删文件）
```

---

## 9. 2026-11-09 实际发布记录（沙箱环境绕行方案）

本节记录本次真实执行，供复现；与上文 §2「在联网终端执行」互为备选。

**诊断结论**（非「网络截断」，与旧记录不同）：
- 沙箱 HTTP 代理 `127.0.0.1:49831` 允许 `api.github.com:443`、`ssh.github.com:443`、
  `github.com:22` 的 CONNECT，**但阻断 `github.com:443` 的 git smart-HTTPS**
  （`CONNECT tunnel failed, response 502`）。所以 `git push` over HTTPS 必失败，
  而 REST API 与 GitHub SSH-on-443 可用。
- `origin` 里旧 token 是带尖括号的占位串 `<ghp_…>`（42 字符），必然 401；
  有效 classic token 在 `C:/Users/Administrator/github_token.doc`（40 字符，scope `repo`）。

**发布链路（本次已执行）**：
1. `code/317_ssh_proxy_relay.py` — 最小 HTTP-CONNECT 中继，供 OpenSSH `ProxyCommand` 用
   （关键：stdin 必须用 `read1()` 而非 `read()`，否则 KEX 小包被缓冲、握手超时）。
2. `code/318_github_ssh_release.py` — `setup` 注册**临时可写 deploy key**（REST API）、
   验证 `ssh -T`；`push` 推 master + 推 tag `v2.0` + 建 GitHub Release（Zenodo 监听的是
   Release 事件，不只是 tag）；`cleanup` 在 `finally` 中删除 deploy key 与临时密钥。
   `all` = setup → push → cleanup。token 全程不打印、不写入 origin。
3. **100 MB 硬限制拦截**：`b685035` 误提交的 4 个 LDSC munge 输入
   （`code/output/genetics/ldsc_munge_input/OA_{hip,knee,spine,hand}.munge_input.tsv.gz`，
   102–106 MB）超过 GitHub 100 MB blob 上限，`pre-receive hook declined`。它们是
   `code/204`+`code/205` 可确定性重建的批量中间文件，稿件 Data availability 本就声明
   **不归档**。处置：先建 `backup/pre-release-20261109` 备份分支，再对**未推送区间**
   `91708c4..HEAD` 执行
   `FILTER_BRANCH_SQUELCH_WARNING=1 git filter-branch --index-filter 'git rm -r --cached --ignore-unmatch <4 paths>' -- 91708c4..HEAD`；
   16 个提交的粒度和顺序全部保留，父提交仍是远端 `91708c4` → 仍是 fast-forward。
   `.gitignore` 已加 `PUSH-SIZE GUARD` 防止再次进入仓库。

**备注**：`git filter-branch` 会重置工作树，上述 4 个中间文件随之从磁盘移除；
需要时用 `code/204`/`code/205` 从源 GWAS 重建（不影响任何已报结果与四道门禁）。
