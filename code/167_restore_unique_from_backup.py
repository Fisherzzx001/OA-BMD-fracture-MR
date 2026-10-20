# -*- coding: utf-8 -*-
"""删除 _backup_pristine.git 之前, 把其中「本地已消失且不可再生」的小体量文件恢复到本地磁盘。

跳过(可再生/已被取代):
  - PLINK *_clump_input.txt (由 data/ 原始 GWAS 经 code/ 脚本再生)
  - build/ 暂存区
  - 已被新版本取代的旧打包 zip
  - PLINK 日志 / .nosex

恢复:
  - TIFF 主图/补充图 (期刊高清投稿资产)
  - code/output/genetics/harmonised/* (分析中间文件)
  - raw extracts (fracture_raw_extract.csv, cohort_ukb.csv)

每个文件写盘后用 sha256 与备份中的 blob 内容比对, 失败即中止。
"""
import hashlib, os, subprocess, sys

PRI = "../_backup_pristine.git"
SKIP = ("_clump_input.txt", "_clump_in.txt")
SKIP_PREFIX = ("build/", "archive/")
SKIP_SUFFIX = (".zip", ".log", ".nosex")


def paths(repo, ref):
    out = subprocess.check_output(["git", "-C", repo, "ls-tree", "-r", ref], text=True)
    m = {}
    for ln in out.splitlines():
        if not ln.strip():
            continue
        meta, path = ln.split("\t", 1)
        m[path] = meta.split()[2]
    return m


pri = paths(PRI, "master")
cur = paths(".", "HEAD")

targets = []
for p, sha in pri.items():
    if p in cur:
        continue                      # 仍在新历史里
    if os.path.exists(p):
        continue                      # 本地还在
    if any(t in p for t in SKIP) or p.startswith(SKIP_PREFIX) or p.endswith(SKIP_SUFFIX):
        continue                      # 可再生 / 已取代
    targets.append((p, sha))

print("待恢复 %d 个文件" % len(targets))
tot = 0
fail = 0
for p, sha in sorted(targets):
    os.makedirs(os.path.dirname(p) or ".", exist_ok=True)
    blob = subprocess.check_output(["git", "-C", PRI, "cat-file", "blob", sha])
    want = hashlib.sha256(blob).hexdigest()
    with open(p, "wb") as f:
        f.write(blob)
    got = hashlib.sha256(open(p, "rb").read()).hexdigest()
    ok = (want == got)
    tot += len(blob)
    fail += 0 if ok else 1
    flag = "OK " if ok else "FAIL"
    print("  [%s] %8.3f MB  %s" % (flag, len(blob) / 1048576, p))
    if not ok:
        print("        >> 哈希不匹配, 中止")
        sys.exit(1)

print()
print("恢复完成: %d 个文件, %.2f MB, 失败 %d" % (len(targets), tot / 1048576, fail))

print("\n=== 复核: Figures_Original_v14 是否齐全 ===")
in_pri = [p for p in pri if p.startswith("Figures_v17_Nature/Figures_Original_v14/")]
miss_orig = [p for p in in_pri if not os.path.exists(p)]
print("  备份中该目录: %d 个 | 磁盘缺失: %d 个" % (len(in_pri), len(miss_orig)))
for p in miss_orig[:5]:
    print("     missing:", p)
