# -*- coding: utf-8 -*-
"""删除 _backup_pristine.git 前的唯一性审计。

比较「清理前镜像」与「当前仓库 + 本地磁盘」:
找出哪些被历史清理移除的路径, 在本地磁盘上也已经不存在 ——
那才是删除备份后会真正丢失的东西。
"""
import subprocess, os, collections

PRI = "../_backup_pristine.git"


def paths(repo, ref):
    out = subprocess.check_output(["git", "-C", repo, "ls-tree", "-r", ref], text=True)
    m = {}
    for ln in out.splitlines():
        if not ln.strip():
            continue
        meta, path = ln.split("\t", 1)
        parts = meta.split()
        m[path] = parts[2]  # blob sha
    return m


def sizes(repo, shas):
    """批量取 blob 大小"""
    if not shas:
        return {}
    q = "\n".join(shas) + "\n"
    out = subprocess.run(["git", "-C", repo, "cat-file", "--batch-check"],
                         input=q, capture_output=True, text=True).stdout
    d = {}
    for ln in out.splitlines():
        p = ln.split()
        if len(p) >= 3 and p[1] == "blob":
            d[p[0]] = int(p[2])
    return d


pri_p = paths(PRI, "master")
cur_p = paths(".", "HEAD")
pri_s = sizes(PRI, list(pri_p.values()))

lost = {p: pri_s.get(sha, 0) for p, sha in pri_p.items() if p not in cur_p}
print("被历史清理移除的路径: %d 个, %.1f MB" % (len(lost), sum(lost.values()) / 1048576))

on_disk, missing = [], []
for p, s in lost.items():
    (on_disk if os.path.exists(p) else missing).append((p, s))
print("  本地磁盘仍在   : %d 个, %.1f MB" % (len(on_disk), sum(s for _, s in on_disk) / 1048576))
print("  本地已不存在   : %d 个, %.1f MB" % (len(missing), sum(s for _, s in missing) / 1048576))


def classify(p):
    if "_clump_input.txt" in p or "_clump_in.txt" in p:
        return "PLINK clump 输入 (由 data/ 经脚本再生)"
    if p.startswith("build/"):
        return "build/ 暂存区 (可再生)"
    if p.endswith(".zip"):
        return "已被取代的旧打包 zip (v16/v17)"
    if p.endswith(".tif"):
        return "TIFF 原图 (已恢复到磁盘)"
    if p.startswith("archive/"):
        return "archive/ 旧产物"
    if "/harmonised/" in p:
        return "harmonised 中间文件"
    if p.endswith("cohort_ukb.csv") or "fracture_raw_extract" in p:
        return "原始抽取结果"
    if p.endswith((".log", ".nosex")):
        return "PLINK 日志 / .nosex"
    return "其他: " + p


g = collections.OrderedDict()
for p, s in sorted(missing, key=lambda x: -x[1]):
    k = classify(p)
    g.setdefault(k, [0, 0, []])
    g[k][0] += 1
    g[k][1] += s
    if len(g[k][2]) < 2:
        g[k][2].append(p)

print("\n=== 删除备份后会真正丢失的内容 (按体积) ===")
for k, (c, s, ex) in sorted(g.items(), key=lambda x: -x[1][1]):
    print("  %-40s %4d 个  %8.1f MB" % (k, c, s / 1048576))
    for e in ex:
        print("        e.g. %s" % e)

# 小体积且不可再生的 -> 建议就地恢复后删备份
small_likely_unique = [(p, s) for p, s in missing
                       if not any(t in p for t in ("_clump_input.txt", "_clump_in.txt"))
                       and not p.startswith("build/")
                       and not p.endswith(".zip")
                       and not p.startswith("archive/")
                       and not p.endswith((".log", ".nosex"))
                       and s <= 20 * 1048576]
print("\n建议恢复的小体量(<=20MB)且非打包/非 stage 文件: %d 个, %.2f MB"
      % (len(small_likely_unique), sum(s for _, s in small_likely_unique) / 1048576))
for p, s in small_likely_unique[:25]:
    print("   %8.3f MB  %s" % (s / 1048576, p))
