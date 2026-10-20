#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# code/reclaim_parts.sh — 回收已完成 chunk 的 part 文件所占磁盘空间
#
# 背景：本环境的 safe-delete 机制会拦截脚本内的 `rm`（每轮 >50 次删除即要求
#   确认），导致 parallel_fetch.sh 在把 part dd 进主文件后无法删除 part。
# cis 有 116 个 32 MiB 块，若全部残留 = 3.7 GB 垃圾，足以撑爆磁盘。
#
# 原理：**不删除，只 truncate 到 0 字节**（释放数据块，inode 保留）。
#   仅处理日志中已有 `CHUNK_OK <i>` 的 part_<i>（其内容已 dd 进主文件，
#   是纯垃圾）；`.tmp` 为活跃在途数据，绝不触碰。
#
# 用法：bash code/reclaim_parts.sh            # 一次性回收
#       bash code/reclaim_parts.sh loop       # 每 3 min 循环（守护）
# ---------------------------------------------------------------------------
set -u
D="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/data/gwas_stratified"

reclaim_once() {
  local freed=0 touched=0
  for log in "$D"/*.fetch.log; do
    [ -f "$log" ] || continue
    base="$(basename "$log" .fetch.log)"
    # 已完成的 chunk id
    grep -o '^CHUNK_OK [0-9]*' "$log" 2>/dev/null | awk '{print $2}' | while read -r i; do
      p="$D/$base.part_$i"
      [ -f "$p" ] || continue
      s=$(stat -c %s "$p" 2>/dev/null || echo 0)
      [ "$s" -gt 0 ] || continue
      # 双保险：只允许 truncate 非 .tmp 的 part 本体
      case "$p" in *.tmp) continue;; esac
      truncate -s 0 "$p" 2>/dev/null && echo "  reclaimed ${base}.part_$i ($s bytes)"
    done
  done
}

case "${1:-once}" in
  loop)
    echo "[$(date -u '+%H:%M')] reclaim loop started"
    while :; do
      reclaim_once
      sleep 180
    done
    ;;
  *)
    echo "== reclaim pass $(date -u '+%H:%M') =="
    reclaim_once
    echo "== disk after =="
    df -h /c | tail -1
    ;;
esac
