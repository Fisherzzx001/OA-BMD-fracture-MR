#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# code/auto_ready_chain.sh — 后台监控：数据就绪即自动执行 F8/F4 分析链
#
# 每 5 min 检查一次；三重校验（存在 + 字节数精确匹配 + gzip -t 通过）后触发：
#   F8: KneeOA + AllOA        -> code/161_tier1_sex_mr.py
#   F4: cis + AF              -> code/166 -> code/165
# 各自用 flag 文件保证只跑一次；全部日志落 data/gwas_stratified/auto_chain.log。
# 超过 MAX_HOURS 小时未就绪则退出（避免永久后台进程）。
# ---------------------------------------------------------------------------
set -u
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA="$ROOT/data/gwas_stratified"
LOG="$DATA/auto_chain.log"
MAX_HOURS="${MAX_HOURS:-26}"

log() { echo "[$(date -u '+%Y-%m-%d %H:%M')] $*" >> "$LOG"; }

ready() { # path expected_bytes
  local p="$1" want="$2"
  [ -f "$p" ] || return 1
  [ "$(stat -L -c %s "$p" 2>/dev/null || echo 0)" -eq "$want" ] || return 1
  gzip -t "$p" 2>/dev/null || return 1
  return 0
}

log "auto_ready_chain started (pid $$, max ${MAX_HOURS}h)"
t0=$(date +%s)
f8_done=0; f4_done=0
while :; do
  # ---- F8: T1 性别分层 ----
  if [ "$f8_done" -eq 0 ] && ready "$DATA/KneeOA_sexspecific.tsv.gz" 744237546 \
                          && ready "$DATA/AllOA_sexspecific.tsv.gz"  918627766; then
    log "F8 inputs READY — launching 161"
    bash "$ROOT/code/run_ready_chain.sh" f8 >> "$LOG" 2>&1
    rc=$?
    log "F8 chain exit=$rc"
    # 释放 T1 raw（公开可重新下载）以回收磁盘；但 T1 现已迁移到 G: 并通过 C: 符号链接访问，
    # 若直接 truncate 符号链接会跟随到 G: 把唯一副本清零（已发生过的事故）。故：仅当是真实文件时才 truncate，
    # 符号链接一律跳过——符号链接本身只占几十字节，无需回收，且 G: 余量充足。
    for f in "$DATA/KneeOA_sexspecific.tsv.gz" "$DATA/AllOA_sexspecific.tsv.gz"; do
      [ -L "$f" ] || truncate -s 0 "$f" 2>/dev/null
    done
    log "freed T1 raw files (truncated real files only; symlinks skipped to protect G: copy)"
    f8_done=1   # 标记已尝试，避免每 5min 重跑；若 rc!=0 由人工按日志重下 T1 复跑 161
  fi
  # ---- F4: eQTLGen coloc ----
  if [ "$f4_done" -eq 0 ] && ready "$DATA/eQTLGen_SNP_AF.txt.gz"    240045342 \
                          && ready "$DATA/eQTLGen_cis_full.tsv.gz" 3880187026; then
    log "F4 inputs READY — launching 166 -> 165"
    bash "$ROOT/code/run_ready_chain.sh" f4 >> "$LOG" 2>&1
    log "F4 chain exit=$?"
    f4_done=1
  fi
  # ---- F8 跑完后：磁盘压力已解除 -> 自动恢复 cis 续传，以便后续跑 F4 ----
  # （策略：串行调度。cis 与 T1 争用带宽会让两边都慢到不可用；先保 F8 产出，
  #   再让 cis 独占 molgenis 主机完成 F4。cis 从 MiB 对齐断点续传，已下字节不丢。）
  if [ "$f8_done" -eq 1 ] && [ "$f4_done" -eq 0 ]; then
    cis_running=$(ps -ef | grep "[p]arallel_fetch" | grep -c "eQTLGen_cis_full" || true)
    if [ "${cis_running:-0}" -eq 0 ]; then
      log "F8 done & disk relieved — resuming cis download"
      nohup bash "$ROOT/code/parallel_fetch.sh" \
        "https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz" \
        "$DATA/eQTLGen_cis_full.tsv.gz" 6 32 >> "$LOG" 2>&1 &
      disown
    fi
  fi
  [ "$f8_done" -eq 1 ] && [ "$f4_done" -eq 1 ] && { log "both chains done — monitor exit"; break; }
  el=$(( ($(date +%s) - t0) / 3600 ))
  if [ "$el" -ge "$MAX_HOURS" ]; then
    log "timeout ${MAX_HOURS}h reached (f8=$f8_done f4=$f4_done) — monitor exit"
    break
  fi
  sleep 300
done
