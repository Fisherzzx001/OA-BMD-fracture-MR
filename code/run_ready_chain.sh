#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# code/run_ready_chain.sh — 数据就绪后的一键执行链（F8 性别分层 / F4 eQTLGen coloc）
#
# 用法（在 OA_B_PROJECT/ 下）：
#   bash code/run_ready_chain.sh f8     # T1: KneeOA+AllOA -> 161 -> 回填要点
#   bash code/run_ready_chain.sh f4     # eQTLGen: cis+AF -> 166 -> 165 -> 回填要点
#
# 每个输入先做"存在 + 大小匹配 + gzip -t"三重校验，缺一不跑（不虚构）。
# 完成后打印 docs/118 回填所需的全部关键数值（人工核对后写入）。
# ---------------------------------------------------------------------------
set -u
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA="$ROOT/data/gwas_stratified"
OUT="$ROOT/code/output/genetics"
PY="C:/Users/Administrator/.workbuddy/binaries/python/envs/default/Scripts/python.exe"

need_file() { # path expected_bytes
  local p="$1" want="$2"
  if [ ! -f "$p" ]; then echo "  [MISS] $p"; return 1; fi
  local got; got=$(stat -L -c %s "$p" 2>/dev/null || echo 0)
  if [ "$got" -ne "$want" ]; then echo "  [SIZE] $p ($got != $want) — still downloading?"; return 1; fi
  if ! gzip -t "$p" 2>/dev/null; then echo "  [GZIP] $p fails integrity — rerun parallel_fetch.sh"; return 1; fi
  echo "  [OK]   $p ($got bytes, gzip valid)"
  return 0
}

case "${1:-}" in
# ---------------------------------------------------------------- F8: T1 性别分层
f8)
  echo "== F8 前置校验 =="
  ok=1
  need_file "$DATA/KneeOA_sexspecific.tsv.gz" 744237546 || ok=0
  need_file "$DATA/AllOA_sexspecific.tsv.gz"  918627766 || ok=0
  [ "$ok" = 1 ] || { echo "F8 输入未就绪，退出。"; exit 1; }
  echo "== 运行 code/161_tier1_sex_mr.py =="
  "$PY" "$ROOT/code/161_tier1_sex_mr.py" || { echo "161 运行失败"; exit 1; }
  echo
  echo "== F8 回填要点（docs/118 F8 行用）=="
  "$PY" - <<'PYEOF'
import json, os
p = os.path.join("code","output","genetics","tier1_sex_mr_summary.json")
d = json.load(open(p))
for oc in ("KneeOA","AllOA"):
    if oc not in d: continue
    r = d[oc]; m = r.get("meta",{})
    F, M = m.get("Female",{}), m.get("Male",{})
    fvm = m.get("F_vs_M",{})
    print(f"[{oc}] Female meta beta={F.get('meta_beta')} (SE {F.get('meta_se')}, RE p={F.get('meta_re_p')}, n_site={F.get('n_sites')})")
    print(f"[{oc}] Male   meta beta={M.get('meta_beta')} (SE {M.get('meta_se')}, RE p={M.get('meta_re_p')}, n_site={M.get('n_sites')})")
    print(f"[{oc}] F vs M diff={fvm.get('diff')} p={fvm.get('p_two_sided')} | site-het Q={m.get('site_het_Q')} p={m.get('site_het_P')}")
PYEOF
  ;;
# ---------------------------------------------------------------- F4: eQTLGen coloc
f4)
  echo "== F4 前置校验 =="
  ok=1
  need_file "$DATA/eQTLGen_SNP_AF.txt.gz"       240045342 || ok=0
  need_file "$DATA/eQTLGen_cis_full.tsv.gz"    3880187026 || ok=0
  [ "$ok" = 1 ] || { echo "F4 输入未就绪，退出。"; exit 1; }
  echo "== 1/2 code/166_eqtlgen_loci_extract.py（全量 cis 流式扫描，约 20-40 min）=="
  "$PY" "$ROOT/code/166_eqtlgen_loci_extract.py" || { echo "166 运行失败"; exit 1; }
  echo "== 2/2 code/165_eqtlgen_blood_coloc.py =="
  "$PY" "$ROOT/code/165_eqtlgen_blood_coloc.py" || { echo "165 运行失败"; exit 1; }
  echo
  echo "== F4 回填要点（docs/118 F4 行用）=="
  "$PY" - <<'PYEOF'
import pandas as pd, os
s = os.path.join("code","output","genetics","eqtlgen_blood_coloc_locus_summary.csv")
if not os.path.exists(s):
    print("!! summary 不存在"); raise SystemExit
d = pd.read_csv(s)
ok = d["PP4_blood_eqtl_OA"].notna()
print(f"loci tested: {len(d)}; with PP4: {int(ok.sum())}")
for lab, col in (("OA_knee","PP4_blood_eqtl_OA"), ("BMD_total","PP4_blood_eqtl_BMD")):
    v = d[col].dropna()
    if len(v):
        print(f"PP.H4(blood_eQTL, {lab}): median={v.median():.3f} mean={v.mean():.3f} "
              f"max={v.max():.3f} n>=0.8={int((v>=0.8).sum())}/{len(v)}")
both = d[ok & d["PP4_blood_eqtl_BMD"].notna()]
if len(both):
    from scipy.stats import wilcoxon
    try:
        st, p = wilcoxon(both["PP4_blood_eqtl_BMD"], both["PP4_blood_eqtl_OA"])
        print(f"paired Wilcoxon BMD vs OA: W={st:.1f} p={p:.4g} (n={len(both)})")
    except Exception as e:
        print("wilcoxon skipped:", e)
print("\nTop loci by PP4(OA):")
print(d.sort_values("PP4_blood_eqtl_OA", ascending=False).head(8)[
    ["index_snp","chrom","pos","best_gene_symbol","n_snp","PP4_blood_eqtl_OA","PP4_blood_eqtl_BMD"]].to_string(index=False))
PYEOF
  ;;
*)
  echo "用法: bash code/run_ready_chain.sh {f8|f4}"; exit 1;;
esac
