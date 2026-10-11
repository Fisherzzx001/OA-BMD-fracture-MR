"""
02_harmonise.py — 原始数据 → 协调后分析表（自动字段探测 + 核对报告）
================================================================================
功能：
  - 读入 CHARLS 原始数据（.dta / .csv）与 UKB 导出（tab 分隔，多列）
  - 按 fields.py 的候选模式自动匹配字段
  - 构造 16 病节点的 0/1 矩阵（CHARLS 用自报，UKB 优先 ICD-10 + 自报辅证）
  - 生成《字段核对报告》—— 标红找不到的字段，供你手动补全
  - 输出 output/data/charls.csv / ukb.csv（与 simulate_data.py 输出一致）

运行：
  python 02_harmonise.py --charls path/to/charls.dta
  python 02_harmonise.py --ukb path/to/ukb_export.tsv
================================================================================
"""

import os
import re
import sys
import argparse
import warnings
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import DKEYS, DLABEL
from fields import (CHARLS_FIELDS, UKB_FIELDS, ICD10_SETS, UKB_SELFREPORT,
                    CHARLS_MULTISITE_SITES)

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)


# ---------------------------------------------------------------- 字段探测
def detect_column(columns, patterns, ignore_case=True, multi=False):
    """
    按候选正则模式顺序尝试匹配，返回：
      single → str（首个匹配）
      multi  → list（所有匹配）
    找不到返回 None。
    """
    cols = [c for c in columns
            if not ignore_case or True]   # 占位
    cols_lc = {c.lower(): c for c in columns}
    found = []
    for pat in patterns:
        rgx = re.compile(pat, re.I if ignore_case else 0)
        for lc, orig in cols_lc.items():
            if rgx.match(lc) and (multi or not found):
                found.append(orig)
                if not multi:
                    return orig
    return found if multi else (found[0] if found else None)


def detect_all(arm, columns, fields_dict):
    report, matched = {}, {}
    for concept, patterns in fields_dict.items():
        hit = detect_column(columns, patterns, multi=concept == "icd10_all"
                                                       or concept == "self_illness"
                                                       or concept == "medication"
                                                       or concept == "pain_types")
        matched[concept] = hit
        status = "✓" if hit else "⚠️ 缺失"
        report[concept] = (status, patterns if not hit else
                            (hit if not isinstance(hit, list) else f"{len(hit)} 列"))
    return matched, report


# ---------------------------------------------------------------- 疾病构造
def build_charls(df, matched):
    out = pd.DataFrame()
    out["pid"] = df[matched["pid"]] if matched["pid"] else df.index.astype(str)

    if matched["age"] and matched["age"] in df.columns:
        out["age"] = df[matched["age"]]
    elif "ba002_1" in df.columns:
        out["age"] = 2011 - df["ba002_1"]
    else:
        out["age"] = np.nan
    out["sex"] = df[matched["sex"]] if matched["sex"] else np.nan

    if matched["bmi"]:
        out["bmi"] = df[matched["bmi"]]
    elif matched["weight_kg"] and matched["weight_cm"]:
        h_m = df[matched["weight_cm"]] / 100
        out["bmi"] = df[matched["weight_kg"]] / (h_m ** 2)
    else:
        out["bmi"] = np.nan

    out["stratum"] = (df[matched["stratum"]] if matched["stratum"]
                      else np.zeros(len(df), dtype=int))
    out["weight"] = (df[matched["samp_weight"]] if matched["samp_weight"]
                     else np.ones(len(df), dtype=float))

    # 16 病节点
    disease_map = {
        "hypertension": "hypertension", "diabetes": "diabetes",
        "dyslipidaemia": "dyslipidaemia", "heart_disease": "heart_disease",
        "stroke": "stroke", "asthma": "asthma", "copd": "copd",
        "cancer": "cancer", "liver_disease": "liver_disease",
        "kidney_disease": "kidney_disease", "digestive": "digestive",
        "depression": "depression", "dementia": "dementia",
        "back_pain": "pain_back", "neck_pain": "pain_neck",
        "multisite_pain": None,            # 由多个 pain 部位衍生
    }
    for out_key, src in disease_map.items():
        if src is None or not matched.get(src):
            out[out_key] = np.nan
        else:
            col = matched[src]
            if col in df.columns:
                out[out_key] = (df[col]
                                .astype(str).str.lower()
                                .isin(["1", "1.0", "yes", "true", "y", "t"])
                                .astype("Int8"))
            else:
                out[out_key] = np.nan

    # 多部位关节痛 = 多个疼痛部位中 ≥2 个为真
    site_cols = [matched[s] for s in CHARLS_MULTISITE_SITES
                 if matched.get(s) and matched[s] in df.columns]
    if site_cols:
        site_mat = df[site_cols].apply(
            lambda c: c.astype(str).str.lower().isin(["1", "1.0", "yes", "true"]))
        out["multisite_pain"] = (site_mat.sum(axis=1) >= 2).astype("Int8")
    else:
        out["multisite_pain"] = np.nan

    # 主结局：T1 any OA
    oa_col = matched.get("oa_any")
    if oa_col and oa_col in df.columns:
        out["oa"] = (df[oa_col].astype(str).str.lower()
                     .isin(["1", "1.0", "yes", "true", "y", "t"]).astype("Int8"))
    else:
        out["oa"] = np.nan
    return out


def _icd10_match(values, codes):
    """Series 元素是否以前缀匹配任一 ICD-10 编码。"""
    if not isinstance(values, pd.Series):
        return pd.Series([False] * len(values), index=values.index if hasattr(values, "index") else None)
    s = values.astype(str).str.upper().str.replace(".", "", regex=False).str.strip()
    out = pd.Series(False, index=values.index)
    for c in codes:
        out |= s.str.startswith(c.upper())
    return out


def build_ukb(df, matched):
    critical = ["pid", "age", "sex", "bmi", "icd10_all", "centre"]
    missing_crit = [c for c in critical if not matched.get(c)]
    if missing_crit:
        print(f"  ⚠️ 关键 UKB 字段缺失 {missing_crit} ——"
              f" 跳过本表生成。请确认 UKB 导出文件包含 eid、Field 31、21001、21022、41270、54 等列。")
        return None
    out = pd.DataFrame()
    out["pid"] = df[matched["pid"]]
    out["age"] = df[matched["age"]]
    out["sex"] = df[matched["sex"]]
    out["bmi"] = df[matched["bmi"]]
    out["stratum"] = (df[matched["centre"]]
                      if matched["centre"] else np.zeros(len(df), dtype=int))
    out["weight"] = np.ones(len(df), dtype=float)

    # ICD-10 数组（多列）：合并为一个 Series
    def _stack(prefix):
        cols = [c for c in df.columns if c.startswith(prefix)]
        if not cols:
            return pd.Series([None] * len(df), index=df.index)
        return df[cols].stack(dropna=True) if False else \
            pd.Series(np.concatenate([df[c].astype("object").fillna("").values
                                      for c in cols]),
                      index=df.index.repeat(len(cols)))

    # 实际：直接对每列做 prefix-match，再按行 any
    def _by_row(prefix, codes):
        cols = [c for c in df.columns if c.startswith(prefix)]
        if not cols:
            return pd.Series(False, index=df.index)
        # 取每个 cell 的字符串，匹配任一 code 前缀
        block = df[cols].astype(str)
        mask = pd.DataFrame(False, index=df.index, columns=cols)
        for c in codes:
            mask |= block.apply(lambda s: s.str.upper().str.startswith(c.upper()))
        return mask.any(axis=1)

    for key, codes in ICD10_SETS.items():
        if key.startswith("oa") or key == "inflammatory_arthritis":
            continue
        out[key] = _by_row("41270", codes).astype("Int8")
        if out[key].sum() == 0:                  # 没有 HES 时退回自报
            self_codes = UKB_SELFREPORT.get(key, [])
            if self_codes and matched.get("self_illness"):
                cols = [c for c in df.columns if c.startswith("20002")]
                if cols:
                    block = df[cols].astype(str)
                    out[key] = block.apply(
                        lambda r: any(str(int(float(v))) in [str(s) for s in self_codes]
                                      if v not in ("", "nan", "None") else False
                                      for v in r),
                        axis=1).astype("Int8")

    # 多部位关节痛：pain_types 3773（多选字段，逗号/分号分隔或多数组）
    pt_cols = [c for c in df.columns if c.startswith("3773")]
    if pt_cols:
        block = df[pt_cols].astype(str)
        # UKB 多选疼痛部位的 code：3 颈/肩 4 背 6 髋 7 膝 1 头痛 ...（⚠️ 待核）
        pain_set = {"3", "4", "6", "7"}
        cnt = block.apply(lambda r: sum(1 for v in r
                                        if isinstance(v, str) and v.strip() in pain_set), axis=1)
        out["multisite_pain"] = (cnt >= 2).astype("Int8")
    else:
        out["multisite_pain"] = np.nan

    # OA 主表型
    out["oa"] = _by_row("41270", ICD10_SETS["oa"]).astype("Int8")
    return out


# ---------------------------------------------------------------- 主流程
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--charls", help="CHARLS 原始数据路径（.dta 或 .csv）")
    ap.add_argument("--ukb", help="UKB 导出路径（.tsv 或 .csv）")
    ap.add_argument("--out", default="output/data")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    lines = ["# 字段核对报告\n", "| 库 | 概念 | 状态 | 匹配字段 | 候选模式 |", "|---|---|---|---|---|"]

    if args.charls:
        print(f"[CHARLS] 读取 {args.charls}")
        df = (pd.read_stata(args.charls) if args.charls.endswith(".dta")
              else pd.read_csv(args.charls, low_memory=False))
        print(f"  {df.shape[0]:,} 行 × {df.shape[1]} 列")
        m, r = detect_all("CHARLS", df.columns, CHARLS_FIELDS)
        miss = []
        for c, (s, info) in r.items():
            lines.append(f"| CHARLS | {c} | {s} | "
                         f"{info if isinstance(info, str) else '|'.join(info)} | "
                         f"{','.join(CHARLS_FIELDS[c][:2])} |")
            if s.startswith("⚠️") and c in {"age", "sex", "oa_any",
                                              "hypertension", "diabetes",
                                              "dyslipidaemia", "heart_disease",
                                              "stroke", "asthma", "copd",
                                              "back_pain", "samp_weight"}:
                miss.append(c)
        if miss:
            print(f"  ⚠️ 关键字段未匹配：{miss}")
        out = build_charls(df, m)
        out.to_csv(os.path.join(args.out, "charls.csv"), index=False)
        print(f"  → output/data/charls.csv ({out.shape[0]:,} 行)")
    else:
        print("[CHARLS] 未提供原始数据，跳过。可用 --charls 重新运行。")

    if args.ukb:
        print(f"[UKB] 读取 {args.ukb}")
        df = pd.read_csv(args.ukb, sep="\t", low_memory=False, nrows=2000) \
            if args.ukb.endswith(".tsv") else pd.read_csv(args.ukb, low_memory=False)
        print(f"  {df.shape[0]:,} 行 × {df.shape[1]} 列")
        m, r = detect_all("UKB", df.columns, UKB_FIELDS)
        for c, (s, info) in r.items():
            lines.append(f"| UKB | {c} | {s} | "
                         f"{info if isinstance(info, str) else '|'.join(info)} | "
                         f"{','.join(UKB_FIELDS[c][:2])} |")
        out = build_ukb(df, m)
        if out is not None:
            out.to_csv(os.path.join(args.out, "ukb.csv"), index=False)
            print(f"  → output/data/ukb.csv ({out.shape[0]:,} 行)")
    else:
        print("[UKB] 未提供原始数据，跳过。")

    rep = os.path.join(args.out, "field_matching_report.md")
    with open(rep, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"\n📋 字段核对报告：{rep}")


if __name__ == "__main__":
    main()
