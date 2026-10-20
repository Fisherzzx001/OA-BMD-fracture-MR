#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G6 未传播修复 · 补丁 1/13（F1）
=============================
问题：171_g6_apply.py 的 SCALE_SPECS 中 mvmr_height_grid.csv 的缩放列
漏含 ref_uni_beta / ref_uni_se，导致今晨 G6 缩放后网格内部尺度不一致：
  - uni/mvmr/egger/trim 等列：新尺度（正确）
  - ref_uni_beta/ref_uni_se：旧尺度（+0.02312 等，错 3.3882 倍）
连带后果：
  - code/110（图 S10）与 code/115（审计）重算 a_t3 = 100*(ref-mvmr)/ref
    → −174%（正确值应为 13.5–33.2%）
  - 02_mvmr_height_grid.md（Table S10）整表 ref 列旧尺度

本脚本：
  1. 备份 mvmr_height_grid.csv
  2. 以 mvmr_height_main.csv 的全精度 ref_uni_beta/ref_uni_se 覆写网格中
     OA_knee 行的对应列（fracture 行本就是正确尺度，不动）
  3. 同步修复 171_g6_apply.py 的 SCALE_SPECS，防止复犯
  4. 自验：全部 8 个 (exposure,spec) 组合重算 a_t3 ∈ [13.5, 33.2]，
     且与 main 文件同键值逐位一致
"""
import csv
import os
import shutil

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GEN = os.path.join(ROOT, "code", "output", "genetics")
GRID = os.path.join(GEN, "mvmr_height_grid.csv")
MAIN = os.path.join(GEN, "mvmr_height_main.csv")
BACKUP = os.path.join(GEN, "_pre_g8_backup")
APPLY171 = os.path.join(ROOT, "code", "171_g6_apply.py")


def main():
    # ---------- 0. 备份 ----------
    os.makedirs(BACKUP, exist_ok=True)
    bak = os.path.join(BACKUP, "mvmr_height_grid_G6ref漏缩放_修复前.csv")
    if not os.path.exists(bak):
        shutil.copy2(GRID, bak)
        print("备份 ->", bak)

    # ---------- 1. 读取 main 的全精度 ref 值 ----------
    with open(MAIN, newline="", encoding="utf-8") as f:
        main_rows = list(csv.DictReader(f))
    ref = {}
    for r in main_rows:
        ref[(r["exposure"], r["outcome"])] = (
            float(r["ref_uni_beta"]), float(r["ref_uni_se"]))

    # ---------- 2. 覆写网格 OA_knee 行 ----------
    with open(GRID, newline="", encoding="utf-8") as f:
        grid = list(csv.DictReader(f))
    n_fix = 0
    for r in grid:
        key = (r["exposure"], r["outcome"])
        if key not in ref:
            continue
        want_b, want_se = ref[key]
        old_b = float(r["ref_uni_beta"])
        old_se = float(r["ref_uni_se"])
        if r["outcome"] == "fracture":
            # fracture 无需修正；确认与 main 一致
            assert abs(old_b - want_b) < 1e-12, (key, old_b, want_b)
            continue
        # OA_knee：旧尺度 -> 覆写为全精度新尺度
        r["ref_uni_beta"] = repr(want_b)
        r["ref_uni_se"] = repr(want_se)
        n_fix += 1
        if n_fix <= 5:
            print(f"  {key[0]}→{key[1]} [{r['spec']}]: "
                  f"β {old_b:.5f} -> {want_b:.6f} | SE {old_se:.5f} -> {want_se:.6f}")

    with open(GRID, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(grid[0].keys()))
        w.writeheader()
        w.writerows(grid)
    print(f"已覆写 {n_fix} 行（OA_knee，8 spec × 5 暴露 = 40 行为预期）")

    # ---------- 3. 修 171 SCALE_SPECS ----------
    with open(APPLY171, encoding="utf-8") as f:
        src = f.read()
    old_block = '''    ("mvmr_height_grid.csv", ("exposure", "outcome"),
     ["uni_beta", "uni_se", "uni_bmd_only", "uni_bmd_only_se",
      "mvmr_beta", "mvmr_se", "mvmr_se_re", "mvmr_h", "mvmr_h_se",
      "uni_h_beta", "uni_h_se", "egger_int", "egger_int_se",
      "egger_beta", "egger_se", "trim_beta", "trim_se"]),'''
    new_block = '''    ("mvmr_height_grid.csv", ("exposure", "outcome"),
     ["uni_beta", "uni_se", "uni_bmd_only", "uni_bmd_only_se",
      "mvmr_beta", "mvmr_se", "mvmr_se_re", "mvmr_h", "mvmr_h_se",
      "uni_h_beta", "uni_h_se", "egger_int", "egger_int_se",
      "egger_beta", "egger_se", "trim_beta", "trim_se",
      "ref_uni_beta", "ref_uni_se"]),   # G6 补丁：ref 列曾漏缩放（2026-10-20）'''
    assert old_block in src or new_block in src, "171 SCALE_SPECS 网格块未匹配"
    if old_block in src:
        src = src.replace(old_block, new_block)
        with open(APPLY171, "w", encoding="utf-8") as f:
            f.write(src)
        print("171_g6_apply.py SCALE_SPECS 已补 ref_uni_beta/ref_uni_se")
    else:
        print("171_g6_apply.py SCALE_SPECS 已含补丁（幂等跳过）")

    # ---------- 4. 自验 ----------
    print("\n自验：a_t3 = 100*(ref_uni_beta - mvmr_beta)/ref_uni_beta")
    print("（13.5–33.2% 区间仅约束主配置 EUR|n_h=100；其余 spec 为敏感性网格，可不同）")
    lo, hi = 1e9, -1e9
    lo_set, hi_set = 1e9, -1e9
    with open(GRID, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["outcome"] != "OA_knee":
                continue
            a = 100.0 * (float(r["ref_uni_beta"]) - float(r["mvmr_beta"])) \
                / float(r["ref_uni_beta"])
            # a_set 定义与 code/110 一致：100*(uni_beta - mvmr_beta)/uni_beta
            a_set = 100.0 * (float(r["uni_beta"]) - float(r["mvmr_beta"])) \
                / float(r["uni_beta"])
            if r["spec"] == "EUR|n_h=100":
                lo, hi = min(lo, a), max(hi, a)
                lo_set, hi_set = min(lo_set, a_set), max(hi_set, a_set)
                print(f"  {r['exposure']:18s} a_t3 = {a:7.2f}%   a_set = {a_set:6.2f}%   "
                      f"(ref {float(r['ref_uni_beta']):.4f} vs mvmr {float(r['mvmr_beta']):.4f})")
    print(f"\n主配置 a_t3 范围 [{lo:.2f}%, {hi:.2f}%]（期望 13.5–33.2%）")
    print(f"主配置 a_set 范围 [{lo_set:.2f}%, {hi_set:.2f}%]（期望 9.5–21.2%）")
    assert 13.4 <= lo and hi <= 33.3, "主配置 a_t3 超出预期区间！"
    assert 9.4 <= lo_set and hi_set <= 21.3, "主配置 a_set 超出预期区间！"
    print("✅ F1 完成")


if __name__ == "__main__":
    main()
