"""
run_all.py — 一键跑通整条管线
================================================================================
用法：
    python run_all.py            # 模拟数据，验证管线（默认）
    python run_all.py --real     # 真实数据（需先放好 data/charls.csv, data/ukb.csv）

阶段：
    00 模拟/载入数据 → 01 队列与匹配 → 03 描述·网络·推断 → 04 出图
================================================================================
"""

import os
import sys
import time
import importlib

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
os.chdir(HERE)

from config import OUT_DIR


def main():
    t0 = time.time()
    real = "--real" in sys.argv
    if not real:
        print("=" * 74)
        print("⚠️  当前为【模拟数据】运行——仅用于验证管线，结果不得写入论文")
        print("=" * 74)

    os.makedirs(OUT_DIR, exist_ok=True)

    c01 = importlib.import_module("01_build_cohort")
    c03 = importlib.import_module("03_network_analysis")
    c04 = importlib.import_module("04_make_outputs")

    print("\n【阶段 1/3】队列构建与 1:1 匹配")
    c01.main(use_real=real)

    print("\n【阶段 2/3】描述 → 疾病对 → 共病网络 → 推断")
    built, k, built_k = c03.main()

    print("\n【阶段 3/3】图件输出")
    c04.main(built, k, built_k)

    print(f"\n✅ 全部完成，用时 {time.time() - t0:.1f}s；产出目录：{os.path.abspath(OUT_DIR)}")


if __name__ == "__main__":
    main()
