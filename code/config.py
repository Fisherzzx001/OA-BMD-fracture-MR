"""
config.py — B 档项目全局配置
================================================================================
修改本文件即可切换分析口径，无需改动分析代码。

【真实数据上线流程】
1. 将协调后的数据放到 data/ 目录：
     data/charls.csv, data/ukb.csv
   必需列见 REQUIRED_COLUMNS；疾病列名为 DISEASES 的 key。
2. 命令行：python run_all.py --real
   不带 --real 时使用模拟数据（仅用于验证管线，绝不可用于论文结果）。
================================================================================
"""

import numpy as np

# ---------------------------------------------------------------- 种子与路径
SEED = 20260908
DATA_DIR = "data"
OUT_DIR = "output"

# ---------------------------------------------------------------- 共同节点集
# (key, 英文标签, 模块)
DISEASES = [
    ("hypertension",    "Hypertension",           "metabolic"),
    ("diabetes",        "Diabetes",               "metabolic"),
    ("dyslipidaemia",   "Dyslipidaemia",          "metabolic"),
    ("heart_disease",   "Heart disease",          "cardiovascular"),
    ("stroke",          "Stroke",                 "cardiovascular"),
    ("asthma",          "Asthma",                 "respiratory"),
    ("copd",            "Chronic lung disease",   "respiratory"),
    ("liver_disease",   "Liver disease",          "organ_other"),
    ("kidney_disease",  "Kidney disease",         "organ_other"),
    ("digestive",       "Digestive disease",      "organ_other"),
    ("cancer",          "Cancer",                 "organ_other"),
    ("depression",      "Depression",             "mental"),
    ("dementia",        "Dementia",               "mental"),
    ("back_pain",       "Chronic back pain",      "musculoskeletal"),
    ("neck_pain",       "Neck/shoulder pain",     "musculoskeletal"),
    ("multisite_pain",  "Multisite joint pain",   "musculoskeletal"),
]

DKEYS = [d[0] for d in DISEASES]
DLABEL = {d[0]: d[1] for d in DISEASES}
DMODULE = {d[0]: d[2] for d in DISEASES}
MODULES = ["metabolic", "cardiovascular", "respiratory",
           "organ_other", "mental", "musculoskeletal"]

# 模块配色（Figure 2 用；色盲友好）
MODULE_COLOR = {
    "metabolic":        "#C0392B",   # 红 —— 代谢
    "cardiovascular":   "#8E44AD",   # 紫 —— 心血管
    "respiratory":      "#2980B9",   # 蓝 —— 呼吸
    "organ_other":      "#7F8C8D",   # 灰 —— 其他器官
    "mental":           "#16A085",   # 青 —— 精神神经
    "musculoskeletal":  "#D35400",   # 橙 —— 骨-关节
}

# ---------------------------------------------------------------- 分析参数
PARAMS = {
    # --- 纳入 / 匹配 ---
    "min_age": 45,
    "match_age_window": 1,        # ±1 岁
    "match_ratio": 1,             # 1:1
    "match_on": ["sex", "stratum"],  # stratum = CHARLS 波次 / UKB 招募中心

    # --- 网络 ---
    "min_prevalence": 0.01,       # 患病率 > 1% 入网
    "oe_threshold": 1.0,          # O/E > 1 才可能成边
    "alpha": 0.05,                # FDR 阈值
    "equal_density": True,        # 等密度稀疏化（跨库可比的前提）

    # --- 推断 ---
    "n_permutation": 2000,        # 生产环境请设 10000
    "n_bootstrap": 2000,          # 生产环境请设 2000–5000

    # --- 模拟数据规模 ---
    "sim_n_ukb": 40000,
    "sim_n_charls": 12000,
}

REQUIRED_COLUMNS = ["pid", "age", "sex", "bmi", "oa", "stratum", "weight"] + DKEYS

# ---------------------------------------------------------------- 模拟参数
# 仅用于 generate_synthetic_data()：让模拟数据呈现"中国代谢驱动 / 英国骨关节驱动"
# 的预设结构，以便验证管线能否检出该差异。**不得用于任何论文结果。**
SIM = {
    "ukb": {
        "n": PARAMS["sim_n_ukb"],
        "age_mean": 62, "age_sd": 8,
        "bmi_mean": 27.5, "bmi_sd": 4.5,
        "global_load": 0.42,
        "oa_prob_base": -1.15, "oa_b_age": 0.45, "oa_b_bmi": 0.38, "oa_b_sex": 0.15,
        "factor_load": {"metabolic": 0.85, "musculoskeletal": 1.60,
                        "cardiovascular": 0.85, "respiratory": 0.70,
                        "organ_other": 0.75, "mental": 0.85},
        "oa_extra_by_module": {"metabolic": 0.30, "musculoskeletal": 1.05,
                               "cardiovascular": 0.20, "respiratory": 0.15,
                               "organ_other": 0.10, "mental": 0.35},
        "prev_scale": {"metabolic": 1.00, "musculoskeletal": 1.15,
                       "cardiovascular": 1.00, "respiratory": 1.00,
                       "organ_other": 0.90, "mental": 1.00},
    },
    "charls": {
        "n": PARAMS["sim_n_charls"],
        "age_mean": 60, "age_sd": 9,
        "bmi_mean": 24.5, "bmi_sd": 3.8,
        "global_load": 0.42,
        "oa_prob_base": -1.30, "oa_b_age": 0.50, "oa_b_bmi": 0.30, "oa_b_sex": 0.25,
        "factor_load": {"metabolic": 1.55, "musculoskeletal": 0.65,
                        "cardiovascular": 0.80, "respiratory": 0.65,
                        "organ_other": 0.75, "mental": 0.85},
        "oa_extra_by_module": {"metabolic": 1.00, "musculoskeletal": 0.35,
                               "cardiovascular": 0.45, "respiratory": 0.20,
                               "organ_other": 0.15, "mental": 0.40},
        "prev_scale": {"metabolic": 1.20, "musculoskeletal": 0.80,
                       "cardiovascular": 0.95, "respiratory": 0.85,
                       "organ_other": 1.00, "mental": 0.90},
    },
}

# 各病基础患病率（logit 尺度的截距，模拟用）
SIM_BASE_LOGIT = {
    "hypertension": -0.30, "diabetes": -1.90, "dyslipidaemia": -1.10,
    "heart_disease": -2.20, "stroke": -3.00, "asthma": -2.90, "copd": -2.60,
    "liver_disease": -3.20, "kidney_disease": -3.10, "digestive": -2.40,
    "cancer": -3.40, "depression": -2.30, "dementia": -3.60,
    "back_pain": -1.40, "neck_pain": -1.90, "multisite_pain": -2.00,
}

rng = np.random.default_rng(SEED)
