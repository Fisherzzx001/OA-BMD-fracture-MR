"""
fields.py — 字段候选清单与疾病编码集（唯一需要按真实 Data Dictionary 回填的文件）
================================================================================
设计原则：
  CHARLS 各波字段名不统一，UKB 字段名形如 "41270-0.0"，且我对部分编号的记忆
  不敢打包票。因此本文件不硬编码"唯一正确"的字段名，而是给出**候选模式列表**，
  由 02_harmonise.py 逐个探测，并把匹配结果写进《字段核对报告》。
  凡是自动探测失败或匹配到多个候选的字段，报告里会标红，必须由你人工确认。

⚠️ 标记含义：该条目来自记忆，编号可能错误，取数时须与官方 Data Dictionary 核对。
================================================================================
"""

# ================================================================= CHARLS
# 每个概念给出若干候选（不区分大小写；支持正则）。越靠前优先级越高。
CHARLS_FIELDS = {
    "pid":        ["^id$", "^householdid$", "^pid$"],
    "age":        ["^rage$", "^r[0-9]agey$", "^agey?$", "^ba002_1$", "^age$"],
    "sex":        ["^ragender$", "^gender$", "^ba000_w2_3$", "^ba000$", "^sex$"],
    "bmi":        ["^r[0-9]bmi$", "^bmi$", "^qi001$"],            # 若只有身高体重则另行计算
    "weight_cm":  ["^qi002$", "^height$"],
    "weight_kg":  ["^qi001$", "^weight$"],
    "waist":      ["^qi003.*", "^waist.*"],
    # OA 主表型 T1：医生诊断「关节炎或风湿病」
    "oa_any":     ["^da007_13_?$", "^r[0-9]arthre$", "^arthritis$"],
    # OA 次表型 T2 所需：疼痛部位（背/腰/颈/髋/膝）
    "pain_back":     ["^da042_8_?$", "^r[0-9]bpainback$"],        # ⚠️ 背
    "pain_waist":    ["^da042_9_?$", "^r[0-9]bpainwaist$"],       # ⚠️ 腰
    "pain_neck":     ["^da042_15_?$", "^r[0-9]bpainneck$"],       # ⚠️ 颈
    "pain_hip":      ["^da042_10_?$", "^r[0-9]bpainhip$"],        # ⚠️ 髋/臀
    "pain_knee":     ["^da042_12_?$", "^r[0-9]bpainknee$"],       # ⚠️ 膝
    "pain_shoulder": ["^da042_2_?$", "^r[0-9]bpainshld$"],        # ⚠️ 肩
    # 14 项自报慢病（DA007 系列）
    "hypertension":   ["^da007_1_?$", "^r[0-9]hibpe$"],
    "dyslipidaemia":  ["^da007_2_?$", "^r[0-9]dyslipe$", "^r[0-9]lipide$"],  # ⚠️
    "diabetes":       ["^da007_3_?$", "^r[0-9]diabe$"],
    "cancer":         ["^da007_4_?$", "^r[0-9]cancre$", "^r[0-9]cancr$"],
    "copd":           ["^da007_5_?$", "^r[0-9]lunge$"],
    "liver_disease":  ["^da007_6_?$", "^r[0-9]livere$"],          # ⚠️
    "heart_disease":  ["^da007_7_?$", "^r[0-9]hearte$"],
    "stroke":         ["^da007_8_?$", "^r[0-9]stroke$"],
    "kidney_disease": ["^da007_9_?$", "^r[0-9]kidneye$"],         # ⚠️
    "digestive":      ["^da007_10_?$", "^r[0-9]stomache$", "^r[0-9]digeste$"],  # ⚠️
    "depression":     ["^da007_11_?$", "^r[0-9]psyche$"],
    "dementia":       ["^da007_12_?$", "^r[0-9]memrye$", "^r[0-9]demente$"],    # ⚠️
    "asthma":         ["^da007_14_?$", "^r[0-9]asthmae$"],        # ⚠️
    # 排除用
    "rheumatoid":   ["^r[0-9]rheumae$", "^da007_[0-9]+_?$"],      # ⚠️ CHARLS 是否单列待核
    # 协变量
    "education":    ["^raeducl$", "^r[0-9]edu.*", "^bd001$"],
    "urban":        ["^urban_nbs$", "^urban$", "^r[0-9]urban$"],
    "smoke":        ["^da065$", "^r[0-9]smoke$", "^r[0-9]smoken$"],
    "alcohol":      ["^da069$", "^r[0-9]drink$", "^r[0-9]alcohol$"],
    "stratum":      ["^wave$", "^province$", "^communityid$"],
    # 复杂抽样权重
    "samp_weight":  ["^r[0-9]wtresp$", "^.*wtresp.*", "^weight$"],
    "psu":          ["^communityid$", "^villageid$", "^psu$"],
    "strata":       ["^stratum$", "^strata$", "^province$"],
}

# ================================================================= UK Biobank
# UKB 导出列为 "eid" 与 "<field>-<instance>.<array>"，故用前缀匹配。
UKB_FIELDS = {
    "pid":       ["^eid$"],
    "age":       ["^21022"],          # Age at recruitment
    "sex":       ["^31-"],            # Sex
    "bmi":       ["^21001"],          # BMI
    "waist":     ["^21002"],          # Waist circumference
    "smoke":     ["^20116"],          # Smoking status
    "alcohol":   ["^1558"],           # Alcohol intake frequency
    "education": ["^6138"],           # Qualifications
    "tdi":       ["^189-"],           # Townsend deprivation index
    "centre":    ["^54-"],            # Assessment centre  → 匹配用 stratum
    "ethnic":    ["^21000"],          # Ethnic background
    # 诊断（数组型，多列）
    "icd10_all":   ["^41270"],        # Diagnoses - ICD10 (住院 HES)
    "icd10_main":  ["^41271"],        # Diagnoses - main ICD10
    "icd10_date":  ["^41280"],        # Date of first in-patient diagnosis
    "self_illness": ["^20002"],       # Non-cancer illness code, self-reported
    "medication":   ["^20003"],       # Treatment/medication code
    "pain_types":   ["^3773"],        # ⚠️ Pain types experienced in last month
    "opcs":         ["^41272"],       # ⚠️ OPCS operative procedures（关节置换）
}

# ================================================================= ICD-10 编码集
# 用于 UKB 臂从 ICD-10 判定各慢病。前缀匹配（如 "I10" 命中 I10, I10.0, I100）
ICD10_SETS = {
    "hypertension":   ["I10", "I11", "I12", "I13", "I15"],
    "diabetes":       ["E10", "E11", "E12", "E13", "E14"],
    "dyslipidaemia":  ["E78"],
    "heart_disease":  ["I20", "I21", "I22", "I23", "I24", "I25", "I50"],
    "stroke":         ["I60", "I61", "I62", "I63", "I64"],
    "asthma":         ["J45", "J46"],
    "copd":           ["J40", "J41", "J42", "J43", "J44"],
    "cancer":         ["C"],                                   # C00–C97
    "liver_disease":  ["K70", "K71", "K72", "K73", "K74",
                       "K75", "K76", "K77", "B18"],
    "kidney_disease": ["N00", "N01", "N02", "N03", "N04", "N05", "N06", "N07",
                       "N08", "N10", "N11", "N12", "N13", "N14", "N15", "N16",
                       "N17", "N18", "N19", "N25", "N26", "N27", "N28", "N29"],
    # 消化：K21–K93 但剔除肝（K70–K77）与肾（不在 K 段）
    "digestive":      ["K21", "K22", "K23", "K24", "K25", "K26", "K27", "K28",
                       "K29", "K30", "K31", "K50", "K51", "K52", "K55", "K56",
                       "K57", "K58", "K59", "K60", "K61", "K62", "K63",
                       "K80", "K81", "K82", "K83", "K85", "K86", "K87",
                       "K90", "K91", "K92", "K93"],
    "depression":     ["F32", "F33"],
    "dementia":       ["F00", "F01", "F02", "F03", "G30"],
    "back_pain":      ["M54", "M51"],                          # 背痛 + 胸腰椎间盘
    "neck_pain":      ["M542", "M75"],                         # ⚠️ 颈痛 + 肩部病损
    "multisite_pain": ["M255", "M797"],                        # ⚠️ 关节痛 + 纤维肌痛
    # 骨关节炎（T1 主表型）
    "oa":             ["M15", "M16", "M17", "M18", "M19"],
    "oa_knee":        ["M17"],
    "oa_hip":         ["M16"],
    # 排除用（炎症性关节病 / 痛风）
    "inflammatory_arthritis": ["M05", "M06", "M08", "M45", "L405", "M10"],
    "osteoporosis":   ["M80", "M81", "M82"],                   # ⚠️ CHARLS 侧需另核
}

# ================================================================= UKB 自报编码
# ⚠️ 以下为记忆中的 UKB Field 20002 编码，必须核对。未命中不影响主流程（ICD 为主）。
UKB_SELFREPORT = {
    "oa":             [1465],            # ⚠️ osteoarthritis
    "hypertension":   [1065, 1072],      # ⚠️
    "diabetes":       [1220, 1221, 1222, 1223],
    "dyslipidaemia":  [1473],            # ⚠️
    "heart_disease":  [1071, 1074, 1075],   # ⚠️ heart attack / angina
    "stroke":         [1081, 1491],      # ⚠️
    "asthma":         [1111],
    "copd":           [1112, 1113],      # ⚠️
    "cancer":         [],                # 自报肿瘤另有字段体系，暂不用
    "depression":     [1286, 1290],      # ⚠️
    "dementia":       [1263],
    "osteoporosis":   [1538],            # ⚠️
    "back_pain":      [1430],            # ⚠️
}

# ================================================================= CHARLS 疼痛部位 → 节点
# 哪些疼痛部位计入"多部位关节痛"（≥2 个部位）
CHARLS_MULTISITE_SITES = ["pain_neck", "pain_shoulder", "pain_back",
                          "pain_waist", "pain_hip", "pain_knee"]
