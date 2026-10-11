"""
06_gwas_qc.py — UKB GWAS 汇总统计量质控管线
============================================
输入：12 个 Neale-lab 格式 GWAS .tsv（chromosome, base_pair_location,
      effect_allele, other_allele, beta, standard_error,
      effect_allele_frequency, p_value, rs_id）
输出：output/genetics/manifest.csv        —— 12 表型质控总表
      output/genetics/QC_<pheno>.txt       —— 逐表型明细
      output/genetics/tophits_<pheno>.csv  —— 各表型 Top20 信号

注意：本文件无样本量(N)列，故 lambda_gc 用全 SNP 中位 chi2 估计；
      MR / LDSC 所需的 N 与 LD 参考需另行补充（见 docs/06_遗传学分析方案.md）。
"""
import os
import glob
import numpy as np
import pandas as pd

SRC_DIR = r"C:\Users\Administrator\Desktop\UBK"
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "genetics")
os.makedirs(OUT_DIR, exist_ok=True)

PHENO_LABEL = {
    "Total_BMD": "Total body BMD",
    "Arms_BMD": "Arms BMD",
    "Femur_neck_BMD": "Femoral neck BMD",
    "Femur_total_BMD": "Femoral total BMD",
    "LS_BMD": "Lumbar spine BMD",
    "Leg_BMD": "Leg BMD",
    "Pelvis_BMD": "Pelvis BMD",
    "Rib_BMD": "Rib BMD",
    "Spine_BMD": "Spine BMD",
    "Trunk_BMD": "Trunk BMD",
    "head_BMD": "Head BMD",
    "fracture_any-type": "Any-type fracture",
}
USECOLS = ["chromosome", "base_pair_location", "effect_allele", "other_allele",
           "beta", "standard_error", "effect_allele_frequency", "p_value", "rs_id"]
DTYPE = {"chromosome": "object", "base_pair_location": "int64",
         "effect_allele": "category", "other_allele": "category",
         "beta": "float64", "standard_error": "float64",
         "effect_allele_frequency": "float64", "p_value": "float64",
         "rs_id": "object"}


def qc_one(path):
    df = pd.read_csv(path, sep="\t", usecols=USECOLS, dtype=DTYPE)
    n = len(df)
    se = df["standard_error"].to_numpy()
    beta = df["beta"].to_numpy()
    p = df["p_value"].to_numpy()
    eaf = df["effect_allele_frequency"].to_numpy()
    valid = np.isfinite(se) & (se > 0) & np.isfinite(beta) & np.isfinite(p)
    chi2 = np.where(valid, (beta / np.where(se > 0, se, np.nan)) ** 2, np.nan)

    # lambda_gc：用常见变异（0.01<eaf<0.99）且 se>0 的中位 chi2
    m = valid & np.isfinite(eaf) & (eaf > 0.01) & (eaf < 0.99)
    lam = float(np.nanmedian(chi2[m]) / 0.456) if m.sum() > 100 else np.nan

    n_gws = int((p < 5e-8).sum())
    n_sug = int((p < 1e-6).sum())
    n_pal = int(((df["effect_allele"].isin(["A", "T"])) &
                 (df["other_allele"].isin(["A", "T"]))).sum() +
                ((df["effect_allele"].isin(["C", "G"])) &
                 (df["other_allele"].isin(["C", "G"]))).sum())
    n_extreme = int((eaf < 0.001).sum() + (eaf > 0.999).sum())
    n_dup = int(df["rs_id"].duplicated().sum())
    n_bad = int((~valid).sum())

    lead_i = int(np.nanargmin(p)) if n_gws > 0 else int(np.nanargmin(p))
    lead = df.iloc[lead_i]

    top = df.loc[p < 1e-5].sort_values("p_value").head(20)
    top_out = top[["chromosome", "base_pair_location", "effect_allele",
                   "other_allele", "beta", "standard_error",
                   "effect_allele_frequency", "p_value", "rs_id"]].copy()

    return dict(n=n, n_gws=n_gws, n_sug=n_sug, lambda_gc=lam,
                n_palindromic=n_pal, n_eaf_extreme=n_extreme,
                n_dup_rsid=n_dup, n_invalid=n_bad,
                min_p=float(np.nanmin(p)),
                lead_rsid=str(lead["rs_id"]),
                lead_chr=str(lead["chromosome"]),
                lead_pos=int(lead["base_pair_location"]),
                lead_beta=float(lead["beta"]), lead_p=float(lead["p_value"]),
                mean_abs_beta=float(np.nanmean(np.abs(beta))),
                top=top_out)


def main():
    files = sorted(glob.glob(os.path.join(SRC_DIR, "*.tsv")))
    rows = []
    print(f"▶ GWAS QC：发现 {len(files)} 个文件\n")
    for path in files:
        base = os.path.splitext(os.path.basename(path))[0]
        label = PHENO_LABEL.get(base, base)
        r = qc_one(path)
        print(f"  [{label:22s}] n={r['n']:,}  GWS={r['n_gws']:,}  "
              f"lambda_gc={r['lambda_gc']:.3f}  minP={r['min_p']:.2e}")
        r2 = {k: v for k, v in r.items() if k != "top"}
        r2["pheno"] = base
        r2["label"] = label
        rows.append(r2)
        r["top"].to_csv(os.path.join(OUT_DIR, f"tophits_{base}.csv"), index=False)
        with open(os.path.join(OUT_DIR, f"QC_{base}.txt"), "w") as f:
            f.write(f"Phenotype: {label} ({base})\n")
            for k, v in r2.items():
                if k in ("pheno", "label", "top"):
                    continue
                f.write(f"  {k:18s}: {v}\n")
    man = pd.DataFrame(rows)[["pheno", "label", "n", "n_gws", "n_sug",
                              "lambda_gc", "min_p", "n_palindromic",
                              "n_eaf_extreme", "n_dup_rsid", "n_invalid",
                              "lead_rsid", "lead_chr", "lead_pos",
                              "lead_beta", "lead_p", "mean_abs_beta"]]
    man.to_csv(os.path.join(OUT_DIR, "manifest.csv"), index=False)
    print(f"\n✅ manifest.csv 已写出 → {OUT_DIR}")
    # 提示：无 N 列，无法直接计算有效样本量
    print("⚠️  注意：源文件无样本量(N)列，lambda_gc 为全 SNP 中位估计；"
          "MR/LDSC 需补充各表型 N 与 LD 参考。")


if __name__ == "__main__":
    main()
