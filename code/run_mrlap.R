#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# run_mrlap.R  --  Sample-overlap (MRlap) correction for the 10 OA/BMD/fracture MR pairs
#
# Real MRlap API (n-mounier/MRlap v0.0.3.3):
#   MRlap(exposure, exposure_name, outcome, outcome_name,
#         ld = LDSC-weights folder, hm3 = SNP/A1/A2 file,
#         do_pruning = TRUE, MR_pruning_LD = 0  -> distance pruning (no plink),
#         MR_threshold = 5e-8)
# Output: res$MRcorrection (observed/corrected effect, p_difference, m_IVs, egger_*),
#         res$LDSC (rg, int_crosstrait=lambda, h2_exp, h2_out)
#
# Inputs are produced by code/prep_mrlap_gwas.py (+ code/fix_mrlap_inputs.py):
#   code/ref/mrlap_gwas/<trait>.tsv  with columns rsid chr pos a1 a2 beta se N
#   * case-control OA_knee is pre-converted to log-OR and labelled "beta" on purpose:
#     MRlap's tidy_inputGWAS has a bug in its "or" branch (dplyr::mutate called without
#     a data argument -> "Can't subset .data outside of a data mask context").
#   * N = manuscript Table-1 values: BMD 31,873 / fracture 264,973 / OA_knee 1,208,397.
#
# run_LDSC always treats traits as continuous (sample.prev=NA): for the case-control
# OA_knee trait the genetic correlation is therefore on the observed (log-OR) scale.
# This is a known MRlap limitation and is documented in the manuscript (Sec 2.8).
# ---------------------------------------------------------------------------

.libPaths("H:/R/R-4.6.1/library")
suppressMessages({library(MRlap); library(data.table)})

BASE    <- "C:/Users/Administrator/WorkBuddy/2026-09-08-22-17-40/OA_B_PROJECT"
GWASDIR <- file.path(BASE, "code/ref/mrlap_gwas")
LD      <- file.path(BASE, "code/ref/ldsc_ref/eur_w_ld_chr")   # GenomicSEM expects <ld>/<chr>.l2.ldscore.gz
HM3     <- file.path(BASE, "code/ref/ldsc_ref/w_hm3.snplist")  # header: SNP A1 A2
OUTDIR  <- file.path(BASE, "code/ref/mrlap_out")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
setwd(OUTDIR)   # MRlap writes temporary .tsv/.sumstats.gz files here

cat(sprintf("MRlap driver\n  GWASDIR=%s\n  LD=%s\n  HM3=%s\n  OUTDIR=%s\n",
            GWASDIR, LD, HM3, OUTDIR))

# 10 MR pairs (exposure -> outcome), identical to code/13_mr_analysis.py :: PAIRS
pairs <- list(
  c("BMD_total",      "fracture"),
  c("BMD_femur_neck", "fracture"),
  c("BMD_total",      "OA_knee"),
  c("BMD_femur_neck", "OA_knee"),
  c("BMD_lumbar",     "OA_knee"),
  c("BMD_arms",       "OA_knee"),
  c("BMD_leg",        "OA_knee"),
  c("OA_knee",        "fracture"),
  c("OA_knee",        "BMD_total"),
  c("OA_knee",        "BMD_femur_neck")
)

# optional: only run a subset, e.g. Rscript run_mrlap.R 3:10
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  idx <- eval(parse(text = args[1]))
  pairs <- pairs[idx]
  cat(sprintf("  -> running subset: %s\n", paste(seq_along(pairs), collapse = ",")))
}

csv_path <- file.path(OUTDIR, "mrlap_overlap_results.csv")
rows <- list()

for (k in seq_along(pairs)) {
  exp <- pairs[[k]][1]; out <- pairs[[k]][2]
  ef  <- file.path(GWASDIR, paste0(exp, ".tsv"))
  of  <- file.path(GWASDIR, paste0(out, ".tsv"))
  cat(sprintf("\n===== [%02d/%02d] %s -> %s =====\n", k, length(pairs), exp, out))

  # Pre-populate ALL result columns so a failed pair still has the right shape.
  row <- data.frame(
    exposure = exp, outcome = out,
    observed_effect = NA_real_, observed_effect_se = NA_real_, observed_effect_p = NA_real_,
    corrected_effect = NA_real_, corrected_effect_se = NA_real_, corrected_effect_p = NA_real_,
    p_difference = NA_real_, m_IVs = NA_integer_, egger_intercept_p = NA_real_,
    rg = NA_real_, int_crosstrait = NA_real_, h2_exp = NA_real_, h2_out = NA_real_,
    status = "", stringsAsFactors = FALSE)

  if (!file.exists(ef) || !file.exists(of)) {
    row$status <- sprintf("MISSING INPUT: exp=%s out=%s", file.exists(ef), file.exists(of))
  } else {
    tryCatch({
      res <- MRlap(ef, exposure_name = exp, of, outcome_name = out,
                   ld = LD, hm3 = HM3,
                   do_pruning = TRUE, MR_pruning_LD = 0, MR_threshold = 5e-8,
                   verbose = TRUE, save_logfiles = FALSE)
      mc       <- res$MRcorrection
      ldsc_res <- res$LDSC          # NB: do NOT name this `ld` (shadows the LD path)
      row$observed_effect      <- mc$observed_effect
      row$observed_effect_se   <- mc$observed_effect_se
      row$observed_effect_p    <- mc$observed_effect_p
      row$corrected_effect     <- mc$corrected_effect
      row$corrected_effect_se  <- mc$corrected_effect_se
      row$corrected_effect_p   <- mc$corrected_effect_p
      row$p_difference         <- mc$p_difference
      row$m_IVs                <- mc$m_IVs
      row$egger_intercept_p    <- mc$egger_intercept_p
      row$rg                   <- ldsc_res$rg
      row$int_crosstrait       <- ldsc_res$int_crosstrait
      row$h2_exp               <- ldsc_res$h2_exp
      row$h2_out               <- ldsc_res$h2_out
      row$status               <- "ok"
      cat(sprintf("  observed=%.4f  corrected=%.4f  p_diff=%.3g  rg=%.3f  h2_exp=%.3f  m_IVs=%d\n",
                  mc$observed_effect, mc$corrected_effect, mc$p_difference,
                  ldsc_res$rg, ldsc_res$h2_exp, mc$m_IVs))
    }, error = function(e) {
      # `<<-` so the assignment lands in the loop environment, not the handler frame
      row$status <<- paste0("ERROR: ", conditionMessage(e))
      cat("  ERROR:", conditionMessage(e), "\n")
    })
  }

  rows[[length(rows) + 1L]] <- row
  # rbindlist(fill=TRUE) tolerates ragged rows -> one bad pair cannot abort the loop
  res_df <- data.table::rbindlist(rows, fill = TRUE)
  write.csv(res_df, csv_path, row.names = FALSE)
  cat(sprintf("  -> saved %d rows to %s\n", nrow(res_df), csv_path))
}

cat("\nALL DONE\n")
print(data.table::rbindlist(rows, fill = TRUE))
