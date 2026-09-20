#!/usr/bin/env Rscript
# =============================================================================
# mrlap_overlap_correction.R
# Sample-overlap correction for the 10 OA/BMD/fracture two-sample MR pairs,
# using the MRlap package (MR-Club; https://github.com/MRCIEU/MRlap).
#
# PURPOSE
#   All 10 exposure-outcome pairs share UK Biobank participants in their GWAS
#   summary statistics. This induces sample overlap that biases the IVW estimate
#   toward the corresponding observational genetic correlation (r_g). The main
#   manuscript (docs/08, S2.8) bounds this bias with the Burgess 2016 formula
#   (IVW bias = overlap fraction f * r_g), giving conservative upper bounds of
#   2.6% (BMD->OA_knee) and 12% (BMD->fracture). MRlap provides an independent,
#   data-driven estimate of the overlap-correlation structure and an
#   overlap-corrected MR effect, serving as a sensitivity analysis.
#
# INPUT
#   code/output/genetics/harmonised/<EXP>__<OUT>.csv  (one file per pair)
#     columns: SNP, CHR, BP, beta_e, se_e, p_e, beta_o, se_o, p_o, eaf, ...
#     beta_e/se_e = exposure; beta_o/se_o = outcome (already allele-harmonised)
#   code/output/genetics/mr_corrected_N.csv  (naive IVW for comparison)
#
# OUTPUT
#   code/output/genetics/mrlap_overlap_results.csv
#     pair, n_snps, ivw_naive_beta, ivw_naive_se, ivw_naive_p,
#     mrlap_beta, mrlap_se, mrlap_p, overlap_rho, pct_bias_ivw_vs_mrlap
#
# REQUIREMENTS (NOT available in the build sandbox; run in a local R >= 4.1):
#   install.packages(c("data.table","dplyr","readr"))
#   remotes::install_github("MRCIEU/MRlap")
#   MRlap also pulls MendelianRandomization, MRInstruments, etc.
# =============================================================================

suppressMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(MRlap)
})

HARMONISED_DIR <- "code/output/genetics/harmonised"
OUT_CSV        <- "code/output/genetics/mrlap_overlap_results.csv"
NAIVE_CSV      <- "code/output/genetics/mr_corrected_N.csv"

# ---- discover the 10 pairs from the harmonised filenames --------------------
files <- list.files(HARMONISED_DIR, pattern = "^.+__.+\\.csv$", full.names = TRUE)
pairs <- data.frame(
  file    = files,
  exp     = sub("^.*/([^/]+)__([^/]+)\\.csv$", "\\1", files),
  out     = sub("^.*/([^/]+)__([^/]+)\\.csv$", "\\2", files),
  stringsAsFactors = FALSE
)
cat(sprintf("Found %d harmonised pair files.\n", nrow(pairs)))

# ---- naive IVW reference (reconstructed-N scale) ---------------------------
naive <- NULL
if (file.exists(NAIVE_CSV)) {
  naive <- fread(NAIVE_CSV)
  # harmonise column names: expect exposure, outcome, beta_ivw, se_ivw, p_ivw
  setnames(naive, tolower(colnames(naive)))
  keycols <- c("exposure", "outcome")
  ivw_cols <- c("beta_ivw", "se_ivw", "p_ivw")
  ivw_cols <- ivw_cols[ivw_cols %in% colnames(naive)]
  if (all(keycols %in% colnames(naive)) && length(ivw_cols) > 0) {
    naive <- naive[, c(keycols, ivw_cols), with = FALSE]
  } else {
    naive <- NULL
  }
}

# ---- per-pair MRlap overlap correction --------------------------------------
results <- list()
for (i in seq_len(nrow(pairs))) {
  pf  <- pairs$file[i]
  exp <- pairs$exp[i]
  out <- pairs$out[i]
  cat(sprintf("\n[%02d/%02d] %s -> %s\n", i, nrow(pairs), exp, out))

  dt <- fread(pf)
  req <- c("SNP", "beta_e", "se_e", "p_e", "beta_o", "se_o", "p_o")
  if (!all(req %in% colnames(dt))) {
    cat("  SKIP: missing required columns\n"); next
  }
  # MRlap format_data expects exposure/outcome frames with MR-base-like columns
  exp_df <- data.frame(
    SNP                   = dt$SNP,
    beta.exposure        = dt$beta_e,
    se.exposure          = dt$se_e,
    pval.exposure        = dt$p_e,
    effect_allele.exposure = "?",
    other_allele.exposure  = "?",
    stringsAsFactors = FALSE
  )
  out_df <- data.frame(
    SNP                  = dt$SNP,
    beta.outcome        = dt$beta_o,
    se.outcome          = dt$se_o,
    pval.outcome        = dt$p_o,
    effect_allele.outcome = "?",
    other_allele.outcome  = "?",
    stringsAsFactors = FALSE
  )

  row <- data.frame(pair = paste(exp, "->", out), n_snps = nrow(dt),
                    ivw_naive_beta = NA_real_, ivw_naive_se = NA_real_,
                    ivw_naive_p = NA_real_, mrlap_beta = NA_real_,
                    mrlap_se = NA_real_, mrlap_p = NA_real_,
                    overlap_rho = NA_real_, pct_bias = NA_real_)

  if (!is.null(naive)) {
    m <- naive[exposure == exp & outcome == out]
    if (nrow(m) == 1) {
      row$ivw_naive_beta <- as.numeric(m$beta_ivw)
      row$ivw_naive_se   <- as.numeric(m$se_ivw)
      row$ivw_naive_p    <- as.numeric(m$p_ivw)
    }
  }

  tryCatch({
    lap <- format_data(exposure = exp_df, outcome = out_df)
    fit <- correlated_exposures(lap, get.beta = TRUE, pval_threshold = 1)
    # correlated_exposures returns a list; extract the overlap-corrected MR beta
    if (!is.null(fit[["beta"]])) {
      row$mrlap_beta <- as.numeric(fit[["beta"]])
      row$mrlap_se   <- if (!is.null(fit[["se"]])) as.numeric(fit[["se"]]) else NA_real_
      row$mrlap_p    <- if (!is.null(fit[["pval"]])) as.numeric(fit[["pval"]]) else NA_real_
    }
    # overlap correlation (rho) between exposure and outcome z-stats, if reported
    if (!is.null(fit[["correlation"]])) {
      row$overlap_rho <- as.numeric(fit[["correlation"]])
    } else if (!is.null(fit[["rho"]])) {
      row$overlap_rho <- as.numeric(fit[["rho"]])
    }
    if (is.finite(row$ivw_naive_beta) && is.finite(row$mrlap_beta) &&
        abs(row$mrlap_beta) > 0) {
      row$pct_bias <- 100 * (row$ivw_naive_beta - row$mrlap_beta) / abs(row$mrlap_beta)
    }
    cat(sprintf("  OK: MRlap beta = %.4f (se %.4f), overlap rho = %.3f, IVW-naive vs MRlap bias = %.1f%%\n",
                row$mrlap_beta, row$mrlap_se, row$overlap_rho,
                ifelse(is.na(row$pct_bias), NaN, row$pct_bias)))
  }, error = function(e) {
    cat(sprintf("  MRlap error: %s\n", conditionMessage(e)))
  })

  results[[length(results) + 1]] <- row
}

out_tbl <- rbindlist(results, fill = TRUE)
fwrite(out_tbl, OUT_CSV)
cat(sprintf("\nWrote %s (%d pairs).\n", OUT_CSV, nrow(out_tbl)))

# ---- summary: does overlap correction move any well-powered estimate? ------
well_powered <- out_tbl[!is.na(ivw_naive_p) & ivw_naive_p < 1e-4, ]
cat(sprintf("\nWell-powered pairs (naive p < 1e-4): %d\n", nrow(well_powered)))
if (nrow(well_powered) > 0) {
  print(well_powered[, .(pair, ivw_naive_beta, mrlap_beta, overlap_rho, pct_bias)])
}
cat("\nInterpretation: if |pct_bias| for all well-powered pairs is << the robust-estimator\n")
cat("absorption already claimed in S2.8, the Burgess 2016 upper bounds are conservative.\n")
