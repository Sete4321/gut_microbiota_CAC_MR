# ============================================================
# 02_run_mr.R
# Two-sample MR: IVW, MR-Egger, Weighted Median, Wald Ratio
# Taxonomy audit → Bonferroni / FDR correction
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(TwoSampleMR)
})
source("R/config.R")

# ── Load prepared data ────────────────────────────────────────
mbg_inst   <- fread(file.path(CONFIG$paths$out_dir, "mbg_instruments.csv"))
dmp_inst   <- fread(file.path(CONFIG$paths$out_dir, "dmp_instruments.csv"))
meta_cac   <- fread(file.path(CONFIG$paths$out_dir, "meta_cac_outcome.csv"))

# ════════════════════════════════════════════════════════════
# 2A. HARMONISE EXPOSURE + OUTCOME
# ════════════════════════════════════════════════════════════

harmonise_exposure_outcome <- function(exp_df, out_df) {
  # TwoSampleMR requires specific column names
  exp_std <- exp_df %>%
    mutate(
      id.exposure  = exposure,
      mr_keep      = TRUE,
      units.exposure = "log CLR abundance"
    )

  out_std <- out_df %>%
    mutate(
      id.outcome  = "meta_CAC",
      outcome     = "meta_CAC",
      mr_keep     = TRUE,
      units.outcome = "log CAC score"
    )

  TwoSampleMR::harmonise_data(
    exposure_dat = exp_std,
    outcome_dat  = out_std,
    action       = 2   # 2 = try to infer strand for palindromic SNPs
  )
}

cat("\n── Harmonising MiBioGen instruments ──\n")
mbg_harm <- harmonise_exposure_outcome(mbg_inst, meta_cac)
cat(sprintf("  Retained %d harmonised SNP-taxon pairs\n", sum(mbg_harm$mr_keep)))

cat("\n── Harmonising DMP instruments ──\n")
dmp_harm <- harmonise_exposure_outcome(dmp_inst, meta_cac)
cat(sprintf("  Retained %d harmonised SNP-taxon pairs\n", sum(dmp_harm$mr_keep)))

# Combine all harmonised data
all_harm <- bind_rows(mbg_harm, dmp_harm) %>% filter(mr_keep)

# ════════════════════════════════════════════════════════════
# 2B. RUN MR PER TAXON
# ════════════════════════════════════════════════════════════

cat("\n── Running MR for", length(unique(all_harm$exposure)), "taxa ──\n")

run_mr_taxon <- function(df, taxon_name) {
  df_t <- df %>% filter(exposure == taxon_name, mr_keep)
  n_iv <- nrow(df_t)

  if (n_iv == 0) return(NULL)

  # IVW (or Wald ratio if n_iv == 1)
  if (n_iv == 1) {
    b_ivw  <- df_t$beta.outcome / df_t$beta.exposure
    se_ivw <- abs(df_t$se.outcome / df_t$beta.exposure)
    p_ivw  <- 2 * pnorm(-abs(b_ivw / se_ivw))
    q_stat <- NA_real_; q_p <- NA_real_
    method <- "Wald ratio"
  } else {
    # Weighted IVW
    w      <- 1 / df_t$se.outcome^2
    b_ivw  <- sum(w * df_t$beta.exposure * df_t$beta.outcome) /
               sum(w * df_t$beta.exposure^2)
    # IVW SE (random-effects if Q > df)
    res_iv <- df_t$beta.outcome - b_ivw * df_t$beta.exposure
    q_stat <- sum(w * res_iv^2)
    q_df   <- n_iv - 1
    q_p    <- pchisq(q_stat, df = q_df, lower.tail = FALSE)
    # SE with over-dispersion correction if heterogeneous
    phi    <- max(1, q_stat / q_df)
    se_ivw <- sqrt(phi / sum(w * df_t$beta.exposure^2))
    p_ivw  <- 2 * pnorm(-abs(b_ivw / se_ivw))
    method <- if (phi > 1) "IVW (random-effects)" else "IVW (fixed-effects)"
  }

  # MR-Egger (needs ≥ 3 instruments)
  egger_b <- egger_se <- egger_p <- egger_int <- egger_int_p <- NA_real_
  if (n_iv >= 3) {
    w       <- 1 / df_t$se.outcome^2
    # Egger regression: beta_outcome ~ beta_exposure (no intercept forced)
    fit_eg  <- lm(beta.outcome ~ beta.exposure, data = df_t, weights = w)
    coef_eg <- summary(fit_eg)$coefficients
    if (nrow(coef_eg) >= 2) {
      egger_b     <- coef_eg["beta.exposure", "Estimate"]
      egger_se    <- coef_eg["beta.exposure", "Std. Error"]
      egger_p     <- coef_eg["beta.exposure", "Pr(>|t|)"]
      egger_int   <- coef_eg["(Intercept)", "Estimate"]
      egger_int_p <- coef_eg["(Intercept)", "Pr(>|t|)"]
    }
  }

  # Weighted Median (needs ≥ 3)
  wm_b <- wm_se <- wm_p <- NA_real_
  if (n_iv >= 3) {
    wr   <- df_t$beta.outcome / df_t$beta.exposure
    w_wm <- abs(df_t$beta.exposure) / df_t$se.outcome
    w_wm <- w_wm / sum(w_wm)
    # Bootstrap SE
    set.seed(42)
    boot_wm <- replicate(1000, {
      bw    <- rexp(n_iv, rate = 1)
      bw    <- abs(df_t$beta.exposure) / df_t$se.outcome * bw
      bw    <- bw / sum(bw)
      matrixStats::weightedMedian(wr, w = bw)
    })
    wm_b  <- matrixStats::weightedMedian(wr, w = w_wm)
    wm_se <- sd(boot_wm)
    wm_p  <- 2 * pnorm(-abs(wm_b / wm_se))
  }

  # Steiger filtering
  r2_exp <- mean((df_t$beta.exposure^2) /
                   (df_t$beta.exposure^2 + df_t$samplesize.exposure * df_t$se.exposure^2),
                 na.rm = TRUE)
  r2_out <- mean((df_t$beta.outcome^2) /
                   (df_t$beta.outcome^2 + (df_t$n_eff + 2) * df_t$se.outcome^2),
                 na.rm = TRUE)
  steiger_correct <- r2_exp > r2_out

  data.frame(
    taxon            = taxon_name,
    n_iv             = n_iv,
    mean_F           = mean(df_t$beta.exposure^2 / df_t$se.exposure^2),
    IVW_b            = b_ivw,
    IVW_se           = se_ivw,
    IVW_p            = p_ivw,
    IVW_Q            = q_stat,
    IVW_Qp           = q_p,
    IVW_method       = method,
    Egger_b          = egger_b,
    Egger_se         = egger_se,
    Egger_p          = egger_p,
    Egger_intercept  = egger_int,
    Egger_int_p      = egger_int_p,
    WM_b             = wm_b,
    WM_se            = wm_se,
    WM_p             = wm_p,
    r2_exposure      = r2_exp,
    r2_outcome       = r2_out,
    steiger_correct  = steiger_correct,
    stringsAsFactors = FALSE
  )
}

# Run for all taxa (with progress indicator)
taxa_list <- unique(all_harm$exposure)
mr_list   <- vector("list", length(taxa_list))

for (i in seq_along(taxa_list)) {
  if (i %% 20 == 0) cat(sprintf("  %d / %d taxa completed\n", i, length(taxa_list)))
  mr_list[[i]] <- tryCatch(
    run_mr_taxon(all_harm, taxa_list[i]),
    error = function(e) { warning(taxa_list[i], ": ", e$message); NULL }
  )
}

mr_results <- bind_rows(mr_list) %>%
  mutate(
    OR     = exp(IVW_b),
    OR_lo  = exp(IVW_b - 1.96 * IVW_se),
    OR_hi  = exp(IVW_b + 1.96 * IVW_se),
    # Multiple testing corrections
    FDR    = p.adjust(IVW_p, method = "BH"),
    sig_Bonf = IVW_p < CONFIG$mr$bonferroni_p,
    sig_FDR  = FDR < CONFIG$mr$fdr_q,
    # Direction consistency across methods
    dir_IVW_Egger = sign(IVW_b) == sign(Egger_b),
    dir_IVW_WM    = sign(IVW_b) == sign(WM_b),
    # Pleiotropy flag
    pleiotropy = !is.na(Egger_int_p) & Egger_int_p < 0.05
  )

cat(sprintf("\n── MR complete: %d taxa ──\n", nrow(mr_results)))
cat(sprintf("  Bonferroni significant (p < %.2e): %d\n",
            CONFIG$mr$bonferroni_p, sum(mr_results$sig_Bonf)))
cat(sprintf("  FDR significant (q < %.2f): %d\n",
            CONFIG$mr$fdr_q, sum(mr_results$sig_FDR)))

fwrite(mr_results, file.path(CONFIG$paths$out_dir, "mr_all_taxa.csv"))

# ════════════════════════════════════════════════════════════
# 2C. TAXONOMY AUDIT
# Remove DMP higher-order taxonomy strings and unknowns
# ════════════════════════════════════════════════════════════

cat("\n── Taxonomy audit ──\n")

is_higher_order <- function(taxon_name) {
  # DMP encodes higher-order groups as full lineage strings
  # Genus-level entries end in "g__GenusName" or "s__SpeciesName"
  # Family-level or above: no g__ at terminal position
  str_detect(taxon_name, "f__[^.]+$") |    # ends at family
  str_detect(taxon_name, "o__[^.]+$") |    # ends at order
  str_detect(taxon_name, "c__[^.]+$") |    # ends at class
  str_detect(taxon_name, "p__[^.]+$") |    # ends at phylum
  str_detect(taxon_name, "unknown|noname", negate = FALSE) |
  str_detect(taxon_name, "unknownid")
}

mr_sig <- mr_results %>% filter(sig_Bonf)
mr_sig <- mr_sig %>%
  mutate(excluded = is_higher_order(taxon))

n_excluded <- sum(mr_sig$excluded)
mr_clean   <- mr_sig %>% filter(!excluded)

cat(sprintf("  Bonferroni-significant before audit: %d\n", nrow(mr_sig)))
cat(sprintf("  Excluded (higher-order taxonomy): %d\n", n_excluded))
cat(sprintf("  Retained (genus-level): %d\n", nrow(mr_clean)))
cat(sprintf("  Risk-increasing (OR>1): %d\n", sum(mr_clean$OR > 1)))
cat(sprintf("  Protective (OR<1): %d\n", sum(mr_clean$OR < 1)))

fwrite(mr_clean,  file.path(CONFIG$paths$out_dir, "mr_clean_genera.csv"))
fwrite(mr_results %>% filter(sig_Bonf, is_higher_order(taxon)),
       file.path(CONFIG$paths$out_dir, "mr_excluded_taxa.csv"))

cat("✅  MR results saved.\n")
