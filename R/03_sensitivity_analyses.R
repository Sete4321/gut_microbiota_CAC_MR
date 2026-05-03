# ============================================================
# 03_sensitivity_analyses.R
# Steiger filtering, MVMR within Lachnospiraceae,
# leave-one-out, funnel plots
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})
source("R/config.R")

mr_clean <- fread(file.path(CONFIG$paths$out_dir, "mr_clean_genera.csv"))
all_harm <- bind_rows(
  fread(file.path(CONFIG$paths$out_dir, "mbg_instruments.csv")),
  fread(file.path(CONFIG$paths$out_dir, "dmp_instruments.csv"))
) %>%
  inner_join(
    fread(file.path(CONFIG$paths$out_dir, "meta_cac_outcome.csv")) %>%
      rename(SNP = SNP),
    by = "SNP"
  ) %>%
  filter(mr_keep.exposure, mr_keep.outcome)

# ════════════════════════════════════════════════════════════
# 3A. STEIGER FILTERING (already done in 02; visualise here)
# ════════════════════════════════════════════════════════════

p_steiger <- mr_clean %>%
  mutate(
    direction = ifelse(steiger_correct, "Correct\n(microbiome→CAC)", "Incorrect"),
    col       = ifelse(OR > 1, "#C0392B", "#2471A3")
  ) %>%
  ggplot(aes(x = r2_exposure, y = r2_outcome, colour = OR > 1)) +
  geom_point(alpha = 0.8, size = 2.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  scale_colour_manual(
    values = c("TRUE" = "#C0392B", "FALSE" = "#2471A3"),
    labels = c("TRUE" = "Risk-increasing (OR>1)", "FALSE" = "Protective (OR<1)")
  ) +
  labs(
    title    = "Steiger Filtering: R² in Exposure vs Outcome",
    subtitle = sprintf("All %d genera: points above diagonal = correct direction", nrow(mr_clean)),
    x        = "R² in microbiome exposure",
    y        = "R² in CAC outcome",
    colour   = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# ════════════════════════════════════════════════════════════
# 3B. MVMR WITHIN LACHNOSPIRACEAE
# ════════════════════════════════════════════════════════════
# Tests whether Blautia and Anaerostipes effects are independent
# or reflect shared Lachnospiraceae family-level variation

cat("\n── MVMR within Lachnospiraceae ──\n")

lach_taxa <- c("blautia", "anaerostipes")   # adjust taxon names as they appear in data

run_mvmr <- function(harm_data, taxa, out_data) {
  # Extract instruments for both genera
  harm_mv <- harm_data %>%
    filter(tolower(exposure) %in% tolower(taxa)) %>%
    select(SNP, exposure, beta.exposure, se.exposure, beta.outcome, se.outcome, n_eff)

  # Pivot to wide: one row per SNP, one column per taxon's beta
  wide <- harm_mv %>%
    select(SNP, exposure, beta.exposure, se.exposure, beta.outcome, se.outcome) %>%
    pivot_wider(
      id_cols   = c(SNP, beta.outcome, se.outcome),
      names_from  = exposure,
      values_from = c(beta.exposure, se.exposure)
    ) %>%
    drop_na()

  n_snps <- nrow(wide)
  cat(sprintf("  MVMR: %d combined instruments for %s\n", n_snps, paste(taxa, collapse=" + ")))

  if (n_snps < length(taxa) + 1) {
    cat("  Insufficient instruments for MVMR\n")
    return(NULL)
  }

  # Build design matrix X (betas for each exposure)
  X_cols <- grep("^beta.exposure_", colnames(wide), value = TRUE)
  X  <- as.matrix(wide[, ..X_cols])
  y  <- wide$beta.outcome
  w  <- 1 / wide$se.outcome^2

  # WLS: theta = (X'WX)^-1 X'Wy
  XtWX   <- t(X) %*% diag(w) %*% X
  XtWy   <- t(X) %*% (w * y)
  theta   <- solve(XtWX, XtWy)
  # Variance
  resid   <- y - X %*% theta
  sigma2  <- sum(w * resid^2) / (n_snps - length(taxa))
  V_theta <- sigma2 * solve(XtWX)
  se_theta <- sqrt(diag(V_theta))
  z_theta  <- theta / se_theta
  p_theta  <- 2 * pnorm(-abs(z_theta))

  data.frame(
    taxon     = gsub("beta.exposure_", "", X_cols),
    MV_beta   = as.numeric(theta),
    MV_se     = se_theta,
    MV_p      = p_theta,
    MV_OR     = exp(as.numeric(theta)),
    n_instruments = n_snps
  )
}

mvmr_res <- run_mvmr(all_harm, lach_taxa, meta_cac)

if (!is.null(mvmr_res)) {
  cat("\n  MVMR results:\n")
  print(mvmr_res)
  # Compare with univariable results
  uv_res <- mr_clean %>%
    filter(tolower(taxon) %in% tolower(lach_taxa)) %>%
    select(taxon, UV_OR = OR, UV_b = IVW_b, UV_se = IVW_se, UV_p = IVW_p)

  mvmr_comparison <- mvmr_res %>%
    left_join(uv_res, by = "taxon")
  cat("\n  Univariable vs Multivariable comparison:\n")
  print(mvmr_comparison %>% select(taxon, UV_OR, MV_OR, UV_p, MV_p))
  fwrite(mvmr_comparison, file.path(CONFIG$paths$out_dir, "mvmr_lachnospiraceae.csv"))
}

# ════════════════════════════════════════════════════════════
# 3C. LEAVE-ONE-OUT ANALYSIS (top 10 genera by |beta|)
# ════════════════════════════════════════════════════════════

cat("\n── Leave-one-out analysis ──\n")

loo_ivw <- function(harm_df, taxon_name) {
  df <- harm_df %>% filter(exposure == taxon_name, mr_keep.exposure)
  if (nrow(df) < 3) return(NULL)

  loo_list <- lapply(seq_len(nrow(df)), function(i) {
    df_loo <- df[-i, ]
    w <- 1 / df_loo$se.outcome^2
    b <- sum(w * df_loo$beta.exposure * df_loo$beta.outcome) /
         sum(w * df_loo$beta.exposure^2)
    se <- sqrt(1 / sum(w * df_loo$beta.exposure^2))
    data.frame(
      SNP_omitted = df$SNP[i],
      taxon       = taxon_name,
      LOO_b       = b,
      LOO_se      = se,
      LOO_OR      = exp(b),
      LOO_OR_lo   = exp(b - 1.96 * se),
      LOO_OR_hi   = exp(b + 1.96 * se)
    )
  })
  bind_rows(loo_list)
}

# Run LOO for top 10 genera
top10 <- mr_clean %>%
  arrange(desc(abs(IVW_b))) %>%
  slice_head(n = 10) %>%
  pull(taxon)

loo_results <- lapply(top10, function(tx) {
  loo_ivw(all_harm, tx)
}) %>%
  bind_rows()

fwrite(loo_results, file.path(CONFIG$paths$out_dir, "loo_top10.csv"))
cat(sprintf("  LOO complete for %d genera\n", length(top10)))

# ════════════════════════════════════════════════════════════
# 3D. FUNNEL PLOT (all instruments, top 12 genera)
# ════════════════════════════════════════════════════════════

top12_taxa <- mr_clean %>%
  arrange(desc(abs(IVW_b))) %>%
  slice_head(n = 12) %>%
  pull(taxon)

funnel_data <- all_harm %>%
  filter(exposure %in% top12_taxa) %>%
  mutate(
    wald_ratio = beta.outcome / beta.exposure,
    precision  = 1 / (se.outcome / abs(beta.exposure))
  )

p_funnel <- ggplot(funnel_data, aes(x = wald_ratio, y = precision)) +
  geom_point(alpha = 0.6, size = 1.8, colour = "#4393C3") +
  geom_vline(data = mr_clean %>% filter(taxon %in% top12_taxa),
             aes(xintercept = IVW_b), colour = "#C0392B",
             linetype = "dashed", linewidth = 0.8) +
  facet_wrap(~ exposure, scales = "free", ncol = 4) +
  labs(
    title = "Funnel Plots: Top 12 Genera by |IVW Beta|",
    x     = "Wald Ratio (per-instrument causal estimate)",
    y     = "Precision (1 / SE_ratio)"
  ) +
  theme_bw(base_size = 10)

# ════════════════════════════════════════════════════════════
# 3E. SENSITIVITY SUMMARY FIGURE
# ════════════════════════════════════════════════════════════

# Panel 1: IVW vs Weighted Median
p_ivw_wm <- mr_clean %>%
  filter(!is.na(WM_b)) %>%
  ggplot(aes(x = IVW_b, y = WM_b, colour = OR > 1)) +
  geom_point(alpha = 0.75, size = 2.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = 0, colour = "grey70", linetype = "dotted") +
  geom_vline(xintercept = 0, colour = "grey70", linetype = "dotted") +
  scale_colour_manual(
    values = c("TRUE" = "#C0392B", "FALSE" = "#2471A3"),
    labels = c("TRUE" = "OR > 1", "FALSE" = "OR < 1"),
    guide  = "none"
  ) +
  labs(title = "A  IVW vs Weighted Median",
       x = "IVW beta", y = "Weighted Median beta") +
  theme_bw(base_size = 12)

# Panel 2: IVW vs MR-Egger
p_ivw_eg <- mr_clean %>%
  filter(!is.na(Egger_b)) %>%
  ggplot(aes(x = IVW_b, y = Egger_b, colour = pleiotropy)) +
  geom_point(alpha = 0.75, size = 2.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  scale_colour_manual(
    values = c("TRUE" = "#E67E22", "FALSE" = "#7F8C8D"),
    labels = c("TRUE" = "Pleiotropy flagged", "FALSE" = "No pleiotropy")
  ) +
  labs(title = "B  IVW vs MR-Egger",
       x = "IVW beta", y = "MR-Egger beta",
       colour = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# Panel 3: Steiger R² scatter
p_stg <- p_steiger + labs(title = "C  Steiger R² Ratio")

fig_sensitivity <- (p_ivw_wm | p_ivw_eg | p_stg) +
  plot_annotation(
    title   = "Sensitivity Analyses",
    subtitle = sprintf("%d Bonferroni-significant genera | Two-sample MR | Meta-CAC N=49,309", nrow(mr_clean))
  )

ggsave(
  file.path(CONFIG$paths$fig_dir, "fig_sensitivity.png"),
  fig_sensitivity,
  width = 15, height = 5,
  dpi   = CONFIG$fig$dpi
)

cat("✅  Sensitivity analyses complete. Figures saved.\n")
