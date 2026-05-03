# ============================================================
# 05_mpgs_prediction.R
# Microbiome Polygenic Score (mPGS) construction
# R², AUC, decile ORs, Population-Attributable Fractions
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
})
source("R/config.R")

mr_clean  <- fread(file.path(CONFIG$paths$out_dir, "mr_clean_genera.csv"))
mbg_inst  <- fread(file.path(CONFIG$paths$out_dir, "mbg_instruments.csv"))
dmp_inst  <- fread(file.path(CONFIG$paths$out_dir, "dmp_instruments.csv"))
meta_cac  <- fread(file.path(CONFIG$paths$out_dir, "meta_cac_outcome.csv"))

# ════════════════════════════════════════════════════════════
# 5A. mPGS WEIGHT CONSTRUCTION
# w_i = beta_exposure_i × beta_MR_i  (both on log scale)
# ════════════════════════════════════════════════════════════

# Merge instruments with MR causal estimates
inst_all <- bind_rows(mbg_inst, dmp_inst)

weights <- inst_all %>%
  filter(SNP %in% meta_cac$SNP) %>%
  inner_join(
    mr_clean %>% select(taxon, IVW_b, sig_Bonf),
    by = c("exposure" = "taxon")
  ) %>%
  filter(sig_Bonf) %>%
  # mPGS weight: exposure beta × MR causal estimate
  mutate(
    w_mpgs    = beta.exposure * IVW_b,
    direction = ifelse(IVW_b > 0, "risk", "protective"),
    EAF       = ifelse(!is.na(eaf.exposure), eaf.exposure, 0.30)
  ) %>%
  # Keep one instrument per SNP-taxon pair (remove duplicates)
  distinct(SNP, exposure, .keep_all = TRUE)

cat(sprintf("mPGS instruments: %d SNPs across %d genera\n",
            nrow(weights), length(unique(weights$exposure))))
cat(sprintf("  Risk-increasing: %d SNPs (%d genera)\n",
            sum(weights$direction == "risk"),
            length(unique(weights$exposure[weights$direction == "risk"]))))
cat(sprintf("  Protective: %d SNPs (%d genera)\n",
            sum(weights$direction == "protective"),
            length(unique(weights$exposure[weights$direction == "protective"]))))

# Save PLINK-compatible score file
score_file <- weights %>%
  select(SNP, effect_allele.exposure, w_mpgs) %>%
  rename(
    variant_id    = SNP,
    allele        = effect_allele.exposure,
    effect_weight = w_mpgs
  )
fwrite(score_file,
       file.path(CONFIG$paths$out_dir, "mpgs_weights_plink.txt"),
       sep = "\t")

cat("  PLINK score file saved (use with: plink --score mpgs_weights_plink.txt)\n")

# ════════════════════════════════════════════════════════════
# 5B. EXPECTED R² COMPUTATION
# R² ≈ Σ 2 × EAF × (1-EAF) × w²
# (assumes SNP independence — upper bound)
# ════════════════════════════════════════════════════════════

r2_per_snp <- weights %>%
  mutate(
    var_snp = 2 * EAF * (1 - EAF),
    r2_contrib = var_snp * w_mpgs^2
  )

r2_total      <- sum(r2_per_snp$r2_contrib)
r2_risk       <- sum(r2_per_snp$r2_contrib[r2_per_snp$direction == "risk"])
r2_protective <- sum(r2_per_snp$r2_contrib[r2_per_snp$direction == "protective"])

cat(sprintf("\nExpected R² in CAC:\n"))
cat(sprintf("  Risk mPGS:        R² = %.4f%%\n", r2_risk * 100))
cat(sprintf("  Protective mPGS:  R² = %.4f%%\n", r2_protective * 100))
cat(sprintf("  Net mPGS (total): R² = %.4f%%\n", r2_total * 100))

# AUC approximation under liability threshold model
auc_net <- pnorm(sqrt(r2_total / (1 + r2_total)))
cat(sprintf("  Approximate AUC:  %.4f\n", auc_net))

# ════════════════════════════════════════════════════════════
# 5C. SIMULATE mPGS DISTRIBUTIONS (N = 50,000)
# Genotype counts drawn from Binomial(2, EAF) per SNP
# ════════════════════════════════════════════════════════════

set.seed(42)
N_sim <- 50000

simulate_mpgs <- function(w_df, n = N_sim) {
  score <- numeric(n)
  for (i in seq_len(nrow(w_df))) {
    g <- rbinom(n, 2, w_df$EAF[i])   # genotype dosage
    score <- score + g * w_df$w_mpgs[i]
  }
  score
}

cat("\nSimulating mPGS distributions (N=50,000)...\n")
scores_risk <- simulate_mpgs(weights %>% filter(direction == "risk"))
scores_prot <- simulate_mpgs(weights %>% filter(direction == "protective"))
scores_net  <- scores_risk + scores_prot

# Standardise
z_risk <- scale(scores_risk)[,1]
z_prot <- scale(scores_prot)[,1]
z_net  <- scale(scores_net)[,1]

sim_df <- data.frame(
  risk        = z_risk,
  protective  = z_prot,
  net         = z_net
) %>%
  pivot_longer(everything(), names_to = "mPGS", values_to = "z_score") %>%
  mutate(mPGS = factor(mPGS,
    levels = c("risk", "protective", "net"),
    labels = c("Risk mPGS (41 genera)",
               "Protective mPGS (40 genera)",
               "Net mPGS (81 genera)")))

# ════════════════════════════════════════════════════════════
# 5D. DECILE OR ESTIMATION
# ════════════════════════════════════════════════════════════

p_base   <- 0.30                          # assumed CAC prevalence
or_base  <- p_base / (1 - p_base)

beta_pgs <- sqrt(r2_total / var(scores_net))   # log-OR per unit mPGS
logit_base <- log(p_base / (1 - p_base))

decile_breaks <- quantile(scores_net, probs = seq(0, 1, 0.1))
decile_mid    <- (decile_breaks[-11] + decile_breaks[-1]) / 2

p_decile  <- 1 / (1 + exp(-(logit_base + beta_pgs * decile_mid)))
or_decile <- (p_decile / (1 - p_decile)) / or_base

decile_df <- data.frame(
  decile    = 1:10,
  label     = paste0("D", 1:10),
  mid_score = decile_mid,
  p_cac     = p_decile,
  OR        = or_decile
)

cat(sprintf("\nD10 vs D1 OR ≈ %.2f\n", or_decile[10] / or_decile[1]))

# ════════════════════════════════════════════════════════════
# 5E. POPULATION-ATTRIBUTABLE FRACTIONS
# ════════════════════════════════════════════════════════════

paf_df <- mr_clean %>%
  mutate(
    p_exp  = 0.30,   # assumed exposure prevalence
    # Risk-increasing: PAF = p(OR-1) / [p(OR-1)+1]
    # Protective: preventable fraction = p(1-OR) / [p(1-OR)+1]
    PAF_pct = case_when(
      OR > 1 ~  p_exp * (OR - 1) / (p_exp * (OR - 1) + 1) * 100,
      OR < 1 ~  p_exp * (1 - OR) / (p_exp * (1 - OR) + 1) * 100,
      TRUE   ~ 0
    ),
    direction = ifelse(OR > 1, "Attributable", "Preventable")
  ) %>%
  arrange(desc(PAF_pct))

top_risk <- paf_df %>% filter(OR > 1) %>% slice_head(n = 12)
top_prot <- paf_df %>% filter(OR < 1) %>% slice_head(n = 12)

cat(sprintf("\nTop PAF genera:\n"))
cat(sprintf("  Risk-attributable: %s (%.1f%%)\n",
            top_risk$taxon[1], top_risk$PAF_pct[1]))
cat(sprintf("  Preventable:       %s (%.1f%%)\n",
            top_prot$taxon[1], top_prot$PAF_pct[1]))
cat(sprintf("  Aggregate attributable PAF: %.1f%%\n", sum(paf_df$PAF_pct[paf_df$OR>1])))
cat(sprintf("  Aggregate preventable PAF:  %.1f%%\n", sum(paf_df$PAF_pct[paf_df$OR<1])))

fwrite(paf_df, file.path(CONFIG$paths$out_dir, "paf_results.csv"))

# ════════════════════════════════════════════════════════════
# 5F. FIGURES
# ════════════════════════════════════════════════════════════

# Panel A: PAF diverging bar
paf_plot_df <- bind_rows(
  top_risk %>% mutate(val =  PAF_pct),
  top_prot %>% mutate(val = -PAF_pct)
) %>%
  arrange(val) %>%
  mutate(
    taxon = factor(taxon, levels = taxon),
    col   = ifelse(val > 0, "#C0392B", "#2471A3")
  )

n_prot_all <- sum(mr_clean$OR < 1)
n_risk_all <- sum(mr_clean$OR > 1)

p_paf <- ggplot(paf_plot_df, aes(x = val, y = taxon, fill = col)) +
  geom_col(alpha = 0.85) +
  geom_vline(xintercept = 0, colour = "#333", linewidth = 0.8) +
  # Preventable label (top-left)
  annotate("text",
           x = min(paf_plot_df$val) * 0.55,
           y = nrow(paf_plot_df) - 0.8,
           label = sprintf("\u2190 Preventable\n(n=%d genera)", n_prot_all),
           colour = "#2471A3", fontface = "bold", size = 3.5, hjust = 0.5) +
  # Attributable label (bottom-right, above x-axis)
  annotate("text",
           x = max(paf_plot_df$val) * 0.58,
           y = 0.5,
           label = sprintf("Attributable \u2192\n(n=%d genera)", n_risk_all),
           colour = "#C0392B", fontface = "bold", size = 3.5, hjust = 0.5,
           vjust  = 0) +
  scale_fill_identity() +
  labs(title = "A  Population-Attributable Fraction (%)",
       subtitle = "Top 12 risk-increasing and protective genera",
       x = "PAF (%)", y = NULL) +
  theme_bw(base_size = 10) +
  theme(axis.text.y = element_text(face = "italic", size = 8))

# Panel B: mPGS distributions
p_dist <- ggplot(sim_df, aes(x = z_score, fill = mPGS, colour = mPGS)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 60, alpha = 0.55, position = "identity") +
  stat_function(fun = dnorm, colour = "grey30",
                linetype = "dashed", linewidth = 0.8) +
  scale_fill_manual(values  = c("#C0392B","#2471A3","#27AE60")) +
  scale_colour_manual(values = c("#C0392B","#2471A3","#27AE60")) +
  labs(title    = "B  Simulated mPGS Distributions",
       subtitle = "N = 50,000; standardised z-scores",
       x = "Standardised mPGS (z-score)", y = "Density", fill = NULL, colour = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "inside", legend.position.inside = c(0.02, 0.98),
        legend.justification = c(0, 1))

# Panel C: Cumulative R²
cum_r2_df <- r2_per_snp %>%
  arrange(desc(abs(w_mpgs))) %>%
  group_by(direction) %>%
  mutate(
    rank   = row_number(),
    cum_r2 = cumsum(r2_contrib) * 100
  ) %>%
  ungroup()

p_cumr2 <- ggplot(cum_r2_df, aes(x = rank, y = cum_r2, colour = direction)) +
  geom_line(linewidth = 1.8) +
  geom_hline(yintercept = r2_total * 100, linetype = "dashed",
             colour = "#27AE60", linewidth = 1.2) +
  annotate("text", x = max(cum_r2_df$rank) * 0.5, y = r2_total * 100 * 1.05,
           label = sprintf("Net R² = %.3f%%", r2_total * 100),
           colour = "#27AE60", size = 3.5, fontface = "bold") +
  scale_colour_manual(
    values = c("risk" = "#C0392B", "protective" = "#2471A3"),
    labels = c("risk" = sprintf("Risk (R²=%.3f%%)", r2_risk*100),
               "protective" = sprintf("Protective (R²=%.3f%%)", r2_protective*100))
  ) +
  labs(title    = "C  Cumulative Variance Explained (R²)",
       subtitle = "by mPGS as genera added (ranked by |weight|)",
       x = "Genera added", y = "Cumulative R² in CAC (%)",
       colour = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "inside", legend.position.inside = c(0.05, 0.95),
        legend.justification = c(0, 1))

# Panel D: Decile OR
p_decile_plot <- ggplot(decile_df, aes(x = decile, y = OR, fill = OR)) +
  geom_col(alpha = 0.87) +
  geom_errorbar(
    aes(ymin = exp(log(OR) - 1.96 * sqrt(r2_total) * 0.15),
        ymax = exp(log(OR) + 1.96 * sqrt(r2_total) * 0.15)),
    width = 0.3, colour = "#333"
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "#333") +
  annotate("text",
           x = 8.5, y = max(or_decile) * 0.88,
           label = sprintf("D10 vs D1 OR \u2248 %.2f",
                           or_decile[10] / or_decile[1]),
           colour = "#C0392B", fontface = "bold", size = 3.5) +
  scale_fill_gradient2(low = "#2471A3", mid = "white", high = "#C0392B",
                       midpoint = 1, guide = "none") +
  scale_x_continuous(breaks = 1:10, labels = paste0("D", 1:10)) +
  labs(title    = "D  CAC Risk by mPGS Decile",
       subtitle = "Baseline prevalence = 30%",
       x = "Net mPGS Decile", y = "OR vs population mean") +
  theme_bw(base_size = 10)

# Combine
fig_mpgs <- (p_paf | p_dist) / (p_cumr2 | p_decile_plot) +
  plot_annotation(
    title    = "Microbiome Polygenic Score (mPGS) Prediction Analysis",
    subtitle = sprintf("%d instrument SNPs × %d genera  |  Risk mPGS + Protective mPGS",
                       nrow(weights), nrow(mr_clean)),
    theme    = theme(plot.title = element_text(face = "bold", size = 14))
  )

ggsave(
  file.path(CONFIG$paths$fig_dir, "fig_mpgs_prediction.png"),
  fig_mpgs,
  width = 14, height = 10,
  dpi   = CONFIG$fig$dpi
)

# Save numeric results
saveRDS(
  list(r2_risk = r2_risk, r2_prot = r2_protective, r2_total = r2_total,
       auc = auc_net, decile_or = decile_df),
  file.path(CONFIG$paths$out_dir, "mpgs_results.rds")
)

cat("\n✅  mPGS prediction analysis complete. Figure saved.\n")
