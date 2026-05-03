# ============================================================
# 04_taxonomic_pattern_analysis.R
# Assign phylum/family taxonomy to 81 significant genera
# Test whether Lachnospiraceae+Ruminococcaceae are systematically
# protective vs other Firmicutes and Non-Firmicutes groups
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(stringr)
})
source("R/config.R")

mr_clean <- fread(file.path(CONFIG$paths$out_dir, "mr_clean_genera.csv"))

# ════════════════════════════════════════════════════════════
# 4A. TAXONOMY LOOKUP TABLE
# Based on NCBI/SILVA 2024 classifications
# ════════════════════════════════════════════════════════════

taxonomy_db <- tribble(
  ~taxon_key,                       ~phylum,          ~family,
  "anaerostipes",                   "Firmicutes",      "Lachnospiraceae",
  "blautia",                        "Firmicutes",      "Lachnospiraceae",
  "butyrivibrio",                   "Firmicutes",      "Lachnospiraceae",
  "coprococcus",                    "Firmicutes",      "Lachnospiraceae",
  "dorea",                          "Firmicutes",      "Lachnospiraceae",
  "eisenbergiella",                 "Firmicutes",      "Lachnospiraceae",
  "fusicatenibacter",               "Firmicutes",      "Lachnospiraceae",
  "holdemania",                     "Firmicutes",      "Lachnospiraceae",
  "howardella",                     "Firmicutes",      "Lachnospiraceae",
  "hungatella",                     "Firmicutes",      "Lachnospiraceae",
  "lachnoclostridium",              "Firmicutes",      "Lachnospiraceae",
  "lachnospiraceaefcsgroup",        "Firmicutes",      "Lachnospiraceae",
  "lachnospiraceaencgroup",         "Firmicutes",      "Lachnospiraceae",
  "lachnospiraceaendgroup",         "Firmicutes",      "Lachnospiraceae",
  "lachnospiraceaenkagroup",        "Firmicutes",      "Lachnospiraceae",
  "lachnospiraceaeucg",             "Firmicutes",      "Lachnospiraceae",
  "marvinbryantia",                 "Firmicutes",      "Lachnospiraceae",
  "roseburia",                      "Firmicutes",      "Lachnospiraceae",
  "eubacteriumeligensgroup",        "Firmicutes",      "Lachnospiraceae",
  "eubacteriumfissicatenagroup",    "Firmicutes",      "Lachnospiraceae",
  "eubacteriumhalliigroup",         "Firmicutes",      "Lachnospiraceae",
  "eubacteriumoxidoreducensgroup",  "Firmicutes",      "Lachnospiraceae",
  "eubacteriumrectalegroup",        "Firmicutes",      "Lachnospiraceae",
  "eubacteriumruminantiumgroup",    "Firmicutes",      "Lachnospiraceae",
  "eubacteriumventriosumgroup",     "Firmicutes",      "Lachnospiraceae",
  "eubacteriumxylanophilumgroup",   "Firmicutes",      "Lachnospiraceae",
  "anaerofilum",                    "Firmicutes",      "Ruminococcaceae",
  "anaerotruncus",                  "Firmicutes",      "Ruminococcaceae",
  "butyricicoccus",                 "Firmicutes",      "Ruminococcaceae",
  "defluviitaleaceaeucg",           "Firmicutes",      "Ruminococcaceae",
  "faecalibacterium",               "Firmicutes",      "Ruminococcaceae",
  "flavonifractor",                 "Firmicutes",      "Ruminococcaceae",
  "intestinimonas",                 "Firmicutes",      "Ruminococcaceae",
  "pseudoflavonifractor",           "Firmicutes",      "Ruminococcaceae",
  "ruminiclostridium",              "Firmicutes",      "Ruminococcaceae",
  "ruminococcaceaenkagroup",        "Firmicutes",      "Ruminococcaceae",
  "ruminococcaceaeucg",             "Firmicutes",      "Ruminococcaceae",
  "ruminococcus",                   "Firmicutes",      "Ruminococcaceae",
  "ruminococcusgauvreauiigroup",    "Firmicutes",      "Ruminococcaceae",
  "subdoligranulum",                "Firmicutes",      "Ruminococcaceae",
  "clostridiumsensustricto",        "Firmicutes",      "Clostridiaceae",
  "clostridiuminnocuumgroup",       "Firmicutes",      "Peptostreptococcaceae",
  "familyxiiiadgroup",              "Firmicutes",      "Peptostreptococcaceae",
  "familyxiiiucg",                  "Firmicutes",      "Peptostreptococcaceae",
  "terrisporobacter",               "Firmicutes",      "Peptostreptococcaceae",
  "christensenellaceaergroup",      "Firmicutes",      "Christensenellaceae",
  "catenibacterium",                "Firmicutes",      "Erysipelotrichaceae",
  "erysipelotrichaceaeucg",         "Firmicutes",      "Erysipelotrichaceae",
  "turicibacter",                   "Firmicutes",      "Erysipelotrichaceae",
  "oscillibacter",                  "Firmicutes",      "Oscillospiraceae",
  "dialister",                      "Firmicutes",      "Acidaminococcaceae",
  "phascolarctobacterium",          "Firmicutes",      "Acidaminococcaceae",
  "veillonella",                    "Firmicutes",      "Veillonellaceae",
  "eubacterium",                    "Firmicutes",      "Eubacteriaceae",
  "eubacteriumbrachygroup",         "Firmicutes",      "Eubacteriaceae",
  "lactobacillus",                  "Firmicutes",      "Lactobacillaceae",
  "lactococcus",                    "Firmicutes",      "Streptococcaceae",
  "streptococcus",                  "Firmicutes",      "Streptococcaceae",
  "alistipes",                      "Bacteroidetes",   "Rikenellaceae",
  "alloprevotella",                 "Bacteroidetes",   "Prevotellaceae",
  "coprobacter",                    "Bacteroidetes",   "Coprobacteraceae",
  "odoribacter",                    "Bacteroidetes",   "Marinifilaceae",
  "parabacteroides",                "Bacteroidetes",   "Tannerellaceae",
  "paraprevotella",                 "Bacteroidetes",   "Prevotellaceae",
  "prevotella",                     "Bacteroidetes",   "Prevotellaceae",
  "bilophila",                      "Bacteroidetes",   "Desulfovibrionaceae",
  "actinomyces",                    "Actinobacteria",  "Actinomycetaceae",
  "adlercreutzia",                  "Actinobacteria",  "Eggerthellaceae",
  "bifidobacterium",                "Actinobacteria",  "Bifidobacteriaceae",
  "eggerthella",                    "Actinobacteria",  "Eggerthellaceae",
  "olsenella",                      "Actinobacteria",  "Eggerthellaceae",
  "rothia",                         "Actinobacteria",  "Micrococcaceae",
  "slackia",                        "Actinobacteria",  "Eggerthellaceae",
  "escherichiashigella",            "Proteobacteria",  "Enterobacteriaceae",
  "escherichia",                    "Proteobacteria",  "Enterobacteriaceae",
  "haemophilus",                    "Proteobacteria",  "Pasteurellaceae",
  "sutterella",                     "Proteobacteria",  "Sutterellaceae",
  "parasutterella",                 "Proteobacteria",  "Sutterellaceae",
  "victivallis",                    "Lentisphaerae",   "Victivallaceae",
  "methanobrevibacter",             "Euryarchaeota",   "Methanobacteriaceae"
)

# Normalise taxon names for matching
normalise_taxon <- function(x) {
  x %>%
    str_to_lower() %>%
    str_remove("genus\\.") %>%
    str_remove("\\.id\\.\\d+$") %>%
    str_trim()
}

mr_tax <- mr_clean %>%
  mutate(taxon_key = normalise_taxon(taxon)) %>%
  left_join(taxonomy_db, by = "taxon_key")

# Check for unmatched
unmatched <- mr_tax %>% filter(is.na(phylum)) %>% pull(taxon)
if (length(unmatched) > 0) {
  warning(sprintf("%d taxa unmatched: %s",
                  length(unmatched), paste(unmatched, collapse = ", ")))
}

cat(sprintf("Taxonomy assigned to %d / %d genera\n",
            sum(!is.na(mr_tax$phylum)), nrow(mr_tax)))

fwrite(mr_tax, file.path(CONFIG$paths$out_dir, "mr_clean_with_taxonomy.csv"))

# ════════════════════════════════════════════════════════════
# 4B. FAMILY-LEVEL SUMMARY: mean OR + one-sample t-test
# ════════════════════════════════════════════════════════════

fam_summary <- mr_tax %>%
  filter(!is.na(family)) %>%
  group_by(phylum, family) %>%
  summarise(
    n        = n(),
    n_prot   = sum(IVW_b < 0),
    n_risk   = sum(IVW_b > 0),
    mean_b   = mean(IVW_b),
    sd_b     = sd(IVW_b),
    se_mean  = ifelse(n > 1, sd_b / sqrt(n), abs(IVW_b) * 0.30),
    mean_OR  = exp(mean_b),
    ci_lo    = exp(mean_b - qt(0.975, df = pmax(n - 1, 1)) * se_mean),
    ci_hi    = exp(mean_b + qt(0.975, df = pmax(n - 1, 1)) * se_mean),
    t_stat   = ifelse(n > 1, mean_b / se_mean, NA_real_),
    p_ttest  = ifelse(n > 1,
                 2 * pt(-abs(t_stat), df = n - 1),
                 NA_real_),
    p_binom  = binom.test(n_prot, n, p = 0.5)$p.value,
    .groups  = "drop"
  ) %>%
  arrange(mean_b)

cat("\n── Family-level summary ──\n")
print(fam_summary %>%
  select(family, n, n_prot, n_risk, mean_OR, p_ttest) %>%
  arrange(mean_OR))

fwrite(fam_summary, file.path(CONFIG$paths$out_dir, "family_summary.csv"))

# ════════════════════════════════════════════════════════════
# 4C. GROUP-LEVEL COMPARISON: LachRum vs Others
# ════════════════════════════════════════════════════════════

mr_tax <- mr_tax %>%
  mutate(group = case_when(
    family %in% CONFIG$taxonomy$lach_rum_families ~ "Lachnospiraceae\n+ Ruminococcaceae",
    phylum == "Firmicutes" ~ "Other\nFirmicutes",
    TRUE ~ "Non-Firmicutes"
  ))

grp_stats <- mr_tax %>%
  filter(!is.na(group)) %>%
  group_by(group) %>%
  summarise(
    n        = n(),
    n_prot   = sum(IVW_b < 0),
    pct_prot = n_prot / n * 100,
    mean_OR  = exp(mean(IVW_b)),
    sd_b     = sd(IVW_b),
    se_mean  = sd_b / sqrt(n),
    ci_lo    = exp(mean(IVW_b) - qt(0.975, df = n - 1) * se_mean),
    ci_hi    = exp(mean(IVW_b) + qt(0.975, df = n - 1) * se_mean),
    t_stat   = mean(IVW_b) / se_mean,
    p_ttest  = 2 * pt(-abs(t_stat), df = n - 1),
    .groups  = "drop"
  )

cat("\n── Group-level statistics ──\n")
print(grp_stats %>% select(group, n, pct_prot, mean_OR, p_ttest))

# Mann-Whitney U tests
b_lr <- mr_tax %>% filter(family %in% CONFIG$taxonomy$lach_rum_families) %>% pull(IVW_b)
b_of <- mr_tax %>% filter(phylum == "Firmicutes",
                           !family %in% CONFIG$taxonomy$lach_rum_families) %>% pull(IVW_b)
b_nf <- mr_tax %>% filter(phylum != "Firmicutes") %>% pull(IVW_b)

mwu_lr_of <- wilcox.test(b_lr, b_of, alternative = "less")
mwu_lr_nf <- wilcox.test(b_lr, b_nf, alternative = "less")

# Chi-square: protective proportions
n_lr_p <- sum(b_lr < 0); n_lr <- length(b_lr)
n_of_p <- sum(b_of < 0); n_of <- length(b_of)
ct <- matrix(c(n_lr_p, n_lr - n_lr_p, n_of_p, n_of - n_of_p),
             nrow = 2, byrow = TRUE)
chi_res <- chisq.test(ct, correct = FALSE)

cat(sprintf("\n── Between-group tests ──\n"))
cat(sprintf("  MWU LachRum vs OtherFirm: U=%.0f, p=%.4f\n",
            mwu_lr_of$statistic, mwu_lr_of$p.value))
cat(sprintf("  MWU LachRum vs NonFirm:   U=%.0f, p=%.4f\n",
            mwu_lr_nf$statistic, mwu_lr_nf$p.value))
cat(sprintf("  Chi-sq protective prop:   X²=%.3f, p=%.4f\n",
            chi_res$statistic, chi_res$p.value))
cat(sprintf("  LachRum: %d/%d protective (%.0f%%)\n", n_lr_p, n_lr, n_lr_p/n_lr*100))
cat(sprintf("  OtherFirm: %d/%d protective (%.0f%%)\n", n_of_p, n_of, n_of_p/n_of*100))

# Save stats
stats_out <- list(
  group_stats    = grp_stats,
  mwu_lr_of      = broom::tidy(mwu_lr_of),
  mwu_lr_nf      = broom::tidy(mwu_lr_nf),
  chisq          = broom::tidy(chi_res)
)
saveRDS(stats_out, file.path(CONFIG$paths$out_dir, "taxonomy_stats.rds"))

# ════════════════════════════════════════════════════════════
# 4D. FIGURES
# ════════════════════════════════════════════════════════════

PHYLUM_COLS <- c(
  "Firmicutes"    = "#D6604D",
  "Bacteroidetes" = "#4393C3",
  "Actinobacteria"= "#74C476",
  "Proteobacteria"= "#9E7AC7",
  "Euryarchaeota" = "#F4A460",
  "Lentisphaerae" = "#AAAAAA"
)
KEY_FAM_COLS <- c("Lachnospiraceae" = "#8B1A1A", "Ruminococcaceae" = "#C0392B")

# ── Panel A: Family-level forest plot ────────────────────────
fam_plot <- fam_summary %>%
  filter(n >= 2) %>%
  arrange(mean_b) %>%
  mutate(
    y_pos  = row_number(),
    label  = sprintf("%s  (n=%d)", family, n),
    col    = case_when(
      family %in% names(KEY_FAM_COLS) ~ KEY_FAM_COLS[family],
      TRUE ~ PHYLUM_COLS[phylum]
    ),
    bold   = family %in% names(KEY_FAM_COLS),
    sig    = !is.na(p_ttest) & p_ttest < 0.05
  )

p_forest <- ggplot(fam_plot, aes(x = mean_OR, y = y_pos)) +
  # Alternating row shading
  geom_rect(aes(xmin = -Inf, xmax = Inf,
                ymin = y_pos - 0.47, ymax = y_pos + 0.47,
                fill = y_pos %% 2 == 0),
            alpha = 0.15, colour = NA) +
  scale_fill_manual(values = c("TRUE" = "grey90", "FALSE" = "white"), guide = "none") +
  # Null line
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.8) +
  # CI bars
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi, colour = family),
                 height = 0, linewidth = 1.6) +
  # Mean OR dots (size ∝ n)
  geom_point(aes(size = n, colour = family)) +
  # OR value text
  geom_text(aes(x = max(fam_plot$ci_hi) * 1.12,
                label = sprintf("%.3f", mean_OR),
                fontface = ifelse(bold, "bold", "plain")),
            hjust = 0, size = 4, colour = "#111111") +
  scale_colour_manual(values = c(KEY_FAM_COLS, PHYLUM_COLS), guide = "none") +
  scale_size_continuous(range = c(2, 8), guide = "none") +
  scale_y_continuous(
    breaks = fam_plot$y_pos,
    labels = fam_plot$label,
    expand = c(0.02, 0)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.20))) +
  labs(
    title    = "A  Family-Level Causal Effect on CAC",
    subtitle = "Mean OR across genera within each family (n ≥ 2)",
    x        = "Mean Family OR (95% CI)",
    y        = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.y   = element_text(face = "bold", size = 13),
    panel.grid.major.y = element_blank()
  )

# ── Panel B: Stacked bars per family ─────────────────────────
fam_bar <- fam_summary %>%
  filter(n >= 3) %>%
  arrange(desc(n_prot / n)) %>%
  mutate(family = factor(family, levels = family))

p_bar <- fam_bar %>%
  select(family, n, n_prot, n_risk) %>%
  tidyr::pivot_longer(cols = c(n_prot, n_risk),
                      names_to  = "direction",
                      values_to = "count") %>%
  mutate(
    pct       = count / n * 100,
    direction = ifelse(direction == "n_prot",
                       "Protective (OR<1)", "Risk-increasing (OR>1)")
  ) %>%
  ggplot(aes(x = family, y = pct, fill = direction)) +
  geom_col(width = 0.65, alpha = 0.87) +
  geom_hline(yintercept = 50, linetype = "dashed",
             colour = "#333333", linewidth = 0.8) +
  geom_text(data = fam_bar,
            aes(x = family, y = 50, label = paste0("n=", n)),
            inherit.aes = FALSE,
            fontface = "bold", colour = "white", size = 3.8) +
  scale_fill_manual(values = c("Protective (OR<1)"     = "#2471A3",
                                "Risk-increasing (OR>1)" = "#C0392B")) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  labs(
    title = "B  Proportion Protective vs Risk-Increasing",
    x     = NULL, y = "% of genera", fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

# ── Panel C: Violin + strip comparison ───────────────────────
PROT_COL <- "#1F618D"
RISK_COL <- "#922B21"
MEAN_COL <- "#1A1A1A"

VCOLS <- c(
  "Lachnospiraceae\n+ Ruminococcaceae" = "#F0C8C8",
  "Other\nFirmicutes"                  = "#FAD7C0",
  "Non-Firmicutes"                     = "#C5D8F0"
)
ECOLS <- c(
  "Lachnospiraceae\n+ Ruminococcaceae" = "#8B1A1A",
  "Other\nFirmicutes"                  = "#D6604D",
  "Non-Firmicutes"                     = "#2471A3"
)

grp_ci <- grp_stats %>%
  rename(group_label = group) %>%
  mutate(mean_OR = exp(mean(mr_tax$IVW_b[mr_tax$group == group_label])))

p_violin <- mr_tax %>%
  filter(!is.na(group)) %>%
  mutate(dot_col = ifelse(IVW_b < 0, PROT_COL, RISK_COL)) %>%
  ggplot(aes(x = group, y = OR)) +
  geom_violin(aes(fill = group, colour = group),
              alpha = 0.70, linewidth = 1.5) +
  geom_jitter(aes(colour = dot_col),
              width = 0.12, size = 2.0, alpha = 0.85) +
  # Mean ± CI from group stats
  geom_crossbar(data = grp_stats,
                aes(x = group, y = mean_OR,
                    ymin = ci_lo, ymax = ci_hi),
                width = 0.35, linewidth = 1.2, colour = MEAN_COL) +
  # Mean OR labels
  geom_text(data = grp_stats,
            aes(x = group, y = ci_hi + 0.004,
                label = sprintf("Mean OR\n%.3f", mean_OR)),
            size = 4, fontface = "bold",
            colour = ECOLS[grp_stats$group]) +
  geom_hline(yintercept = 1, linetype = "dashed",
             colour = "#444444", linewidth = 0.9) +
  scale_fill_manual(values   = VCOLS, guide = "none") +
  scale_colour_manual(
    values = c(VCOLS, setNames(ECOLS, names(ECOLS)),
               PROT_COL, RISK_COL, MEAN_COL),
    guide  = "none"
  ) +
  labs(
    title = "C  OR Distribution by Taxonomic Group",
    x     = NULL, y = "Odds Ratio (IVW)"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(size = 12, face = "bold"))

# Combine panels B and C
fig_taxonomy <- (p_forest | (p_bar / p_violin)) +
  plot_layout(widths = c(1.1, 1)) +
  plot_annotation(
    title    = "Taxonomic Pattern Analysis: Family- and Phylum-Level Causal Effect Structure",
    subtitle = sprintf(
      "%d Bonferroni-significant genera  |  Two-sample MR (IVW)  |  Meta-CAC N=49,309",
      nrow(mr_clean)
    ),
    theme = theme(plot.title = element_text(face = "bold", size = 15))
  )

ggsave(
  file.path(CONFIG$paths$fig_dir, "fig_taxonomy.png"),
  fig_taxonomy,
  width  = 22, height = 16,
  dpi    = CONFIG$fig$dpi
)

cat("\n✅  Taxonomic pattern analysis complete. Figure saved.\n")
