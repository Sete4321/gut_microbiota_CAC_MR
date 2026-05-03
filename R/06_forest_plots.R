# ============================================================
# 06_forest_plots.R
# Publication-quality forest plots:
# Figure 1: 41 risk-increasing genera
# Figure 2: 40 protective genera
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(stringr)
  library(forcats)
})
source("R/config.R")

mr_tax <- fread(file.path(CONFIG$paths$out_dir, "mr_clean_with_taxonomy.csv"))

# ════════════════════════════════════════════════════════════
# DISPLAY NAMES
# Maps raw taxon identifiers to readable genus names
# ════════════════════════════════════════════════════════════

NAME_MAP <- c(
  "eubacteriumhalliigroup"       = "E. hallii group",
  "ruminococcusgauvreauiigroup"  = "R. gauvreauii group",
  "lachnospiraceaendgroup"       = "[Lachno.] ND group",
  "clostridiumsensustricto"      = "C. sensu stricto 1",
  "escherichiashigella"          = "Escherichia/Shigella",
  "ruminococcaceaeucg"           = "[Rumin.] UCG",
  "defluviitaleaceaeucg"         = "Defluviitaleaceae UCG",
  "familyxiiiucg"               = "[Family XIII] UCG",
  "lachnospiraceaefcsgroup"      = "[Lachno.] FCS group",
  "lachnospiraceaenkagroup"      = "[Lachno.] NKA group",
  "christensenellaceaergroup"    = "Christensenellaceae R",
  "clostridiuminnocuumgroup"     = "C. innocuum group",
  "eubacteriumrectalegroup"      = "E. rectale group",
  "eubacteriumbrachygroup"       = "E. brachy group",
  "familyxiiiadgroup"            = "[Family XIII] AD group",
  "ruminococcaceaenkagroup"      = "[Rumin.] NKA group",
  "lachnospiraceaencgroup"       = "[Lachno.] NC group",
  "lachnospiraceaeucg"           = "[Lachno.] UCG",
  "eubacteriumeligensgroup"      = "E. eligens group",
  "eubacteriumfissicatenagroup"  = "E. fissicatena group",
  "eubacteriumoxidoreducensgroup"= "E. oxidoreducens group",
  "eubacteriumruminantiumgroup"  = "E. ruminantium group",
  "eubacteriumventriosumgroup"   = "E. ventriosum group",
  "eubacteriumxylanophilumgroup" = "E. xylanophilum group"
)

PHYLUM_COLS <- c(
  "Firmicutes"     = "#D6604D",
  "Bacteroidetes"  = "#4393C3",
  "Actinobacteria" = "#74C476",
  "Proteobacteria" = "#9E7AC7",
  "Euryarchaeota"  = "#F4A460",
  "Lentisphaerae"  = "#AAAAAA"
)

format_name <- function(taxon) {
  tl <- str_to_lower(str_remove(taxon, "genus\\."))
  tl <- str_remove(tl, "\\.id\\.\\d+$")
  ifelse(tl %in% names(NAME_MAP), NAME_MAP[tl],
         str_to_sentence(tl))
}

# ════════════════════════════════════════════════════════════
# FOREST PLOT FUNCTION
# ════════════════════════════════════════════════════════════

make_forest_plot <- function(df, direction = "risk", title_label = "A") {

  col_scale <- if (direction == "risk") {
    scale_colour_gradient(low = "#F5B7B1", high = "#922B21",
                          name = "|beta|", guide = "none")
  } else {
    scale_colour_gradient(low = "#AED6F1", high = "#1A5276",
                          name = "|beta|", guide = "none")
  }

  # Shape by data source
  df <- df %>%
    mutate(
      display_name = format_name(taxon),
      display_name = factor(display_name, levels = rev(display_name)),
      source_shape = case_when(
        str_detect(tolower(source), "mbg|mibiogen") ~ 16,   # circle
        str_detect(tolower(source), "dmp")          ~ 18,   # diamond
        TRUE                                        ~ 15    # square (pooled)
      ),
      abs_beta = abs(IVW_b),
      label_text = sprintf("%.3f [%.3f, %.3f]", OR, OR_lo, OR_hi)
    ) %>%
    arrange(OR)

  ggplot(df, aes(x = OR, y = display_name, colour = abs_beta)) +

    # Alternating row bands
    geom_rect(aes(xmin = -Inf, xmax = Inf,
                  ymin = as.numeric(display_name) - 0.47,
                  ymax = as.numeric(display_name) + 0.47,
                  fill = as.numeric(display_name) %% 2 == 0),
              inherit.aes = FALSE, alpha = 0.12, colour = NA) +
    scale_fill_manual(values = c("TRUE" = "grey88", "FALSE" = "white"),
                      guide  = "none") +

    # Null line
    geom_vline(xintercept = 1, linetype = "dashed",
               colour = "grey50", linewidth = 0.7) +

    # CI bars
    geom_errorbarh(aes(xmin = OR_lo, xmax = OR_hi),
                   height = 0, linewidth = 1.3, alpha = 0.85) +

    # Mean OR points (size ∝ n instruments)
    geom_point(aes(size = n_iv, shape = factor(source_shape)),
               alpha = 0.90) +
    scale_shape_manual(
      values = c("16" = 16, "18" = 18, "15" = 15),
      labels = c("16" = "MiBioGen", "18" = "DMP", "15" = "Pooled"),
      name   = "Source"
    ) +

    # OR label (right of CI)
    geom_text(aes(x = max(df$OR_hi) * 1.02,
                  label = label_text),
              hjust = 0, size = 2.6, colour = "#222222") +

    col_scale +
    scale_size_continuous(range = c(2, 5), guide = "none") +
    scale_x_continuous(
      expand = expansion(mult = c(0.02, 0.35))
    ) +

    labs(
      title    = sprintf("%s  %s genera (OR %s 1)",
                         title_label,
                         if (direction == "risk") "Risk-increasing" else "Protective",
                         if (direction == "risk") ">" else "<"),
      subtitle = "Bonferroni-significant (p < 3.03×10⁻⁴) | Meta-CAC N = 49,309",
      x        = "Odds Ratio (95% CI) per SD CLR-transformed abundance",
      y        = NULL
    ) +

    theme_bw(base_size = 11) +
    theme(
      axis.text.y    = element_text(face = "bold.italic", size = 10.5),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      legend.position    = "bottom"
    )
}

# ════════════════════════════════════════════════════════════
# GENERATE AND SAVE
# ════════════════════════════════════════════════════════════

risk_df <- mr_tax %>% filter(OR > 1) %>% arrange(OR)
prot_df <- mr_tax %>% filter(OR < 1) %>% arrange(desc(OR))

p_risk <- make_forest_plot(risk_df, direction = "risk",  title_label = "Figure 1")
p_prot <- make_forest_plot(prot_df, direction = "protective", title_label = "Figure 2")

# Risk forest plot
ggsave(
  file.path(CONFIG$paths$fig_dir, "fig1_forest_risk.png"),
  p_risk,
  width  = 14, height = 0.45 * nrow(risk_df) + 3,
  dpi    = CONFIG$fig$dpi,
  limitsize = FALSE
)

# Protective forest plot
ggsave(
  file.path(CONFIG$paths$fig_dir, "fig2_forest_protective.png"),
  p_prot,
  width  = 14, height = 0.45 * nrow(prot_df) + 3,
  dpi    = CONFIG$fig$dpi,
  limitsize = FALSE
)

cat(sprintf("✅  Forest plots saved:\n"))
cat(sprintf("   fig1_forest_risk.png (%d genera)\n", nrow(risk_df)))
cat(sprintf("   fig2_forest_protective.png (%d genera)\n", nrow(prot_df)))
