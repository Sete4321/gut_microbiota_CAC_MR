# ============================================================
# config.R
# Central configuration: paths, thresholds, parameters
# Edit paths to match your local setup before running
# ============================================================

# ── Directory structure ─────────────────────────────────────
CONFIG <- list(

  # Input data paths (edit these)
  paths = list(
    mbg_hits    = "data/raw/MBG_allHits_p1e4.txt",          # MiBioGen top hits
    dmp_hits    = "data/raw/dmp_summary_stats_taxa.csv",     # DMP top hits
    cac1_gwas   = "data/raw/GCST90278456_h_tsv.gz",         # CAC1 (Kavousi 2023)
    cac2_gwas   = "data/raw/GCST90503074_h_tsv.gz",         # CAC2/SIS (Gummesson 2025)
    out_dir     = "results/",
    fig_dir     = "figures/"
  ),

  # ── MiBioGen instrument selection ──────────────────────────
  mbg = list(
    p_threshold  = 1e-4,    # MiBioGen provides p < 1e-4 hits
    f_min        = 10,      # minimum F-statistic
    n_mbg        = 18340    # MiBioGen sample size
  ),

  # ── DMP instrument selection ────────────────────────────────
  dmp = list(
    p_threshold  = 1e-4,
    f_min        = 10,
    n_dmp        = 7738
  ),

  # ── Clumping parameters ─────────────────────────────────────
  clump = list(
    window_kb    = 500,     # distance-based window (no LD reference used)
    # For LD-based clumping, use TwoSampleMR::clump_data() with
    # pop = "EUR" and r2 = 0.001
    r2           = 0.001,
    use_ld       = FALSE    # set TRUE to use LD-based clumping via IEU API
  ),

  # ── MR analysis ─────────────────────────────────────────────
  mr = list(
    bonferroni_n   = 165,               # number of taxa tested (before taxonomy audit)
    bonferroni_p   = 0.05 / 165,        # 3.03e-4
    fdr_q          = 0.05,
    n_cac1         = 26909,
    n_cac2         = 24811
  ),

  # ── Taxonomy audit ──────────────────────────────────────────
  taxonomy = list(
    # Strings indicating DMP higher-order taxonomy (not genus-level)
    dmp_exclude_patterns = c(
      "^p__", "^c__", "^o__", "^f__",  # phylum/class/order/family prefixes
      "Proteobacteria$", "Firmicutes$",  # phylum names
      "unknown", "noname"
    ),
    # Known genus-level families for taxonomic pattern analysis
    lach_rum_families = c("Lachnospiraceae", "Ruminococcaceae")
  ),

  # ── Figure aesthetics ───────────────────────────────────────
  fig = list(
    dpi    = 300,
    width  = 14,   # inches
    height = 10,
    format = "png" # "pdf", "svg", "tiff" also supported
  )
)

# Create output directories if they don't exist
for (d in c(CONFIG$paths$out_dir, CONFIG$paths$fig_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

cat("✅  Config loaded.\n")
cat("   Bonferroni threshold: p <", formatC(CONFIG$mr$bonferroni_p, format = "e", digits = 2), "\n")
