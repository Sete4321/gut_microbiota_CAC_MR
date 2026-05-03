# ============================================================
# run_pipeline.R
# Master script: runs the complete analysis pipeline in order
#
# Usage:
#   Rscript run_pipeline.R
#   or source("run_pipeline.R") from RStudio
#
# Steps:
#   01  Prepare data (load instruments + meta-analyse CAC GWAS)
#   02  Run MR (IVW, Egger, WM) + taxonomy audit
#   03  Sensitivity analyses (Steiger, MVMR, LOO, funnel)
#   04  Taxonomic pattern analysis (family/group-level tests)
#   05  mPGS construction and prediction analysis
#   06  Publication-quality forest plots
# ============================================================

cat("═══════════════════════════════════════════════════════\n")
cat(" Gut Microbiota → CAC: MR Meta-Analysis Pipeline\n")
cat("═══════════════════════════════════════════════════════\n\n")

start_time <- Sys.time()

steps <- list(
  list(script = "R/01_prepare_data.R",             label = "01  Prepare data"),
  list(script = "R/02_run_mr.R",                   label = "02  Run MR + taxonomy audit"),
  list(script = "R/03_sensitivity_analyses.R",     label = "03  Sensitivity analyses"),
  list(script = "R/04_taxonomic_pattern_analysis.R",label = "04  Taxonomic pattern analysis"),
  list(script = "R/05_mpgs_prediction.R",          label = "05  mPGS prediction"),
  list(script = "R/06_forest_plots.R",             label = "06  Forest plots")
)

for (step in steps) {
  cat(sprintf("\n──── %s ────\n", step$label))
  t0 <- Sys.time()
  tryCatch(
    source(step$script, local = new.env()),
    error = function(e) {
      cat(sprintf("❌  ERROR in %s:\n    %s\n", step$script, e$message))
      stop(sprintf("Pipeline stopped at: %s", step$script))
    }
  )
  cat(sprintf("    Done in %.1f seconds\n", as.numeric(Sys.time() - t0, units = "secs")))
}

total_time <- round(as.numeric(Sys.time() - start_time, units = "mins"), 1)
cat(sprintf("\n✅  Pipeline complete in %.1f minutes.\n", total_time))
cat(sprintf("   Results: %s\n", file.path(getwd(), "results/")))
cat(sprintf("   Figures: %s\n", file.path(getwd(), "figures/")))
