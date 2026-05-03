# ============================================================
# 01_prepare_data.R
# Step 1A: Load and meta-analyse CAC GWAS outcome datasets
# Step 1B: Load and pool microbiome exposure instruments
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})
source("R/config.R")

# ════════════════════════════════════════════════════════════
# 1A. LOAD EXPOSURE INSTRUMENTS
# ════════════════════════════════════════════════════════════

cat("\n── Loading MiBioGen instruments ──\n")

mbg_raw <- fread(CONFIG$paths$mbg_hits, quote = '"')

# Clean column names (MiBioGen wraps names in quotes)
colnames(mbg_raw) <- gsub('"', '', colnames(mbg_raw))

# Rename to standard TwoSampleMR format
mbg <- mbg_raw %>%
  rename(
    taxon        = bac,
    SNP          = rsID,
    chr.exposure = chr,
    pos.exposure = bp,
    other_allele.exposure  = ref.allele,
    effect_allele.exposure = eff.allele,
    beta.exposure  = beta,
    se.exposure    = SE,
    pval.exposure  = P.weightedSumZ,
    samplesize.exposure = N
  ) %>%
  mutate(
    taxon          = gsub('"', '', taxon),
    SNP            = gsub('"', '', SNP),
    other_allele.exposure  = gsub('"', '', other_allele.exposure),
    effect_allele.exposure = gsub('"', '', effect_allele.exposure),
    # F-statistic
    F_stat         = (beta.exposure / se.exposure)^2,
    eaf.exposure   = NA_real_,   # not provided in MiBioGen hits file
    exposure       = taxon,
    id.exposure    = taxon,
    mr_keep.exposure = TRUE
  ) %>%
  filter(
    pval.exposure <= CONFIG$mbg$p_threshold,
    F_stat >= CONFIG$mbg$f_min
  )

cat(sprintf("  MiBioGen: %d SNPs across %d taxa\n",
            nrow(mbg), length(unique(mbg$taxon))))

# ── Load DMP instruments ─────────────────────────────────────
cat("\n── Loading DMP instruments ──\n")

dmp_raw <- fread(CONFIG$paths$dmp_hits)

dmp <- dmp_raw %>%
  rename(
    SNP                    = id,
    chr.exposure           = chr,
    pos.exposure           = pos,
    other_allele.exposure  = ref,
    effect_allele.exposure = alt,
    eaf.exposure           = AF_Allele2,
    beta.exposure          = beta,
    se.exposure            = SE,
    pval.exposure          = pval,
    samplesize.exposure    = num
  ) %>%
  mutate(
    taxon   = short,               # genus-level short name
    F_stat  = (beta.exposure / se.exposure)^2,
    exposure    = taxon,
    id.exposure = taxon,
    mr_keep.exposure = TRUE
  ) %>%
  filter(
    pval.exposure <= CONFIG$dmp$p_threshold,
    F_stat >= CONFIG$dmp$f_min,
    !is.na(SNP)
  )

cat(sprintf("  DMP: %d SNPs across %d taxa\n",
            nrow(dmp), length(unique(dmp$taxon))))

# ════════════════════════════════════════════════════════════
# 1B. DISTANCE-BASED CLUMPING PER TAXON
# (replaces LD clumping when reference panel unavailable)
# For LD-based clumping see 01b_ld_clumping.R
# ════════════════════════════════════════════════════════════

distance_clump <- function(df, window_kb = 500) {
  # Within each taxon, keep the most significant SNP in each window
  # Sort by p-value; greedily retain SNPs >window_kb from any retained SNP
  df <- df %>% arrange(pval.exposure)
  keep <- logical(nrow(df))
  retained_pos <- numeric(0)
  retained_chr <- character(0)

  for (i in seq_len(nrow(df))) {
    chr_i <- as.character(df$chr.exposure[i])
    pos_i <- df$pos.exposure[i]
    # Check distance to all retained SNPs on same chromosome
    same_chr <- retained_chr == chr_i
    if (!any(same_chr)) {
      keep[i] <- TRUE
    } else {
      dists <- abs(retained_pos[same_chr] - pos_i)
      if (all(dists > window_kb * 1000)) keep[i] <- TRUE
    }
    if (keep[i]) {
      retained_pos <- c(retained_pos, pos_i)
      retained_chr <- c(retained_chr, chr_i)
    }
  }
  df[keep, ]
}

cat("\n── Distance-based clumping (", CONFIG$clump$window_kb, "kb) ──\n")

mbg_clumped <- mbg %>%
  group_by(taxon) %>%
  group_modify(~ distance_clump(.x, CONFIG$clump$window_kb)) %>%
  ungroup()

dmp_clumped <- dmp %>%
  group_by(taxon) %>%
  group_modify(~ distance_clump(.x, CONFIG$clump$window_kb)) %>%
  ungroup()

cat(sprintf("  MiBioGen after clumping: %d SNPs\n", nrow(mbg_clumped)))
cat(sprintf("  DMP after clumping: %d SNPs\n", nrow(dmp_clumped)))

# ── Save instruments ─────────────────────────────────────────
fwrite(mbg_clumped, file.path(CONFIG$paths$out_dir, "mbg_instruments.csv"))
fwrite(dmp_clumped, file.path(CONFIG$paths$out_dir, "dmp_instruments.csv"))

cat("\n✅  Instruments saved.\n")

# ════════════════════════════════════════════════════════════
# 1C. LOAD AND META-ANALYSE CAC OUTCOME GWAS
# ════════════════════════════════════════════════════════════

cat("\n── Loading CAC GWAS datasets ──\n")

# Collect all unique SNPs needed
all_snps <- unique(c(mbg_clumped$SNP, dmp_clumped$SNP))
cat(sprintf("  Looking up %d unique instrument SNPs in outcome GWAS\n", length(all_snps)))

load_cac_gwas <- function(path, snp_set, dataset_name, n_default) {
  cat(sprintf("  Reading %s ...\n", dataset_name))
  # fread handles .gz automatically
  dt <- fread(path, showProgress = FALSE)

  # Harmonise column names across the two datasets
  # GCST90278456 columns: chromosome, base_pair_location, effect_allele,
  #   other_allele, beta, standard_error, effect_allele_frequency, p_value,
  #   variant_id, rsid, n
  # GCST90503074: same structure, rsid column

  rsid_col <- intersect(c("rsid", "SNP", "variant_id"), colnames(dt))[1]
  dt <- dt %>%
    rename(
      SNP                   = !!rsid_col,
      beta.outcome          = beta,
      se.outcome            = standard_error,
      effect_allele.outcome = effect_allele,
      other_allele.outcome  = other_allele,
      chr.outcome           = chromosome,
      pos.outcome           = base_pair_location,
      pval.outcome          = p_value
    ) %>%
    mutate(
      eaf.outcome         = if ("effect_allele_frequency" %in% colnames(.))
                              effect_allele_frequency else NA_real_,
      samplesize.outcome  = if ("n" %in% colnames(.)) n else n_default,
      outcome             = dataset_name,
      id.outcome          = dataset_name
    )

  # Filter to instrument SNPs only
  dt_filt <- dt %>% filter(SNP %in% snp_set)
  cat(sprintf("    Found %d / %d instrument SNPs\n", nrow(dt_filt), length(snp_set)))
  dt_filt
}

cac1 <- load_cac_gwas(CONFIG$paths$cac1_gwas, all_snps, "CAC1_Kavousi2023",  CONFIG$mr$n_cac1)
cac2 <- load_cac_gwas(CONFIG$paths$cac2_gwas, all_snps, "CAC2_Gummesson2025", CONFIG$mr$n_cac2)

# ── Fixed-effects IVW meta-analysis at each SNP position ─────
cat("\n── Meta-analysing CAC outcomes ──\n")

harmonise_alleles <- function(b2, ea1, oa1, ea2, oa2) {
  # Flip beta2 if alleles are swapped
  need_flip <- (ea1 == oa2) & (oa1 == ea2)
  b2[need_flip] <- -b2[need_flip]
  b2
}

meta_cac <- cac1 %>%
  select(SNP, chr.outcome, pos.outcome,
         effect_allele.outcome, other_allele.outcome,
         beta.outcome, se.outcome, eaf.outcome,
         samplesize.outcome, pval.outcome) %>%
  rename(b1 = beta.outcome, se1 = se.outcome,
         ea1 = effect_allele.outcome, oa1 = other_allele.outcome,
         n1  = samplesize.outcome) %>%
  inner_join(
    cac2 %>%
      select(SNP, beta.outcome, se.outcome,
             effect_allele.outcome, other_allele.outcome, samplesize.outcome) %>%
      rename(b2 = beta.outcome, se2 = se.outcome,
             ea2 = effect_allele.outcome, oa2 = other_allele.outcome,
             n2  = samplesize.outcome),
    by = "SNP"
  ) %>%
  mutate(
    # Align alleles
    b2_aligned = harmonise_alleles(b2, ea1, oa1, ea2, oa2),
    # Fixed-effects IVW
    w1       = 1 / se1^2,
    w2       = 1 / se2^2,
    beta_meta = (w1 * b1 + w2 * b2_aligned) / (w1 + w2),
    se_meta   = sqrt(1 / (w1 + w2)),
    # Effective N (harmonic mean)
    n_eff     = (n1 * n2) / (n1 + n2) * 2,
    # I² heterogeneity
    Q         = (b1 - beta_meta)^2 / se1^2 + (b2_aligned - beta_meta)^2 / se2^2,
    I2        = pmax(0, (Q - 1) / Q),
    # p-value from z-score
    z_meta    = beta_meta / se_meta,
    pval_meta = 2 * pnorm(-abs(z_meta))
  ) %>%
  select(SNP, chr.outcome, pos.outcome,
         effect_allele.outcome = ea1,
         other_allele.outcome  = oa1,
         eaf.outcome, n_eff,
         beta.outcome  = beta_meta,
         se.outcome    = se_meta,
         pval.outcome  = pval_meta,
         Q_meta = Q, I2_meta = I2)

# Add SNPs only in one dataset
cac1_only <- cac1 %>%
  filter(!SNP %in% cac2$SNP) %>%
  rename(beta.outcome = beta.outcome, se.outcome = se.outcome,
         pval.outcome = pval.outcome) %>%
  mutate(n_eff = samplesize.outcome, Q_meta = NA_real_, I2_meta = NA_real_) %>%
  select(SNP, chr.outcome, pos.outcome, effect_allele.outcome, other_allele.outcome,
         eaf.outcome, n_eff, beta.outcome, se.outcome, pval.outcome, Q_meta, I2_meta)

cac2_only <- cac2 %>%
  filter(!SNP %in% cac1$SNP) %>%
  rename(beta.outcome = beta.outcome, se.outcome = se.outcome,
         pval.outcome = pval.outcome) %>%
  mutate(n_eff = samplesize.outcome, Q_meta = NA_real_, I2_meta = NA_real_) %>%
  select(SNP, chr.outcome, pos.outcome, effect_allele.outcome, other_allele.outcome,
         eaf.outcome, n_eff, beta.outcome, se.outcome, pval.outcome, Q_meta, I2_meta)

meta_cac_full <- bind_rows(meta_cac, cac1_only, cac2_only) %>%
  mutate(outcome    = "meta_CAC",
         id.outcome = "meta_CAC",
         mr_keep.outcome = TRUE)

cat(sprintf("  Meta-CAC: %d SNPs (median N_eff = %d)\n",
            nrow(meta_cac_full),
            as.integer(median(meta_cac_full$n_eff, na.rm = TRUE))))
cat(sprintf("  GWS in meta-CAC (p<5e-8): %d\n",
            sum(meta_cac_full$pval.outcome < 5e-8, na.rm = TRUE)))

fwrite(meta_cac_full, file.path(CONFIG$paths$out_dir, "meta_cac_outcome.csv"))
cat("✅  Meta-CAC outcome saved.\n")
