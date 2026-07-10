#!/usr/bin/env Rscript
"
Extract observed summary statistics from empirical IS6110 data.

Matches the Python niche_model_numba.py summary statistics exactly.
Used for ABC tolerance & posterior predictive checks.
"

library(dplyr)
library(tidyr)


#' Compute Gini coefficient
#'
#' @param arr numeric vector
#' @return Gini coefficient (0 = perfect equality, 1 = perfect inequality)
#'
#' @examples
#' gini(c(1, 1, 1, 1))  # 0
#' gini(c(0, 0, 0, 100))  # ~0.75
#'
gini <- function(arr) {
  arr <- sort(as.numeric(arr))
  n <- length(arr)
  if (n == 0 || sum(arr) == 0) return(0.0)
  
  idx <- seq_len(n)
  (2.0 * sum(idx * arr)) / (n * sum(arr)) - (n + 1.0) / n
}


#' Compute site frequency spectrum statistics
#'
#' @param occupancy numeric vector of per-site occupancy counts
#' @param n_tips integer, number of tips in tree
#'
#' @return list with SFS metrics
#'
get_site_frequency_spectrum <- function(occupancy, n_tips) {
  
  if (length(occupancy) == 0 || sum(occupancy) == 0) {
    return(list(
      n_singletons = 0L,
      n_doubletons = 0L,
      n_rare = 0L,
      n_common = 0L,
      singleton_prop = 0.0,
      tajimas_d = 0.0
    ))
  }
  
  n_singletons <- sum(occupancy == 1)
  singleton_prop <- n_singletons / length(occupancy)
  
  # Tajima's D proxy: deviation from Ewens expectation
  tajimas_d <- singleton_prop - (1.0 / log(n_tips) - 0.5)
  if (!is.finite(tajimas_d)) tajimas_d <- 0.0
  
  list(
    n_singletons = as.integer(n_singletons),
    n_doubletons = as.integer(sum(occupancy == 2)),
    n_rare = as.integer(sum(occupancy <= 5)),
    n_common = as.integer(sum(occupancy > 0.1 * n_tips)),
    singleton_prop = as.numeric(singleton_prop),
    tajimas_d = as.numeric(tajimas_d)
  )
}


#' Extract all summary statistics from observed IS6110 data
#'
#' Computes the full set of 16 summary statistics used in ABC, matching
#' the Python niche_model_numba implementation.
#'
#' @param copy_number_matrix matrix or data.frame, binary (0/1) with strains in rows,
#'                           IS sites (insertion regions) in columns.
#' @param metadata data.frame with columns:
#'   - 'GNUMBER': strain/tip names (must match rownames of copy_number_matrix)
#'   - 'LINEAGE_x': lineage assignment for each strain
#'
#' @return data.frame with single row containing all 16 summary statistics:
#'
#'   **Copy number (per-strain IS6110 count):**
#'   - mean_cn, std_cn, min_cn, q25_cn, median_cn, q75_cn, max_cn
#'
#'   **Lineage variance (between-lineage copy number heterogeneity):**
#'   - lineage_var: variance of per-lineage mean copy numbers
#'   - lineage_sd: standard deviation of per-lineage means
#'
#'   **Hotspot occupancy (among compatible sites):**
#'   - gini_occupancy: Gini coefficient (inequality in per-site occupancy)
#'   - n_sites_freq_05: number of sites occupied in ≥5% of tips
#'   - n_sites_freq_10: number of sites occupied in ≥10% of tips
#'   - max_occupancy: maximum per-site occupancy count
#'   - mean_occupancy: mean per-site occupancy
#'
#'   **Site frequency spectrum:**
#'   - n_singletons: insertions at exactly 1 tip
#'   - n_doubletons: insertions at exactly 2 tips
#'   - n_rare: insertions at ≤5 tips
#'   - n_common: insertions at >10% of tips
#'   - singleton_prop: fraction of sites occupied only once
#'   - tajimas_d: deviation from neutral expectation
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # Load empirical data
#'   mat <- read.csv('workflow/results/matrix5.tsv', sep='\t', row.names=1)
#'   metadata <- read.csv('workflow/results/metadata.tsv', sep='\t')
#'
#'   # Extract observed statistics
#'   obs <- get_all_summary_stats(mat, metadata)
#'
#'   # Use in ABC
#'   abc_result <- abc(
#'     target = as.numeric(obs),
#'     param = abc_params,
#'     sumstat = abc_summaries,
#'     tol = 0.01,
#'     method = "rejection"
#'   )
#' }
#'
get_all_summary_stats <- function(copy_number_matrix, metadata) {
  
  mat <- as.matrix(copy_number_matrix)
  n_strains <- nrow(mat)
  n_sites <- ncol(mat)
  
  # ========================================================================
  # COPY NUMBER STATISTICS
  # ========================================================================
  
  # Per-strain IS6110 count
  cn_arr <- rowSums(mat == 1)
  
  mean_cn <- mean(cn_arr)
  std_cn <- sd(cn_arr)
  min_cn <- min(cn_arr)
  q25_cn <- as.numeric(quantile(cn_arr, 0.25))
  median_cn <- as.numeric(quantile(cn_arr, 0.50))
  q75_cn <- as.numeric(quantile(cn_arr, 0.75))
  max_cn <- max(cn_arr)
  
  # ========================================================================
  # LINEAGE VARIANCE
  # ========================================================================
  
  # Map strains to lineages
  strain_names <- rownames(mat)
  if (is.null(strain_names)) {
    strain_names <- seq_len(n_strains)
  }
  
  lineage_vec <- metadata$LINEAGE_x[match(strain_names, metadata$GNUMBER)]
  
  # Per-lineage mean copy number
  lineage_stats <- tibble(cn = cn_arr, lineage = lineage_vec) %>%
    filter(!is.na(lineage)) %>%
    group_by(lineage) %>%
    summarise(mean_cn = mean(cn), .groups = "drop")
  
  lineage_var <- var(lineage_stats$mean_cn, na.rm = TRUE)
  lineage_sd <- sd(lineage_stats$mean_cn, na.rm = TRUE)
  
  # Handle single-lineage case
  if (is.na(lineage_var)) lineage_var <- 0.0
  if (is.na(lineage_sd)) lineage_sd <- 0.0
  
  # ========================================================================
  # HOTSPOT OCCUPANCY STATISTICS
  # ========================================================================
  
  # Per-site occupancy (how many strains have insertion at this site)
  occupancy <- colSums(mat == 1)
  
  # Compatible sites: those with at least one insertion or compatible genotype
  compat <- colSums(mat >= 0) > 0
  occ_compat <- occupancy[compat]
  
  if (length(occ_compat) == 0) {
    gini_occupancy <- 0.0
    n_sites_freq_05 <- 0L
    n_sites_freq_10 <- 0L
    max_occupancy <- 0L
    mean_occupancy <- 0.0
  } else {
    freq_compat <- occ_compat / n_strains
    
    gini_occupancy <- gini(occ_compat)
    n_sites_freq_05 <- sum(freq_compat >= 0.05)
    n_sites_freq_10 <- sum(freq_compat >= 0.10)
    max_occupancy <- max(occ_compat)
    mean_occupancy <- mean(occ_compat)
  }
  
  # ========================================================================
  # SITE FREQUENCY SPECTRUM
  # ========================================================================
  
  sfs <- get_site_frequency_spectrum(occ_compat, n_strains)
  
  # ========================================================================
  # ASSEMBLE OUTPUT
  # ========================================================================
  
  data.frame(
    # Copy number
    mean_cn = mean_cn,
    std_cn = std_cn,
    min_cn = min_cn,
    q25_cn = q25_cn,
    median_cn = median_cn,
    q75_cn = q75_cn,
    max_cn = max_cn,
    
    # Lineage variance
    lineage_var = lineage_var,
    lineage_sd = lineage_sd,
    
    # Hotspot occupancy
    gini_occupancy = gini_occupancy,
    n_sites_freq_05 = n_sites_freq_05,
    n_sites_freq_10 = n_sites_freq_10,
    max_occupancy = max_occupancy,
    mean_occupancy = mean_occupancy,
    
    # Site frequency spectrum
    n_singletons = sfs$n_singletons,
    n_doubletons = sfs$n_doubletons,
    n_rare = sfs$n_rare,
    n_common = sfs$n_common,
    singleton_prop = sfs$singleton_prop,
    tajimas_d = sfs$tajimas_d,
    
    check.names = FALSE,
    row.names = NULL
  )
}


#' Wrapper: Extract observed stats from standard MTBC IS6110 workflow files
#'
#' @param matrix_path path to matrix file (TSV or CSV with strains as rows)
#' @param metadata_path path to metadata (TSV with GNUMBER, LINEAGE_x columns)
#' @param side which matrix to use: 5 or 3 (default: 5)
#'
#' @return data.frame with observed summary statistics
#'
#' @export
get_observed_stats_from_files <- function(matrix_path, metadata_path, side = 5) {
  
  # Load matrix (strain × site binary matrix)
  mat <- read.csv(matrix_path, sep = '\t', row.names = 1)
  
  # Load metadata
  metadata <- read.csv(metadata_path, sep = '\t')
  
  # Extract statistics
  obs <- get_all_summary_stats(mat, metadata)
  
  return(obs)
}


# ============================================================================
# EXAMPLE USAGE (run standalone)
# ============================================================================

if (!interactive()) {
  
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 2) {
    cat("Usage: extract_observed_stats.R <matrix.tsv> <metadata.tsv> [output.tsv]\n")
    cat("\nExample:\n")
    cat("  Rscript extract_observed_stats.R \\
      workflow/results/matrix5.tsv \\
      workflow/results/metadata.tsv \\
      workflow/results/niche_model/observed_stats.tsv\n")
    quit(status = 1)
  }
  
  matrix_path <- args[1]
  metadata_path <- args[2]
  output_path <- args[3]
  
  if (is.na(output_path)) {
    output_path <- "observed_stats.tsv"
  }
  
  cat("Loading matrix from:", matrix_path, "\n")
  cat("Loading metadata from:", metadata_path, "\n\n")
  
  # Extract
  obs <- get_observed_stats_from_files(matrix_path, metadata_path)
  
  cat("Observed summary statistics:\n")
  print(obs)
  
  # Save
  write.csv(obs, output_path, row.names = FALSE)
  cat("\nSaved to:", output_path, "\n")
}

