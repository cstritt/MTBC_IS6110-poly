# Site Frequency Spectrum functions for IS6110 analysis

#' Calculate Gini coefficient
#'
#' @param x numeric vector of non-negative values
#'
#' @return Gini coefficient (0 = perfect equality, 1 = perfect inequality)
#'
#' @export
gini <- function(x) {
  x <- sort(as.numeric(x))
  n <- length(x)
  if (n == 0 || sum(x) == 0) {
    return(0.0)
  }
  index <- seq_len(n)
  (2 * sum(index * x)) / (n * sum(x)) - (n + 1) / n
}


#' Calculate site frequency spectrum statistics
#'
#' Given a presence-absence matrix where rows are strains and columns are IS sites,
#' compute summary statistics of the site frequency spectrum (SFS).
#'
#' @param mat matrix or data.frame, binary (0/1) with strains in rows, sites in columns.
#'            Can also be -1 (lethal), 0 (absent), 1 (present) from simulations.
#'
#' @return list with components:
#'   - n_singletons: number of sites present in exactly 1 strain
#'   - n_doubletons: number of sites present in exactly 2 strains
#'   - n_rare: number of sites present in 1-5 strains (allele freq < ~5%)
#'   - n_common: number of sites present in >10% of strains
#'   - singleton_prop: proportion of occupied sites that are singletons
#'   - tajimas_d: singleton proportion minus expected under neutral model
#'
#' @export
get_site_frequency_spectrum <- function(mat) {
  
  mat <- as.matrix(mat)
  n_tips <- nrow(mat)
  
  # Count occupancy per site
  occupancy <- colSums(mat == 1)
  
  # Filter to compatible sites (present in at least one strain OR marked compatible in simulation)
  compat <- colSums(mat >= 0) > 0
  occ <- occupancy[compat]
  
  if (length(occ) == 0) {
    return(list(
      n_singletons = 0,
      n_doubletons = 0,
      n_rare = 0,
      n_common = 0,
      singleton_prop = 0.0,
      tajimas_d = 0.0
    ))
  }
  
  n_singletons <- sum(occ == 1)
  n_doubletons <- sum(occ == 2)
  n_rare <- sum(occ <= 5)
  n_common <- sum(occ > 0.1 * n_tips)
  singleton_prop <- n_singletons / length(occ)
  
  # Tajima's D approximation
  # Under neutral model, E[singleton_prop] ≈ 1/log(n_tips)
  expected_singleton <- if (n_tips > 1) 1.0 / log(n_tips) else 0.5
  tajimas_d <- singleton_prop - expected_singleton
  
  list(
    n_singletons = n_singletons,
    n_doubletons = n_doubletons,
    n_rare = n_rare,
    n_common = n_common,
    singleton_prop = singleton_prop,
    tajimas_d = tajimas_d
  )
}


#' Extract all summary statistics from observed IS6110 data
#'
#' Computes the full set of 20 summary statistics used in ABC, matching
#' the Python implementation.
#'
#' @param copy_number_matrix matrix or data.frame, binary (0/1) with strains in rows,
#'                           IS sites in columns.
#' @param metadata data.frame with columns 'GNUMBER' (strain names) and 'LINEAGE_x'
#'                 (lineage assignment). Rows must match columns of copy_number_matrix
#'                 or be provided in the same order.
#'
#' @return data.frame with single row containing all 20 summary statistics:
#'   - cn_mean, cn_std, cn_min, cn_25, cn_50, cn_75, cn_max
#'   - cladevar, cladesd
#'   - hs_gini, hs_n_top5, hs_n_top10, hs_occ_max, hs_occ_mean
#'   - sfs_n_singletons, sfs_n_doubletons, sfs_n_rare, sfs_n_common,
#'     sfs_singleton_prop, sfs_tajimas_d
#'
#' @export
get_all_summary_stats <- function(copy_number_matrix, metadata) {
  
  library(dplyr)
  
  mat <- as.matrix(copy_number_matrix)
  n_strains <- nrow(mat)
  
  # Copy numbers per strain
  cn_arr <- rowSums(mat == 1)
  
  # --- CN statistics ---
  cn_mean <- mean(cn_arr)
  cn_std <- sd(cn_arr)
  cn_min <- min(cn_arr)
  cn_25 <- unname(quantile(cn_arr, 0.25))
  cn_50 <- unname(quantile(cn_arr, 0.50))
  cn_75 <- unname(quantile(cn_arr, 0.75))
  cn_max <- max(cn_arr)
  
  # --- Between-lineage variance ---
  # Match strains to lineages
  strain_names <- rownames(mat)
  if (is.null(strain_names)) {
    strain_names <- metadata$GNUMBER[seq_len(n_strains)]
  }
  
  lineage_vec <- metadata$LINEAGE_x[match(strain_names, metadata$GNUMBER)]
  
  lineage_means <- tibble(cn = cn_arr, lineage = lineage_vec) %>%
    filter(!is.na(lineage)) %>%
    group_by(lineage) %>%
    summarise(mean_cn = mean(cn), .groups = "drop") %>%
    pull(mean_cn)
  
  cladevar <- var(lineage_means, na.rm = TRUE)
  cladesd <- sd(lineage_means, na.rm = TRUE)
  
  # Handle single-lineage case
  if (is.na(cladevar)) cladevar <- 0.0
  if (is.na(cladesd)) cladesd <- 0.0
  
  # --- Hotspot statistics ---
  occupancy <- colSums(mat == 1)
  compat <- colSums(mat >= 0) > 0
  occ_compat <- occupancy[compat]
  freq_compat <- occ_compat / n_strains
  
  hs_gini <- gini(occ_compat)
  hs_n_top5 <- sum(freq_compat >= 0.05)
  hs_n_top10 <- sum(freq_compat >= 0.10)
  hs_occ_max <- max(occ_compat)
  hs_occ_mean <- mean(occ_compat)
  
  # --- Site frequency spectrum ---
  sfs <- get_site_frequency_spectrum(mat)
  
  # Return as single-row data frame
  data.frame(
    cn_mean = cn_mean,
    cn_std = cn_std,
    cn_min = cn_min,
    cn_25 = cn_25,
    cn_50 = cn_50,
    cn_75 = cn_75,
    cn_max = cn_max,
    cladevar = cladevar,
    cladesd = cladesd,
    hs_gini = hs_gini,
    hs_n_top5 = hs_n_top5,
    hs_n_top10 = hs_n_top10,
    hs_occ_max = hs_occ_max,
    hs_occ_mean = hs_occ_mean,
    sfs_n_singletons = sfs$n_singletons,
    sfs_n_doubletons = sfs$n_doubletons,
    sfs_n_rare = sfs$n_rare,
    sfs_n_common = sfs$n_common,
    sfs_singleton_prop = sfs$singleton_prop,
    sfs_tajimas_d = sfs$tajimas_d,
    check.names = FALSE
  )
}


#' Visualize site frequency spectrum
#'
#' @param mat matrix, binary (0/1) with strains in rows, sites in columns
#' @param title character, plot title
#'
#' @return ggplot object
#'
#' @export
plot_sfs <- function(mat, title = "Site Frequency Spectrum") {
  
  library(ggplot2)
  
  mat <- as.matrix(mat)
  n_tips <- nrow(mat)
  occupancy <- colSums(mat == 1)
  compat <- colSums(mat >= 0) > 0
  occ <- occupancy[compat]
  
  freq <- occ / n_tips
  
  freq_df <- data.frame(
    frequency = freq,
    count = table(factor(occ, levels = seq_len(n_tips)))
  )
  
  ggplot(freq_df, aes(x = frequency, y = count)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    theme_bw() +
    labs(
      title = title,
      x = "Site frequency (proportion of strains)",
      y = "Number of sites"
    ) +
    scale_x_continuous(limits = c(0, 1))
}
