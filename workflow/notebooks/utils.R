" 
Library of functions that would clutter the Rmarkdown scripts 
"

# 1_benchmarking ---------------------------------------------------------


# 2_copynumbers ---------------------------------------------------------

infer_ancestral_copy_numbers <- function(){}


# 3_birthdeath ------------------------------------------------------

get_descendants <- function(tree, node) {
  " Get all descendant nodes of a given node in a phylogenetic tree
  "
  children <- tree$edge[tree$edge[, 1] == node, 2]
  if (length(children) == 0) return(NULL)
  descendants <- children
  for (child in children) {
    descendants <- c(descendants, get_descendants(tree, child))
  }
  return(sort(descendants))
}


# 4_genomic_niches -----------------------------------------------------





# 4_chromosome_organization --------------------------------------------

add_blast_info <- function(
    matrix_metadata, 
    query5='workflow/results/detettore/ALL_presence-absence.5prime.fasta', 
    query3='workflow/results/detettore/ALL_presence-absence.3prime.fasta', 
    target='data/MTBC0/MTBC0_CDS_intergenic.fasta'
    ){
  
  blast_five = blastn(query5, target)
  blast_three = blastn(query3, target)
  
  blast_summary = data.frame()
  
  for (i in 1:nrow(matrix_metadata)){
    
    row = matrix_metadata[i,]  
    query = paste(row$side, row$site_id, sep='_')
    
    if (row$side == 5){
      hits = blast_five[blast_five$qseqid==query,]
      query_simple = gsub('5_','', query)
    }
    else {
      hits = blast_three[blast_three$qseqid==query,]
      query_simple = gsub('3_','', query)
    }
    
    if (nrow(hits)==0){next}
    
    besthit = hits[1,]
    
    # Assign category
    pref = substr(besthit$sseqid,1,2)
    if (pref=='mt'){categ = get_alternative_gene_name(besthit$sseqid, genemap)}
    else if (pref=='ig'){categ='intergenic'}
    else if (pref=='id'){categ='repeat'}
    else if (pref=='ge'){categ='pseudogene'}
    
    df = data.frame(
      'site_id' = query_simple,
      'blasthit' = besthit$sseqid,
      'query_len' = besthit$qlen,
      'query_cov' = round(besthit$length/besthit$qlen,2),
      'qstart' = besthit$qstart,
      'qend' = besthit$qend
    )
    blast_summary = rbind(blast_summary, df)
  }
  
  # Add information to matrix_metadata
  matrix_metadata <- merge(matrix_metadata, blast_summary, by = "site_id", all.x = TRUE)
  return(matrix_metadata)
}


augment_genemap <- function(genemap, dejesus_path, resistance_path){
  
  # Add DeJesus essentiality
  dejesus <- read.csv(dejesus_path, sep= '\t')
  genemap$essentiality <- as.factor(dejesus$Final.Call[match(genemap$h37rv_id, dejesus$ORF.ID)])
  
  # Distance to oriC
  genemap$dist_ori <- sapply(genemap$h37rv_start, dist_to_oriC)
  
  
  # Add resistance association
  resist = read.delim(resistance_path, sep='\t')
  
  resistance_genes = c()
  for (locus in resist$variant){
    gene = strsplit(locus, '_') %>% unlist()
    resistance_genes = c(resistance_genes, gene[1])
  }
  resistance_genes = unique(resistance_genes)
  resistance_genes
  
  genemap$resistance = rep('no', nrow(genemap))
  genemap$resistance[which(genemap$gene_id %in% resistance_genes | genemap$h37rv_id %in% resistance_genes)] <- 'yes'  
  
  return(genemap)
}


circlize_plot <- function(genome, bed_list, colors){
  
  circos.clear()
  
  circos.par(start.degree=90)

  circos.genomicInitialize(
    genome, 
    plotType = c("axis"),
    major.by = 250000,
    axis.labels.cex = 0.5,
    sector.names = NULL,
    labels.cex = 1
  )
  
  #circos.genomicDensity(bed_list, track.height = 0.05, col="grey", window.size = 2e4)
  
  circos.genomicTrack(
    bed_list,
    stack=TRUE,
    track.height = 0.7,
    bg.border = NA,
    panel.fun = function(region, value, ...) {
      i = getI(...)
      circos.genomicRect(region, value, col = colors[i], border = colori[i], ...)
    }
  )

  # Show origin of replication
  #circos.genomicTrackPlotRegion(
  #  data.frame(
  #    chr = 'MTBC0',
  #    start = 1349,
  #    end = 2161
  #  ),
  #  ylim = c(0, 1),
  #  track.height = 0.03,
  #  bg.border = NA,
  #  panel.fun = function(region, value, ...) {
  #    circos.genomicRect(region, value, col = "red", border = NA, ...)
  #  }
  #)

  # Add symmetry axis
  #circos.segments(
  #  x0 = 99000, x1 = 99000,
  #  y0 = 0,   y1 = 1,
  #  col = "black",
  #  lwd = 2,
  #  lty = 2,
  #  sector.index = 'MTBC0'
  #)
}

fisher_test <- function(genic_observed, intergenic_observed, chromosome_length, gff_path){
  
  annot <- read.csv(gff_path, sep='\t', comment.char = '#', header=F)
  annot = subset(annot, V3=='gene')
  annot$length <- annot$V5 - annot$V4
  
  genic_total = sum(annot$length)
  intergenic_total = chromosome_length - genic_total
  
  contingency <- matrix(c(intergenic_observed,
                          genic_observed,
                          intergenic_total - intergenic_observed,
                          genic_total - genic_observed),
                        nrow = 2,
                        byrow = TRUE)
  
  rownames(contingency) <- c("Intergenic", "Genic")
  colnames(contingency) <- c("Observed", "NotObserved")
  
  ft <- fisher.test(contingency)
  
  return(list(contingency, ft))
  
}




















"
Symmetry statistic: sum of squared differences between mirrored bins

"

angles_from_positions <- function(positions, genome_length) {
  (2 * pi * (positions %% genome_length) / genome_length) %% (2*pi)
}


symmetry_stat_for_axis <- function(angles, weights = NULL, axis, m = 72) {
  if (is.null(weights)) weights <- rep(1, length(angles))
  
  breaks <- seq(0, 2*pi, length.out = m+1)
  
  # Count in bins
  counts1 <- hist(angles, breaks = breaks, plot = FALSE, include.lowest = TRUE,
                  right = FALSE, weights = weights)$counts
  # Count in mirrored bins
  counts2 <- hist((2*axis - angles) %% (2*pi), breaks = breaks, plot = FALSE,
                  include.lowest = TRUE, right = FALSE, weights = weights)$counts
  
  sum((counts1 - counts2)^2)
}

find_best_axis <- function(angles, weights = NULL, n_axes = 180, m = 72) {
  axes <- seq(0, 2*pi, length.out = n_axes)
  stats <- sapply(axes, function(phi) symmetry_stat_for_axis(angles, weights, phi, m))
  idx <- which.min(stats)
  list(best_axis = axes[idx], best_stat = stats[idx],
       axes = axes, stats = stats)
}

permutation_test_best_axis <- function(angles, weights = NULL,
                                       n_perm = 1000, n_axes = 180, m = 72) {
  obs <- find_best_axis(angles, weights, n_axes = n_axes, m = m)
  obs_stat <- obs$best_stat
  
  perm_stats <- numeric(n_perm)
  for (i in 1:n_perm) {
    rot <- runif(1, 0, 2*pi)
    angles_rot <- (angles + rot) %% (2*pi)
    perm_stats[i] <- find_best_axis(angles_rot, weights, n_axes = n_axes, m = m)$best_stat
  }
  
  pval <- (sum(perm_stats <= obs_stat) + 1) / (n_perm + 1)
  
  list(best_axis = obs$best_axis,
       best_stat = obs_stat,
       axes = obs$axes,
       stats = obs$stats,
       perm_stats = perm_stats,
       pval = pval)
}


plot_axis_scan <- function(res, genome_length = NULL) {
  df <- data.frame(
    axis_angle = res$axes,
    stat = res$stats
  )
  
  p <- ggplot(df, aes(x = axis_angle, y = stat)) +
    geom_line(color = "steelblue", size = 1) +
    geom_vline(xintercept = res$best_axis, color = "red", linetype = "dashed") +
    labs(x = "Axis angle (radians)", y = "Symmetry statistic",
         title = "Symmetry statistic across axes") +
    theme_minimal()
  
  if (!is.null(genome_length)) {
    # add secondary axis in base pairs
    p <- p + scale_x_continuous(
      sec.axis = sec_axis(~ . * genome_length / (2*pi), name = "Axis position (bp)")
    )
  }
  
  p
}

# ----------------------------------
# Polar histogram of weighted counts
# ----------------------------------
plot_circular_histogram <- function(angles, weights = NULL, m = 72,
                                    best_axis = NULL, genome_length = NULL) {
  if (is.null(weights)) weights <- rep(1, length(angles))
  
  breaks <- seq(0, 2*pi, length.out = m+1)
  counts <- hist(angles, breaks = breaks, plot = FALSE, include.lowest = TRUE,
                 right = FALSE, weights = weights)$counts
  
  mids <- head(breaks, -1) + diff(breaks)/2
  df <- data.frame(angle = mids, count = counts)
  
  p <- ggplot(df, aes(x = angle, y = count)) +
    geom_col(width = 2*pi/m, fill = "steelblue", color = "white") +
    coord_polar(start = 0, direction = -1) +
    labs(x = NULL, y = "Weighted count", title = "Circular histogram of insertions") +
    theme_void()
  
  if (!is.null(best_axis)) {
    # add axis line
    df_axis <- data.frame(angle = c(best_axis, best_axis + pi), 
                          r = c(max(df$count), max(df$count)))
    p <- p + geom_line(data = df_axis, aes(x = angle, y = r), inherit.aes = FALSE,
                       color = "red", size = 1)
  }
  p
}


plot_paired_bins <- function(angles, weights = NULL, best_axis, m = 72) {
  if (is.null(weights)) weights <- rep(1, length(angles))
  
  breaks <- seq(0, 2*pi, length.out = m+1)
  mids   <- head(breaks, -1) + diff(breaks)/2
  
  counts <- hist(angles, breaks = breaks, plot = FALSE, include.lowest = TRUE,
                 right = FALSE, weights = weights)$counts
  counts_ref <- hist((2*best_axis - angles) %% (2*pi), breaks = breaks,
                     plot = FALSE, include.lowest = TRUE, right = FALSE,
                     weights = weights)$counts
  
  df <- data.frame(bin = 1:m, mid_angle = mids,
                   observed = counts, mirrored = counts_ref)
  
  library(ggplot2)
  ggplot(df, aes(x = bin)) +
    geom_col(aes(y = observed, fill = "Observed"), position = "dodge") +
    geom_col(aes(y = mirrored, fill = "Mirrored"), position = "dodge") +
    labs(x = "Bin index", y = "Weighted counts",
         title = "Observed vs mirror bins") +
    scale_fill_manual(values = c("Observed" = "steelblue", "Mirrored" = "red")) +
    theme_minimal()
}




# 5_natural_mutagenesis ---------------------------------------------------






# 6_parallelism ------------------------------------------------------------


# 7_niche_model ------------------------------------------------------------








# X_other_stuff ------------------------------------------------------------













# C: birth and death



# Phylogenetic heritability of CN -----------------------------------------

calc_heritability <- function(tree, trait){
  library(ape)
  library(MCMCglmm)
  
  tree <- read.tree("your_tree.newick")

  trait <- read.csv("trait_data.csv", row.names=1)
  trait <- trait[tree$tip.label, , drop=FALSE]
  
  invA <- inverseA(tree)$Ainv
  model <- MCMCglmm(trait ~ 1, random=~animal, ginverse=list(animal=invA), data=trait, pedigree=tree)
  
  heritability <- model$VCV[,"animal"] / rowSums(model$VCV)
  return(heritability)

}



# GWAS utils --------------------------------------------------------------

hit_on_tree <- function(gwas_results, tree, copy_numbers, metadata,variant){
  " Get circular tree, showing in which strains a variant is 
  present and IS copy number as a heatmap
  
  "
  library(ggtree)
  library(viridis)
  library(ggtreeExtra)
  
  gene_burden_df <- data.frame(
    'strain' = copy_numbers$strain,
    'copy_number' = copy_numbers$CN,
    'lineage' = copy_numbers$lineage
  )
  
  variant_df = subset(gwas_results, variant == variant)
  strains_with_variant <- strsplit(variant_df$k.samples, ',') %>% unlist()
  
  gene_burden_df[[variant]] = rep('0', nrow(gene_burden_df))
  gene_burden_df[[variant]][which(gene_burden_df$strain %in% strains_with_variant)] <- '1'
  
  
  p0 <- ggtree(tree, layout='circular', size=0.05) %<+% copy_numbers
  p1 <- p0 + geom_tippoint(aes(col=CN), size=0.7, alpha=0.4, stroke=NA) + 
    scale_color_viridis(option = 'turbo') 
  
  p2 <- p1 + 
    geom_fruit(
      data=gene_burden_df,
      geom=geom_tile,
      mapping=aes_string(y="strain", x=variant, fill=variant),
      offset=0.03,  
      pwidth=0.05
    ) +
    scale_fill_manual(
      values=c("#1B9E77", "#D95F02"),
      guide=guide_legend(keywidth=0.5, keyheight=0.5, order=3)
    )
  
  # Boxplot CN and burden
  #bp <- ggplot(gene_burden_df, aes_string(y=copy_number, fill=variant, x=lineage))  + 
  #  geom_boxplot(position='dodge') + theme_bw() + scale_fill_brewer(palette="Dark2") 
  
  t <- table(gene_burden_df$lineage,gene_burden_df[[variant]])
  
  return(list(
    tree = p2,
    #boxplot = bp,
    table = t
  ))
  
}

# Blast
blastn <- function(query, target){
  
  blast_output <- tempfile(fileext = ".txt")
  
  command <- paste(
    "blastn",
    "-query", query,
    "-subject", target,
    "-out", blast_output,
    "-outfmt", "\"6 qseqid sseqid pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore\"",
    "-max_target_seqs", "1"
  )
  system(command)
  
  blast_df <- read.table(blast_output, header = FALSE, sep = "\t")
  
  colnames(blast_df) <- c(
    "qseqid", "sseqid", "pident", "qlen","length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore"
  )
  
  unlink(blast_output)
  return(blast_df)
}


# Infer birth-death rates per branch --------------------------------------

write_beast_nexus <- function(tree, anc_states, metadata, trait_id, outpath){
  " Create an annotated nexus file for exploration with figtree or other
  interactive tree viewers. To all internal and terminal nodes, the likelihood 
  of the presence of the IS is added. 
  "
  p <- ggtree(tree, size=0.0001) %<+% metadata
  td <- as.treedata(p)
  td@data$state <- anc_states[[trait_id]][match(td@data$node, anc_states$node)]
  
  path = paste(outpath, trait_id, sep='/')
  
  write.beast(td, paste(path, 'anc_states.nex', sep='.'))
}



infer_ancestral_likelihoods <- function(tree, matrix){
  " Returns a data frame where the first colum is the node on the tree, 
  starting with terminal nodes followed by internal. Each additional column 
  gives the likelihood of IS presence at a site (column) in the matrix.
  
  "
  
  library(ape)
  library(castor)
  library(magrittr)
    
  site_states <- data.frame(
    'node' = 1:(length(tree$tip.label) + tree$Nnode)
  )
  
  # Loop over insertion sites and extract likelihood of IS presence
  for (i in 1:ncol(matrix)) {
    
    if (i%%500 == 0){
      print(i)
    }
    
    site_vec <- matrix[, i]
    site_id = colnames(matrix)[i]
    
    # Skip constant or NA-only sites
    if (any(is.na(site_vec)) || length(unique(site_vec)) == 1) next
    
    # Convert to numeric starting with 1
    site_vec_mapped <- map_to_state_space(site_vec)
    
    # Reconstruct ancestral states
    asr_result <- asr_max_parsimony(
      tree, 
      tip_states = site_vec_mapped$mapped_states
    )
    
    # For each node, store likelihood of presence
    lh_presence <- asr_result$ancestral_likelihoods[,2]
    
    # Combine tip and internal node states
    full_states <- c(site_vec, lh_presence)
    full_states <- as.data.frame(full_states)
    names(full_states) <- site_id
    site_states <- cbind(site_states, full_states)
  }
  return(site_states)
}

get_AS_summaries <- function(tree, matrix, anc_likelihoods){
  "
  After inferring ancestral state likelihoods (infer_ancestral_likelihoods), 
  summaries IS births and deaths per branch and per site.
  
  "
  edge_ids <- paste(tree$edge[,1], tree$edge[,2], sep = "-")
  birth_per_branch <- setNames(rep(0, nrow(tree$edge)), edge_ids)
  death_per_branch <- setNames(rep(0, nrow(tree$edge)), edge_ids)
  
  site_ids <- colnames(matrix)
  birth_per_site <- setNames(rep(0,ncol(matrix)), site_ids)
  death_per_site <- setNames(rep(0,ncol(matrix)), site_ids)
  
  for (i in 2:ncol(anc_likelihoods)){
    
    if (i%%100==0){print(i)}
    
    site_id = site_ids[i-1]
    site_states = anc_likelihoods[,i]
    
    for (j in 1:nrow(tree$edge)) {
      
      parent = tree$edge[j, 1]
      p1_parent = site_states[parent]
      p0_parent = 1 - p1_parent
      
      child  = tree$edge[j, 2]
      p1_child = site_states[child]
      p0_child = 1 - p1_child
      
      if (p1_child == p1_parent) next
      
      p_birth = p0_parent * p1_child
      p_death = p1_parent * p0_child
      
      edge_key <- paste(parent, child, sep = "-")
      
      birth_per_branch[edge_key] = birth_per_branch[edge_key] + p_birth
      death_per_branch[edge_key] = death_per_branch[edge_key] + p_death
      
      birth_per_site[site_id] = birth_per_site[site_id] + p_birth
      death_per_site[site_id] = death_per_site[site_id] + p_death
    }
  }
  
  split_edges <- strsplit(edge_ids, "-")
  split_edges <- do.call(rbind, split_edges)
  
  branch_df = data.frame(
    parent = split_edges[,1],
    child = split_edges[,2],
    births = unname(birth_per_branch),
    deaths = unname(death_per_branch)
  )
  
  site_df = data.frame(
    site_id = site_ids,
    births = unname(birth_per_site),
    deaths = unname(death_per_site)
  )
  
  return(list(
    branch_summary = branch_df,
    site_summary = site_df
  ))
}


# Fast binary correlation (phi coefficient)
# Used to test if absence at one site correlates with presence at others
binary_cor <- function(X) {
  # Convert to numeric matrix
  X <- as.matrix(X)
  
  # Compute cross-products
  n <- nrow(X)
  A <- t(X) %*% X          # joint presence
  row_sums <- colSums(X)   # trait totals
  
  B <- outer(row_sums, row_sums, function(a, b) a * b / n)
  C <- sqrt((row_sums * (n - row_sums)) %*% t(row_sums * (n - row_sums)))
  
  phi <- (A - B) / C
  diag(phi) <- NA
  return(phi)
}



is_birth_death <- function(tree, matrix){
  # Not used anymore, use more modular approach: first get ancestral likelihoods
  # for each trait, then parse this data frame for summaries
  library(castor)
  
  # Initialize per-branch insertion/deletion counters
  edge_ids <- paste(tree$edge[,1], tree$edge[,2], sep = "-")
  ins_per_branch <- setNames(rep(0, nrow(tree$edge)), edge_ids)
  del_per_branch <- setNames(rep(0, nrow(tree$edge)), edge_ids)
  
  site_states = c()
  site_ids <- colnames(matrix)
  ins_per_site <- setNames(rep(0,ncol(matrix)), site_ids)
  del_per_site <- setNames(rep(0,ncol(matrix)), site_ids)
  
  # Loop over each insertion site (column)
  for (i in 1:ncol(matrix)) {
    
    site_vec <- matrix[, i]
    site_id = colnames(matrix)[i]
    
    # Skip constant or NA-only sites
    if (any(is.na(site_vec)) || length(unique(site_vec)) == 1) next
    
    # Ensure it's numeric binary. Add 1 because a trait can' be zero with asr_max_parsimony
    site_vec <- as.integer(site_vec)
    site_vec <- as.integer(site_vec) + 1
    
    # Reconstruct ancestral states
    asr_result <- asr_max_parsimony(tree, tip_states = site_vec, Nstates = 2)
    states <- asr_result$ancestral_states  # matrix: Nnodes x Nstates
    assigned_states <- states - 1  # Reconvert to 0/1
    
    # Combine tip and internal node states
    full_states <- c(site_vec-1, assigned_states)
    names(full_states) <- 1:(length(tree$tip.label) + tree$Nnode)
    
    site_states_x = data.frame()
    
    # Count 0->1, 1->0 transitions per edge
    for (j in 1:nrow(tree$edge)) {
      
      parent <- tree$edge[j, 1]
      child  <- tree$edge[j, 2]
      state <- paste(full_states[parent], full_states[child], sep = '_')
      
      row = data.frame(
        'parent' = parent, 
        'child' = child,
        'state' = state
      )
      site_states_x = rbind(site_states_x, row)
      
      # birth of element
      if (full_states[parent] == 0 && full_states[child] == 1) {
        edge_key <- paste(parent, child, sep = "-")
        ins_per_branch[edge_key] <- ins_per_branch[edge_key] + 1
        ins_per_site[site_id] <- ins_per_site[site_id] + 1
      }
      # death of element
      else if (full_states[parent] == 1 && full_states[child] == 0) {
        edge_key <- paste(parent, child, sep = "-")
        del_per_branch[edge_key] <- del_per_branch[edge_key] + 1
        del_per_site[site_id] <- del_per_site[site_id] + 1
      }
    }
    
    # Add transitions for single site to list
    site_states_x$id <- paste(site_states_x$parent, site_states_x$child, sep='-')
  }
  
  branch_summary <- data.frame(
    parent = tree$edge[, 1],
    child = tree$edge[, 2],
    edge_id = paste(lineage,tree$edge[,1], tree$edge[,2], sep = "-"),
    insertions = ins_per_branch,
    deletions = del_per_branch,
    branch_length = tree$edge.length,
    birth_rate = ins_per_branch / tree$edge.length,
    death_rate = del_per_branch / tree$edge.length
  )
  return(list(
    branch_summary = branch_summary,
    site_summary = site_states))
}



# Subset matrix and tree ----------------------------------------------------

subset_tree_matrix <- function(tree, matrix, metadata, l, level){
  " Get a sub tree and matrix for a specified lineage or sublineage
  "
  if (level=='lineage'){
    md_l <- subset(metadata, lineage==l) %>% droplevels()
  }
  else if (level=='sublineage'){
    md_l <- subset(metadata, sublineage==l) %>% droplevels()
  }
  
  gnrs <- md_l$GNUMBER
  
  # Subset tree
  tree_l <- keep.tip(tree, gnrs) 

  # Lineage-specific matrix: remove sites with no insertion, convert to factor
  m_l <- matrix[gnrs,]
  m_l_names <- rownames(m_l)
  m_l <- m_l[,colSums(m_l) != 0]
  
  m_l <- as.matrix(m_l)
  rownames(m_l) <- m_l_names

  return(list(tree_l, m_l))
  
}








# Gene annotation tools ---------------------------------------------------
" Genes in the MTBC0 reference lack references to their Rv orthologs.
Add this information to the reference insertions table. 

"

get_alternative_gene_name <- function(genename, genemap){
  " Replace mtbc0 gene names with gene or Rv names from
  H37Rv. Only use mtbc0 name when the gene is not present in 
  H37Rv.
  "
  
  pref = substr(genename, 1,5)
  
  if (pref=='mtbc0'){
    gene_id = genemap$gene_id[which(genemap$mtbc0_gene_id==genename)]
    if (nchar(gene_id) == 0){
      alt_name = genemap$h37rv_id[which(genemap$mtbc0_gene_id==genename)]
    }
    else {alt_name = gene_id}
  }
  else {
    alt_name = genename
    #rv = genename
    #rv = genemap$h37rv_id[which(genemap$gene_id==genename)]
  }
  
  #if (length(rv) > 1){
  #  rv = rv[[1]]
  #  print(paste('Two genemap entries for ', genename, sep=" "))
  #}
  if (length(alt_name) > 0){return(alt_name)}
  else if (length(alt_name) == 0){return('None')}
}


add_alternative_gene_names <- function(detettore_context, genemap){
  " Replace MTBCO gene names with more useful H37Rv locus tags (modus='rv')
  or gene names (modus='gene').
  
  Complication: detettore output contains mix of mtbco locus tags and
  gene names.
  
  "
  rv_names = c()

  for (context in detettore_context){
    
    # intergenic
    if (grepl(';', context)){
      
      context_splt = strsplit(context, ';') %>% unlist()
      
      geneA = context_splt[1]
      geneB = context_splt[2]
      
      geneA_new = get_alternative_gene_name(geneA, genemap)
      geneB_new = get_alternative_gene_name(geneB, genemap)
      context_new = paste(geneA_new, geneB_new, sep = ';')
      rv_names = c(rv_names, context_new)
    }
    
    # genic
    else {
      gene_new = get_alternative_gene_name(context, genemap)
      rv_names = c(rv_names, gene_new)
    }
  }
  return(rv_names)
}

add_resistance_info <- function(gene_names, resistance_genes){

}


add_gene_and_resistance_info <- function(context_list, genemap, resistance_genes){
  " Go through the gene context provided in the detettore output. Replace
  mtbc0 gene names with Rv gene names, if they exist. Also check if the 
  gene is present in the WHO catalogue of resistance mutations

  "
  context_cat = c()
  rv_names = c()
  resistance = c()
  
  for (context in context_list){
    
    if (grepl(';', context)){
      
      context_cat = c(context_cat, 'intergenic')
      
      context_splt = strsplit(context, ';') %>% unlist()
      geneA = context_splt[1]
      geneB = context_splt[2]
      
      geneA_new = get_alternative_gene_name(geneA, genemap)
      geneB_new = get_alternative_gene_name(geneB, genemap)
      context_new = paste(geneA_new, geneB_new, sep = ';')
      rv_names = c(rv_names, context_new)
      
      if (geneA_new %in% resistance_genes | geneB_new %in% resistance_genes){
        resistance_cat = 'Yes'
      }
      else {resistance_cat = 'No'}
    }
    
    else {
      context_cat = c(context_cat, 'genic')
      gene_new = get_alternative_gene_name(context, genemap)
      rv_names = c(rv_names, gene_new)
      
      if (gene_new %in% resistance_genes){
        resistance_cat = 'Yes'
      }
      else {resistance_cat = 'No'}
    }
    resistance = c(resistance, resistance_cat)
  }
  return(list(context_cat, rv_names, resistance))
}



# Circos plot function ----------------------------------------------------




# Intergenic overrepresentation -------------------------------------------





# Chromosome organization -------------------------------------------------



dist_to_oriC <- function(position, genome_size=4435783, oriC_start=1349, oriC_end=2161){
  "
  Get the distance to the origin of replication in H37Rv.
  oriC coordinates were obtained from Qin et al. 1999.
  "
  if (is.na(position)){distance = NA}
  else if (position < oriC_start) {distance <- genome_size - oriC_end + position} 
  else if (position > oriC_end) {distance <- position - oriC_start} 
  else {distance <- 0}
  return(min(distance, genome_size - distance))
}
