#!/usr/bin/Rscript

args = commandArgs(trailingOnly = TRUE)

library(ape)
library(castor)
library(magrittr)

tree <- read.tree(args[1])

matrix <- read.csv(args[2], sep='\t')
rownames(matrix) <- matrix[,1]
matrix <- matrix[,-1]

outpath <- args[3]
prefix <- args[4]
outpath = paste(outpath, prefix, sep='/')


infer_ancestral_likelihoods <- function(tree, matrix){
  " Returns a data frame where the first colum is the node on the tree, 
  starting with terminal nodes followed by internal. Each additional column 
  gives the likelihood of IS presence at a site (column) in the matrix.
  
  "
  
  site_states <- data.frame(
    'node' = 1:(length(tree$tip.label) + tree$Nnode)
  )

  # Reorder matrix rows so they match the tree!
  phyl_order <- match(tree$tip.label, rownames(matrix))
  matrix <- matrix[phyl_order,]

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

ancestral_states <- infer_ancestral_likelihoods(tree, matrix)
summaries <- get_AS_summaries(tree, matrix, ancestral_states)


write.table(ancestral_states, paste(outpath, 'ancestral_likelihoods.tsv', sep='.'), quote=F, sep='\t',row.names=F)
write.table(summaries$site_summary, paste(outpath, 'site_summary.tsv', sep='.'), quote=F, sep='\t',row.names=F)
write.table(summaries$branch_summary, paste(outpath, 'branch_summary.tsv', sep='.'), quote=F, sep='\t',row.names=F)