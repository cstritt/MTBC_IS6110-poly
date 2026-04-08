#!/usr/bin/Rscript

args = commandArgs(trailingOnly = TRUE)

library(ape)
library(castor)
library(magrittr)

tree <- read.tree(args[1])
refins <- read.delim(args[2], sep='\t', header=TRUE)
outpath <- args[3]


# Create strainXregion matrix --------------------------

strains = tree$tip.label
regions = unique(refins$gene)

# Empty matrix with strains in rows and regions in columns
regionmat = matrix(
  nrow=length(strains), 
  ncol=length(regions),
  dimnames=list(strains, regions)
  )

npos = list()  # Number of positions per region

# Fill matrix
for (i in 1:nrow(refins)){

  if (i%%1000==0){print(i)}

  row = refins[i,]
  strain = row$strain
  region = row$gene
  position = row$position
  strand = row$strand

  regionmat[strain,region] = 1

  # Collect different positions per region
  if (!(region %in% names(npos))){
    npos[[region]] = c(position)
    }
  else {
    npos[[region]] = c(npos[[region]], position)
    }
}

# Convert NA to 0
regionmat[is.na(regionmat)] <- 0


# Loop through regions and infer number of independent insertions --------------------------

hotspot_summary = data.frame()

for (i in 1:ncol(regionmat)){

  if (i%%100==0){print(i)}
  
  site_vec <- regionmat[, i]
  site_id = colnames(regionmat)[i]

  # Skip constant or NA-only sites
  if (any(is.na(site_vec)) || length(unique(site_vec)) == 1) next
      
  # Convert to numeric starting with 1
  site_vec_mapped <- map_to_state_space(site_vec)
  index2name = names(site_vec_mapped$name2index)
  index2name <- setNames(index2name, unname(site_vec_mapped$name2index))

  # Reconstruct ancestral states
  asr_result <- asr_max_parsimony(tree, tip_states = site_vec_mapped$mapped_states)

  # Go through tree and record transitions
  ntips = length(tree$tip.label)
  # Nr. of strains that contain an insertion in the region
  nstrains = sum(site_vec_mapped$mapped_states!=1) 
  # Number of different uFS
  #nstates = site_vec_mapped$Nstates - 1

  # For each node, store likelihood of presence
  lh_presence <- asr_result$ancestral_likelihoods[,2]
  # Combine tip and internal node states
  full_states <- c(unname(site_vec), lh_presence)

  # Count number of independent birth/insertion events
  nbirths = 0
  ndeaths = 0

  for (j in 1:nrow(tree$edge)) {
        
    parent = tree$edge[j, 1]
    child  = tree$edge[j, 2]

    parent_state = full_states[parent]
    child_state = full_states[child]

    # No change
    if (parent_state == child_state) next
    # Birth
    if (parent_state==0 & child_state==1){nbirths = nbirths+1}
    # Death
    else if (parent_state==1 & child_state==0){ndeaths = ndeaths+1}
  }

  row = data.frame(
    'region' = site_id,
    'nstrains' = nstrains,
    'nbirths' = nbirths,
    'ndeaths' = ndeaths
  )
  print(row)
  hotspot_summary = rbind(hotspot_summary, row)
}

write.table(hotspot_summary, outpath, sep='\t', quote=F, row.names = F)