#!/usr/bin/Rscript

args = commandArgs(trailingOnly = TRUE)

library(ape)
library(castor)
library(magrittr)

# Load phylogeny
tree <- read.tree(args[1])

# Five prime presence-absence matrix
matrix5 <- read.csv(args[2], sep='\t')
rownames(matrix5) <- matrix5[,1]
matrix5 <- matrix5[,-1]
colnames(matrix5) <- gsub('X', '', colnames(matrix5))

# Matrix metadata
matrix_metadata <- read.csv(args[3], sep='\t')

# Output path
outpath <- args[4]


# Disrupted regions
regions = unique(matrix_metadata$context) %>% na.omit()

# Reorder matrix rows so they match the tree!
phyl_order <- match(tree$tip.label, rownames(matrix5))
matrix5 <- matrix5[phyl_order,]
all(rownames(matrix5) == tree$tip.label)

hotspot_summary = data.frame()

for (region in regions){

  ufs = subset(matrix_metadata, side==5 & context == region)$site_id
  if (length(ufs)==0){next}

  regionmat = matrix5[,ufs]

  if (length(ufs) == 1){regionmat_vec = regionmat}
  else {regionmat_vec = as.integer(rowSums(regionmat) > 0)}
  nstrains = sum(regionmat_vec)

  # Convert to numeric starting with 1
  regionmat_vec_mapped <- map_to_state_space(regionmat_vec)

  # Reconstruct ancestral states
  asr_result <- asr_max_parsimony(tree, tip_states = regionmat_vec_mapped$mapped_states)

  # Which column in ancestral_likelihoods corresponds to presence (1)?
  presence_state_index <- regionmat_vec_mapped$name2index[["1"]]
  # For each node, store likelihood of presence
  lh_presence <- asr_result$ancestral_likelihoods[, presence_state_index]

  # Combine tip and internal node states
  # Use mapped tip states (1=absent, 2=present) converted to 0/1
  tip_states_01 <- ifelse(regionmat_vec_mapped$mapped_states == 1, 0, 1)
  full_states <- c(tip_states_01, lh_presence)
  full_states_binary <- ifelse(full_states > 0.5, 1, 0)

  # Count number of independent birth/insertion events
  nbirths = 0
  ndeaths = 0

  for (j in 1:nrow(tree$edge)) {
        
    parent = tree$edge[j, 1]
    child  = tree$edge[j, 2]

    parent_state = full_states_binary[parent]
    child_state = full_states_binary[child]

    # No change
    if (parent_state == child_state) next
    # Birth
    if (parent_state==0 & child_state==1){nbirths = nbirths+1}
    # Death
    else if (parent_state==1 & child_state==0){ndeaths = ndeaths+1}
  }

  row = data.frame(
    'region' = region,
    'nstrains' = nstrains,
    'nbirths' = nbirths,
    'ndeaths' = ndeaths
  )
  #print(row)
  hotspot_summary = rbind(hotspot_summary, row)
}

write.table(hotspot_summary, outpath, sep='\t', quote=F, row.names = F)