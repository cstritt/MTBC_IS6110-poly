#!/usr/bin/Rscript

library(ape)
library(ggtree)
library(castor)

args = commandArgs(trailingOnly = TRUE)

tree <- read.tree(args[1])
metadata <- read.csv(args[2], sep='\t')
copy_numbers <- read.csv(args[3], sep='\t')
outdir <- args[4]

# Reorder matrix rows so they match the tree!
phyl_order <- match(tree$tip.label, copy_numbers$strain)
copy_numbers_ordered <- as.integer(copy_numbers[phyl_order,]$CN)
copy_numbers_ordered <- copy_numbers_ordered + 1

# Convert to numeric starting with 1
site_vec_mapped <- map_to_state_space(copy_numbers_ordered)

# Reconstruct ancestral states
asr_result <- asr_max_parsimony(tree, tip_states = copy_numbers_ordered)

# Get the state with the highest likelihood
ancestral_states = c()
for (i in nrow(asr_result$ancestral_likelihoods)){
  lhs = asr_result$ancestral_likelihoods[i,]
  max_i = which(lhs==max(lhs))
  if (length(max_i) > 1){anc_state = NA}  # ambiguous states
  else {anc_state = max_i -1 }
  ancestral_states = c(ancestral_states, anc_state)
}
ancestral_states

# Create treedata object, with a copy number for each internal node
p <- ggtree(tree, layout='circular',size=0.1, alpha=0.5)  %<+% metadata
p$data$copy_number <- c(copy_numbers_ordered, asr_result$ancestral_states)

# Write to tsv and beast format (for visualisation with figtree)
write.table(p$data, paste(outdir, 'treedata.CN_ASR.tsv', sep='/'), sep='\t', quote=F, row.names=F)

#td <- as.treedata(p)
#write.beast(td, 'workflow/results/phylogeny/snp_alignment.rooted.treefile.annotated') 