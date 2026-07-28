#!/usr/bin/Rscript

library(ape)

args = commandArgs(trailingOnly = TRUE)

tree = read.tree(args[1])
metadata = read.csv(args[2], sep='\t')
ntips = as.numeric(args[3])
outfile = args[4]

include = c('L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'La1', 'La3')
n_per_lineage = ntips / length(include)

# Randomly select strains
simstrains = c()

for (l in include){
  gnrs_lineage = subset(metadata, LINEAGE_x==l)$GNUMBER
  gnrs_sample = sample(gnrs_lineage, n_per_lineage)
  simstrains = c(simstrains, gnrs_sample)
}

simtree = keep.tip(tree, simstrains)

# Rescale branch length such that they are the same as used when estimating birth/death rates!

# First, account for non-variable positions
genome_size = 4411532  
excluded_positions = 673325
var_positions = 548052 # from workflow/results/phylogeny/snp_alignment.uniqueseq.phy
scaling_factor <- (genome_size - excluded_positions)/ var_positions
simtree$edge.length = simtree$edge.length * scaling_factor

# Scale branch length by callable genome size (from 0_prologue.Rmd)
callable_genome_size = 4411532  - 673325
simtree$edge.length = simtree$edge.length * callable_genome_size


write.tree(simtree, outfile)

