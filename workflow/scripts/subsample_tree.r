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
write.tree(simtree, outfile)

