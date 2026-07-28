#!/usr/bin/Rscript

args = commandArgs(trailingOnly = TRUE)

library(ape)

tree <- read.tree(args[1])
outfile <- args[3]

# Remove '' from tip labels
tree$tip.label <- gsub("\\'","" , tree$tip.label)

# Root tree
tree_r <- root(tree, outgroup=args[2], resolve.root = TRUE)

# Drop the outgroup
tree_r <- drop.tip(tree_r, args[2])

# Write to file
write.tree(tree_r, file=outfile)
