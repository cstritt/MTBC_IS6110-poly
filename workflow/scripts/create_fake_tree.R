#!/usr/bin/Rscript

library(ape)
library(phytools)
library(ggtree)
library(ggplot2)


# 1. Simulate a random phylogenetic tree with 20 samples
set.seed(123)  # for reproducibility
n_taxa <- 50
tree <- rtree(n_taxa)

# Plot the tree
plot(tree, main = "Simulated Phylogenetic Tree")
tree <- root(tree, 't1')

# Simulate 10 traits evolving on the tree
Q <- matrix(c(-0.1, 0.1,
              0.1, -0.1), 
            nrow = 2, byrow = TRUE)

rownames(Q) <- colnames(Q) <- c("0", "1")

n_traits <- 40
trait_matrix <- matrix(NA, nrow = n_taxa, ncol = n_traits)

for (i in 1:n_traits) {
  trait <- sim.Mk(tree, Q = Q, anc = "0")
  trait_matrix[, i] <- as.numeric(trait)
}
rownames(trait_matrix) <- tree$tip.label
colnames(trait_matrix) <- paste0("Trait_", 1:n_traits)

trait_matrix <- apply(trait_matrix, c(1, 2), as.character)

# Plot
p <- ggtree(tree)

pm <- gheatmap(p, trait_matrix, colnames = FALSE) +  
  scale_fill_manual(values=c('2'='black', '1'='white')) +
  theme(legend.position = "none") +
  scale_x_reverse()

pm