# estimate transposition rate

library(ape)
library(phytools)

tree <- read.tree("phylogeny.nwk")  # Load phylogenetic tree
data <- read.table("presence_absence.tsv", header=TRUE, row.names=1)  # TE matrix

# ASR
fit_mk <- fitMk(tree, data, model="ARD")  # Asymmetric gain/loss model
print(fit_mk)  # Shows estimated λ (insertion rate) and μ (loss rate)

#stochastic character mapping for visualization
simmap <- make.simmap(tree, data, model="ARD", nsim=100)
plotSimmap(simmap)




# Birth-death model

library(ape)
library(diversitree)

# Load the tree
tree <- read.tree("phylogeny.nwk")

# Load presence-absence data
data <- read.table("presence_absence.tsv", header=TRUE, row.names=1)

# Convert presence-absence matrix into a named vector
te_states <- as.numeric(data[,1])  # First column of the matrix
names(te_states) <- rownames(data)

# Create a birth-death likelihood function
lik_bd <- make.bd(tree)

# Fit the model, estimating λ (insertion) and μ (loss)
fit_bd <- find.mle(lik_bd, c(lambda=0.1, mu=0.05))  # Initial guesses
print(fit_bd)


# You can compare models with AIC to determine whether rates vary across lineages.


#Time-Dependent Birth-Death Model
library(TESS)

# Define time grid (e.g., split time into intervals for rate variation)
time.grid <- seq(0, max(branching.times(tree)), length.out=10)

# Fit a time-dependent birth-death model
bd_model <- tess.analysis(tree, empirical.hyperprior = TRUE, samplingProbability=1.0)
summary(bd_model)


# Visualize rates on tree
library(phytools)

# Assume branch-wise rates are stored in `fit_bd$par`
rate_vector <- fit_bd$par["lambda"] - fit_bd$par["mu"]  # Net diversification

# Map rates to tree
contMap(tree, rate_vector, fsize=0.8, legend=TRUE)



