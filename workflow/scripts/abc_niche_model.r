#!/usr/bin/Rscript

library(processx)

args = commandArgs(trailingOnly = TRUE)
tree = args[1]
nsim = as.numeric(args[2])
metadata = read.csv(args[3])
output = args[4]

tree='workflow/results/niche_model/subsampled_tree.rooted.nex'
nsim = 100
metadata = read.csv('metadata/10k/metadata_reduced.tsv', sep='\t')
output = 'abc_niche_model.RData'

max_niches = 40

get_summary_stats <- function(copy_numbers, metadata){
  copy_numbers$lineage = metadata$LINEAGE_x[match(copy_numbers[,1], metadata$GNUMBER)]
  summdata = summary(copy_numbers[,2]) %>% as.list()

  cladevar = copy_numbers %>% 
    group_by(lineage) %>%
    summarize(Mean = mean(.[[2]]), SD = sd(.[[2]]))

  summdata$cladevar = var(cladevar$Mean)
  summdata$cladesd = sd(cladevar$Mean)

  return(unlist(summdata))
}

run_simulation <- function(parameters){
  p <- processx::run(
    "python",
    c(
      "workflow/scripts/constrained_niche_simulation.py", 
      "--tree", parameters['tree'],
      "--p_lethal", parameters['p_lethal'],
      "--r_birth", parameters['r_birth'],
      "--r_death", parameters['r_death'],
      "--max_niches", parameters['max_niches']
      ),
    error_on_status = TRUE)

  out = read.table(text = p$stdout, header = TRUE, stringsAsFactors = FALSE) 
  return(out)
}

params_mat <- matrix(NA, nrow=nsim, ncol=3)
colnames(params_mat) <- c('p_lethal', 'r_birth', 'r_death')

summaries_mat <- matrix(NA, nrow=nsim, ncol=8) # 8 summary stats

for(i in seq_len(nsim)){

  # Sample from prior distributions
  p_lethal <- runif(1, 0.1,0.99)
  r_birth <- runif(1, 100, 10000)
  r_death <- runif(1, 100, 10000)

  params_mat[i, ] <- c(p_lethal, r_birth, r_death)
  
  params = c(
    "tree" = tree, 
    "p_lethal" = p_lethal,
    "r_birth" = r_birth,
    "r_death" = r_death, 
    "max_niches" = max_niches)
  
  cn_sim <- run_simulation(params)
  cn_stats <- get_summary_stats(cn_sim, metadata)
  
  summaries_mat[i, ] <- cn_stats
}

save.image(output)