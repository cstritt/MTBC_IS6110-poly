#!/usr/bin/Rscript

#' Collect and harmonize workflow outputs for downstream analyses.
#'
#' This script reads the main inputs produced by the phylogeny and detettore
#' workflow steps, prepares them into a shared R environment, and saves the
#' resulting objects as an RData image. The goal is to avoid repeated data
#' loading and to ensure that downstream notebooks use the same preprocessed
#' data.
#'
#' The script is designed to be called by Snakemake, but it also supports
#' interactive use by falling back to local test paths when the Snakemake
#' object is not available.

library(ape)
library(ggplot2)
library(reshape2)
library(magrittr)


# ---- Configure input and output paths ----

if (exists("snakemake")) {
  paths <- list(
      # INPUT
      tree             = snakemake@input$tree,
      matrix5          = snakemake@input$matrix5,
      matrix3          = snakemake@input$matrix3,
      matrix_metadata  = snakemake@input$matrix_metadata,
      refins           = snakemake@input$refins,
      anchormap        = snakemake@input$anchormap,
      copy_numbers     = snakemake@input$copy_numbers,
      # PARAMS
      metadata         = snakemake@params$metadata,
      genemap          = snakemake@params$genemap,
      dejesus          = snakemake@params$dejesus,
      resistance       = snakemake@params$resistance,
      # OUTPUT
      outfile          = snakemake@output[[1]]
    )
} else {
  message("No snakemake object found — using default paths for interactive testing")
  paths <- list(
      # INPUT
      metadata = 'metadata/10k/metadata_reduced.tsv',
      tree = 'workflow/results/phylogeny/snp_alignment.rooted.treefile',
      matrix5 = 'workflow/results/detettore/ALL_presence-absence.5prime.tsv',
      matrix3 = 'workflow/results/detettore/ALL_presence-absence.3prime.tsv',
      matrix_metadata = 'workflow/results/detettore/ALL_presence-absence.metadata.tsv',
      refins = 'workflow/results/detettore/ALL_reference_insertions.tsv',
      anchormap = 'workflow/results/detettore/anchor_map.tsv',
      copy_numbers = 'workflow/results/detettore/ALL_copy_numbers.tsv',
      # PARAMS
      genemap = 'data/MTBC0/mtbc0_to_h37rv.genemap.tsv',
      dejesus = 'data/DeJesus2017/orfs.csv',
      resistance = 'data/who_resistance/WHO-UCN-TB-2023.7-eng_genomic_coordinates.txt',
      # OUTPUT
      outfile = 'workflow/results/is6110.RData'
    )
}

# ---- Load results from disk ----

# Read sample metadata.
metadata <- read.csv(paths$metadata, sep='\t')

# Read the rooted phylogenetic tree.
tree <- read.tree(paths$tree)

# Read copy-number estimates for the insertion sites.
copy_numbers <- read.csv(paths$copy_numbers, sep='\t')

# Read the 5' and 3' presence/absence matrices and set the first column as row names.
matrix5 <- read.csv(paths$matrix5, sep='\t')
rownames(matrix5) <- matrix5[,1]
matrix5 <- matrix5[,-1]
print(paste('Number of 5 prime insertion sites: ', ncol(matrix5)))

matrix3 <- read.csv(paths$matrix3, sep='\t')
rownames(matrix3) <- matrix3[,1]
matrix3 <- matrix3[,-1]
print(paste('Number of 3 prime insertion sites: ', ncol(matrix3)))

# Matrix metadata
matrix_metadata <- read.csv(paths$matrix_metadata, sep='\t')

# Anchor-reference map
anchor_map <- read.csv(paths$anchormap, sep='\t', header=F)

# Reference insertions
refins <- read.csv(paths$refins, sep='\t', header=TRUE)
print(paste('Number of reference insertion sites: ', length(unique(refins$position))))



# ---- Define plotting defaults and color palettes ----

# Load LaTeX fonts once if they are available locally.
#library(extrafont)
#loadfonts()
# Run the following command once to install the fonts on the system.
#font_import()

# Set the plotting theme used by downstream figures.

plot_theme = theme_bw(base_family = "CMU Sans Serif", base_size=12)
theme_set(plot_theme)


# Color scales, fonts etc. 

mtb_col_scale = c(
  "L1" = "#ff00ff", 
  "L2" = "#0000ff", 
  "L3" = "#a000cc", 
  "L4" = "#ff0000", 
  "L5" = "#663200", 
  "L6" = "#00cc33", 
  "L7" = "#ede72e", 
  "L8" ="#FF7F00" ,
  "L9" = "#006400",
  "L10" = "black",
  "La1"="#fbf390", 
  "La2"="#c1f9dd",
  "La3"="#dfe1f3",
  "La4"="#f7cfe7"
)

mtb_col_scale_animals_black <- c(
  "A1"="#fbf390", 
  "A2"="#c1f9dd", 
  "A3"="#dfe1f3", 
  "A4"="#f7cfe7", 
  "L1" = "#ff00ff", 
  "L2" = "#0000ff", 
  "L3" = "#a000cc", 
  "L4" = "#ff0000", 
  "L5" = "#663200", 
  "L6" = "#00cc33", 
  "L7" = "#ede72e", 
  "L8" ="#FF7F00" ,
  "L9" = "#006400",
  "L10" = "black",
  "La1"="black",
  "La2" ="black",
  "La3"="black",
  "La4"="black"
)
  
lineage_order = c('L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8', 'L9', 'L10','La1', 'La2', 'La3', 'LaX', 'Dassie bacillus')


# ---- Load and augment the gene map ----

#' Compute the shortest circular distance from a genomic coordinate to the
#' origin of replication in H37Rv.
#'
#' The oriC coordinates are taken from the reference used in the project.
#'
#' @param position A genomic coordinate.
#' @param genome_size Total genome size in base pairs.
#' @param oriC_start Start coordinate of the oriC region.
#' @param oriC_end End coordinate of the oriC region.
#' @return The shortest circular distance to oriC.
dist_to_oriC <- function(position, genome_size = 4435783, oriC_start = 1349, oriC_end = 2161) {
  if (is.na(position)) {
    distance <- NA
  } else if (position < oriC_start) {
    distance <- genome_size - oriC_end + position
  } else if (position > oriC_end) {
    distance <- position - oriC_start
  } else {
    distance <- 0
  }
  return(min(distance, genome_size - distance))
}

#' Add functional annotations to the gene map.
#'
#' This helper merges the gene map with information about essentiality,
#' distance to oriC, and resistance-associated genes.
#'
#' @param genemap A data frame containing the MTBC0-to-H37Rv gene map.
#' @param dejesus_path Path to the DeJesus essentiality table.
#' @param resistance_path Path to the WHO resistance coordinate table.
#' @return A gene map data frame with added annotation columns.
augment_genemap <- function(genemap, dejesus_path, resistance_path) {
  # Add DeJesus essentiality annotations.
  dejesus <- read.csv(dejesus_path, sep = '\t')
  genemap$essentiality <- as.factor(dejesus$Final.Call[match(genemap$h37rv_id, dejesus$ORF.ID)])

  # Add distance to oriC.
  genemap$dist_ori <- sapply(genemap$h37rv_start, dist_to_oriC)

  # Add resistance association.
  resist <- read.delim(resistance_path, sep = '\t')

  resistance_genes <- c()
  for (locus in resist$variant) {
    gene <- strsplit(locus, '_') %>% unlist()
    resistance_genes <- c(resistance_genes, gene[1])
  }
  resistance_genes <- unique(resistance_genes)

  genemap$resistance <- rep('no', nrow(genemap))
  genemap$resistance[which(genemap$gene_id %in% resistance_genes | genemap$h37rv_id %in% resistance_genes)] <- 'yes'

  return(genemap)
}


# Add information to gene metadata: essentiality, distance to origin, resistance.
genemap <- read.delim(paths$genemap, sep='\t')

# Add essentiality and resistance information
genemap = augment_genemap(genemap, paths$dejesus, paths$resistance)

# Add blast hits to matrix metadata ---------------------------
#matrix_metadata <- add_blast_info(matrix_metadata)


# ---- Normalize lineage labels ----
metadata$lineage <- metadata$LINEAGE_x
metadata$sublineage <- metadata$SUBLINEAGE_x
metadata$lineage[which(metadata$lineage == 'C')] <- 'LaX'
metadata$lineage[which(metadata$lineage == 'D')] <- 'Dassie bacillus'


# ---- Rescale branch lengths to account for non-variable positions ----
# Information from https://github.com/cstritt/large_variable_alignment/

genome_size = 4411532  
excluded_positions = 673325
var_positions = 548052 # from workflow/results/phylogeny/snp_alignment.uniqueseq.phy

scaling_factor <- (genome_size - excluded_positions)/ var_positions
tree$edge.length = tree$edge.length * scaling_factor


# ---- Infer MRCA nodes for each lineage ----

# Check monophyly and retrieve the MRCA node for each lineage.
mrcas <- list()
nomono <- list()

for (l in lineage_order){
  md_l <- subset(metadata, lineage==l)
  gnrs <- md_l$GNUMBER[which(md_l$GNUMBER %in% tree$tip.label)]
  
  is_mono <- is.monophyletic(tree, gnrs)
  mrca_node <- getMRCA(tree, gnrs)
  
  if (is_mono){
    mrcas[[l]] <- mrca_node
  }
  else {
    print(paste("Not monophyletic: ", l))
    all_descendants <- extract.clade(tree, mrca_node)$tip.label
    subtree <- keep.tip(tree, all_descendants)
    nomono[[l]] <- subtree
  }
}
mrcas
mrcas = as.data.frame(mrcas) %>% melt()
names(mrcas) <- c("lineage", "node")
mrcas


# Lineages 8, 10, and D are missing as they have a single sample
table(metadata$lineage)
subset(metadata, lineage == 'L8')
subset(metadata, lineage == 'L9')
subset(metadata, lineage == 'L10')
subset(metadata, lineage == 'Dassie bacillus')

# Use L9  mrca to label L9, 10 and Dassie bacillus, which are too close anyway but not monophyletic
singles <- data.frame(
  'lineage' = c('L8', 'L10', 'Dassie bacillus', 'L9/L10/Dassie'),
  'node' = c(
    which(tree$tip.label=='G26689'), 
    which(tree$tip.label=='G50409'),
    which(tree$tip.label=='G01494'),
    18206)
)

mrcas <- rbind(mrcas, singles)

# ---- Clean up and save the final R image ----
outfile = paths$outfile
rm(singles, md_l, paths)
save.image(outfile)
