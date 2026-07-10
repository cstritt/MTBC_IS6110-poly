#!/usr/bin/Rscript

" Collect results, parse them and save as an RData image. 
To avoid repetitions and ensure that all the following analyses use 
exactly the same data imported from the same R image. 
"

library(ape)
library(ggplot2)
library(reshape2)


args = commandArgs(trailingOnly = TRUE)

tree             = args[1]
metadata         = args[2]
matrix5          = args[3]
matrix3          = args[4]
matrix_metadata  = args[5]
refins           = args[6]
anchormap        = args[7]
copy_numbers     = args[8]


# Load results -----------------------------------------------------


# Metadata
metadata <- read.csv(paths$metadata, sep='\t')

# Tree
tree <- read.tree(paths$tree)

# Copy numbers
copy_numbers <- read.csv(paths$copy_numbers, sep='\t')

# Presence/absence matrices
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



# Color schemes ----------------------------------------------

# Load latex fonts (only has to be done once)
#library(extrafont)
#loadfonts()
# run the following command once
#font_import()

# Set plotting theme (not useful, ggsave onyl renders correct font when specified in the plot itself)

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


# Load and augment gene map matching MTBC0 to H37Rv gene names ----------------

augment_genemap <- function(genemap, dejesus_path, resistance_path){
  
  # Add DeJesus essentiality
  dejesus <- read.csv(dejesus_path, sep= '\t')
  genemap$essentiality <- as.factor(dejesus$Final.Call[match(genemap$h37rv_id, dejesus$ORF.ID)])
  
  # Distance to oriC
  genemap$dist_ori <- sapply(genemap$h37rv_start, dist_to_oriC)
  
  
  # Add resistance association
  resist = read.delim(resistance_path, sep='\t')
  
  resistance_genes = c()
  for (locus in resist$variant){
    gene = strsplit(locus, '_') %>% unlist()
    resistance_genes = c(resistance_genes, gene[1])
  }
  resistance_genes = unique(resistance_genes)
  resistance_genes
  
  genemap$resistance = rep('no', nrow(genemap))
  genemap$resistance[which(genemap$gene_id %in% resistance_genes | genemap$h37rv_id %in% resistance_genes)] <- 'yes'  
  
  return(genemap)
}


# Add information to gene metadata: essentiality, distance to origin, resistance
genemap = read.delim("data/MTBC0/mtbc0_to_h37rv.genemap.tsv", sep='\t')

# Add essentiality and resistance information
genemap = augment_genemap(
  genemap,
  'manuscript/literatura/DeJesus2017/orfs.csv',
  'data/who_resistance/WHO-UCN-TB-2023.7-eng_genomic_coordinates.txt'
)


# Add blast hits to matrix metadata ---------------------------
matrix_metadata <- add_blast_info(matrix_metadata)



# Add alternative lineage names -------------------------------
metadata$lineage = metadata$LINEAGE_x
metadata$sublineage = metadata$SUBLINEAGE_x
metadata$lineage[which(metadata$lineage=='C')] <- 'LaX'
metadata$lineage[which(metadata$lineage=='D')] <- 'Dassie bacillus'


# Rescale branch lengths to account for non-variable positions -----------------
#Information from https://github.com/cstritt/large_variable_alignment/

genome_size = 4411532  
excluded_positions = 673325
var_positions = 548052 # from workflow/results/phylogeny/snp_alignment.uniqueseq.phy

scaling_factor <- (genome_size - excluded_positions)/ var_positions
tree$edge.length = tree$edge.length * scaling_factor


# Get MRCA of lineages -----------------------------------------


# Check monophyly and get mrca
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


# Clean up data and save image ---------------------------
rm(singles, row, md_l, ua_asr, ua_meta)
rm(paths)
save.image('is6110.RData')
