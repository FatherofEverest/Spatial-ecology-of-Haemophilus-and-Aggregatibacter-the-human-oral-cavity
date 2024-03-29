
# load libraries
library(dendextend)
library(ape)
library(dplyr)
library(phytools)
library(phylogram)
library(gplots)

# load trees
PANtree<-ape::read.tree("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/09_PANGENOME/09_PANGENOME/P_0622_Haemophilus_Aggregatibacter-98ANI/gene_cluster_frequencies_newick")

Bac_71_MLtree<-ape::read.tree("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/14_PHYLOGENOMICS/REDO/Dereped_98ANI_fastani_Bac71_fasta_outgroup_removed.clean.fa.contree")

# convert to dendrogram objects
dend_Bac_71_MLtree <- ape::chronos(Bac_71_MLtree)
dend_Bac_71_MLtree_2 <- as.dendrogram.phylo(dend_Bac_71_MLtree)
dend_PANtree <- ape::chronos(PANtree)
dend_PANtree_2 <- rev(as.dendrogram.phylo(dend_PANtree))

# Extract labels from pangenome dendrogram 
labels <- dend_PANtree_2 %>% set("labels_to_char") %>% labels 

# create a color labels data frame for the pangenome
colourslabels <- as.data.frame(labels)

# save pangenome labels 
write.csv(colourslabels, "/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/DATA/tanglegram_pan_lables_colors.csv", row.names = FALSE)

# Extract labels from bac71 phylo dendrogram 
bac71_phylo_labels <- dend_Bac_71_MLtree_2 %>% set("labels_to_char") %>% labels 

# create a color labels data frame for the bac71 phylo
bac71_phylo_labels <- as.data.frame(bac71_phylo_labels)

# save bac71 phylo labels 
write.csv(bac71_phylo_labels, "/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/DATA/tanglegram_bac71_phylo_labels_colors.csv", row.names = FALSE)


# manually add colours to the pangenome labels data frame and then re-load into R
pan_labels_colors <- read.csv("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/DATA/tanglegram_pan_lables_colors.csv", header = TRUE)

# manually add colours to the bac71 phylo labels data frame
bac71_phylo_labels_colors <- read.csv("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/DATA/tanglegram_bac71_phylo_labels_colors.csv", header = TRUE)


# set variables for plotting
cols_pan <- col2hex(pan_labels_colors$Colours)
bac71_cols_phy <- col2hex(bac71_phylo_labels_colors$Colours)

# make dendrogram list; second dedrogram will be fixed
bac71_SCG_PAN <- dendlist(dend_Bac_71_MLtree_2 %>% 
    set("labels_col", value = bac71_cols_phy),
    dend_PANtree_2 %>% 
      set("labels_col", value = cols_pan)) 


# test untangle methods
# SCG_PAN %>% dendextend::untangle(method="random") %>% entanglement() #0.07311872
# SCG_PAN %>% dendextend::untangle(method="labels") %>% entanglement() #0.2584771
# SCG_PAN %>% dendextend::untangle(method="ladderize") %>% entanglement() #0.08805879
# SCG_PAN %>% dendextend::untangle(method="step2side") %>% entanglement() #0.02238234 ###WINNER...Usually is
# SCG_PAN %>% dendextend::untangle(method="step1side") %>% entanglement() #0.03265482
# SCG_PAN %>% dendextend::untangle(method="DendSer") %>% entanglement() #0.3936369

# set height for pdf
rows <- as.data.frame(Bac_71_MLtree$tip.label)
height=nrow(rows) *0.0825

# build tanglegram
bac71_SCG_PAN_TANGLEGRAM <- bac71_SCG_PAN %>% dendextend::untangle(method="step1side") 

# set colors of lines; needs to be based on lefthand side dendrogram
bac71_PAN_TANGLEGRAM_labels <- bac71_SCG_PAN_TANGLEGRAM[[1]] %>%  labels
bac71_PAN_TANGLEGRAM_labels <- as.data.frame(bac71_PAN_TANGLEGRAM_labels)
line_colors <- merge(bac71_PAN_TANGLEGRAM_labels, bac71_phylo_labels_colors, by.x="bac71_PAN_TANGLEGRAM_labels", by.y="bac71_phylo_labels", sort=F)
line_colors2 <- as.character(line_colors$Colours)
line_colors3 <- col2hex(line_colors2)

# plot tanglegram
pdf(file = "/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/Bac_71_Phylogeny_vs_Pan_tanglegram_colors.pdf", width = 15, height = height) 

bac71_SCG_PAN_TANGLEGRAM %>% plot(common_subtrees_color_lines=FALSE,
                                  highlight_distinct_edges=FALSE,
                                  common_subtrees_color_branches=FALSE,
                                  highlight_branches_lwd=FALSE,
                                  lwd=3, 
                                  lab.cex = 0.65, edge.lwd = 2, 
                                  margin_inner = 20, 
                                  columns_width = c(1, 0.5, 1), 
                                  axes=FALSE,
                                  color_lines = line_colors3,
                                  main = paste("entanglement =",
                                               round(entanglement(bac71_SCG_PAN_TANGLEGRAM), 4)))



dev.off()






