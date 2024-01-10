## ---------------------------
##
## Script name: BioF_gene_depths_of_coverage.R
##
## Purpose of script: Plot mean gene depth of coverage, mean genome depth of coverage and species mean depth of coverage for a selected gene. 
##
## Author: Dr. Jonathan J. Giacomini
##
## Date Created: 07-25-2023
##
## Copyright (c) Jonathan J. Giacomini, 2023
## Email: jonjgiacomini@gmail.com
##
## ---------------------------
##
## load up the packages and functions we will need: 

library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(scales)
library(vegan)
library(janitor)
library(egg)   
library(ggdendro)
library(forcats)
library(RColorBrewer)
library(grid)


# source("functions/packages.R")      

## ---------------------------

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### The followiong was run in bash environment to create the gene caller ID files for the BioF gene 
# Need to add an NA if there is no gene caller ID so that R can read the file

for site in PP BM TD; do
for genome in $(cat /Volumes/JJG_Ext_SSD/metapangenomics/genome_list.txt); do
grep_result=$(grep '7-keto-8-aminopelargonate synthetase or related enzyme (BioF) (PDB:1BS0)' "/Volumes/JJG_Ext_SSD/metapangenomics/P_0622_Haemophilus_Aggregatibacter_$site-profile/bin_by_bin/$genome/$genome-gene_calls.txt")
if [ -z "$grep_result" ]; then
echo "NA" > "/Volumes/JJG_Ext_SSD/metapangenomics/P_0622_Haemophilus_Aggregatibacter_$site-profile/bin_by_bin/$genome/$genome-BioF-gene_caller_ID.txt"
else
  echo "$grep_result" | awk -F"\t" '{print $1}' > "/Volumes/JJG_Ext_SSD/metapangenomics/P_0622_Haemophilus_Aggregatibacter_$site-profile/bin_by_bin/$genome/$genome-BioF-gene_caller_ID.txt"
fi
done
done
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####  




# Set up variables for this script


# set up oral site list
site_list <- list("TD", "PP", "BM")

gene_caller_ID_BioF <- "176362"


# focal genome ID
focal_genome_ID <- "H_parainfluenzae_str_CCUG_58848_id_GCA_001679405_1"

# make list of genomes in species group
genome_list <- list("H_influenzae_str_159_HINF_id_GCA_001055565_1",
                    "H_influenzae_str_841_HINF_id_GCA_001058575_1",
                    "H_parainfluenzae_str_1128_HPAR_id_GCA_001053915_1",
                    "H_parainfluenzae_str_1209_HPAR_id_GCA_001053035_1",
                    "H_parainfluenzae_str_137_HINF_id_GCA_001053535_1",
                    "H_parainfluenzae_str_146_HPAR_id_GCA_001053575_1",
                    "H_parainfluenzae_str_155_HPAR_id_GCA_001054475_1",
                    "H_parainfluenzae_str_209_HPAR_id_GCA_001055885_1",
                    "H_parainfluenzae_str_432_HPAR_id_GCA_001055095_1",
                    "H_parainfluenzae_str_488_HPAR_id_GCA_001057005_1",
                    "H_parainfluenzae_str_60884_B_Hi_2_id_GCA_001949895_1",
                    "H_parainfluenzae_str_65114_B_Hi_3_id_GCA_001949885_1",
                    "H_parainfluenzae_str_901_HPAR_id_GCA_001059815_1",
                    "H_parainfluenzae_str_ATCC_9796_id_GCA_001680775_1",
                    "H_parainfluenzae_str_C2004000280_id_GCA_003252725_1",
                    "H_parainfluenzae_str_C2004002727_id_GCA_003252885_1",
                    "H_parainfluenzae_str_C2004002729_id_GCA_003252795_1",
                    "H_parainfluenzae_str_C2005004058_id_GCA_003252915_1",
                    "H_parainfluenzae_str_C2006002596_id_GCA_003252755_1",
                    "H_parainfluenzae_str_C2008001710_id_GCA_003252775_1",
                    "H_parainfluenzae_str_CCUG_58848_id_GCA_001679405_1",
                    "H_parainfluenzae_str_COPD_014_E1_O_id_GCA_009914785_1",
                    "H_parainfluenzae_str_Haemophilus_parainfluenzae_BgEED17_id_GCA_901873375_1",
                    "H_parainfluenzae_str_HK262_id_GCA_000259485_1",
                    "H_parainfluenzae_str_LC_1315_18_id_GCA_008868695_1",
                    "H_parainfluenzae_str_M1C111_2_id_GCA_014982385_1",
                    "H_parainfluenzae_str_M1C113_1_id_GCA_014931475_1",
                    "H_parainfluenzae_str_M1C116_1_id_GCA_014982375_1",
                    "H_parainfluenzae_str_M1C120_2_id_GCA_014931455_1",
                    "H_parainfluenzae_str_M1C125_4_id_GCA_014931435_1",
                    "H_parainfluenzae_str_M1C130_2_id_GCA_014931415_1",
                    "H_parainfluenzae_str_M1C137_2_id_GCA_014931395_1",
                    "H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1",
                    "H_parainfluenzae_str_M1C146_1_id_GCA_014931355_1",
                    "H_parainfluenzae_str_M1C147_1_id_GCA_014931335_1",
                    "H_parainfluenzae_str_M1C149_1_id_GCA_014931315_1",
                    "H_parainfluenzae_str_M1C152_1_id_GCA_014931295_1",
                    "H_parainfluenzae_str_M1C160_1_id_GCA_014931275_1",
                    "H_parainfluenzae_str_M27794_id_GCA_003390455_1",
                    "H_parainfluenzae_str_NCTC_7857_id_GCA_900450845_1",
                    "H_parainfluenzae_str_NCTC10665_id_GCA_900638025_1",
                    "H_parainfluenzae_str_NCTC10672_id_GCA_900450995_1",
                    "H_parainfluenzae_str_T3T1_id_GCA_000210895_1",
                    "H_parainfluenzae_str_UMB0748_id_GCA_002884755_1",
                    "H_sp_str_CCUG_60358_id_GCA_001679485_1",
                    "H_sp_str_HMSC061E01_id_GCA_001810345_1",
                    "H_sp_str_HMSC066D03_id_GCA_001811025_1",
                    "H_sp_str_HMSC068C11_id_GCA_001815355_1",
                    "H_sp_str_HMSC073C03_id_GCA_001814055_1",
                    "H_sp_str_HMSC61B11_id_GCA_001838615_1",
                    "H_sp_str_HMSC71H05_id_GCA_001838635_1")

#write.table(genome_list, file = "/Volumes/JJG_Ext_SSD/metapangenomics/genome_list.txt", quote = FALSE, row.names = FALSE)

## ---------------------------

# per oral site (TD, BM or SUPP), load the three sets of data (gene coverage, genome coverage, species coverage), then merge the data frames into one. 

# For each sample there will be a column value for site, gene depth of coverage, focal genome depth of coverage, species depth of coverage 



# make final df

final_df <- data.frame(sample = character(),
                       site = character(),
                       gene_coverage = numeric(),
                       genome_coverage = numeric(),
                       species_coverage = numeric(),
                       species_median_coverage = numeric())



for (oral_site in site_list) {
  
  
  ##### 1.   load gene coverage data   ##### 
  gene_coverage_df <- read.table(paste0("/Volumes/JJG_Ext_SSD/metapangenomics/P_0622_Haemophilus_Aggregatibacter_",oral_site,"-profile/bin_by_bin/",focal_genome_ID,"/",focal_genome_ID,"-gene_coverages.txt"), header = TRUE)
  
  gene_coverage_df_filtered_1 <- gene_coverage_df %>% 
    dplyr::filter(gene_callers_id == gene_caller_ID_BioF) %>% # select focal gene of interest
    dplyr::select(-gene_callers_id) #drop gene caller ID col
  
  gene_coverage_df_filtered2 <- reshape2::melt(gene_coverage_df_filtered_1) %>%  # melt from wide to long
    dplyr::rename(sample = variable, gene_coverage = value) %>% #rename cols to match main df
    dplyr::mutate(site = oral_site) # make site column
  
  ##### 2.  load genome coverage data (for focal first)  ##### 
  focal_taxa_mean_coverage_df <- read.table(paste0("/Volumes/JJG_Ext_SSD/metapangenomics/P_0622_Haemophilus_Aggregatibacter_",oral_site,"-profile/bin_by_bin/",focal_genome_ID,"/",focal_genome_ID,"-mean_coverage_Q2Q3.txt"), header = TRUE)
  
  focal_taxa_mean_coverage_df_1 <- focal_taxa_mean_coverage_df %>% 
    dplyr::select(-bin) #drop bin col
  
  focal_taxa_mean_coverage_df_2 <- reshape2::melt(focal_taxa_mean_coverage_df_1) %>%  # melt from wide to long
    dplyr::rename(sample = variable, genome_coverage = value) %>% #rename cols to match main df
    dplyr::mutate(site = oral_site) # make site column
  
  
  ##### 3. load species coverage for all of the genomes for the focal taxa   ##### 
  
  # make data frame
  species_coverage <- data.frame(sample = focal_taxa_mean_coverage_df_2$sample,
                                 site = focal_taxa_mean_coverage_df_2$site)
  
  
  for (genome in genome_list) {
    
    get_mean_coverage <- read.table(paste0("/Volumes/JJG_Ext_SSD/metapangenomics/P_0622_Haemophilus_Aggregatibacter_",oral_site,"-profile/bin_by_bin/",genome,"/",genome,"-mean_coverage_Q2Q3.txt"), header = TRUE)
    
    get_mean_coverage_1 <- get_mean_coverage %>% 
      dplyr::select(-bin) #drop bin col
    
    get_mean_coverage_2 <- reshape2::melt(get_mean_coverage_1) %>%  # melt from wide to long
      dplyr::rename(sample = variable, genome_coverage = value) %>% #rename cols to match main df
      dplyr::mutate(site = oral_site) %>% # make site column
      dplyr::select(sample, site, genome_coverage)
    
    colnames(get_mean_coverage_2) <- c('sample','site', genome)
    
    get_mean_coverage_3 <- get_mean_coverage_2 %>% 
      dplyr::select(genome)
    
    # Append to species_coverage data frame
    species_coverage <- cbind(species_coverage, get_mean_coverage_3)
    
  }
  
  
  
  # Calculate median coverage function
  median_coverage <- function(x) {
    median(x, na.rm = TRUE)
  }
  
  # Take average 
  species_coverage <- species_coverage %>% 
    dplyr::mutate(species_coverage = rowMeans(.[, -(1:2)])) %>% 
    dplyr::mutate(species_median_coverage = apply(.[, -(1:2)], 1, median_coverage)) %>% 
    dplyr::select(sample, site, species_coverage, species_median_coverage)
  
  ##### 4. merge with final df    ##### 
  
  
  site_df <- data.frame(sample = gene_coverage_df_filtered2$sample,
                        site = gene_coverage_df_filtered2$site,
                        gene_coverage = gene_coverage_df_filtered2$gene_coverage,
                        genome_coverage = focal_taxa_mean_coverage_df_2$genome_coverage,
                        species_coverage = species_coverage$species_coverage,
                        species_median_coverage = species_coverage$species_median_coverage)
  
  final_df <- rbind(final_df, site_df)
  
}




#####  Gene coverage averaged across all H. para genomes  #####  


for (oral_site in site_list) {
  
  if ("TD" %in% oral_site) {
    # make final data frame for each site
    site_df <- final_df %>% 
      filter(site == "TD")
    
    genomes_gene_coverage_TD <- data.frame(sample = site_df$sample,
                                           site = site_df$site)
  }
  
  if ("BM" %in% oral_site) {
    # make final data frame for each site
    site_df <- final_df %>% 
      filter(site == "BM")
    
    genomes_gene_coverage_BM <- data.frame(sample = site_df$sample,
                                           site = site_df$site)
  }
  
  if ("PP" %in% oral_site) {
    # make final data frame for each site
    site_df <- final_df %>% 
      filter(site == "PP")
    
    genomes_gene_coverage_PP <- data.frame(sample = site_df$sample,
                                           site = site_df$site)
  }
  
  for (genome in genome_list) {
    
    ##### 1. gene coverage  ##### 
    
    # load gene coverage data
    gene_coverage_df <- read.table(paste0("/Volumes/JJG_Ext_SSD/metapangenomics/P_0622_Haemophilus_Aggregatibacter_",oral_site,"-profile/bin_by_bin/",genome,"/",genome,"-gene_coverages.txt"), header = TRUE)
    
    # load gene caller ID for BioF gene for each genome
    gene_caller_ID_df <- read.table(paste0("/Volumes/JJG_Ext_SSD/metapangenomics/P_0622_Haemophilus_Aggregatibacter_",oral_site,"-profile/bin_by_bin/",genome,"/",genome,"-BioF-gene_caller_ID.txt"), header = FALSE)       
    
    gene_caller_ID <- gene_caller_ID_df[1,1]
    
    # get gene coverage for focal gene (e.g. BioF)
    gene_coverage_df_filtered_1 <- gene_coverage_df %>% 
      dplyr::filter(gene_callers_id == gene_caller_ID) %>% # select focal gene of interest
      dplyr::select(-gene_callers_id) #drop gene caller ID col
    
    if (NA %in% gene_caller_ID) {
      # Create a new row of 0's with 183 columns
      new_row <- data.frame(matrix(0, nrow = 1, ncol = ncol(gene_coverage_df_filtered_1)))
      colnames(new_row) <- colnames(gene_coverage_df_filtered_1)
      
      # Append the new row to the data frame
      gene_coverage_df_filtered_1 <- rbind(gene_coverage_df_filtered_1, new_row)
    }
    
    gene_coverage_df_filtered2 <- reshape2::melt(gene_coverage_df_filtered_1) %>%  # melt from wide to long
      dplyr::rename(sample = variable, gene_coverage = value) %>% #rename cols to match main df
      dplyr::mutate(site = oral_site) %>% # make site column
      dplyr::select(sample, site, gene_coverage)
    
    colnames(gene_coverage_df_filtered2) <- c('sample','site', genome)
    
    gene_coverage_df_filtered3 <- gene_coverage_df_filtered2 %>% 
      dplyr::select(genome)
    
    if ("TD" %in% oral_site) {
      # Append to species_coverage data frame
      genomes_gene_coverage_TD <- cbind(genomes_gene_coverage_TD, gene_coverage_df_filtered3)
    }
    
    if ("BM" %in% oral_site) {
      # Append to species_coverage data frame
      genomes_gene_coverage_BM <- cbind(genomes_gene_coverage_BM, gene_coverage_df_filtered3)
    }
    
    if ("PP" %in% oral_site) {
      # Append to species_coverage data frame
      genomes_gene_coverage_PP <- cbind(genomes_gene_coverage_PP, gene_coverage_df_filtered3)
    }
  }
  
}


#Combine sites

genomes_gene_coverages_combo <- rbind(genomes_gene_coverage_TD, genomes_gene_coverage_BM, genomes_gene_coverage_PP)

# Take average 
genomes_gene_coverages_combo_average <- genomes_gene_coverages_combo %>% 
  dplyr::mutate(mean_gene_coverage = rowMeans(.[, -(1:2)])) %>% 
  dplyr::select(sample, site, mean_gene_coverage)

BioF_average_cov <- genomes_gene_coverages_combo_average %>% 
  dplyr::select(mean_gene_coverage) %>% 
  dplyr::rename(BioF_average_cov = mean_gene_coverage)

final_df <- cbind(final_df, BioF_average_cov)



# Set order of sites
final_df$site <- factor(final_df$site, levels = c("TD", "PP", "BM"))


df_tidy <- final_df %>%
  select(sample, site, species_coverage, BioF_average_cov) %>% 
  gather(key = "Category", value = "Value", species_coverage, BioF_average_cov)

# Calculate the order based on gene_coverage values
sample_order <- df_tidy %>%
  filter(Category == "species_coverage") %>%
  arrange(desc(Value)) %>%
  pull(sample)

# Reorder Sample_ID based on the calculated order
df_tidy$sample <- factor(df_tidy$sample, levels = sample_order)

# install.packages("stringr")          # Install stringr package
# library("stringr")                   # Load stringr

Gene_and_species_cov_plot <- ggplot(df_tidy, aes(x = sample, y = Value, fill = Category)) +
  facet_grid(. ~ site, scales = "free", space = "free", switch = "x") + 
  geom_bar(stat = "identity") +
  scale_y_continuous(expand = c(0, 0), limits=c(0, 40), breaks=seq(0,40,by=10), oob=squish, 
                     labels = scales::comma) +
  scale_fill_manual(values = c("blue", "red"))+
  labs(y = paste0("Mean coverage of", "\n", "H. parainfluenzae & BioF gene"),
       x = NULL) +
  theme_classic() +
  theme(text = element_text(size = 12),
        axis.ticks.y = element_line(linewidth = 1),
        axis.text.x = element_blank(),
        axis.text.y =  element_text(color = "black"),
        #axis.line.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.5, color = "black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        panel.background = element_blank(),
        strip.text.x = element_text(color = "black"),
        strip.background = element_rect(colour=NA, fill=NA),
        plot.margin=unit(c(1,0,0.1,0), "cm")) 


#########  #########  #########  #########  #########  #########  #########  #########  #########

df_tidy_focal <- final_df %>%
  select(sample, site, genome_coverage, gene_coverage) %>% 
  gather(key = "Category", value = "Value", genome_coverage, gene_coverage)

# Calculate the order based on gene_coverage values
sample_order_genome <- df_tidy_focal %>%
  filter(Category == "genome_coverage") %>%
  arrange(desc(Value)) %>%
  pull(sample)

# Reorder Sample_ID based on the species coverage
df_tidy_focal$sample <- factor(df_tidy_focal$sample, levels = sample_order_genome)


Gene_and_focal_genome_cov_plot <- ggplot(df_tidy_focal, aes(x = sample, y = Value, fill = Category)) +
  facet_grid(. ~ site, scales = "free", space = "free", switch = "x") + 
  geom_bar(stat = "identity") +
  scale_y_continuous(expand = c(0, 0), limits=c(0, 100), breaks=seq(0,100,by=10), oob=squish, labels = scales::comma) +
  scale_fill_manual(values = c("blue", "red"))+
  xlab(NULL)+
  labs(y = paste0("Mean coverage of", "\n", "H. parainfluenzae SUPP specialist genome", "\n", 
                  "&", "\n", "BioF gene"),
       x = NULL) +
  theme_classic() +
  theme(text = element_text(size = 12),
        axis.ticks.y = element_line(linewidth = 1),
        axis.text.x = element_blank(),
        axis.text.y =  element_text(color = "black"),
        #axis.line.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.5, color = "black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        panel.background = element_blank(),
        strip.text.x = element_text(color = "black"),
        strip.background = element_rect(colour=NA, fill=NA),
        plot.margin=unit(c(1,0,0.1,0), "cm")) 



# combine total mapped reads plot and relab tile plot
final_plot_version_2 <- egg::ggarrange(Gene_and_species_cov_plot,
                                       Gene_and_focal_genome_cov_plot,
                                       ncol = 1,
                                       heights = c(1,1))

# save final plot
ggsave("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/FIGURES/Gene_level/BioF_gene_genome_species_coverage_plot.pdf",final_plot_version_2, width = 10, height = 6)

# Save data
write.csv(df_tidy_focal, file = "/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/FIGURES/Gene_level/BioF_focal_genome_coverage.csv", quote = FALSE, row.names = FALSE)
write.csv(df_tidy, file = "/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/FIGURES/Gene_level/BioF_species_coverage.csv", quote = FALSE, row.names = FALSE)
