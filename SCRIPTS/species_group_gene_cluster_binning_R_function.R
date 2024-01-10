library(janitor)
library(dplyr)

species_group_gene_cluster_binning <- function(group) {
  # load data
  pan_summary_df <- read.table("/Users/home/Veillonella/06_PANGENOME/SUMMARY/P_1221_Veillonella_gene_clusters_summary_subset.txt", header = TRUE, sep = "\t", fill = TRUE)
  pan_summary_df$genome_name <- as.factor(pan_summary_df$genome_name)
  group_df = read.table(paste0("/Users/home/Veillonella/06_PANGENOME/",group,"_collection_list.txt"), header = TRUE, sep = "\t")
  
  # count number of rows in group collection list
  n_group = nrow(group_df)
  
  # filter data by genome group
  df1 <- pan_summary_df %>% 
    filter(genome_name %in% group_df$Genome_ID) %>% 
    droplevels()
  df1$gene_callers_id <- as.factor(df1$gene_callers_id)
  df1$gene_cluster_id <- as.factor(df1$gene_cluster_id)
  
  # Create wide df that tallys the count of gene-cluster_ID for each genome
  df2 <- df1 %>%
    # count number of rows for each combination of genome_name and gene_cluster_id
    group_by(genome_name, gene_cluster_id) %>%
    tally() %>%
    # pivot the protocol names over the columns
    pivot_wider(names_from=gene_cluster_id, values_from=n) %>%
    # replace NA values in all columns with 0
    mutate(across(everything(), .fns=~replace_na(., 0)))
  
  # remove first column
  df4<-df2[,-1]
  
  # convert matrix to binary
  df5 <- as.matrix((df4>0)+0)
  
  # convert matrix to data frame
  df5 <- as.data.frame(df5)
  
  # add back first column
  df5 <- cbind(df2$genome_name, df5)
  
  df6 <-df5 %>% 
    rename(genome_name = 'df2$genome_name')
  
  # Add new row that is the sum of each column 
  df7 <- df6 %>%
    adorn_totals("row")
  
  # drop column 1
  df8 <- df7 %>%
    select(-c(genome_name))
  
  #drop rows (only keep last row  that contains sum total for each col)
  df9 <- df8[-c(1:n_group), ]
  
  #transpose
  df10 <- df9 %>% 
    melt() %>% 
    rename(gene_cluster_id = variable, sum = value)
  
  # add col for bin
  df11 <- df10 %>% 
    mutate(Bin = case_when(sum == n_group ~ paste0(group,"_Core"),
                           sum > 1 & sum < n_group ~ paste0(group,"_Acc"),
                           sum == 1 ~ paste0(group,"_Single")))
  
  # merge with group_pan_summary_df (df1) by gene-cluster-ID
  df12 <- merge(df11, df1, by = "gene_cluster_id")
  
  # write data frame
  write.table(df12, paste0("/Users/home/Veillonella/06_PANGENOME/",group,"_pan_summary.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  
}

species_group_gene_caller_binning <- function(group,focal_genome) {
  
  # load data frame
  focal_df <- read.table(paste0("/Users/home/Veillonella/06_PANGENOME/",group,"_pan_summary.txt"), header = TRUE, sep = "\t")
  
  # For focal genome, filter by the genome ID, subset columns so that we have Bin, bin_name and gene_callers_id
  # rename Bin to group_bin, bin_name to genus_bin 
  # re-order so that gene_callers_id is first column
  focal_df_2 <- focal_df %>% 
    filter(genome_name == focal_genome) %>% 
    select(c(gene_callers_id, Bin, bin_name)) %>% 
    rename(group_bin = Bin, genus_bin = bin_name)
  
  
  # re-code genus_bin vvalues to reflect their status as genus
  focal_df_3 <- focal_df_2 %>% 
    mutate(genus_bin=recode(genus_bin, 'Genus_core'='Genus_core', 'Accessory'='Genus_Acc', 'Singletons' = 'Genus_Single'))
  
  # write data frame
  write.table(focal_df_3, paste0("/Users/home/Veillonella/06_PANGENOME/",focal_genome,"gene_caller_ID_binning.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  
}
