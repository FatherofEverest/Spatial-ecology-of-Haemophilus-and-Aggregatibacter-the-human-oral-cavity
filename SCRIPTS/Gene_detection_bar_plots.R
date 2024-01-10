## ---------------------------
##
## Script name: Gene_detection_bar_plots
##
## Purpose of script: Generate gene detection bar plots for Haemophilus and Aggregatibacter metapangenomics study
##
## Author: Dr. Jonathan Giacomini
##
## Date Created: 2023-05-05
##
## Copyright (c) Jonathan Giacomini, 2023
## Email: jonjgiacomini@gmail.com
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

## set working directory for Mac and PC

setwd("/Users/home")      # Jon's working directory (mac)

## ---------------------------

## load up the packages we will need: 

library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)


## ---------------------------



genomes = list("H_parahaemolyticus_str_NCTC_8479_id_GCA_900450675_1",
  "H_parahaemolyticus_str_NCTC10794_id_GCA_900450875_1",
  "H_parahaemolyticus_str_C2006000788_id_GCA_003252925_1",
  "H_parahaemolyticus_str_G321_id_GCA_000826045_1",
  "H_parahaemolyticus_str_C2010039593_id_GCA_003253075_1",
  "H_sp_str_Marseille_Q0026_id_GCA_916856725_1",
  "H_paraphrohaemolyticus_str_NCTC10671_id_GCA_900451065_1",
  "H_paraphrohaemolyticus_str_C2014016342_id_GCA_003253115_1",
  "H_sputorum_str_CCUG_13788_id_GCA_000238795_2",
  "H_sputorum_str_C2015005679_id_GCA_003252875_1",
  "H_sputorum_str_C2002001239_id_GCA_003252715_1",
  "H_ducreyi_str_VAN1_id_GCA_001647655_1",
  "H_ducreyi_str_CCUG_4438_id_GCA_002015155_1",
  "A_actinomycetemcomitans_str_SA15_id_GCA_021309365_1",
  "A_actinomycetemcomitans_str_HK_974_id_GCA_008085185_1", 
  "A_actinomycetemcomitans_str_NCTC_9710_id_GCA_008085305_1", 
  "A_actinomycetemcomitans_str_ANH9776_id_GCA_001596315_1",
  "H_massiliensis_str_NA_id_GCA_000752475_1",
  "A_aphrophilus_str_PN_517_id_GCA_003130125_1",
  "A_aphrophilus_str_ATCC_7901_id_GCA_001680805_1", 
  "A_aphrophilus_str_HK83_id_GCA_003703745_1", 
  "A_aphrophilus_str_HK_319_id_GCA_003130365_1",
  "A_aphrophilus_str_NJ8700_id_GCA_000022985_1", 
  "A_aphrophilus_str_FDAARGOS_248_id_GCA_002083125_2", 
  "A_aphrophilus_str_NCTC5906_id_GCA_900636915_1",
  "A_kilianii_str_PN_491_id_GCA_003129965_1", 
  "A_aphrophilus_str_C2008003249_id_GCA_003252995_1", 
  "A_kilianii_str_PN_528_id_GCA_003130255_1", 
  "A_aphrophilus_str_C2010020251_id_GCA_003252835_1", 
  "A_kilianii_str_PN_649_id_GCA_003130015_1", 
  "A_segnis_str_HK_296_id_GCA_003130075_1", 
  "A_segnis_str_C2001002503_id_GCA_003252685_1", 
  "A_segnis_str_C2000002669_id_GCA_003252645_1", 
  "A_sp_str_Marseille_P9115_id_GCA_920939465_1", 
  "A_segnis_str_NCTC10977_id_GCA_900476035_1", 
  "A_segnis_str_PN_651_id_GCA_003130325_1", 
  "A_segnis_str_PN_450_id_GCA_003130345_1", 
  "A_sp_str_W11186_id_GCA_018128305_1", 
  "A_segnis_str_933_AAPH_id_GCA_001059425_1", 
  "A_sp_str_2125159857_id_GCA_017798005_1", 
  "A_sp_HMT_458_str_W10330_id_GCA_000466335_1", 
  "H_parainfluenzae_str_NCTC10672_id_GCA_900450995_1", 
  "H_parainfluenzae_str_65114_B_Hi_3_id_GCA_001949885_1", 
  "H_parainfluenzae_str_60884_B_Hi_2_id_GCA_001949895_1", 
  "H_parainfluenzae_str_Haemophilus_parainfluenzae_BgEED17_id_GCA_901873375_1", 
  "H_parainfluenzae_str_155_HPAR_id_GCA_001054475_1", 
  "H_parainfluenzae_str_1128_HPAR_id_GCA_001053915_1", 
  "H_influenzae_str_159_HINF_id_GCA_001055565_1", 
  "H_sp_str_HMSC073C03_id_GCA_001814055_1", 
  "H_pittmaniae_str_NCTC13334_id_GCA_900186995_1", 
  "H_parainfluenzae_str_C2006002596_id_GCA_003252755_1", 
  "H_influenzae_str_841_HINF_id_GCA_001058575_1", 
  "H_parainfluenzae_str_NCTC10665_id_GCA_900638025_1", 
  "H_parainfluenzae_str_1209_HPAR_id_GCA_001053035_1", 
  "H_parainfluenzae_str_LC_1315_18_id_GCA_008868695_1", 
  "H_parainfluenzae_str_M27794_id_GCA_003390455_1", 
  "H_parainfluenzae_str_COPD_014_E1_O_id_GCA_009914785_1", 
  "H_parainfluenzae_str_M1C116_1_id_GCA_014982375_1", 
  "H_parainfluenzae_str_C2005004058_id_GCA_003252915_1", 
  "H_parainfluenzae_str_C2004002727_id_GCA_003252885_1", 
  "H_parainfluenzae_str_ATCC_9796_id_GCA_001680775_1",
  "H_parainfluenzae_str_137_HINF_id_GCA_001053535_1",
  "H_parainfluenzae_str_UMB0748_id_GCA_002884755_1",
  "H_sp_str_HMSC61B11_id_GCA_001838615_1",
  "H_parainfluenzae_str_HK262_id_GCA_000259485_1",
  "H_parainfluenzae_str_CCUG_58848_id_GCA_001679405_1",
  "H_parainfluenzae_str_NCTC_7857_id_GCA_900450845_1", 
  "H_parainfluenzae_str_488_HPAR_id_GCA_001057005_1",
  "H_parainfluenzae_str_M1C160_1_id_GCA_014931275_1", 
  "H_parainfluenzae_str_M1C152_1_id_GCA_014931295_1", 
  "H_parainfluenzae_str_T3T1_id_GCA_000210895_1", 
  "H_parainfluenzae_str_M1C137_2_id_GCA_014931395_1",
  "H_parainfluenzae_str_M1C120_2_id_GCA_014931455_1", 
  "H_parainfluenzae_str_M1C130_2_id_GCA_014931415_1", 
  "H_parainfluenzae_str_M1C125_4_id_GCA_014931435_1", 
  "H_parainfluenzae_str_209_HPAR_id_GCA_001055885_1", 
  "H_parainfluenzae_str_C2004000280_id_GCA_003252725_1",
  "H_sp_str_HMSC068C11_id_GCA_001815355_1", 
  "H_parainfluenzae_str_C2004002729_id_GCA_003252795_1",
  "H_parainfluenzae_str_M1C146_1_id_GCA_014931355_1",
  "H_sp_str_HMSC061E01_id_GCA_001810345_1", 
  "H_parainfluenzae_str_901_HPAR_id_GCA_001059815_1", 
  "H_parainfluenzae_str_432_HPAR_id_GCA_001055095_1", 
  "H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1", 
  "H_parainfluenzae_str_M1C113_1_id_GCA_014931475_1", 
  "H_parainfluenzae_str_M1C149_1_id_GCA_014931315_1", 
  "H_parainfluenzae_str_M1C147_1_id_GCA_014931335_1", 
  "H_parainfluenzae_str_146_HPAR_id_GCA_001053575_1", 
  "H_parainfluenzae_str_M1C111_2_id_GCA_014982385_1", 
  "H_sp_str_HMSC066D03_id_GCA_001811025_1", 
  "H_sp_str_CCUG_60358_id_GCA_001679485_1",
  "H_sp_str_HMSC71H05_id_GCA_001838635_1", 
  "H_parainfluenzae_str_C2008001710_id_GCA_003252775_1",
  "H_haemolyticus_str_M26167_id_GCA_003493125_1",
  "H_sp_str_C1_id_GCA_001276515_1",
  "H_haemolyticus_str_CCUG30047_id_GCA_013401795_1", 
  "H_haemolyticus_str_M26161_id_GCA_003493685_1", 
  "H_haemolyticus_str_S32F2_id_GCA_013401805_1", 
  "H_haemolyticus_str_M19079_id_GCA_003490595_1",
  "H_haemolyticus_str_M19122_id_GCA_003493365_1", 
  "H_haemolyticus_str_M19080_id_GCA_003494285_1", 
  "H_haemolyticus_str_M19155_id_GCA_003490305_1", 
  "H_haemolyticus_str_M19140_id_GCA_003490235_1", 
  "H_haemolyticus_str_M19071_id_GCA_003493965_1", 
  "H_haemolyticus_str_M26166_id_GCA_003492745_1", 
  "H_haemolyticus_str_M19197_id_GCA_003492365_1", 
  "H_haemolyticus_str_M19164_id_GCA_003491025_1", 
  "H_haemolyticus_str_M19501_id_GCA_000222025_2", 
  "H_haemolyticus_str_M19187_id_GCA_003493605_1", 
  "H_haemolyticus_str_M19107_id_GCA_000222005_2",
  "H_haemolyticus_str_11P18_id_GCA_001008205_1", 
  "H_haemolyticus_str_M19066_id_GCA_003494695_1", 
  "H_haemolyticus_str_M19345_id_GCA_003351405_1", 
  "H_haemolyticus_str_M11818_id_GCA_003494635_1", 
  "H_haemolyticus_str_M19099_id_GCA_003490935_1", 
  "H_haemolyticus_str_CCUG_39154_id_GCA_001679445_1", 
  "H_haemolyticus_str_CCUG_24149_id_GCA_001679135_1", 
  "H_haemolyticus_str_M26157_id_GCA_003494105_1", 
  "H_haemolyticus_str_M25342_id_GCA_003493465_1", 
  "H_haemolyticus_str_27P25_id_GCA_001008275_1", 
  "H_haemolyticus_str_M26160_id_GCA_003494485_1", 
  "H_haemolyticus_str_3P5_id_GCA_001008225_1", 
  "H_haemolyticus_str_ATCC_33390_id_GCA_004368535_1",
  "H_haemolyticus_str_16_549009_id_GCA_004362455_1",
  "H_haemolyticus_str_M26174_id_GCA_003493245_1", 
  "H_haemolyticus_str_M21127_id_GCA_000222045_2", 
  "H_haemolyticus_str_NCTC10839_id_GCA_900477945_1", 
  "H_haemolyticus_str_1P26_id_GCA_001008215_1", 
  "H_haemolyticus_str_2019_19_id_GCA_019973675_1", 
  "H_sp_str_F0397_id_GCA_000242295_1", 
  "H_haemolyticus_str_M28908_id_GCA_003494545_1", 
  "H_haemolyticus_str_M26156_id_GCA_003490655_1", 
  "H_haemolyticus_str_M26164_id_GCA_003494525_1", 
  "H_haemolyticus_str_HI2028_id_GCA_004368395_1", 
  "H_haemolyticus_str_M21621_id_GCA_000222065_2", 
  "H_sp_str_SZY_H36_id_GCA_018829625_1", 
  "H_haemolyticus_str_60971_B_Hi_3_id_GCA_006439235_1", 
  "H_sp_HMT_036_str_F0629_id_GCA_002998595_1", 
  "H_haemolyticus_str_Haemophilus_haemolyticus_BgEED18_id_GCA_901873495_1", 
  "H_sp_str_SZY_H8_id_GCA_018829755_1", 
  "H_influenzae_str_839_HINF_id_GCA_001059185_1", 
  "H_haemolyticus_str_PN24_id_GCA_013401825_1", 
  "H_sp_HMT_036_str_ccug_66565_id_GCA_001679495_1", 
  "H_haemolyticus_str_CCUG_31732_id_GCA_006439125_1", 
  "H_haemolyticus_str_60982_B_Hi_1_id_GCA_006439215_1", 
  "H_haemolyticus_str_60819_B_Hi_1_id_GCA_006439245_1", 
  "H_haemolyticus_str_60824_B_Hi_4_id_GCA_006439315_1", 
  "H_haemolyticus_str_65117_B_Hi_3_id_GCA_006439285_1", 
  "H_haemolyticus_str_CCUG_15949_id_GCA_006439275_1", 
  "H_seminalis_str_SZY_H1_id_GCA_006384255_1", 
  "H_haemolyticus_str_CCUG_30218_id_GCA_006439145_1", 
  "H_influenzae_str_NCTC11931_id_GCA_900475535_1", 
  "H_influenzae_str_PTHi_3421_id_GCA_900407925_1", 
  "H_influenzae_str_PTHi_5709_id_GCA_900408095_1", 
  "H_influenzae_str_P603_4482_id_GCA_003415355_2", 
  "H_influenzae_str_HI1408_id_GCA_001184705_1",
  "H_influenzae_str_1_id_GCA_003203015_1",
  "H_influenzae_str_124P40H1_id_GCA_002991045_1", 
  "H_influenzae_str_1P16H5_id_GCA_002991295_1",
  "H_influenzae_str_WAPHL1_id_GCA_002237715_1", 
  "H_influenzae_str_M21384_id_GCA_003351425_1",
  "H_influenzae_str_PN134_id_GCA_013402015_1", 
  "H_influenzae_str_149P2H1_id_GCA_002991135_1", 
  "H_influenzae_str_39P18H1_id_GCA_002990315_1", 
  "H_influenzae_str_39P25H_id_GCA_002991365_1", 
  "H_influenzae_str_19P81H1_id_GCA_002985065_1", 
  "H_aegyptius_str_CCUG_628_id_GCA_001679305_1", 
  "H_aegyptius_str_ATCC_11116_id_GCA_000195005_1", 
  "H_influenzae_str_F3039_id_GCA_005888205_1", 
  "H_influenzae_str_P608_8895_id_GCA_003415305_1", 
  "H_influenzae_str_51P4H_id_GCA_002985485_1", 
  "H_influenzae_str_P645_8193_id_GCA_003414225_1",
  "H_influenzae_str_M10489_id_GCA_003494325_1", 
  "H_influenzae_str_M03842_id_GCA_003491005_1", 
  "H_influenzae_str_NCTC8143_id_GCA_001457655_1", 
  "H_influenzae_str_492_HINF_id_GCA_001076835_1", 
  "H_influenzae_str_22_4_21_id_GCA_000169855_1", 
  "H_influenzae_str_CHBN_IV_1_id_GCA_019703735_1", 
  "H_influenzae_str_M14666_id_GCA_003494725_1", 
  "H_influenzae_str_PTHi_1402_id_GCA_900407805_1", 
  "H_influenzae_str_60373_BAL_Hi1_id_GCA_004802295_1", 
  "H_influenzae_str_PittHH_id_GCA_000169815_1", 
  "H_influenzae_str_RHH_3_id_GCA_015694745_1", 
  "H_influenzae_str_P612_8066_id_GCA_003415255_1", 
  "H_influenzae_str_148P4H1_id_GCA_002985005_1", 
  "H_influenzae_str_12P37H2_id_GCA_002984945_1", 
  "H_influenzae_str_M25147_id_GCA_003496925_1", 
  "H_influenzae_str_PTHi_13064_id_GCA_900407625_1", 
  "H_influenzae_str_P596_8591_id_GCA_003415525_1", 
  "H_influenzae_str_40P44H1_id_GCA_002985305_1", 
  "H_influenzae_str_PTHi_2060_id_GCA_900407905_1", 
  "H_influenzae_str_40P92H1_id_GCA_002985345_1", 
  "H_influenzae_str_CGSHiCZ412602_id_GCA_000698365_1", 
  "H_influenzae_str_24_id_GCA_003203275_1", 
  "H_influenzae_str_60370_BAL_Hi2_id_GCA_004802355_1", 
  "H_influenzae_str_AS1_id_GCA_021391215_1", 
  "H_influenzae_str_P658_8889_id_GCA_003414605_1", 
  "H_influenzae_str_P677_2580_id_GCA_003414465_1", 
  "H_influenzae_str_AS012764_id_GCA_010604085_1", 
  "H_influenzae_str_60295_NP_Hi1_id_GCA_004802075_1", 
  "H_influenzae_str_60294_NP_Hi1_id_GCA_004802405_1",
  "H_influenzae_str_PTHi_8273_id_GCA_900408165_1")

for (genome in genomes) {
  
  # load data frame
  df <- read.table(paste0("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/07_GENE_LEVEL/COMBO_30_DATA/",genome,"-gene_detection.txt-COMBO-30"), sep = "\t", header = TRUE)
  
  
  # make gene_callers_id a factor
  df$gene_callers_id <- as.factor(df$gene_callers_id)
  # transpose df1 wide to long
  melted_df <- melt(df)
  
  melted_df$variable <- as.factor(melted_df$variable)
  melted_df$value<- as.numeric(melted_df$value)
  
  # # Filter to only inlcude PP metagenomes
  # PP_filtered_df <- melted_PP_df %>% 
  #   filter(grepl("PP",variable))
  
  # remove suffixes at end of variable
  melted_df$variable <- gsub(".PP", "", melted_df$variable)
  melted_df$variable <- gsub(".TD", "", melted_df$variable)
  melted_df$variable <- gsub(".BM", "", melted_df$variable)
  
  # create new variable called prop_detected (>= 0.9 == "yes"; <0.9 == "no")
  melted_df_detection <- melted_df %>% 
    mutate(detected = case_when(value >= 0.900 ~ 'yes', value < 0.900 ~ 'no'))
  
  # add QC_total_reads variable from mapping_stats.df
  mapping_stats.df$Sample_ID <- as.factor(mapping_stats.df$Sample_ID)
  melted_df_detection$QC_total_reads <- mapping_stats.df$QC_total_reads[match(melted_df_detection$variable, mapping_stats.df$Sample_ID)]
  
  
  
  # sort variable by QC_total_reads
  melted_df_detection_sorted <- melted_df_detection %>%
    mutate(sorted_variable = fct_reorder(variable, desc(QC_total_reads))) 
  
  # make site variable
  # make oral site variable for facet grid
  melted_df_detection_sorted <- melted_df_detection_sorted %>% mutate(site = case_when(grepl('TD_', variable) ~ 'TD', 
                                                             grepl('BM_', variable) ~ 'BM', 
                                                             grepl('PP_', variable) ~ 'SUPP'))
  
  
  # specify order for site terms
  melted_df_detection_sorted$site <- factor(melted_df_detection_sorted$site, levels = c("SUPP", "BM", "TD"))
  
  #Plot
  p <-ggplot(melted_df_detection_sorted, aes(x=variable, fill = detected)) + 
    geom_bar(position="fill") + 
    facet_grid(. ~ site, scales = "free", space = "free", switch = "x") +
    scale_y_continuous(breaks = seq(0, 1, .2), 
                       expand = c(0, 0)) + 
    scale_fill_manual(values = c("red","steelblue")) +
    labs(y = "Proportion of genes detected/un-detected", 
         fill = "Gene detection",
         x = NULL) +
    theme_classic() + 
    theme(text = element_text(size = 12),
          axis.ticks.y = element_line(linewidth = 1),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y =  element_text(size=12, colour = "black"),
          axis.line.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_rect(fill = "black")) + 
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          strip.text=element_text(size=12))
  
  ggsave(paste0("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/07_GENE_LEVEL/Gene_level_bar_plots_TD_BM_SUPP/",genome,".pdf"), plot = p, width = 10, height = 5)
  
  # Save data
  melted_df_detection_sorted$detected <- as.factor( melted_df_detection_sorted$detected)
  total_count <- nlevels(melted_df_detection_sorted$gene_callers_id)
  
  summary_detected <- melted_df_detection_sorted %>% 
    mutate(detected_binary = case_when(detected == "yes" ~ 1,
                                       detected == "no" ~ 0)) %>% 
    group_by(sorted_variable, site) %>% 
    plyr::count("sorted_variable","detected_binary") %>% 
    mutate(Total = total_count) %>% 
    mutate(prop_detected = freq/Total)
  
  # add QC_total_reads variable from mapping_stats.df
  summary_detected$QC_total_reads <- mapping_stats.df$QC_total_reads[match(summary_detected$sorted_variable, mapping_stats.df$Sample_ID)]
  summary_detected <- summary_detected %>% 
    arrange(desc(QC_total_reads))
  
  # add back oral site
  summary_detected <- summary_detected %>% mutate(site = case_when(grepl('TD_', sorted_variable) ~ 'TD', 
                                                                                       grepl('BM_', sorted_variable) ~ 'BM', 
                                                                                       grepl('PP_', sorted_variable) ~ 'SUPP'))
  
  write.csv(summary_detected,
            paste0("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/07_GENE_LEVEL/Gene_level_bar_plots_TD_BM_SUPP/",genome,".csv"),
            row.names = FALSE)
}

