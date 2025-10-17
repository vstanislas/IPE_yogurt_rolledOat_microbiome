library(readxl)
library(dplyr)

# contributing_species_signif_genes.tsv:
# - Indicates for each genes (from table SignificantGenes_946461.xlsx) if they were found to originate from Streptococcus thermophilus or Lactobacillus delbrueckii.
# - 1 indicates presence and 0 indicates absence.
# - The columns Streptococcus_thermophilus, Lactobacillus_delbrueckii and Streptococcus_thermophilus_CAG_236 are independent from each others, if one is "1" the others two are "0"
# - Genes are duplicated if they are find to originated from both species


###### IMPORT ######

# Import list of 3730 significant genes (SignificantGenes.xlsx)
significantGenesFile <- "ResultsFunctional/On_946461genes_fromCPM/SignificantGenes.xlsx"
all_sig_increase <- read_excel(significantGenesFile, sheet = "all_sig_increase")
all_sig_decrease <- read_excel(significantGenesFile, sheet = "all_sig_decrease")
onlyY_sig_increase <- read_excel(significantGenesFile, sheet = "onlyY_sig_increase")
onlyYO_sig_increase <- read_excel(significantGenesFile, sheet = "onlyYO_sig_increase")
onlyY_sig_decrease <- read_excel(significantGenesFile, sheet = "onlyY_sig_decrease")
onlyYO_sig_decrease <- read_excel(significantGenesFile, sheet = "onlyYO_sig_decrease")

# Import list of genes from yogurt contributing species (contributing_species_signif_genes.tsv)
genes_yogspec <- read.delim("ResultsFunctional/On_946461genes_fromCPM/contributing_species_signif_genes.tsv")
colnames(genes_yogspec)[6:8] <- c("Streptococcus_thermophilus", "Lactobacillus_delbrueckii", "Streptococcus_thermophilus_CAG_236")
# [1] "UnirefTarget": genes names corresponding to the ones listed in all_sig 
# [2] "Type": either all_sig_increase, all_sig_decrease...
# [3] "Uniref_ID": correspond to Uniref90 + one id number (1308 for Streptococcus_thermophilus, 1584 for Lactobacillus_delbrueckii and 1263110 for Streptococcus_thermophilus_CAG_236)
# [4] "Uniref90" 
# [5] "AA_Length"
# [6] "Streptococcus_thermophilus": 1 if presence, 0 if abscence
# [7] "Lactobacillus_delbrueckii": 1 if presence, 0 if abscence
# [8] "Streptococcus_thermophilus_CAG_236": 1 if presence, 0 if abscence

###### PREPARE DATA ######

### SIGNIFICANT GENES FILE
## Merge significant genes in one file
all_sig_increase$Type <- "all_sig_increase"
all_sig_decrease$Type <- "all_sig_decrease"
onlyY_sig_increase$Type <- "onlyY_sig_increase"
onlyYO_sig_increase$Type <- "onlyYO_sig_increase"
onlyY_sig_decrease$Type <- "onlyY_sig_decrease"
onlyYO_sig_decrease$Type <- "onlyYO_sig_decrease"

all_sig <- bind_rows(list(all_sig_increase, all_sig_decrease, onlyY_sig_increase, onlyYO_sig_increase, onlyY_sig_decrease, onlyYO_sig_decrease))
# dim(all_sig) # [1] 3730   10
# length(unique(all_sig$feature)) # [1] 3730 


### Yogurt species FILE
# dim(genes_yogspec) # [1] 4769    8
# sum(duplicated(genes_yogspec$UnirefTarget)) # 2078
# length(unique(genes_yogspec$UnirefTarget)) # 2691 
# sum(duplicated(genes_yogspec$Uniref_ID)) # 0
# length(unique(genes_yogspec$Uniref_ID)) # 4769 
# sum(duplicated(genes_yogspec$Uniref90)) # 1423
# length(unique(genes_yogspec$Uniref90)) # 3346 

# sum(genes_yogspec$UnirefTarget %in% all_sig$feature) # 4769
# all genes in genes_yogspec$UnirefTarget come from all_sig
# sum(all_sig$feature %in% genes_yogspec$UnirefTarget) # 2691
# 2691 of the 3730 significant genes are in genes_yogspec, (3730-2691) 1039 are missing

# separate genes_yogspec by yogurt species
genes_yogspec_strep <- filter(genes_yogspec, Streptococcus_thermophilus==1) # [1] 2841    8
genes_yogspec_lact <- filter(genes_yogspec, Lactobacillus_delbrueckii==1) # [1] 493   8
genes_yogspec_strep236 <- filter(genes_yogspec, Streptococcus_thermophilus_CAG_236==1) # [1] 1435    8
# 2841 + 493+ 1435 = 4769
sum(genes_yogspec_strep$Lactobacillus_delbrueckii) # 0
sum(genes_yogspec_strep$Streptococcus_thermophilus_CAG_236) # 0
sum(genes_yogspec_lact$Streptococcus_thermophilus) # 0
sum(genes_yogspec_lact$Streptococcus_thermophilus_CAG_236) # 0
sum(genes_yogspec_strep236$Lactobacillus_delbrueckii) # 0
sum(genes_yogspec_strep236$Streptococcus_thermophilus) # 0
# the columns Streptococcus_thermophilus, Lactobacillus_delbrueckii and Streptococcus_thermophilus_CAG_236 are independent from each others, if one is "1" the others two are "0"

dim(genes_yogspec_strep) # [1] 2841    8
length(unique(genes_yogspec_strep$UnirefTarget)) # [1] 2415
length(unique(genes_yogspec_strep$Uniref_ID)) # [1] 2841
length(unique(genes_yogspec_strep$Uniref90)) # [1] 2841
sum(all_sig$feature %in% genes_yogspec_strep$UnirefTarget) # 2415

dim(genes_yogspec_lact) # [1] 493    8
length(unique(genes_yogspec_lact$UnirefTarget)) # [1] 404
length(unique(genes_yogspec_lact$Uniref_ID)) # [1] 493
length(unique(genes_yogspec_lact$Uniref90)) # [1] 493
sum(all_sig$feature %in% genes_yogspec_lact$UnirefTarget) # 404

dim(genes_yogspec_strep236) # [1] 1435    8
length(unique(genes_yogspec_strep236$UnirefTarget)) # [1] 1429
length(unique(genes_yogspec_strep236$Uniref_ID)) # [1] 1435
length(unique(genes_yogspec_strep236$Uniref90)) # [1] 1435
sum(all_sig$feature %in% genes_yogspec_strep236$UnirefTarget) # 1429


sum(unique(genes_yogspec_lact$UnirefTarget) %in% unique(genes_yogspec_strep$UnirefTarget)) # 133
sum(unique(genes_yogspec_strep$UnirefTarget) %in% unique(genes_yogspec_lact$UnirefTarget)) # 133




#### Venn diagrams

library(ggvenn)

# With Streptococcus Thermophilus and Lactobacillus Delbrueckii only
df <- data.frame(genes= all_sig$feature)
df <- df %>% mutate(`Streptococcus \nThermophilus` = ifelse(genes %in% genes_yogspec_strep$UnirefTarget, T, F))
df <- df %>% mutate(`Lactobacillus \nDelbrueckii` = ifelse(genes %in% genes_yogspec_lact$UnirefTarget, T, F))
gg <- ggvenn(df,  show_outside="always",
       fill_color = c("#700054", "#02a0c8"),
       stroke_size=0.7,
       set_name_size =5) + update_geom_defaults("text", list(fontface = "italic"))
gg
ggsave("VennDiagram2.png", gg)



