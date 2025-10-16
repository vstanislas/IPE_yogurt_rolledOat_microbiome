library(readxl)
library(dplyr)


###### IMPORT ######
# path <- "~/Work/IPE-Freiburg/Manuscript Yogurt/DataAnalysis/Edited figures and additional results 09-24/"
significantGenesFile <- "Q:/IPE-P-Joghurtstudie/Joghurt und Haferflocken/Analysis/DataAnalysis/Results_Complete_Analysis/ResultsFunctional/On_946461genes_fromCPM/ResultsBetadiv/NBZIMM_asinsqrt_prevalence_10/SignificantGenes.xlsx"

# Import list of 3730 significant genes
all_sig_increase <- read_excel(significantGenesFile, sheet = "all_sig_increase")
all_sig_decrease <- read_excel(significantGenesFile, sheet = "all_sig_decrease")
onlyY_sig_increase <- read_excel(significantGenesFile, sheet = "onlyY_sig_increase")
onlyYO_sig_increase <- read_excel(significantGenesFile, sheet = "onlyYO_sig_increase")
onlyY_sig_decrease <- read_excel(significantGenesFile, sheet = "onlyY_sig_decrease")
onlyYO_sig_decrease <- read_excel(significantGenesFile, sheet = "onlyYO_sig_decrease")

# Import Kun file (list of genes from yogurt contributing species)
genes_yogspec <- read.delim("contributing_species_signif_genes.tsv")
colnames(genes_yogspec)[6:8] <- c("Streptococcus_thermophilus", "Lactobacillus_delbrueckii", "Streptococcus_thermophilus_CAG_236")
# [1] "UnirefTarget": genes names corresponding to the ones listed in all_sig (file send to Kun)
# [2] "Type": either all_sig_increase, all_sig_decrease...
# [3] "Uniref_ID": correspond to Uniref90 + one id number (probably one for each species? 1308 for Streptococcus_thermophilus, 1584 for Lactobacillus_delbrueckii and 1263110 for Streptococcus_thermophilus_CAG_236)
# [4] "Uniref90": 
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


### KUN FILE
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

###### ANALYSIS ######
all_sig %>% group_by(Type) %>% summarise(n = n()) # number of genes by "Types"
# Type                      n
# <chr>                 <int>
# 1 all_sig_decrease        3
# 2 all_sig_increase     3017
# 3 onlyYO_sig_decrease    99
# 4 onlyYO_sig_increase   128
# 5 onlyY_sig_decrease    136
# 6 onlyY_sig_increase    347

genes_yogspec %>% group_by(Type) %>% summarise(n = n()) 
# Type                      n
# <chr>                 <int>
# 1 all_sig_increase     4375
# 2 onlyYO_sig_increase    41
# 3 onlyY_sig_decrease      2
# 4 onlyY_sig_increase    351
# all_sig_decrease, onlyYO_sig_decrease are not in the file from Kun


genes_yogspec_strep %>% group_by(Type) %>% summarise(n = n())
# Type                    n
# <chr>               <int>
# 1 all_sig_increase     2782
# 2 onlyYO_sig_increase    23
# 3 onlyY_sig_decrease      2
# 4 onlyY_sig_increase     34

genes_yogspec_lact %>% group_by(Type) %>% summarise(n = n())
# Type                    n
# <chr>               <int>
# 1 all_sig_increase      186
# 2 onlyYO_sig_increase     5
# 3 onlyY_sig_increase    302

genes_yogspec_strep236 %>% group_by(Type) %>% summarise(n = n())
# Type                    n
# <chr>               <int>
# 1 all_sig_increase     1407
# 2 onlyYO_sig_increase    13
# 3 onlyY_sig_increase     15


#### Venn diagrams
library(ggvenn)
x = list(
  Strept=genes_yogspec_strep$UnirefTarget ,
  Lact=genes_yogspec_lact$UnirefTarget ,
  Strept_236=genes_yogspec_strep236$UnirefTarget
)
ggvenn(x,  c("Strept" , "Lact" , "Strept_236"),show_percentage = F,show_outside="always")  


# With Streptococcus Thermophilus, Lactobacillus Delbrueckii and Streptococcus Thermophilus CAG 236
df <- data.frame(genes= all_sig$feature)
df <- df %>% mutate(`Streptococcus \nThermophilus` = ifelse(genes %in% genes_yogspec_strep$UnirefTarget, T, F))
df <- df %>% mutate(`Lactobacillus \nDelbrueckii` = ifelse(genes %in% genes_yogspec_lact$UnirefTarget, T, F))
df <- df %>% mutate(`Streptococcus \nThermophilus CAG 236` = ifelse(genes %in% genes_yogspec_strep236$UnirefTarget, T, F))
ggvenn(df,  show_outside="always",
       fill_color = c("#700054", "#02a0c8", "#fde725ff"),
       stroke_size=0.7,
       set_name_size =5)
       #padding=0.5)  
#+  theme(plot.title = element_text(face = "italic"))
# save manually with width=1400 height=700


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



# c("#700054", '#02a0c8', '#fde725ff')
# c("#440154ff", '#21908dff', '#fde725ff')
# c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF")



