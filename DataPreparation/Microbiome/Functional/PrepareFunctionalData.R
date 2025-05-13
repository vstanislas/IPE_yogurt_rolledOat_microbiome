# with read.table, the # on the first line need to be remove manually and quote should be set to quote="" to allow to read quoting character (') correctly, exemple:
# coverage <- read.table("pathway_coverage_all_704_samples2.tsv", header=T, sep = "\t", quote = "")

# with read.delim we can turn off comments by setting comment.char = "", exemple:
# coverage <- read.delim("pathway_coverage_all_704_samples.tsv",  comment.char = "")

library(tidyr)

################################
#     KEGG table               #
################################
# Only yoghurt samples and JH144 already present
KEGG_JH <- read.delim("ko_merged_440samples_1cpm_10prev.tsv", comment.char = "")
colnames(KEGG_JH)[1] <- "KEGG"

dim(KEGG_JH) # 967 441   
rs <- rowSums(KEGG_JH[, 3:ncol(KEGG_JH)], na.rm =T) 

## remove empty rows (no empty rows)
which(rs==0) # integer(0)

save(KEGG_JH, file="JH_KEGG.Rdata")


################################
#     GO table               #
################################
# Only yoghurt samples and JH144 already present
GO_JH <- read.delim("go_merged_440samples_1cpm_10prev.tsv", comment.char = "")
colnames(GO_JH)[1] <- "GO"

dim(GO_JH) # 2402  441   
rs <- rowSums(GO_JH[, 3:ncol(GO_JH)], na.rm =T) 

## remove empty rows (no empty rows)
which(rs==0) # integer(0)

save(GO_JH, file="JH_GO.Rdata")


######################################
#     KEGG and GO annotation         #
######################################
library(stringr)
library(dplyr)
library(tidyr)
kegg_go <- read.delim("kegg_go_annotation.tsv")
#kegg_go0 <- kegg_go
kegg_go <- kegg_go0

length(table(kegg_go$ko_id)) # 4296
sum(kegg_go$ko_id=="unknown") # 914448/946460*100 = 97% missing

length(table(kegg_go$go_id)) # 15017
sum(kegg_go$go_id=="unknown") # 453355/946460*100 = 48% missing



### Counts the number of associated KEGG/GO id and names for each genes
# genes can be associated to several KEGG or GO terms, they are separated by ";"
var_to_deconstruct <- c("ko_id", "ko_name", "go_id", "go_name")
for(var in var_to_deconstruct){
  kegg_go$x <- str_count(kegg_go[, var], ";") + 1
  kegg_go[which(kegg_go[, var] == "unknown"), "x"] <- 0
  name_nb <- paste(var, "_number", sep="" )
  kegg_go <- rename(kegg_go, !!name_nb:="x")
}
# verify if there is the same number of id and names for each gene
kegg_go[which(kegg_go$ko_id_number != kegg_go$ko_name_number), c("ko_id", "ko_name")] #  43  2
kegg_go[which(kegg_go$go_id_number != kegg_go$go_name_number), c("go_id", "go_name")] #   7  2
# correspond to KEGG or GO with an id but no description 

table(kegg_go[, c("ko_id_number")])
#      0      1      2      3 
# 914448  31815    161     36 

table(kegg_go[, c("go_id_number")])
#      0      1      2      3      4      5      6      7      8      9     10     11     12     13 
# 453355 183543 126033  92732  46296  24541  10398   5060   3045    844    420     83     37     13 
#     14     15     16     17     18     20     21     22     23     26     33     39     43     76 
#     20     11      5      3      4      3      3      2      1      1      2      1      1      1 
#     88 
#      2
# 2 genes have 88 GO description
which(kegg_go$go_id_number==88)



### create new data frame go_data and ko_data with one row for each genes and unique GO/KEGG id
deconstruct_var <- function(name_var){
  data <- kegg_go
  
  number_item_max <- max(data[, paste(name_var, "_number", sep="" )])
  item <- seq(number_item_max)
  into_col <- as.character(item)
  name_order <- "order"
  name_val <- paste(name_var, "_value", sep="")
  
  
  sub_data <- separate(data[, c("uniref", name_var)], name_var, into_col, sep=";", remove=F)
  sub_data <- reshape2::melt(sub_data, id.vars=c("uniref", name_var),
                             measure.vars=into_col,
                             variable.name = name_order,
                             value.name = name_val)
  sub_data <- subset(sub_data, !is.na(get(name_val)))
  sub_data <- sub_data[order(sub_data$uniref),]
  
  nb_new_row <- nrow(sub_data) - nrow(data)
  print(paste("number of new rows:", nb_new_row))
  
  return(sub_data)
}
ko_id_data <- deconstruct_var("ko_id") #  "number of new rows: 233"
ko_name_data <- deconstruct_var("ko_name") #  "number of new rows: 233"
go_id_data <- deconstruct_var("go_id") #  "number of new rows: 665343"
go_name_data <- deconstruct_var("go_name") #  "number of new rows: 665343"

go_data <- merge(go_id_data, go_name_data, by=c("uniref", "order"), all=T)
ko_data <- merge(ko_id_data, ko_name_data, by=c("uniref", "order"), all=T)

length(table(go_data$go_id_value)) # 4852 (4851 GO terms + 1 "unknown" category)
length(table(go_data$go_name_value)) # 4824

length(table(ko_data$ko_id_value)) # 4280 (4279 KEGG terms + 1 "unknown" category)
length(table(ko_data$ko_name_value)) # 3844



### verify that each GO term have one unique description
all_go <- unique(go_data$go_id_value)
info_go <- c()
for(i in seq(length(all_go))){
  sub_go_i <- subset(go_data, go_id_value==all_go[i])
  nb_name_i <- length(table(sub_go_i$go_name_value, useNA ="ifany"))
  name_i <- unique(sub_go_i$go_name_value)  
  nb_genes_i <- nrow(sub_ko_i)
  
  info_go_i <- c(GO=all_go[i], nb_name=nb_name_i, name=name_i, nb_genes=nb_genes_i)
  info_go <- rbind(info_go, info_go_i)
}
rownames(info_go) <- NULL
info_go <- as.data.frame(info_go)
table(info_go$nb_name) # if 2 correspond to NA
#    1 
# 4852 
subset(info_go , name=="unknown")   # 29   4 (29 GO terms do not have a description)
info_go[duplicated(info_go$name),]  # 28   4 (no GO terms share the same description, 28 correspond to the number of repeated "unknown")


### verify that each KEGG term have one unique description
all_ko <- unique(ko_data$ko_id_value)
info_ko <- c()
for(i in seq(length(all_ko))){
  sub_ko_i <- subset(ko_data, ko_id_value==all_ko[i])
  nb_name_i <- length(table(sub_ko_i$ko_name_value, useNA ="ifany"))
  name_i <- unique(sub_ko_i$ko_name_value)  
  nb_genes_i <- nrow(sub_ko_i)
  
  info_ko_i <- c(KEGG=all_ko[i], nb_name=nb_name_i, name=name_i, nb_genes=nb_genes_i)
  info_ko <- rbind(info_ko, info_ko_i)
}
rownames(info_ko) <- NULL
info_ko <- as.data.frame(info_ko) 
table(info_ko$nb_name) # if 2 correspond to NA
#    1 
# 4280 
subset(info_ko , name=="unknown")   #  11   4 (11 KEGG do not have a description)
info_ko[duplicated(info_ko$name),] #  436   4 (some KEGGs share the same description)


save(kegg_go0, go_data, ko_data, info_ko, info_go, file="KEGG_GO_annotation.Rdata")

######################################
#     Gene annotation                #
######################################
## REDUCED RPK TABLE: genefamilies_exp86_440_relab_unstratified_10prev_0001relab.tsv
# Only yoghurt samples and JH144 already present
genes_JH <- read.delim("genefamilies_exp86_440_relab_unstratified_10prev_0001relab.tsv", comment.char = "")
colnames(genes_JH)[1] <- "Genes"

## separate genes info 
genes_JH <- separate(genes_JH, c("Genes", "Descr"), col="Genes", sep=":")
dim(genes_JH) # 946461   
rs <- rowSums(genes_JH[, 4:ncol(genes_JH)], na.rm =T) 

## remove empty rows (no empty rows)
which(rs==0) # integer(0)

save(genes_JH, file="JH_946461genes_RPK.Rdata")



## REDUCED CPM TABLE: genefamilies_exp86_440_cpm_unstratified_10prev.tsv
# Only yoghurt samples and JH144 already present
genes_JH <- read.delim("genefamilies_exp86_440_cpm_unstratified_10prev.tsv", comment.char = "")
colnames(genes_JH)[1] <- "Genes"

## separate genes info 
genes_JH <- separate_wider_delim(genes_JH,  cols="Genes", delim=": ", names=c("Genes", "Descr"), too_few = "align_start", too_many = "merge")
genes_JH <- as.data.frame(genes_JH)
dim(genes_JH) # 946461   
rs <- rowSums(genes_JH[, 3:ncol(genes_JH)], na.rm =T) 

## remove empty rows (no empty rows)
which(rs==0) # integer(0)

save(genes_JH, file="JH_946461genes_CPM.Rdata")



## REDUCED CPM TABLE (with additionaly CPM>1): genefamilies_exp86_440_cpm_unstratified_10prev_1cpm.tsv
# Only yoghurt samples and JH144 already present
genes_JH <- read.delim("genefamilies_exp86_440_cpm_unstratified_10prev_1cpm.tsv", comment.char = "")
colnames(genes_JH)[1] <- "Genes"

## separate genes info 
genes_JH <- separate_wider_delim(genes_JH,  cols="Genes", delim=": ", names=c("Genes", "Descr"), too_few = "align_start", too_many = "merge") # will keep in same column further iteration of ": "
genes_JH <- as.data.frame(genes_JH)
dim(genes_JH) # 191601   
rs <- rowSums(genes_JH[, 3:ncol(genes_JH)], na.rm =T) 

## remove empty rows (no empty rows)
which(rs==0) # integer(0)

save(genes_JH, file="JH_191601genes_CPM.Rdata")

######################################
#             Coverage               #
######################################
coverage <- read.delim("Merged_data_Joghurt_Covid/pathway_coverage_all_704_samples.tsv",  comment.char = "")
colnames(coverage)[1] <- "Pathway"

## YOGHURT STUDY
id_JH <- grepl("JH", colnames(coverage), fixed = TRUE)
id_JH[1] <- TRUE
coverage_JH <- coverage[, id_JH]
coverage_JH <- dplyr::rename(coverage_JH, JH277_Coverage=Ext_86_JH277_Coverage)

## add the missing sample JH144
coverage_JH144 <- read.delim("JH144/JH144_pathcoverage.tsv",  comment.char = "")
colnames(coverage_JH144)[1] <- "Pathway"
coverage_JH_full <- merge(coverage_JH, coverage_JH144, by="Pathway", sort=T, all=T)
dim(coverage_JH_full) - dim(coverage_JH) # 0  1 (0 new pathway, 1 new sample)

## reorder samples
n <- ncol(coverage_JH_full)
coverage_JH_full <- coverage_JH_full[, order(colnames(coverage_JH_full))]
coverage_JH_full <- coverage_JH_full[, c(n, 1:(n-1))]

# CHECK everything correct
# rownames(coverage_JH) <- coverage_JH$Pathway
# coverage_JH <- coverage_JH[coverage_JH_full$Pathway,]
# rownames(coverage_JH) <- rownames(coverage_JH_full) <- NULL
# identical(coverage_JH, coverage_JH_full[,-145])
# 
# coverage_JH144_2 <- coverage_JH_full[,c(1, 145), drop = FALSE]
# rownames(coverage_JH144_2) <- coverage_JH144_2$Pathway
# coverage_JH144_2 <- coverage_JH144_2[coverage_JH144$Pathway, , drop = FALSE]
# rownames(coverage_JH144) <- rownames(coverage_JH144_2) <-NULL
# identical(coverage_JH144, coverage_JH144_2)

## separate pathway info in categories
coverage_JH_full <- separate(coverage_JH_full, c("Pathway", "Taxonomy"), col="Pathway", sep="\\|")
coverage_JH_full <- separate(coverage_JH_full, c("Pathway", "Descr"), col="Pathway", sep=":")
coverage_JH <- coverage_JH_full

## remove empty rows
dim(coverage_JH) # 27630   
rs <- rowSums(coverage_JH[, 4:ncol(coverage_JH)], na.rm =T) 
# ess <- coverage_JH[which(rs==0), ]
# dim(ess) # 5185  
# sum(ess[, 4:ncol(ess)], na.rm=T) # 0
# sum(coverage_JH[, 4:ncol(coverage_JH)], na.rm=T) #  462035.7
coverage_JH <- coverage_JH[which(rs>0), ] 
# sum(coverage_JH[, 4:ncol(coverage_JH)], na.rm=T) # 462035.7
# dim(coverage_JH) # 22445      

save(coverage_JH, file="JH_coverage.Rdata")


#######################
#             Pathway                #
######################################
pathway <- read.delim("Merged_data_Joghurt_Covid/pathways_all_704_samples_cpm.tsv",  comment.char = "")
colnames(pathway)[1] <- "Pathway"

## YOGHURT STUDY
id_JH <- grepl("JH", colnames(pathway), fixed = TRUE)
id_JH[1] <- TRUE
pathway_JH <- pathway[, id_JH]
pathway_JH <- dplyr::rename(pathway_JH, JH277_Abundance.CPM=Ext_86_JH277_Abundance.CPM)

## add the missing sample JH144
pathway_JH144 <- read.delim("JH144/JH144_pathabundance_cpm.tsv",  comment.char = "")
colnames(pathway_JH144)[1] <- "Pathway"
pathway_JH_full <- merge(pathway_JH, pathway_JH144, by="Pathway", sort=F, all=T)
dim(pathway_JH_full) - dim(pathway_JH) # 0 1 ( 0 new pathway, 1 new sample)

## reorder samples
n <- ncol(pathway_JH_full)
pathway_JH_full <- pathway_JH_full[, order(colnames(pathway_JH_full))]
pathway_JH_full <- pathway_JH_full[, c(n, 1:(n-1))]


# CHECK everything correct
# rownames(pathway_JH) <- pathway_JH$Pathway
# pathway_JH <- pathway_JH[pathway_JH_full$Pathway,]
# rownames(pathway_JH) <- rownames(pathway_JH_full) <- NULL
# identical(pathway_JH, pathway_JH_full[,-145])
# 
# pathway_JH144_2 <- pathway_JH_full[,c(1, 145), drop = FALSE]
# rownames(pathway_JH144_2) <- pathway_JH144_2$Pathway
# pathway_JH144_2 <- pathway_JH144_2[pathway_JH144$Pathway, , drop = FALSE]
# rownames(pathway_JH144) <- rownames(pathway_JH144_2) <-NULL
# identical(pathway_JH144, pathway_JH144_2)

## separate pathway info in categories
pathway_JH_full <- separate(pathway_JH_full, c("Pathway", "Taxonomy"), col="Pathway", sep="\\|")
pathway_JH_full <- separate(pathway_JH_full, c("Pathway", "Descr"), col="Pathway", sep=":")
pathway_JH <- pathway_JH_full


## remove empty rows
dim(pathway_JH) # 27630   443
rs <- rowSums(pathway_JH[, 4:ncol(pathway_JH)], na.rm =T) 
# ess <- pathway_JH[which(rs==0), ]
# dim(ess) # 4309  443
# sum(ess[, 4:ncol(ess)], na.rm=T) # 0
# sum(pathway_JH[, 4:ncol(pathway_JH)], na.rm=T) #  734221599
pathway_JH <- pathway_JH[which(rs>0), ] 
# sum(pathway_JH[, 4:ncol(pathway_JH)], na.rm=T) # 734221599
# dim(pathway_JH) # 23321   443

save(pathway_JH, file="JH_pathway.Rdata")




######################################
#        Unpacked Pathway            #
######################################
unpck_pathway <- read.delim("Merged_data_Joghurt_Covid/unpacked_pathways_all_704_samples_cpm.tsv",  comment.char = "")
colnames(unpck_pathway)[1] <- "Pathway"

## YOGHURT STUDY
id_JH <- grepl("JH", colnames(unpck_pathway), fixed = TRUE)
id_JH[1] <- TRUE
unpck_pathway_JH <- unpck_pathway[, id_JH]
unpck_pathway_JH <- dplyr::rename(unpck_pathway_JH, JH277_Abundance.CPM=Ext_86_JH277_Abundance.CPM)

## add the missing sample JH144
unpck_pathway_JH144 <- read.delim("JH144/JH144_pathabundance_cpm_unpacked.tsv",  comment.char = "")
colnames(unpck_pathway_JH144)[1] <- "Pathway"
unpck_pathway_JH_full <- merge(unpck_pathway_JH, unpck_pathway_JH144, by="Pathway", sort=T, all=T)
dim(unpck_pathway_JH_full) - dim(unpck_pathway_JH) # 39  1 (39 new pathway, 1 new sample)

## reorder samples
n <- ncol(unpck_pathway_JH_full)
unpck_pathway_JH_full <- unpck_pathway_JH_full[, order(colnames(unpck_pathway_JH_full))]
unpck_pathway_JH_full <- unpck_pathway_JH_full[, c(n, 1:(n-1))]

# CHECK everything correct
JH144_pathways <- unpck_pathway_JH_full$Pathway[!(unpck_pathway_JH_full$Pathway %in% unpck_pathway_JH$Pathway)]
unpck_pathway_JH_full_toCheck <- subset(unpck_pathway_JH_full, !(Pathway %in% JH144_pathways))
rownames(unpck_pathway_JH) <- unpck_pathway_JH$Pathway
unpck_pathway_JH <- unpck_pathway_JH[unpck_pathway_JH_full_toCheck$Pathway,]
rownames(unpck_pathway_JH) <- rownames(unpck_pathway_JH_full_toCheck) <- NULL
identical(unpck_pathway_JH, unpck_pathway_JH_full_toCheck[,-145])

unpck_pathway_JH144_2 <- unpck_pathway_JH_full[,c(1, 145), drop = FALSE]
rownames(unpck_pathway_JH144_2) <- unpck_pathway_JH144_2$Pathway
unpck_pathway_JH144_2 <- unpck_pathway_JH144_2[unpck_pathway_JH144$Pathway, , drop = FALSE]
rownames(unpck_pathway_JH144) <- rownames(unpck_pathway_JH144_2) <-NULL
identical(unpck_pathway_JH144, unpck_pathway_JH144_2)

## separate pathway info in categories
unpck_pathway_JH_full <- separate(unpck_pathway_JH_full, c("Pathway", "Taxonomy", "Uniref"), 
                                  col="Pathway", sep="\\|")
unpck_pathway_JH_full <- separate(unpck_pathway_JH_full, c("Pathway", "Descr"), 
                                  col="Pathway", sep=":")
unpck_pathway_JH <- unpck_pathway_JH_full

## remove empty rows
dim(unpck_pathway_JH) #  1355931     444
rs <- rowSums(unpck_pathway_JH[, 5:ncol(unpck_pathway_JH)], na.rm =T) 
ess <- unpck_pathway_JH[which(rs==0), ]
dim(ess) #  130259    444
sum(ess[, 5:ncol(ess)], na.rm=T) # 0
sum(unpck_pathway_JH[, 5:ncol(unpck_pathway_JH)], na.rm=T) #  772533119
unpck_pathway_JH <- unpck_pathway_JH[which(rs>0), ] 
sum(unpck_pathway_JH[, 5:ncol(unpck_pathway_JH)], na.rm=T) # 772533119
dim(unpck_pathway_JH) # 1225672     444

save(unpck_pathway_JH , file="JH_unpck_pathway.Rdata")

