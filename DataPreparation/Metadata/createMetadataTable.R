library(readxl)
library(xlsx)
library(dplyr)
library(tidyr)
library(ggplot2)

#########################################
# Import and filter Participant tables  #
#########################################
PartiInfo <- read_excel("Joghurt-Haferflocken-Studie_Probandeninformationen.xlsx", sheet="Probandeninformationen", range = "A1:Q120" )
PartiInfo <- PartiInfo[, c(1, 3:9, 17)]
colnames(PartiInfo) <- c("Studienphase", "Participant_ID", "Geschlecht", "Alter", 
                         "Groesse", "Gewicht1", "BMI", "Gewicht2", "Intervention1")

# Filter participants who didn't finish the study "Ausgeschieden"
table(PartiInfo$Studienphase, useNA = "a")
# Ausgeschieden Studie abgeschlossen                 <NA> 
#             9                  110                    0 
PartiInfo <- subset(PartiInfo, Studienphase=="Studie abgeschlossen")

# Create variable Gruppe
PartiInfo$Gruppe <- case_when(PartiInfo$Intervention1=="J" ~ "A", 
                              PartiInfo$Intervention1=="J+H"  ~ "B", 
                              TRUE ~ NA_character_)
table(PartiInfo$Intervention1, PartiInfo$Gruppe, useNA = "always")
#         A  B <NA>
#   J    53  0    0
# J+H     0 57    0
# <NA>    0  0    0

PartiInfo <- PartiInfo[,-1]


####### Weight change
ggplot(PartiInfo, aes(x=Gewicht1, y=Gewicht2)) + geom_point() 
#some people had an important weight variation during the study
PartiInfo$BMI_2 <- PartiInfo$Gewicht2 / (PartiInfo$Groesse^2)
PartiInfo$BMI_var <- PartiInfo$BMI_2 - PartiInfo$BMI
PartiInfo[order(abs(PartiInfo$BMI_var), decreasing=T), ] 
# Participant_ID Geschlecht Alter Groesse Gewicht1      BMI Gewicht2 Intervention1 Gruppe    BMI_2   BMI_var
# 1          ID037          w    39   1.750     79.0 25.79592     60.1           J+H      B 19.62449 -6.171429
# 2          ID094          m    54   1.778     86.0 27.20414     69.7           J+H      B 22.04800 -5.156133
# 3          ID102          m    46   1.710     88.0 30.09473     94.2             J      A 32.21504  2.120311
# 4          ID083          m    50   1.699     60.5 20.95891     66.6             J      A 23.07212  2.113212

# save(PartiInfo, file="Info_participant.Rdata")
PartiInfo <- PartiInfo[, c("Participant_ID", "Gruppe", "Geschlecht", "Alter", "BMI", "BMI_var")]



###############################
# Import and filter Stool ID  #
###############################
SamplesInfo_strow <- read_excel("Probanden-Pseudonymisierung-Stuhl.xlsx", sheet="Strowig", range = "A1:C441")
colnames(SamplesInfo_strow) <- c("Sample_ID", "Participant_ID", "Sample_order")
SamplesInfo_strow <- separate(SamplesInfo_strow, "Sample_order", c("Timepoint", "Number"))
SamplesInfo_strow$Timepoint <- paste("T", SamplesInfo_strow$Timepoint, sep="")
# table(SamplesInfo_strow$Timepoint, SamplesInfo_strow$Participant_ID, useNA = "always")
# table(SamplesInfo_strow$Number, useNA="a")

SamplesInfo_meta <- read_excel("Probanden-Pseudonymisierung-Stuhl.xlsx", sheet="Metabolon", range = "A1:C463")
colnames(SamplesInfo_meta) <- c("Sample_ID", "Participant_ID", "Sample_order")
SamplesInfo_meta <- separate(SamplesInfo_meta, "Sample_order", c("Timepoint", "Number"))
SamplesInfo_meta$Timepoint <- paste("T", SamplesInfo_meta$Timepoint, sep="")
# table(SamplesInfo_meta$Timepoint, SamplesInfo_meta$Participant_ID, useNA = "always")
# table(SamplesInfo_meta$Number, useNA="a")
#     1    2 <NA>
#   440   22    0
# Number == 2 correspond to replicates
identical(SamplesInfo_strow, SamplesInfo_meta[1:440,]) # TRUE

SamplesInfo_strow <- SamplesInfo_strow[, c("Participant_ID","Sample_ID", "Timepoint")]
# SamplesInfo_meta <- SamplesInfo_meta[, c("Participant_ID","Sample_ID", "Timepoint")]



################################
# Import and filter Boston ID  #
################################
SamplesInfo_Boston_S <- read_excel("Plasma and Serum Samples Study ID 22-207ND_Karin Michels.xlsx", sheet="inter Serum", range = "A2:D464")
colnames(SamplesInfo_Boston_S) <- c("ID_num", "Boston_S_ID", "Sample_ID", "rep")
SamplesInfo_Boston_S$Sample_ID <- gsub("JH-", "ID", SamplesInfo_Boston_S$Sample_ID)
SamplesInfo_Boston_S <- separate(SamplesInfo_Boston_S, Sample_ID, c("Participant_ID", "Timepoint"), sep="-", remove=F)
# generate 2 warning in rows [82, 126] correspond to the samples which do not have timepoints
# table(SamplesInfo_Boston_S$Timepoint, SamplesInfo_Boston_S$Participant_ID, useNA = "always")
SamplesInfo_Boston_S$Timepoint <- paste("T", SamplesInfo_Boston_S$Timepoint, sep="")
SamplesInfo_Boston_S$Sample_ID <- NULL

SamplesInfo_Boston_P <- read_excel("Plasma and Serum Samples Study ID 22-207ND_Karin Michels.xlsx", sheet="intern Plasma", range = "A2:D464")
colnames(SamplesInfo_Boston_P) <- c("ID_num", "Boston_P_ID", "Sample_ID", "rep")
SamplesInfo_Boston_P$Sample_ID <- gsub("JH-", "ID", SamplesInfo_Boston_P$Sample_ID)
SamplesInfo_Boston_P <- separate(SamplesInfo_Boston_P, Sample_ID, c("Participant_ID", "Timepoint"), sep="-", remove=F)
# generate 3 warning in rows [126, 423, 455] correspond to the samples which do not have timepoints
# table(SamplesInfo_Boston_P$Timepoint, SamplesInfo_Boston_P$Participant_ID, useNA = "always")
SamplesInfo_Boston_P$Timepoint <- paste("T", SamplesInfo_Boston_P$Timepoint, sep="")
SamplesInfo_Boston_P$Sample_ID <- NULL


#################################
# Merge with info participants  #
#################################
JoinInfo <- merge(SamplesInfo_strow, PartiInfo, by="Participant_ID", all=TRUE)
JoinInfob <- arrange(JoinInfo, Sample_ID)

# JoinInfo2 <- full_join(SamplesInfo_strow, PartiInfo, by="Participant_ID")
# JoinInfo2 <- as.data.frame(JoinInfo2)
# identical(JoinInfob, JoinInfo2)


JoinInfo_meta <- merge(SamplesInfo_meta, PartiInfo, by="Participant_ID", all=TRUE)
JoinInfo_meta <- arrange(JoinInfo_meta, Participant_ID, desc(Timepoint))

ess <- subset(JoinInfo_meta, Number =="1")[, c(1,2,3,5,6,7,8,9)]
rownames(ess) <- NULL
rownames(JoinInfo) <- NULL
identical(JoinInfo, ess) # TRUE


# JoinInfo_Boston_S <- merge(SamplesInfo_Boston_S, PartiInfo, by="Participant_ID", all=TRUE)
# JoinInfo_Boston_S <- arrange(JoinInfo_Boston_S, Boston_S_ID_num )
# 
# JoinInfo_Boston_P <- merge(SamplesInfo_Boston_P, PartiInfo, by="Participant_ID", all=TRUE)
# JoinInfo_Boston_P <- arrange(JoinInfo_Boston_P, Boston_P_ID_num )


###########
# Export  #
###########
write.xlsx(JoinInfo, "Metadata.xlsx", sheetName = "Sheet1", row.names = F)
write.xlsx(JoinInfo_meta, "Metadata_metabolome.xlsx", sheetName = "Sheet1", row.names = F)

write.xlsx(as.data.frame(SamplesInfo_Boston_S), "Metadata_Boston_S.xlsx", sheetName = "Sheet1", row.names = F)
write.xlsx(as.data.frame(SamplesInfo_Boston_P), "Metadata_Boston_P.xlsx", sheetName = "Sheet1", row.names = F)
