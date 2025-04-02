colorize <- function(x, color) {
  if (knitr::is_latex_output()) {
    sprintf("\\textcolor{%s}{%s}", color, x)
  } else if (knitr::is_html_output()) {
    sprintf("<span style='color: %s;'>%s</span>", color,
            x)
  } else x
}


library(dplyr)
# Evaluate direction of change for a given variable "measure" between treatment or between two time points 
computeDiff <- function(DF, x, measure){
  
  if(!("Before" %in% unique(DF$Treatment))) x = "TreatmentAfter"
  
  if(x=="TreatmentAfter"){
    DF2 <- reshape2::dcast(DF,  Participant_ID  ~ Inter_treat, value.var = measure)

    DF2b <- DF2[,c("Participant_ID", "After.Yogurt", "After.YogurtOatmeal")]
    DF2b$Diff <- DF2b$After.Yogurt - DF2b$After.YogurtOatmeal
    DF2b_melt <- reshape2::melt(DF2b, id.vars=c("Participant_ID", "Diff"))
    DF2 <- DF2b_melt
    
    DF2$Diff <- if_else(DF2$Diff>0, "Increase", if_else(DF2$Diff<0, "Decrease", "Unchange")) 
    colnames(DF2) <- c("Participant_ID", "Diff", "Inter_treat", "value")
    DF  <- merge(DF, DF2, by=c("Participant_ID", "Inter_treat"))
  }
  
  if(x=="Treatment"){
    DF2 <- reshape2::dcast(DF,  Participant_ID  ~ Inter_treat, value.var = measure)
    
    if("Yogurt" %in% unique(DF$Intervention)){
      DF2b <- DF2[,c("Participant_ID", "Before.Yogurt", "After.Yogurt")]
      DF2b$Diff <- DF2b$After.Yogurt - DF2b$Before.Yogurt
      DF2b_melt <- reshape2::melt(DF2b, id.vars=c("Participant_ID", "Diff"))
      if(length(unique(DF$Intervention)) == 1) DF2 <- DF2b_melt
    }
    
    if("YogurtOatmeal" %in% unique(DF$Intervention)){
      DF2c <- DF2[,c("Participant_ID", "Before.YogurtOatmeal", "After.YogurtOatmeal")]
      DF2c$Diff <- DF2$After.YogurtOatmeal  - DF2$Before.YogurtOatmeal 
      DF2c_melt <- reshape2::melt(DF2c, id.vars=c("Participant_ID", "Diff"))
      if(length(unique(DF$Intervention)) == 1) DF2 <- DF2c_melt
    }
    
    if(length(unique(DF$Intervention)) == 2)  DF2 <- rbind(DF2b_melt, DF2c_melt)
    
    DF2$Diff <- if_else(DF2$Diff>0, "Increase", if_else(DF2$Diff<0, "Decrease", "Unchange")) 
    colnames(DF2) <- c("Participant_ID", "Diff", "Inter_treat", "value")
    DF  <- merge(DF, DF2, by=c("Participant_ID", "Inter_treat"))
  }
  
  if(x=="Timepoint"){
    DF2 <- reshape2::dcast(DF,  Participant_ID  ~ Timepoint, value.var = measure)
    
    DF2b <- DF2[,c("Participant_ID", "T1", "T2")]
    DF2b$Diff <- DF2b$T2 - DF2b$T1
    DF2c <- DF2[,c("Participant_ID", "T2", "T3")]
    DF2c$Diff <- DF2$T3  - DF2$T2 
    DF2d <- DF2[,c("Participant_ID", "T3", "T4")]
    DF2d$Diff <- DF2$T4  - DF2$T3 
    
    DF2b_melt <- reshape2::melt(DF2b, id.vars=c("Participant_ID", "Diff"))
    DF2c_melt <- reshape2::melt(DF2c, id.vars=c("Participant_ID", "Diff"))
    DF2d_melt <- reshape2::melt(DF2d, id.vars=c("Participant_ID", "Diff"))
    DF2 <- rbind(subset(DF2b_melt, variable=="T1"), subset(DF2c_melt, variable=="T2"), DF2d_melt)
    DF2$Diff <- if_else(DF2$Diff>0, "Increase", if_else(DF2$Diff<0, "Decrease", "Unchange")) 
    colnames(DF2) <- c("Participant_ID", "Diff", "Timepoint", "value")
    DF  <- merge(DF, DF2, by=c("Participant_ID", "Timepoint"))
  }
  return(DF)
}




subset_DATA <- function(DF, measures){
  
  DF_Y = subset(DF,  Intervention=="Yogurt")
  DF_Y <- droplevels(DF_Y)
  attr(DF_Y, "NameDF") <- "Intervention == Yogurt"
  
  DF_YO = subset(DF,  Intervention=="YogurtOatmeal")
  DF_YO <- droplevels(DF_YO)
  attr(DF_YO, "NameDF") <- "Intervention == YogurtOatmeal"
  
  DF_A <- subset(DF, Treatment=="After")
  DF_A <- droplevels(DF_A)
  DF_A_full <- DF_A
  attr(DF_A, "NameDF") <- "Treatment == After"
  
  for(i in seq(length(measures))){
    
    measure_i <- measures[i]
    DF_B <- reshape2::dcast(subset(DF, Treatment=="Before"),  Participant_ID + Gruppe + Geschlecht + Alter + Treatment   ~  Timepoint, value.var = measure_i)
    DF_B$T3_T1 <- DF_B$T3 - DF_B$T1
    colnames(DF_B)[6:8] <- paste(measure_i, colnames(DF_B)[6:8], sep="_")
    # DF_A_full <- merge(DF_A_full, DF_B[, c(1, 6:8)], by="Participant_ID", all=TRUE, sort=F) # keep baseline values
    DF_A_full <- merge(DF_A_full, DF_B[, c(1, 8)], by="Participant_ID", all=TRUE, sort=F) # only baseline differences
  }
  attr(DF_A_full, "NameDF") <- "Treatment == After + baseline difference"
  
  
  return(list(DF_Y=DF_Y, DF_YO=DF_YO, DF_A=DF_A, DF_A_full=DF_A_full))
}




filter_res_mymodel <- function(res, name_levelsig=NULL){ # qval < 0.25 for maaslin
  
  # we do the comparison in the following direction:
  # YO-Y, Y-base, YO-base
  # change the sign of estimators if the contrast is not in the good direction
  if(unique(res$contrast) == "Yogurt - YogurtOatmeal"){
    res$estimate <- res$estimate * -1
    res$contrast <- "YogurtOatmeal - Yogurt"
  } 
  if(unique(res$contrast) == "Before - After"){
    res$estimate <- res$estimate * -1
    res$contrast <- "After - Before"
  } 
  cat("Contrast: ", unique(res$contrast), "\n")
  #res[which(res$q.value > 0.25), "q.value"] <- NA
  res$levelsig <- -log(res$q.value) * sign(res$estimate)
  res <- res[, c("outcomes", "levelsig", "estimate", "q.value")]
  if(!is.null(name_levelsig)){
    res <- dplyr::rename(res, !!name_levelsig := levelsig)
  }
  return(res)
}

extract_my_model <- function(file, q.value_tresh, name_levelsig){
  res <- read.xlsx(file, sheetName="PRS_all")
  res <- subset(res, q.value < q.value_tresh)
  res <- filter_res_mymodel(res, name_levelsig=name_levelsig)
  colnames(res)[-c(1:2)] <- paste(name_levelsig, c("est", "q.val"), sep=".")
  res <- dplyr::rename(res, feature=outcomes)
  return(res)
}






compare_model <- function(list_model, name, order=F, col_rot=0){
  
  
  res_1 <- list_model[[1]]
  res.qval <- res_1[, c(1,4)]
  res <- res_1[, 1:2]
  
  for(i in seq(length(list_model))[-1]){
    res_i <- list_model[[i]]
    res_qval_i <- res_i[, c(1,4)]
    res_i <- res_i[, 1:2]
    
    res <- merge(res, res_i, by="feature", all=T, sort=F)
    res.qval <- merge(res.qval, res_qval_i, by="feature", all=T, sort=F)
  }

  rownames(res) <- res$feature
  rownames(res.qval) <- res.qval$feature
  if(order) res <- res[order(res$feature),]
  if(order) res.qval <- res.qval[order(res.qval$feature),]
  
  res$feature <-NULL
  res.qval$feature <-NULL
  res.qval[is.na(res.qval)] <- 1
  
  
  cell_fun <- function(j, i, x, y, w, h, fill) {
    if(res.qval[i, j] < 0.001) {
      grid.text("***", x, y)
    } else if(res.qval[i, j] < 0.01) {
      grid.text("**", x, y)
    }else if(res.qval[i, j] < 0.05) {
      grid.text("*", x, y)
    }
  }
  
  HM <- Heatmap(as.matrix(res), name=name, cell_fun=cell_fun, col= col_fun_div, 
                cluster_rows=F, row_names_side = "left", cluster_columns=F, 
                row_names_max_width = unit(9, "cm"), column_names_rot = col_rot, 
                column_names_centered =T)
  
  return(HM)
}


 print_model <- function(method, path_analysis, model, q.value_tresh, name=model, print=T, order=F){
  
  
  path_analysis_model <- paste(path_analysis, model, sep="")
  
  if(method=="my_model"){
    # Use results from pairwise comparison in PRS_all
    res_after <- extract_my_model(paste(path_analysis_model , "Model1_After_results.xlsx", sep=""), q.value_tresh, "YO-Y")
    res_afterwBase <- extract_my_model(paste(path_analysis_model , "Model2_AfterwBase_results.xlsx", sep=""), q.value_tresh, "YO-YwB")
    res_Y <- extract_my_model(paste(path_analysis_model , "Model3_Yogurt_results.xlsx", sep=""), q.value_tresh, "Y-B")
    res_YO <- extract_my_model(paste(path_analysis_model , "Model4_YogurtOatmeal_results.xlsx", sep=""), q.value_tresh, "YO-B")
    
  } else if(method=="maaslin"){
    res_after <- extract_maaslin_res(paste(path_analysis_model, "Maaslin2 After/", sep=""), 
                               compar_var="Intervention", name_levelsig="YO-Y")
    res_Y <- extract_maaslin_res(paste(path_analysis_model, "Maaslin2 Yogurt/", sep=""), 
                               compar_var="Treatment", name_levelsig="Y-B")
    res_YO <- extract_maaslin_res(paste(path_analysis_model, "Maaslin2 YogurtOatmeal/", sep=""), 
                               compar_var="Treatment", name_levelsig="YO-B")
    
  }
  
  res_after_prt <- res_after[, c("feature", "YO-Y.est", "YO-Y.q.val")]
  res_after <- res_after[, 1:2]
  
  if(exists("res_afterwBase")){
    res_afterwBase_prt <- res_afterwBase[, c("feature", "YO-YwB.est", "YO-YwB.q.val")]
    res_afterwBase <- res_afterwBase[, 1:2]
  }
  
  res_Y_prt <- res_Y[, c("feature", "Y-B.est", "Y-B.q.val")]
  res_Y <- res_Y[, 1:2]
  
  res_YO_prt <- res_YO[, c("feature", "YO-B.est", "YO-B.q.val")]
  res_YO <- res_YO[, 1:2]
  
  
  if(exists("res_afterwBase")){
    res <- merge(res_after, res_afterwBase, by="feature", all=T, sort=F)
  } else res <- res_after
  res <- merge(res, res_Y, by="feature", all=T, sort=F)
  res <- merge(res, res_YO, by="feature", all=T, sort=F)
  #res[is.na(res)] <- 0
  rownames(res) <- res$feature
  if(order) res <- res[order(res$feature),]
  res$feature <-NULL
  
  if(exists("res_afterwBase")){
    res_prt <- merge(res_after_prt, res_afterwBase_prt, by="feature", all=T, sort=F)
  } else res_prt <- res_after_prt
  res_prt <- merge(res_prt, res_Y_prt, by="feature", all=T, sort=F)
  res_prt <- merge(res_prt, res_YO_prt, by="feature", all=T, sort=F)
  res_prt$model <- model
  rownames(res_prt) <- res_prt$feature
  if(order) res_prt <- res_prt[order(res_prt$feature),]
  res_prt$feature <-NULL
  if(print) print(res_prt)
  
  if(exists("res_afterwBase")){
    res.qval <- res_prt[, c("YO-Y.q.val", "YO-YwB.q.val", "Y-B.q.val", "YO-B.q.val")]
  } else  res.qval <- res_prt[, c("YO-Y.q.val", "Y-B.q.val", "YO-B.q.val")]
 
  res.qval[is.na(res.qval)] <- 1
  cell_fun <- function(j, i, x, y, w, h, fill) {
    if(res.qval[i, j] < 0.001) {
      grid.text("***", x, y)
    } else if(res.qval[i, j] < 0.01) {
      grid.text("**", x, y)
    }else if(res.qval[i, j] < 0.05) {
      grid.text("*", x, y)
    }
  }
  
  HM <- Heatmap(as.matrix(res), name=name, cell_fun=cell_fun, col= col_fun_div, cluster_rows=F, row_names_side = "left", cluster_columns=F, row_names_max_width = unit(9, "cm"), column_names_rot = 0, column_names_centered =T)
  
  return(HM)
}


filter_res_maaslin <- function(sig_res, compar_var, name_levelsig=NULL){ # qval < 0.25 for maaslin
  sig_res_inter <- subset(sig_res,  metadata==compar_var)[,c(1,2,4,9)]
  sig_res_inter$levelsig <- -log(sig_res_inter$qval) * sign(sig_res_inter$coef)
  if(!is.null(name_levelsig)){
    sig_res_inter <- dplyr::rename(sig_res_inter, !!name_levelsig := levelsig)
  }
  # print(sig_res_inter)
  sig_res_inter <- sig_res_inter[, c(1,5,3,4)]
  return(sig_res_inter)
}


extract_maaslin_res <- function(path, compar_var, name_levelsig){
  res <- read.table(file=paste(path, "/significant_results.tsv", sep="") , sep = '\t', header = TRUE) 
  res <- filter_res_maaslin(res, compar_var=compar_var, name_levelsig=name_levelsig)
  colnames(res)[-c(1:2)] <- paste(name_levelsig, c("est", "q.val"), sep=".")
  res$feature <- gsub("Species.", "", res$feature)
  return(res)
}




extract_Rosner_fitStat <- function(path_models, models_df, compa){
  
  res <- c()
  for(i in seq(nrow(models_df))){
   model_i <- models_df[i, "models"]
   file <- paste(path_models, model_i, "/", compa, "/", model_i, "_", compa, "_results.xlsx", sep="")
   res_i <- read.xlsx(file=file, sheetName="fit_stat")
   res_i$names_model <- models_df[i, "names_models"]
   res_i$data_type <- models_df[i, "data_type"]
   if(ncol(res_i) > 7) res_i <- res_i[, -1]
   colnames(res_i)[1] <- "feature"
   res <- rbind(res, res_i)
  }

  return(res)
}

