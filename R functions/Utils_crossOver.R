
plot_sample_measure_CO <- function(obj, measure, subtitle="", x= "Treatment", colorInter=F, trans =NULL, 
                                   print_table=T, lineMean=F, printMean=F, ind_filter=NULL, line_diff=T, 
                                   linewidth_diff=1.4, ylim=NULL, counts=F, color_lines=NULL, ...){
  
  if(class(obj) == "TreeSummarizedExperiment"){
    DF <- as.data.frame(colData(obj))
  }  
  else DF <- obj
  
  if(measure == "abundance") {
    DF <-  data.table::setnames(DF, "abundance", subtitle)
    measure <- subtitle
    subtitle <- NULL
  }
  
  pp <-   plot_sample_measure(DF, measure, subtitle, x, colorInter, trans, 
                              lineMean, printMean, ind_filter, line_diff, 
                              linewidth_diff, ylim, counts, color_lines, ...)
  
  
  if(x=="Treatment" & "Intervention" %in% colnames(DF)){
    pp <-  pp  + facet_grid( ~ Intervention, scales = "free_x", space = "free_x") 
  } 
  if(x=="Timepoint" & "Gruppe" %in% colnames(DF)){
    pp <-  pp  + facet_grid( ~ Gruppe ) 
  } 
  
  if(print_table){
    if(line_diff) DF <- computeDiff(DF, x, measure)
    pp <- pp + tableGrob(table(DF$Diff, DF$Intervention)/2, theme=ttheme_default(base_size =7)) + 
      plot_layout(widths = c(3, 1))
  }  
  
  return(pp)
}


library(ggedit)
plot_outliers_path <- function(sub, ...){
  # to_label <- group_by(sub, Participant_ID) %>% summarise(rep = n())
  # to_label$id <- sapply(to_label$rep, function(x) sample.int(x, 1))
  sub_label <- group_by(sub, Participant_ID) %>% slice_sample(n=1) %>% ungroup()
  
  p <-  plot_sample_measure_CO(...) 
  p <- remove_geom(p, 'line') # library(ggedit)
  p <- p + geom_line(data= sub, aes(group=Participant_ID, color= Participant_ID),  linewidth=1.2, alpha=0.7)+
    guides(colour="legend") + 
    geom_label(data= sub_label, aes(label=Participant_ID, color= Participant_ID, fontface=2), show.legend =F, label.size=0.2)
  return(p)
}



prepare_DATA <- function(DF, trans="NULL", name_group_var=NULL){
  # inters <- as.character(unique(DF$Intervention))
  # inter1 <- inters[2]
  # inter2 <- inters[1]
  inter1 <- as.character(unique(subset(DF, Gruppe=="A" & Timepoint=="T1")$Intervention))
  inter2 <- as.character(unique(subset(DF, Gruppe=="B" & Timepoint=="T1")$Intervention))
  
  cat(crayon::bold("Create new variables needed for regression model:\n"))
  
  cat(paste("inter1 correspond to : ", inter1, "\n", sep=""))
  cat(paste("inter2 correspond to : ", inter2, "\n", sep=""))
  
  
  ### Create new variables
  DF$B <-  case_when(DF$Treatment %in% c("Before", "Baseline") ~ 1,
                     DF$Treatment == "After" ~ 0,
                     TRUE ~ NA_real_)
  
  DF$P <-  case_when(DF$Timepoint %in% c("T1", "T2") ~ 0,
                     DF$Timepoint %in% c("T3", "T4") ~ 1,
                     TRUE ~ NA_real_)
  
  DF$T <-  case_when(DF$Inter_treat == paste("After.", inter2, sep="") ~ 1,
                     DF$Inter_treat != paste("After.", inter2, sep="") ~ 0,
                     TRUE ~ NA_real_)
  
  DF$C1 <- case_when(DF$Timepoint %in% c("T3", "T4") & DF$Gruppe =="A" ~ 1,
                     TRUE ~ 0)
  
  DF$C2 <- case_when(DF$Timepoint %in% c("T3", "T4") & DF$Gruppe =="B" ~ 1,
                     TRUE ~ 0)  
  
  # DF$B <- as.factor(DF$B)
  # DF$P <- as.factor(DF$P)
  # DF$T <- as.factor(DF$T)
  # DF$C1 <- as.factor(DF$C1)
  # DF$C2 <- as.factor(DF$C2)
  
  DF$B <- factor(DF$B, levels = c("0", "1"), labels=c("After", "Baseline"))
  DF$B <- relevel(DF$B, ref="Baseline")
  DF$P <- factor(DF$P, levels = c("0", "1"), labels=c("Period1", "Period2"))
  DF$T <- factor(DF$T, levels = c("0", "1"), labels=c("Other", paste("After_", inter2, sep="")))
  DF$C1 <- factor(DF$C1, levels = c("0", "1"), labels=c("Other", paste("CO_", inter1, sep="")))
  DF$C2 <- factor(DF$C2, levels = c("0", "1"), labels=c("Other", paste("CO_", inter2, sep="")))
  
  cat(paste("B = 1 for baseline samples (Treatment == \"Before\" or  \"Baseline\")\n", sep=""))
  cat(paste("P = 1 for 2nd Period (Timepoint %in% c(\"T3\", \"T4\"))\n", sep=""))
  cat(paste("T = 1 for ",  paste("After.", inter2, sep=""), " samples\n", sep=""))
  cat(paste("C1 = 1 carry over for Gruppe ==\"A\"\n", sep=""))
  cat(paste("C2 = 1 carry over for Gruppe ==\"B\"\n", sep=""))
  
  
  
  ### Transform dependent variables
  cat(crayon::bold("\nTransformation of the dependent variables:\n"))
  if(!is.null(trans)){
    DF <- transform_data(DF, trans)
  } 
  
  if(is.null(trans)){
    cat(blue(paste("No transformation realized.\n")))
  }
  
  if(!is.null(name_group_var)){
    cat(crayon::bold(paste("\nCreate a new dichotomical variable for each group defined in ", name_group_var, ":\n", sep="")))
    DF <- add_group_variable(DF, name_group_var)
  }
  
  DF_int1 = subset(DF,  Intervention== inter1)
  DF_int1 <- droplevels(DF_int1)  
  
  DF_int2 = subset(DF,  Intervention== inter2)
  DF_int2 <- droplevels(DF_int2)
  
  rownames(DF) <- DF$Sample_ID
  rownames(DF_int1) <- DF_int1$Sample_ID
  rownames(DF_int2) <- DF_int2$Sample_ID
  
  return(list(DF=DF, DF_int1=DF_int1, DF_int2=DF_int2))
}


add_group_variable <- function(DF, name_group_var){
  
  groups <- unique(DF[, name_group_var])
  
  for(i in groups){
    
    DF$Gi <-  case_when(DF[,name_group_var] == i ~ 1,
                        DF[,name_group_var] != i  ~ 0,
                        TRUE ~ NA_real_)
    name_var <- paste("G", i, sep="_")
    DF <- dplyr::rename(DF, !!name_var:=Gi)
    
    cat(paste(name_var, " = 1 for participants in group", i, "\n",sep=""))
  }
  
  return(DF=DF)
}


library(progress) # progress()
# library(xlsx) # write.xlsx() # change to openxlsx because of java errors with large dataset
library(writexl)
library(openxlsx)
library(crayon)
# First version of the code
combined_results <- function(outcomes, groupvars, DF, compar_var, baseline=F, 
                             verbose=T, savingName=NULL, logTrans=F, id_totrans=NULL,
                             plot_res=F, form=NULL, sig_threshold=0.05, path=NULL){
  cat('\n')
  if(!is.null(attr(DF, "NameDF"))) print(paste("DATA subset: ", attr(DF, "NameDF")))
  
  DF0 <- DF
  if(logTrans){
    cat(paste("Log transform",  length(id_totrans), "variables from", 
              names(DF)[id_totrans[1]], "to", names(DF)[tail(id_totrans, 1)]), "\n")
    DF[, id_totrans] <- log_transform(DF[, id_totrans])
  } 
  
  
  res_anov_all <- c()
  emm_all <- c()
  PRS_all <- c()
  res_all <- c()
  id_bas_var <- which(groupvars=="baseline_diff")
  
  nbOut <- length(outcomes)
  pb <- progress_bar$new(format = "[:bar] :current/:total (:percent)", total = nbOut)
  
  
  for(i in seq(nbOut)){
    
    # svMisc::progress(i, max.value=nbOut, progress.bar=T)
    # pb$tick()
    
    outcome_i <- outcomes[i]
    if(baseline==T){
      groupvars[id_bas_var] <- paste(outcomes[i], "_T3_T1", sep="")
    }
    
    
    MM <- result_mixedmod(outcome_i, groupvars=groupvars, DF, compar_var, verbose = verbose, plot=plot_res, form)
    
    res_all <- rbind(res_all, MM$para)
    
    res_anov_i <- round(as.data.frame(MM$ANOVA), 5)
    res_anov_i$names <- rownames(res_anov_i)
    rownames(res_anov_i) <- NULL
    res_anov_i$outcomes <- outcomes[i]
    res_anov_i <- res_anov_i[, c(5,4, 1:3)]
    res_anov_all <- rbind(res_anov_all, res_anov_i)
    
    emm_i <-  as.data.frame(MM$emm)
    emm_i[, 2:6] <- round(emm_i[, 2:6], 2)
    emm_i$outcomes <- outcomes[i]
    emm_i <- emm_i[, c(7, 1:6)]
    emm_all <- rbind(emm_all, emm_i)
    
    PRS_i <-  as.data.frame(MM$PRS)
    PRS_i[, 2:6] <-  round(PRS_i[, 2:6] , 2)
    PRS_i[, 7:8] <-  round(PRS_i[, 7:8] , 4)
    PRS_i$outcomes <- outcomes[i]
    PRS_i$prev <-  round(sum(DF0[,outcome_i] > 0)/nrow(DF),2)
    PRS_i$abun <-  round(mean(DF0[,outcome_i]), 4) 
    PRS_i <- PRS_i[, c(9, 1:6, 10, 11, 7, 8)]#use data before log transformation
    PRS_all <- rbind(PRS_all, PRS_i)
    
    if(plot_res){
      cat(green(paste("Saving residual plots ", path, "/Results_MM/Plots_", savingName, "/Resid_plots \n\n", sep="")))
      ggsave(paste(path, "Results_MM/Plots_", savingName,"/Resid_plots/Resid_panel_",
                   outcome_i, ".png", sep=""),
             MM$res_plot, width=9, height=7)
    }
  }
  
  res_all <- dplyr::rename(res_all, p.value=`Pr(>|t|)`)
  res_anov_all <- dplyr::rename(res_anov_all, p.value=`Pr(>Chisq)`)
  
  res_all$q.value <- round(p.adjust(res_all$p.value, method="fdr"), 4)
  res_anov_all$q.value <- round(p.adjust(res_anov_all$p.value, method="fdr"), 4)
  PRS_all$q.value <- round(p.adjust(PRS_all$p.value, method="fdr"), 4)
  
  res_all <- res_all[order(res_all$q.value ),]
  res_anov_all <- res_anov_all[order(res_anov_all$q.value ),]
  PRS_all <- PRS_all[order(PRS_all$q.value ),]
  
  to_return <- list(res_all=res_all, res_anov_all=res_anov_all, emm_all=emm_all, PRS_all=PRS_all)
  
  
  if(!is.null(savingName)){
    cat(green(paste("Saving results files in ...", path, "/Results_MM/", savingName, "\n", sep="")))
    cat(green("Saving results models as excel files \n"))
    write.xlsx(res_all, file=paste(path, "Results_MM/", savingName, "_results.xlsx", sep=""), sheet="res_all")
    write.xlsx(res_anov_all, file=paste(path, "Results_MM/", savingName, "_results.xlsx", sep=""), sheet="res_anov_all", append = T)
    write.xlsx(emm_all, file=paste(path, "Results_MM/", savingName, "_results.xlsx", sep=""), sheet="emm_all", append = T)
    write.xlsx(PRS_all, file=paste(path, "Results_MM/", savingName, "_results.xlsx", sep=""), sheet="PRS_all", append = T)
    
    if(nrow(res_anov_all) <= 50){
      cat(green("Saving results models as one png \n"))
      pp <- print_result_MM(to_return, sig_threshold)
      ggsave(paste(path, "Results_MM/", savingName, ".png", sep=""), 
             wrap_elements(pp$tg1) + pp$tg2 / pp$tg3 + plot_layout(widths = c(1, 2)),
             height=12, width=19)
    }
    
    sig_res_anov <- subset(res_anov_all, p.value < sig_threshold) # [, c("outcomes", "names", "p.value", "q.value")]
    sig_res_PRS <- subset(PRS_all, p.value < sig_threshold) # [, c("outcomes", "p.value", "q.value")]
    if(nrow(sig_res_anov)>0) {
      pp_anov <- plots_varsignificant(sig_res_anov, DF, savingName="res_anov", 
                                      path=paste(path, "Results_MM/Plots_", savingName, "/", sep=""))
      ggsave_adaptSize(paste(path, "Results_MM/Plots_", savingName, "/sigTaxa_anov.png", sep=""), pp_anov)
    }
    if(nrow(sig_res_PRS)>0) {
      sig_res_PRS$names <- compar_var
      pp_PRS <- plots_varsignificant(sig_res_PRS, DF, savingName="res_PRS", 
                                     path=paste(path, "Results_MM/Plots_", savingName, "/", sep=""))
      ggsave_adaptSize(paste(path, "Results_MM/Plots_", savingName, "/sigTaxa_PRS.png", sep=""), pp_PRS)
      
    }
  }
  
  return(to_return)
}



library(car)
library(ggrepel)
# Adaptation of combined_results to newest model version
# logTrans, if TRUE, will transform the id_totrans variables in log base 2 
run_model <- function(outcomes, groupvars, DF, compar_var, var_check=NULL,
                      verbose=T, savingName=NULL, logTrans=F, id_totrans=NULL,
                      plot_res=F, sig_threshold=0.05, path=NULL, stopCode=NULL){
  
  
  if(is.null(var_check)) var_check <- groupvars[!(groupvars %in% grep("[|]", groupvars, value = TRUE))]
  DF <- droplevels(DF)
  DF0 <- DF
  
  
  res_anov_all <- c()
  emm_all <- c()
  PRS_all <- c()
  res_all <- c()
  inf_obs_all <- c()
  inf_part_all <- c()
  tzi_all <- c()
  tdisp_all <- c()
  
  nbOut <- length(outcomes)
  
  if(is.null(stopCode)) stopCode <- rep(F, nbOut)
  
  for(i in seq(nbOut)){
    DF <- DF0
    outcome_i <- outcomes[i]
    cat(crayon::bold(green(paste(outcome_i, "\n", sep=""))))
    
    id_NA <- unique(which(is.na(DF[, c(outcome_i, var_check)]), arr.ind = T)[,"row"])
    if(length(id_NA) > 0){
      if(verbose) cat(green(paste("NA detected in outcome, groupvars or var_check, ", length(id_NA), " rows removed: \n", sep="")))
      if(verbose) print(DF[id_NA, c("Participant_ID", "Timepoint", outcome_i, var_check)])
      DF =  DF[-id_NA,]
    } 
    
    MM <- result_mixedmod(outcome_i, groupvars=groupvars, DF, compar_var, verbose = verbose, plot=plot_res, stopCode=stopCode[i])
    
    res_all <- rbind(res_all, MM$para)
    
    res_anov_i <- round(as.data.frame(MM$ANOVA), 5)
    res_anov_i$names <- rownames(res_anov_i)
    rownames(res_anov_i) <- NULL
    res_anov_i$outcomes <- outcomes[i]
    res_anov_i <- res_anov_i[, c(5,4, 1:3)]
    res_anov_all <- rbind(res_anov_all, res_anov_i)
    
    emm_i <-  as.data.frame(MM$emm)
    emm_i[, 2:6] <- round(emm_i[, 2:6], 2)
    emm_i$outcomes <- outcomes[i]
    emm_i <- emm_i[, c(7, 1:6)]
    emm_all <- rbind(emm_all, emm_i)
    
    PRS_i <-  as.data.frame(MM$PRS)
    digits = 1
    if(abs(PRS_i$estimate) < 100) digits = 2
    if(abs(PRS_i$estimate) < 1) digits = 3
    if(abs(PRS_i$estimate) < 0.1) digits = 4
    if(abs(PRS_i$estimate) < 0.01) digits = 5
    if(abs(PRS_i$estimate) < 0.001) digits = 6
    if(abs(PRS_i$estimate) < 0.00001) digits = 7
    if(abs(PRS_i$estimate) < 0.0000001) digits = 8
    PRS_i[, 2:6] <-  round(PRS_i[, 2:6] , digits)
    PRS_i[, 7:8] <-  round(PRS_i[, 7:8] , 4)
    PRS_i$outcomes <- outcomes[i]
    PRS_i$prev <-  round(sum(DF[,outcome_i] > 0)/nrow(DF0),2)
    PRS_i$abun <-  round(mean(DF[,outcome_i]), 4) 
    PRS_i <- PRS_i[, c(9, 1:6, 10, 11, 7, 8)]#use data before log transformation
    PRS_all <- rbind(PRS_all, PRS_i)
    
    tzi <- testZeroInflation(MM$simRes, plot=F)
    tzi <- data.frame(tzi$statistic, tzi$p.value, outcome_i)
    tzi_all <- rbind(tzi_all, tzi)
    
    tdisp <- testDispersion(MM$simRes_cond, plot=F)
    tdisp <- data.frame(tdisp$statistic, tdisp$p.value, outcome_i) 
    tdisp_all <- rbind(tdisp_all, tdisp)
    
    
    if(nrow(MM$inf_obs)>0){
      #inf_obs <- cbind(MM$inf_obs, DF[as.numeric(MM$inf_obs$id), c("Participant_ID", "Timepoint", outcome_i)], outcome_i)
      inf_obs <- cbind(MM$inf_obs[, c("Participant_ID", "Timepoint", "cooksd")], outcome_i)
      inf_obs_all <- rbind(inf_obs_all, inf_obs)
    }
    
    if(nrow(MM$inf_part)>0){
      inf_part <- cbind(MM$inf_part[, c("Participant_ID", "cooksd")], outcome_i)
      inf_part_all <- rbind(inf_part_all, inf_part)
    }
    
    if(plot_res){
      
      cat(green(paste("Saving residual plots ", path, "/Plots_", savingName, "/Resid_plots \n\n", sep="")))
      
      # Residual from ggResidpanel: conditional pearson residual from lmer
      ggsave(paste(path, "Plots_", savingName,"/Resid_plots/", outcome_i, "_Resid_panel.png", sep=""),
             MM$res_plot, width=12, height=7)
      
      # Random effects distribution
      DF_RE <- merge(DF, MM$RE, by="Participant_ID")
      ll_plot <- list()
      ll_plot[[1]] <- MM$ranef
      for(i in seq(length(var_check))){
        if(class(DF_RE[,var_check[i]])=="factor"){
          ll_plot[[i+1]] <- ggplot(DF_RE, aes(x=.data[[var_check[i]]], y=RE_int)) + geom_boxplot()+ 
            geom_text_repel(data=DF_RE[is_outlier(DF_RE$RE_int),], aes(label=Sample_ID))
        } 
        else ll_plot[[i+1]] <- ggplot(DF_RE, aes(x=.data[[var_check[i]]], y=RE_int)) + geom_point() + 
            geom_text_repel(data=DF_RE[is_outlier(DF_RE$RE_int),], aes(label=Sample_ID))
      }
      ggsave_adaptSize(paste(path, "Plots_", savingName,"/Resid_plots/", outcome_i, "_Ranef.png", sep=""),
                       ll_plot, 4, 3)
      
      # DHARMa Residuals 
      # id_NA <- unique(which(is.na(DF[, c(outcome_i, var_check)]), arr.ind = T)[,"row"])
      # #if(length(id_NA) > 0) cat(green(paste("NA detected, ", length(id_NA), " rows removed", sep="")))
      # DF_NA =  DF[-id_NA,]
      
      png(paste(path, "Plots_", savingName,"/Resid_plots/", outcome_i, "_DHARMa_main.png", sep=""),
          width=7, height=5, units="in", res=300)
      try(plot(MM$simRes), silent=T)
      dev.off()
      
      if(length(var_check)<=4){
        width=7 
        height=7
      } else {
        width=10 
        height=10
      }
      png(paste(path, "Plots_", savingName,"/Resid_plots/", outcome_i, "_DHARMa_var.png", sep=""),
          width=width, height=height, units="in", res=300)
      par(mfrow = n2mfrow(length(var_check)))
      for(i in seq(length(var_check))){
        plotResiduals(MM$simRes, form = DF[ ,var_check[i]], xlab=var_check[i])
        mtext(var_check[i], 1, line=2)
      }
      dev.off()
      
      # Predicted effects
      png(paste(path, "Plots_", savingName,"/Resid_plots/", outcome_i, "_pred_effect.png", sep=""),
          width=8, height=8, units="in", res=300)
      plot(MM$pred_effects)
      dev.off()
      
      # # Influential points from car package
      # png(paste(path, "Plots_", savingName,"/Resid_plots/", outcome_i, "_inf_points.png", sep=""),
      #     width=7, height=7, units="in", res=300)
      # inf_pts <- influencePlot(MM$model)
      # dev.off()
      # inf_pts <- cbind(inf_pts, DF[rownames(inf_pts), c("Participant_ID", "Timepoint")], outcome_i)
      # inf_pts_all <- rbind(inf_pts_all, inf_pts)
      
      # Influential points from HLMdiag package
      
      if(nrow(MM$inf_obs)>0 & nrow(MM$inf_part)>0){
        ggsave(paste(path, "Plots_", savingName,"/Resid_plots/", outcome_i, "_inf_points.png", sep=""),
               wrap_plots(MM$inf_obs_plot, MM$inf_part_plot), width=9, height=7)
      } else if(nrow(MM$inf_obs)>0){
        ggsave(paste(path, "Plots_", savingName,"/Resid_plots/", outcome_i, "_inf_points.png", sep=""),
               MM$inf_obs_plot, width=4.5, height=7)
      }else if(nrow(MM$inf_part)>0){
        ggsave(paste(path, "Plots_", savingName,"/Resid_plots/", outcome_i, "_inf_points.png", sep=""),
               MM$inf_part_plot, width=4.5, height=7)
      }
      
    }
  }
  
  DF <- DF0
  
  res_all <- dplyr::rename(res_all, p.value=`Pr(>|t|)`)
  res_anov_all <- dplyr::rename(res_anov_all, p.value=`Pr(>Chisq)`)
  
  res_all$q.value <- round(p.adjust(res_all$p.value, method="fdr"), 4)
  res_anov_all$q.value <- round(p.adjust(res_anov_all$p.value, method="fdr"), 4)
  PRS_all$q.value <- round(p.adjust(PRS_all$p.value, method="fdr"), 4)
  
  res_all <- res_all[order(res_all$q.value ),]
  res_anov_all <- res_anov_all[order(res_anov_all$q.value ),]
  PRS_all <- PRS_all[order(PRS_all$q.value ),]
  
  to_return <- list(res_all=res_all, res_anov_all=res_anov_all, 
                    emm_all=emm_all, PRS_all=PRS_all, tzi_all=tzi_all,
                    tdisp_all=tdisp_all, inf_obs_all=inf_obs_all, inf_part_all=inf_part_all)
  
  
  if(!is.null(savingName)){
    cat(green(paste("Saving results files in ...", path, "/", savingName, "\n", sep="")))
    cat(green("Saving results models as excel files \n"))
    write.xlsx(res_all, file=paste(path, savingName, "_results.xlsx", sep=""), sheet="res_all")
    write.xlsx(res_anov_all, file=paste(path, savingName, "_results.xlsx", sep=""), sheet="res_anov_all", append = T)
    write.xlsx(emm_all, file=paste(path, savingName, "_results.xlsx", sep=""), sheet="emm_all", append = T)
    write.xlsx(PRS_all, file=paste(path, savingName, "_results.xlsx", sep=""), sheet="PRS_all", append = T)
    write.xlsx(tzi_all, file=paste(path, savingName, "_results.xlsx", sep=""), sheet="tzi_all", append = T)
    write.xlsx(tdisp_all, file=paste(path, savingName, "_results.xlsx", sep=""), sheet="tdisp_all", append = T)
    if(!is.null(inf_obs_all)) write.xlsx(inf_obs_all, file=paste(path, savingName, "_results.xlsx", sep=""), sheet="inf_obs_all", append = T)
    if(!is.null(inf_part_all)) write.xlsx(inf_part_all, file=paste(path, savingName, "_results.xlsx", sep=""), sheet="inf_part_all", append = T)
    
    if(nrow(res_anov_all) <= 50){
      cat(green("Saving results models as one png \n"))
      pp <- print_result_MM(to_return, sig_threshold)
      ggsave(paste(path, savingName, ".png", sep=""), 
             wrap_elements(pp$tg1) + pp$tg2 / pp$tg3 + plot_layout(widths = c(1, 2)),
             height=12, width=19)
    }
    
    sig_res_anov <- subset(res_anov_all, p.value < sig_threshold) # [, c("outcomes", "names", "p.value", "q.value")]
    sig_res_PRS <- subset(PRS_all, p.value < sig_threshold) # [, c("outcomes", "p.value", "q.value")]
    if(nrow(sig_res_anov)>0) {
      pp_anov <- plots_varsignificant(sig_res_anov, DF, savingName="res_anov", 
                                      path=paste(path, "Plots_", savingName, "/", sep=""))
      ggsave_adaptSize(paste(path, "Plots_", savingName, "/sigTaxa_anov.png", sep=""), pp_anov)
    }
    if(nrow(sig_res_PRS)>0) {
      sig_res_PRS$names <- compar_var
      pp_PRS <- plots_varsignificant(sig_res_PRS, DF, savingName="res_PRS", 
                                     path=paste(path, "Plots_", savingName, "/", sep=""))
      ggsave_adaptSize(paste(path, "Plots_", savingName, "/sigTaxa_PRS.png", sep=""), pp_PRS)
      
    }
  }
  
  return(to_return)
}


try_err_model <- function(mod, outcome){
  if(any(class(mod) %in% "try-error")){
    logging::loginfo(paste("cannot run model for ", outcome, "\n", sep=""))
    # print(cat(red(paste("cannot run model for ", outcome, "\n", sep=""))))
    mod <- NULL
  }
  return(mod)
}


# From https://delladata.fr/introduction-aux-glmm-avec-donnees-de-proportion/ 
# after Ben Bolker proposition https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#overdispersion
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}



run_mixed_model <- function(outcomes, FE, RE, compar_var, DF){
  
  nbOut <- length(outcomes)
  PRS_all <- c()
  res_all <- c()
  res_anov_all <- c()
  for(i in seq(nbOut)){
    outcome_i <- outcomes[i]
    
    random.effects <- paste("(1|", RE, ")", sep = "", collapse = " + ")
    formula <- formula(paste(outcome_i, " ~ ", paste(paste(paste(FE, collapse = " + "), random.effects, sep=" + "))))
    
    mod <- lmer(formula, data=DF, REML = TRUE)
    lm_summary <- coef(summary(mod))
    para <- as.data.frame(lm_summary)[-1, -c(3:4)]
    CI <- confint(mod, parm=rownames(para))
    CI <- as.data.frame(CI)
    
    simRes <- simulateResiduals(fittedModel = mod, plot = F)
    print(plot(simRes, title = paste("DHARMa residual for", outcome_i)))
    
    para <- cbind(para[, 1], CI[,1:2], para[, 2:3])
    colnames(para) <- c("est", "ci_low", "ci_high", "std.e", "p.value")
    para$names <- rownames(lm_summary)[-1]
    para$feature <- outcome_i
    rownames(para) <- NULL
    para <- para[, c(7, 6, 1:5)]
    
    ANOVA <- car::Anova(mod, type="II")
    res_anov_i <- round(as.data.frame(ANOVA), 5)
    res_anov_i$names <- rownames(res_anov_i)
    rownames(res_anov_i) <- NULL
    res_anov_i$feature <- outcome_i
    res_anov_i <- dplyr::rename(res_anov_i, p.value=`Pr(>Chisq)`)
    res_anov_i <- res_anov_i[, c(5,4, 1:3)]
    
    if(!is.null(compar_var)){
      emm <- emmeans(mod, compar_var)
      PRS_i <- as.data.frame(pairs(emm, infer=T))
      # PRS_i[, 2:6] <-  round(PRS_i[, 2:6] , 4)
      # PRS_i[, 7:8] <-  round(PRS_i[, 7:8] , 4)
      PRS_i[, 2:6] <-  PRS_i[, 2:6]
      PRS_i[, 7:8] <-  PRS_i[, 7:8] 
      PRS_i$feature <- outcome_i
      PRS_i <- PRS_i[, c(9, 1:8)]
    } else{
      PRS_i <- data.frame(feature=outcome_i, contrast = NA, estimate=NA, SE=NA, p.value=NA, q.value=NA)
    }
    
    
    res_all <- rbind(res_all, para)
    res_anov_all <- rbind(res_anov_all, res_anov_i)
    PRS_all <- rbind(PRS_all, PRS_i)
  }
  if(nbOut > 2) res_all$q.value <- round(p.adjust(res_all$p.value, method="fdr"), 4)
  if(nbOut > 2) res_anov_all$q.value <- round(p.adjust(res_anov_all$p.value, method="fdr"), 4) 
  if(nbOut > 2) PRS_all$q.value <- round(p.adjust(PRS_all$p.value, method="fdr"), 4)
  
  
  res_all <- res_all[order(res_all$p.value ),]
  if(!is.null(compar_var)) PRS_all <- PRS_all[order(abs(PRS_all$t.ratio), decreasing = T),]
  res_anov_all <- res_anov_all[order(res_anov_all$Chisq, decreasing = T),]
  
  return(list(res_all=res_all, res_anov_all=res_anov_all, PRS_all=PRS_all))
}



library(GLMMadaptive)
library(cplm)
library(emmeans)

run_uniMarker_model <- 
  run_taxa_model <- function(outcomes, package, model=NULL, counts=F, FE, RE, DF, compar_var, var_check=NULL, 
                             savingName=NULL, sig_threshold=0.05, check_mod=F, save_mod=F, path=NULL, verbose=F, 
                             save_final_files=T, trans_plot = "psd_log", Xinputs=NULL,  ...){
    
    if (!file.exists(path)){
      dir.create(path, recursive =T)
    }
    log_file <- file.path(path, "run_uniMarker_model.log")
    logging::basicConfig(level = "FINEST")
    logging::addHandler(logging::writeToFile, file = log_file, level = "DEBUG")
    logging::setLevel(20, logging::getHandler("basic.stdout"))
    logging::logdebug("outcomes: %s", outcomes)
    logging::logdebug("Xinputs: %s", Xinputs)
    logging::logdebug("package: %s", package)
    logging::logdebug("model: %s", model)
    logging::logdebug("counts: %s", counts)
    logging::logdebug("FE: %s", FE)
    logging::logdebug("RE: %s", RE)
    logging::logdebug("compar_var: %s", compar_var)
    logging::logdebug("var_check: %s", var_check)
    logging::logdebug("sig_threshold: %s", sig_threshold)
    logging::logdebug("check_mod: %s", check_mod)
    logging::logdebug("save_mod: %s", save_mod)
    
    
    DF <- droplevels(DF)
    DF0 <- DF
    if(is.null(var_check) &!is.null(FE)) var_check <- unique(unlist(strsplit(FE, "\\*")))
    
    res_anov_all <- c()
    emm_all <- c()
    PRS_all <- c()
    res_all <- c()
    fit_stat <- c()
    inf_obs_all <- c()
    inf_part_all <- c()
    
    nbOut <- length(outcomes)
    
    # For models with interaction effects
    if(!is.null(FE)){
      is_FEinteract  <- grepl("\\*", FE)
      interaction_in_model <- any(is_FEinteract)
      interact_vars <- unlist(strsplit(FE[is_FEinteract ], "\\*"))
      interact_var <- interact_vars[!interact_vars==compar_var]
    }else{
      interaction_in_model <- F
      interact_var <- NULL
    }
    FE0 <- FE
    
    # When we run independent model on multiple inputs and one outcome (taxa as independent variables)
    if(!is.null(Xinputs)){
      nbOut <- length(Xinputs)
      outcome_i <- outcomes
      FE0 <- FE
      interact_var0 <- interact_var
    }
    
    
    for(i in seq(nbOut)){
      DF <- DF0
      
      if(is.null(Xinputs)){
        outcome_i <- outcomes[i]
        cat(crayon::bold(green(paste(outcome_i, ": ", i, "/", nbOut,"\n", sep=""))))
      } else {
        FE <- sub("Xinputs", Xinputs[i], FE0) 
        if(length(FE)==0) FE <- c(Xinputs[i], FE0) # for compatibility with version before interaction when FE=NULL (COVID study)
        interact_var <- sub("Xinputs", Xinputs[i], interact_var0) 
        var_check <- unique(unlist(strsplit(FE, "\\*")))
        cat(crayon::bold(green(paste(Xinputs[i], ": ", i, "/", nbOut,"\n", sep=""))))
      }
      
      
      DF$expr <- DF[, outcome_i]
      
      id_NA <- unique(which(is.na(DF[, c(outcome_i, var_check)]), arr.ind = T)[,"row"])
      if(length(id_NA) > 0){
        cat(green(paste("NA detected in outcome, groupvars or var_check, ", length(id_NA), " rows removed. \n", sep="")))
        if(verbose) cat(green("Removed samples: \n"))
        if(verbose) print(DF[id_NA, c("Participant_ID", "Timepoint", outcome_i, var_check)])
        DF =  DF[-id_NA,]
        DF <- droplevels(DF)
      } 
      
      
      var_to_keep <- c("Participant_ID", "Sample_ID", "Timepoint", outcome_i, var_check, "expr")
      var_to_keep <- var_to_keep[var_to_keep %in% colnames(DF)]
      DF_mod <<- DF[, var_to_keep] 
      # otherwise predictorEffects() does not want to work (added Timepoint otherwise no influential plot)
      
      
      cat(green("Run model \n"))
      
      
      if(package=="cplm"){
        formula <- formula(paste(outcome_i, " ~ ", paste(paste(FE, collapse = " + "))))
        mod <- try(cpglm(formula, data=DF_mod), T)
        mod <- try_err_model(mod, outcome_i)
        mod_name <- "cpglm"
        print(outcome_i)
        if(!is.null(mod)){
          invisible(capture.output(lm_summary <- coef(summary(mod))))
          para <- as.data.frame(lm_summary)[-1, -3]
          CI <- data.frame(ci_low=NA, ci_high=NA)
        }
      }
      
      if(package=="lm"){
        formula_char <- paste(outcome_i, " ~ ", paste(paste(FE, collapse = " + ")))
        print(formula_char)
        formula <- formula(formula_char)
        mod <- try(lm(formula, data=DF_mod), T)
        mod <- try_err_model(mod, outcome_i)
        mod_name <- "lm"
        if(!is.null(mod)){
          invisible(capture.output(lm_summary <- coef(summary(mod))))
          para <- as.data.frame(lm_summary)[-1, -3]
          CI <- as.data.frame(confint(mod, parm=rownames(para)))
        }
      }
      # print(subset(DF_mod, group_info=="Obesity I"))
      # print(summary(mod))
      
      if(package=="glm"){
        formula_char <- paste(outcome_i, " ~ ", paste(paste(FE, collapse = " + ")))
        print(formula_char)
        formula <- formula(formula_char)
        mod <- try(glm(formula, data=DF_mod, family = "binomial"), T)
        mod <- try_err_model(mod, outcome_i)
        mod_name <- "glm"
        if(!is.null(mod)){
          invisible(capture.output(lm_summary <- coef(summary(mod))))
          para <- as.data.frame(lm_summary)[-1, -3]
          CI <- data.frame(ci_low=NA, ci_high=NA)
          #CI <- suppressMessages(as.data.frame(confint(mod)[rownames(para), c(1,2), drop=F]))
        }
      }
      
      if(package=="lme4"){
        random.effects <- paste("(1|", RE, ")", sep = "", collapse = " + ")
        formula <- formula(print(paste(outcome_i, " ~ ", paste(paste(paste(FE, collapse = " + "), random.effects, sep=" + ")))))
        
        if(is.null(model)) model <- "LMM"
        
        if(model == "LMM"){ #  Linear Mixed-effects Model
          mod <- lmer(formula, data=DF_mod, REML = TRUE)
          mod_name <- "lmer"
          lm_summary <- coef(summary(mod))
          para <- as.data.frame(lm_summary)[-1, -c(3:4)]
          CI <- confint(mod, parm=rownames(para))
        }
        
        if(model == "GLMM"){ # Generalized Linear Mixed-effects Model
          mod <- glmer(formula, data=DF_mod, family = binomial)
          mod_name <- "glmer"
          lm_summary <- coef(summary(mod))
          para <- as.data.frame(lm_summary)[-1, -3]
          CI <- confint(mod, parm=rownames(para), method="Wald") #  faster-but-less-accurate Wald confidence intervals
        }
        
      }
      if(package=="glmmTMB"){
        random.effects <- paste("(1|", RE, ")", sep = "", collapse = " + ")
        formula <- formula(paste("expr ~ ", paste(paste(FE, collapse = " + "), random.effects, sep=" + ")))
        if(counts){ # ZINB
          suppressWarnings(
            mod <- try(glmmTMB::glmmTMB(formula,
                                        data = DF_mod, 
                                        family = glmmTMB::nbinom2(link = "log"), 
                                        ziformula = ~ 1), T))
          mod <- try_err_model(mod, outcome_i)
          mod_name <- "ZINB"
        }else{ # CPLM
          suppressWarnings(
            mod <- try(glmmTMB::glmmTMB(formula,
                                        data = DF_mod,
                                        family = glmmTMB::tweedie(link = "log"),
                                        ziformula = ~ 0), T))
          mod <- try_err_model(mod, outcome_i)
          mod_name <- "CPLM"
        }
        if(!is.null(mod)){
          lm_summary <- coef(summary(mod))$cond
          para <- as.data.frame(lm_summary)[-1, -3]
          CI <- confint(mod, parm=rownames(para))
        }
      }
      
      if(package=="GLMMadaptive"){
        form_RE <- formula(paste( "~ 1|", RE, collapse=" + "))
        form_FE <- formula(paste(outcome_i, paste(FE, collapse=" + "), sep=" ~ "))
        if(counts){
          if(model == "ZINB"){
            mod <-  try(mixed_model(form_FE, 
                                    random = form_RE, data = DF_mod,
                                    family = zi.negative.binomial(), zi_fixed = ~ 1), T)
            mod <- try_err_model(mod, outcome_i)
            mod_name <- "ZINB"
            
            # TRY model with a more complex zi random structure
            # if(!is.null(mod)){
            #   mod2 <- try(update(mod, zi_random = form_RE), T)
            #   mod2 <- try_err_model(mod2, outcome_i)
            #   if(!is.null(mod2)){
            #     mod_comp <- anova(mod, mod2)
            #     if(mod_comp$p.value < 0.05){
            #       mod <- mod2
            #       logging::logdebug(paste(outcome_i, ": keep model with zi_random= ", RE, "\n", sep=""))
            #       cat(green(paste(outcome_i, ": keep model with zi_random= ", RE, "\n", sep="")))
            #       mod_name <- "ZINB_zi_random"
            #     }
            #   }
            # }
          }
          if(model == "HuNB"){
            mod <-  try(mixed_model(form_FE, 
                                    random = form_RE, data = DF_mod,
                                    family = hurdle.negative.binomial(), zi_fixed = ~ 1), T)
            mod <- try_err_model(mod, outcome_i)
            mod_name <- "HuNB"
            
            # TRY model with a more complex zi random structure
            # if(!is.null(mod)){
            #   mod2 <- try(update(mod, zi_random = form_RE), T)
            #   mod2 <- try_err_model(mod2, outcome_i)
            #   if(!is.null(mod2)){
            #     mod_comp <- anova(mod, mod2)
            #     if(mod_comp$p.value < 0.05){
            #       mod <- mod2
            #       logging::logdebug(paste(outcome_i, ": keep model with zi_random= ", RE, "\n", sep=""))
            #       cat(green(paste(outcome_i, ": keep model with zi_random=", RE,"\n", sep="")))
            #       mod_name <- "HuNB_zi_random"
            #     }
            #   }
            # }
          }
        }
        else{
          if(model == "HuLN"){
            mod <-  try(mixed_model(form_FE, 
                                    random = form_RE, data = DF_mod,
                                    family = hurdle.lognormal(), zi_fixed = ~ 1), T)
            mod <- try_err_model(mod, outcome_i)
            mod_name <- "HuLN"
            
            # TRY model with a more complex zi random structure
            # if(!is.null(mod)){
            #   mod2 <- try(update(mod, zi_random = form_RE), T)
            #   mod2 <- try_err_model(mod2, outcome_i)
            #   if(!is.null(mod2)){
            #     mod_comp <- anova(mod, mod2)
            #     if(mod_comp$p.value < 0.05){
            #       mod <- mod2
            #       logging::logdebug(paste(outcome_i, ": keep model with zi_random= ", RE, "\n", sep=""))
            #       cat(green(paste(outcome_i, ": keep model with zi_random=", RE,"\n", sep="")))
            #       mod_name <- "HuLN_zi_random"
            #     }
            #   }
            # }
          }
          if(model == "HuB"){
            mod <-  try(mixed_model(form_FE, 
                                    random = form_RE, data = DF_mod,
                                    family = hurdle.beta.fam(), zi_fixed = ~ 1), T)
            mod <- try_err_model(mod, outcome_i)
            mod_name <- "HuB"
            
            # TRY model with a more complex zi random structure
            # if(!is.null(mod)){
            #   mod2 <- try(update(mod, zi_random = form_RE), T)
            #   mod2 <- try_err_model(mod2, outcome_i)
            #   if(!is.null(mod2)){
            #     mod_comp <- anova(mod, mod2)
            #     if(mod_comp$p.value < 0.05){
            #       mod <- mod2
            #       logging::logdebug(paste(outcome_i, ": keep model with zi_random= ", RE, "\n", sep=""))
            #       cat(green(paste(outcome_i, ": keep model with zi_random=", RE,"\n", sep="")))
            #       mod_name <- "HuB_zi_random"
            #     }
            #   }
            # }
          }
        }
        if(!is.null(mod)){
          lm_summary <- coef(summary(mod))
          para <- as.data.frame(lm_summary)[-1, -3]
          CI <- confint(mod)[rownames(para), c(1,3)]
        }
      }
      
      if(package=="NBZIMM"){
        form_RE <- formula(paste( "~ 1|", RE, collapse=" + "))
        form_FE <- formula(paste(outcome_i, paste(FE, collapse=" + "), sep=" ~ "))
        form_FE2 <- formula(paste("out_asinsqrt", paste(FE, collapse=" + "), sep=" ~ "))
        if(counts){ 
          # Zero inflated negative binomial mixed models (ZINBMM)
          mod <- try(glmm.zinb(form_FE, data = DF_mod, 
                               random = form_RE, zi_fixed = ~1, zi_random =NULL), T)
          mod <- try_err_model(mod, outcome_i)
          mod_name <- "ZINB"
        } else{ 
          # Zero inflated Gaussian mixed models (ZIGMM)
          mod <- try(lme.zig(form_FE, data = DF_mod, 
                             random = form_RE, zi_fixed = ~1, zi_random =NULL), T)
          mod <- try_err_model(mod, outcome_i)
          mod_name <- "ZIG"
        }
        
        if(!is.null(mod)){
          lm_summary <- coef(summary(mod))
          para <- as.data.frame(lm_summary)[-1, -c(3,4)]
          CI <- try(intervals(mod)$fixed[rownames(para), c(1,3)], T)
          if(any(class(CI) %in% "try-error")) CI <- data.frame(ci_low=NA, ci_high=NA)
        }
      }
      
      
      
      
      
      if(is.null(mod)){
        para <- data.frame(feature=outcome_i, names=NA, est=NA, ci_low = NA, ci_high=NA, std.e=NA, p.value=NA)
        if(!is.null(Xinputs)){
          para$prev_xinput <-  round(sum(DF[,Xinputs[i]] > 0)/nrow(DF0),2)
          para$abun_xinput <-  round(mean(DF[,Xinputs[i]]), 4) 
          para <- relocate(para, prev_xinput, .before=est)
          para <- relocate(para, abun_xinput, .before=est)
        }
        
        if(mod_name == "lm") res_anov_i <- data.frame(feature=outcome_i, names=NA, Sum.Sq=NA, Df=NA, F.value=NA, p.value = NA)
        else res_anov_i <- data.frame(feature=outcome_i, names=NA, Chisq=NA, Df=NA, p.value = NA)
        emm_i <- data.frame(feature=outcome_i, compar_var=NA, emmean=NA, SE=NA, df = NA, lower.CL = NA, upper.CL=NA)
        if(!is.null(compar_var)) emm_i <- dplyr::rename(emm_i, !!compar_var:=compar_var)
        if(!interaction_in_model) emm_i[,interact_var] <- NULL
        PRS_i <- data.frame(feature=outcome_i, contrast=NA, estimate=NA, SE = NA, df=NA, lower.CL = NA, upper.CL=NA, 
                            prev= round(sum(DF[,outcome_i] > 0)/nrow(DF0),2),
                            abun=round(mean(DF[,outcome_i]), 4),
                            ratio=NA, p.value=NA)
        if(!is.null(Xinputs)){
          PRS_i$prev_xinput <-  round(sum(DF[,Xinputs[i]] > 0)/nrow(DF0),2)
          PRS_i$abun_xinput <-  round(mean(DF[,Xinputs[i]]), 4) 
          PRS_i <- relocate(PRS_i, prev_xinput, .before=ratio)
          PRS_i <- relocate(PRS_i, abun_xinput, .before=ratio)
        }
        
        if(!interaction_in_model)PRS_i[,interact_var] <- NULL
        fit_stat_i <- data.frame(feature=outcome_i, model=mod_name, aic=NA, bic=NA, logLik = NA, 
                                 overdisp_test=NA, tzi=NA, tdisp=NA)
        
        inf_obs_i <- data.frame(Participant_ID=NA, Timepoint=NA, cooksd = NA, outcome_i=outcome_i)
        inf_part_i <- data.frame(Participant_ID=NA, cooksd = NA, outcome_i=outcome_i)
        
      }
      
      if(!is.null(mod)){
        cat(green("Check model \n"))
        
        if(check_mod){
          if(is.null(Xinputs)) saving_name = outcome_i else saving_name = Xinputs[i]
          c_m <- check_model(mod, outcome_i, DF_mod, var_check, print=F, save=T, 
                             path=paste(path, "Plots/", sep=""), verbose=verbose, 
                             resid_panel_pck="none", saving_name=saving_name)
        } 
        if(save_mod){
          if (!file.exists(paste(path, "models/", sep=""))){
            dir.create(paste(path, "models/", sep=""), recursive =T)
          }
          save(mod, file=paste(path, "models/", "mod_", outcome_i, ".Rdata", sep=""))
        } 
        
        aic=AIC(mod)
        bic=tryCatch(BIC(mod), error = function(e) {return(NA)}) 
        logLik = tryCatch(logLik(mod), error = function(e) {return(NA)}) 
        overdisp_test = tryCatch(overdisp_fun(mod)["p"], error = function(e) {return(NA)}) 
        
        if(exists("c_m")){
          tzi <- tryCatch(testZeroInflation(c_m$simRes, plot=F)$p.value, error = function(e) {return(NA)}) 
          tdisp <- tryCatch(testDispersion(c_m$simRes_cond, plot=F)$p.value, error = function(e) {return(NA)}) 
          
          if(nrow(c_m$inf_obs)>0){
            inf_obs_i <- cbind(c_m$inf_obs[, c("Participant_ID", "Timepoint", "cooksd")], outcome_i)
          } else inf_obs_i <- data.frame(Participant_ID=NA, Timepoint=NA, cooksd = NA, outcome_i=outcome_i)
          
          if(nrow(c_m$inf_part)>0){
            inf_part_i <- cbind(c_m$inf_part[, c("Participant_ID", "cooksd")], outcome_i)
          } else inf_part_i <- data.frame(Participant_ID=NA, cooksd = NA, outcome_i=outcome_i)
        }else {
          tzi <- NA
          tdisp <- NA
          inf_obs_i <- data.frame(Participant_ID=NA, Timepoint=NA, cooksd = NA, outcome_i=outcome_i)
          inf_part_i <- data.frame(Participant_ID=NA, cooksd = NA, outcome_i=outcome_i)
        }
        
        
        fit_stat_i <- data.frame(feature=outcome_i, model=mod_name, aic=aic, bic=bic, logLik = logLik, 
                                 overdisp_test=overdisp_test, tzi=tzi, tdisp=tdisp)
        
        
        
        para <- cbind(para[, 1], CI[,1:2], para[, 2:3])
        colnames(para) <- c("est", "ci_low", "ci_high", "std.e", "p.value")
        para$names <- rownames(lm_summary)[-1]
        para$feature <- outcome_i
        rownames(para) <- NULL
        para <- para[, c(7, 6, 1:5)]
        if(!is.null(Xinputs)){
          para$prev_xinput <-  round(sum(DF[,Xinputs[i]] > 0)/nrow(DF0),2)
          para$abun_xinput <-  round(mean(DF[,Xinputs[i]]), 4) 
          para <- relocate(para, prev_xinput, .before=est)
          para <- relocate(para, abun_xinput, .before=est)
        }
        
        ANOVA <- try(car::Anova(mod, type="II"), T) # type="III"
        
        if(!is.null(compar_var)){
          if(interaction_in_model){
            # if(class(DF_mod[, compar_var])=="numeric"){# compar_var continue and interact_var discrete
            #   # slope estimation
            #   emm <- tryCatch(emtrends(mod, interact_var, var= compar_var), error = function(e) {return(NA)})
            # }else{# compar_var discrete
            if(any(class(DF_mod[,interact_var])=="numeric")){
              # interact_var is continuous: slope estimation
              emm <- tryCatch(emtrends(mod, compar_var, var= interact_var), error = function(e) {return(NA)})
            } else{ # interact_var is factor
              emm <- tryCatch(emmeans(mod, compar_var, by= interact_var), error = function(e) {return(NA)})
              # }
            }
          }else{
            emm <- tryCatch(emmeans(mod, compar_var), error = function(e) {return(NA)})
            # emm are provided on the transformed scale (eg log scale when using a model with log link)
            # we should add type = "response" to have result on the response and not on the transformed scale 
            if(class(emm) != "emmGrid"){ # before was if(is.na(emm)) but generates warning message when emmeans doesn't generate error
              # from: https://github.com/rvlenth/emmeans/issues/116
              rg <- qdrg(formula(paste(" ~ ", paste(paste(FE, collapse = " + ")))), 
                         data=DF_mod, coef = coef(mod), vcov = vcov(mod))
              emm <- emmeans(rg, compar_var)
            }
          }
          
          if(compar_var %in% c("B", "T")) reverse = T else reverse=F # to have the comparison After - Before and After_YogurtOatmeal - Other
          
          PRS <- pairs(emm, infer=T, reverse=reverse)
          # PRS2 <- pairs(pairs(emm, infer=T, reverse=T), by=NULL) # contrasts of contrasts possible
          # https://stackoverflow.com/questions/64446605/contrast-of-contrast-with-emmeans-second-differences
          
        }
        
        
        if(!(any(class(ANOVA) %in% "try-error"))){
          res_anov_i <- round(as.data.frame(ANOVA), 5)
          res_anov_i <- res_anov_i[which(rownames(res_anov_i) != "Residuals"),]
          res_anov_i$names <- rownames(res_anov_i)
          rownames(res_anov_i) <- NULL
          res_anov_i$feature <- outcome_i
          if("Pr(>Chisq)" %in% colnames(res_anov_i)) res_anov_i <- dplyr::rename(res_anov_i, p.value=`Pr(>Chisq)`)
          if("Pr(>F)" %in% colnames(res_anov_i)) res_anov_i <- dplyr::rename(res_anov_i, p.value=`Pr(>F)`)
          #res_anov_i <- res_anov_i[, c(5,4, 1:3)]
          res_anov_i <- relocate(res_anov_i, names)
          res_anov_i <- relocate(res_anov_i, feature)
          colnames(res_anov_i) <- gsub(" ", ".", colnames(res_anov_i))
        }else{
          res_anov_i <- data.frame(feature=outcome_i, names=NA, Chisq=NA, Df=NA, p.value = NA)
        }
        
        
        class_out <- class(DF[,outcome_i])
        if(!is.null(compar_var)){
          emm_i <-  as.data.frame(emm, destroy.annotations =T)
          emm_i[, sapply(emm_i, is.numeric)] <- round(emm_i[, sapply(emm_i, is.numeric)], 2)
          emm_i$feature <- outcome_i
          emm_i <- relocate(emm_i, feature)
          id_CL <- which(grepl("CL", colnames(emm_i)))
          colnames(emm_i)[id_CL] <- c("lower.CL", "upper.CL")
          # colnames(emm_i)[6:7] <- c("lower.CL", "upper.CL")
          id_trend <- which(grepl("trend", colnames(emm_i)))
          colnames(emm_i)[id_trend] <- c("trend")
          
          PRS_i <-  as.data.frame(PRS, destroy.annotations =T)
          if(interaction_in_model & any(class(DF_mod[,interact_var])=="numeric")){
            PRS_i$comments <- "slope difference" # equal interaction coefficient
          }
          digits = 1
          if(min(abs(PRS_i$estimate), na.rm = T) < 100) digits = 2
          if(min(abs(PRS_i$estimate), na.rm = T) < 1) digits = 3
          if(min(abs(PRS_i$estimate), na.rm = T) < 0.1) digits = 4
          if(min(abs(PRS_i$estimate), na.rm = T) < 0.01) digits = 5
          if(min(abs(PRS_i$estimate), na.rm = T) < 0.001) digits = 6
          if(min(abs(PRS_i$estimate), na.rm = T) < 0.00001) digits = 7
          if(min(abs(PRS_i$estimate), na.rm = T) < 0.0000001) digits = 8
          PRS_i[, sapply(emm_i, is.numeric)] <-  round(PRS_i[, sapply(emm_i, is.numeric)] , digits)
          id_CL <- which(grepl("CL", colnames(PRS_i)))
          id_ratio <- which(grepl("ratio", colnames(PRS_i)))
          colnames(PRS_i)[id_CL] <- c("lower.CL", "upper.CL")
          colnames(PRS_i)[id_ratio] <- "ratio"
          PRS_i[, c("ratio", "p.value")] <-  round(PRS_i[, c("ratio", "p.value")] , 4)
          PRS_i$feature <- outcome_i
          if(class_out == "factor"){
            PRS_i$prev <- NA
            PRS_i$abun <- NA
          } else{
            PRS_i$prev <-  round(sum(DF[,outcome_i] > 0)/nrow(DF0),2)
            PRS_i$abun <-  round(mean(DF[,outcome_i]), 4) 
          }
          if(!is.null(Xinputs)){
            PRS_i$prev_xinput <-  round(sum(DF[,Xinputs[i]] > 0)/nrow(DF0),2)
            PRS_i$abun_xinput <-  round(mean(DF[,Xinputs[i]]), 4) 
          }
          PRS_i <- relocate(PRS_i, feature)
          PRS_i <- relocate(PRS_i, prev, .before=ratio)
          PRS_i <- relocate(PRS_i, abun, .before=ratio)
          if(!is.null(Xinputs)){
            PRS_i <- relocate(PRS_i, prev_xinput, .before=ratio)
            PRS_i <- relocate(PRS_i, abun_xinput, .before=ratio)
          }
          # PRS_i <- PRS_i[, c(9, 1:6, 10, 11, 7, 8)]#use data before log transformation
          # colnames(PRS_i)[c(6, 7, 10)] <- c("lower.CL", "upper.CL", "ratio")
        } else{
          emm_i <- data.frame(feature=outcome_i, compar_var=NA, emmean=NA, SE=NA, df = NA, lower.CL = NA, upper.CL=NA)
          
          if(class_out == "factor"){
            PRS_i <- data.frame(feature=outcome_i, contrast=NA, estimate=NA, SE = NA, df=NA, lower.CL = NA, upper.CL=NA, 
                                prev= NA, abun= NA, ratio=NA, p.value=NA)
          }else{
            PRS_i <- data.frame(feature=outcome_i, contrast=NA, estimate=NA, SE = NA, df=NA, lower.CL = NA, upper.CL=NA, 
                                prev= round(sum(DF[,outcome_i] > 0)/nrow(DF0),2),
                                abun= round(mean(DF[,outcome_i]), 4),
                                ratio=NA, p.value=NA)
          }
          
          if(!is.null(Xinputs)){
            PRS_i$prev_xinput <-  round(sum(DF[,Xinputs[i]] > 0)/nrow(DF0),2)
            PRS_i$abun_xinput <-  round(mean(DF[,Xinputs[i]]), 4) 
            PRS_i <- relocate(PRS_i, prev_xinput, .before=ratio)
            PRS_i <- relocate(PRS_i, abun_xinput, .before=ratio)
          }
        }
      }
      if(!is.null(Xinputs)){
        para$Xinputs <- Xinputs[i]
        res_anov_i$Xinputs <- Xinputs[i]
        fit_stat_i$Xinputs <- Xinputs[i]
        inf_obs_i$Xinputs <- Xinputs[i]
        inf_part_i$Xinputs <- Xinputs[i]
        emm_i$Xinputs <- Xinputs[i]
        PRS_i$Xinputs <- Xinputs[i]
        
        para <- relocate(para, Xinputs, .after=feature)
        res_anov_i <- relocate(res_anov_i, Xinputs, .after=feature)
        fit_stat_i <- relocate(fit_stat_i, Xinputs, .after=feature)
        
        emm_i <- relocate(emm_i, Xinputs, .after=feature)
        PRS_i <- relocate(PRS_i, Xinputs, .after=feature)
        
      }
      
      res_all <- rbind(res_all, para)
      res_anov_all <- rbind(res_anov_all, res_anov_i)
      fit_stat <- rbind(fit_stat, fit_stat_i)
      inf_obs_all <- rbind(inf_obs_all, inf_obs_i)
      inf_part_all <- rbind(inf_part_all, inf_part_i)
      emm_all <- rbind(emm_all, emm_i)
      PRS_all <- rbind(PRS_all, PRS_i)
    }
    
    
    res_all$q.value <- round(p.adjust(res_all$p.value, method="fdr"), 4)
    res_all <- res_all[order(res_all$p.value ),]
    
    # if(!(any(class(ANOVA) %in% "try-error"))){
    #   res_anov_all$q.value <- round(p.adjust(res_anov_all$p.value, method="fdr"), 4)
    #   res_anov_all <- res_anov_all[order(res_anov_all$Chisq, decreasing = T),]
    # } else {
    #   res_anov_all <- data.frame()
    # }
    res_anov_all$q.value <- round(p.adjust(res_anov_all$p.value, method="fdr"), 4) # never tested outside if(!(any(class(ANOVA) %in% "try-error"))){
    
    if("LR.Chisq" %in% colnames(res_anov_all)) res_anov_all <- res_anov_all[order(res_anov_all$p.value, decreasing = F),]
    if("Chisq" %in% colnames(res_anov_all)) res_anov_all <- res_anov_all[order(res_anov_all$Chisq, decreasing = T),]
    if("F.value" %in% colnames(res_anov_all)) res_anov_all <- res_anov_all[order(res_anov_all$F.value, decreasing = T),]
    
    PRS_all$q.value <- round(p.adjust(PRS_all$p.value, method="fdr"), 4)
    PRS_all <- PRS_all[order(abs(PRS_all$ratio), decreasing = T),]
    
    
    to_return <- list(res_all=res_all, res_anov_all=res_anov_all, 
                      emm_all=emm_all, PRS_all=PRS_all, fit_stat=fit_stat)
    
    if(!save_final_files) print(to_return)
    if(save_final_files){
      ### excel file
      cat(green("Saving results models as excel files \n"))
      # write.xlsx(res_all, file=paste(path, savingName, "_results.xlsx", sep=""), sheet="res_all", row.names = F)
      # #if(!(any(class(ANOVA) %in% "try-error"))) write.xlsx(res_anov_all, file=paste(path, savingName, "_results.xlsx", sep=""), sheet="res_anov_all", append = T, row.names = F)
      # write.xlsx(res_anov_all, file=paste(path, savingName, "_results.xlsx", sep=""), sheet="res_anov_all", append = T, row.names = F)
      # write.xlsx(emm_all, file=paste(path, savingName, "_results.xlsx", sep=""), sheet="emm_all", append = T, row.names = F)
      # write.xlsx(PRS_all, file=paste(path, savingName, "_results.xlsx", sep=""), sheet="PRS_all", append = T, row.names = F)
      # write.xlsx(fit_stat, file=paste(path, savingName, "_results.xlsx", sep=""), sheet="fit_stat", append = T, row.names = F)
      # if(check_mod) write.xlsx(inf_obs_all, file=paste(path, savingName, "_results.xlsx", sep=""), sheet="inf_obs_all", append = T, row.names = F)
      # if(check_mod) write.xlsx(inf_part_all, file=paste(path, savingName, "_results.xlsx", sep=""), sheet="inf_part_all", append = T, row.names = F)
      DF_tosave <- list(res_all=res_all, res_anov_all=res_anov_all, emm_all=emm_all, PRS_all=PRS_all, fit_stat=fit_stat)
      if(check_mod) DF_tosave <- append(DF_tosave, list(inf_obs_all=inf_obs_all, inf_part_all=inf_part_all)) 
      write_xlsx(DF_tosave, paste(path, savingName, "_results.xlsx", sep=""))
      
      
      if(nrow(res_anov_all) > 0 & nrow(res_anov_all) <= 50 ){
        cat(green("Saving results models as one png \n"))
        pp <- print_result_MM(to_return, sig_threshold)
        if(interaction_in_model){
          ggsave(paste(path, savingName, "_tables.png", sep=""),
                 wrap_elements(pp$tg1) + pp$tg2 + pp$tg3 + plot_layout(widths = c(1, 2, 2)),
                 height=12, width=30)
        }else{
          ggsave(paste(path, savingName, "_tables.png", sep=""),
                 wrap_elements(pp$tg1) + pp$tg2 / pp$tg3 + plot_layout(widths = c(1, 2)),
                 height=12, width=19)
        }
      }
      
      plot_results(res_all, res_anov_all, PRS_all, DF0, FE0, sig_threshold, 
                   savingName, path, counts, trans_plot, compar_var, ...)# change DF to DF0 (21/07/23)
    }
  }


plot_results <- function(res_all, res_anov_all, PRS_all, DF, FE, sig_threshold, 
                         savingName, path, counts=F, trans_plot="psd_log", compar_var, ...){
  
  ## Separate interaction from main effects
  if(!is.null(FE)){
    is_FEinteract  <- grepl("\\*", FE)
    interaction_in_model <- any(is_FEinteract)
    interact_vars <- unlist(strsplit(FE[is_FEinteract ], "\\*"))
    interact_var <- interact_vars[!interact_vars==compar_var]
    FE <- unique(unlist(strsplit(FE, "\\*")))
  }else{
    interaction_in_model <- F
    interact_var <- NULL
  }
  
  ## keep only significant results
  sig_res_all <- subset(res_all, p.value < sig_threshold) 
  if("p.value" %in% colnames(res_anov_all)){
    sig_res_anov <- subset(res_anov_all, p.value < sig_threshold) 
  } else sig_res_anov <- data.frame()
  sig_res_PRS <- subset(PRS_all, p.value < sig_threshold)
  
  
  ## Separate variable names from labels in res_all
  if(nrow(sig_res_all)>0){
    if("Xinputs" %in% colnames(sig_res_all)){
      FE <- FE[FE!="Xinputs"]
      all_var <- c(unique(sig_res_all$Xinputs), FE) 
    } else all_var <- FE
    lev <- data.frame()
    for(i in seq(length(all_var))){
      level_var_i <- levels(DF[,all_var[i]])
      if(is.null(level_var_i)) level_var_i <- ""
      lev_i <- cbind(all_var[i], level_var_i, paste(all_var[i], level_var_i, sep=""))
      lev <- rbind(lev, lev_i) #all possible combinations of variable and labels
    }
    colnames(lev) <- c("variable", "label", "names")
    sig_res_all <- merge(sig_res_all, lev, by="names")# will remove all interactions variables as they are not define in lev 
    sig_res_all <- dplyr::rename(sig_res_all, variable_label = names)
    sig_res_all <- sig_res_all[order(sig_res_all$p.value ),]
    
    ### Remove variables that are not straightforward to interpret
    var_to_remove <- NULL
    if(!is.null(compar_var)){
      if(interaction_in_model){
        if(compar_var=="T") var_to_remove <- c("B", "T") 
        if(compar_var=="B") var_to_remove <- "B"
      }else{
        if(compar_var=="T") var_to_remove <- c("B") 
        if(compar_var=="B") var_to_remove <- NULL
      }
    }
    
    if(interaction_in_model){
      if(interact_var == "Xinputs") var_to_remove <- c(sig_res_all$Xinputs, interact_vars, var_to_remove)
      else var_to_remove <- c(interact_vars, var_to_remove)
    }
    
    
    sig_res_all <- subset(sig_res_all, !(variable %in% var_to_remove))
    
    if(nrow(sig_res_all) > 100) sig_res_all <- sig_res_all[1:100,] # to not generate more than 100 plots
  }
  
  
  ### Forest plot for significant res_all
  if(nrow(sig_res_all)>0 & !("Xinputs" %in% colnames(sig_res_all))) {
    cat(green("Create forest plot for significant res_all \n"))
    sig_res_all$color <- case_when(sig_res_all$q.value < 0.05 ~  "q.val<0.05",
                                   sig_res_all$q.value >= 0.05 ~ "q.val>=0.05")
    var <- unique(sig_res_all$variable_label)
    ll_plot <- list()
    for(i in seq(var)){
      sub <- subset(sig_res_all, variable_label==var[i])
      sub$feature <- factor(sub$feature, levels = rev(sub$feature))
      pp <- ggplot(sub, aes(x=est, y=feature)) +  geom_point(size=1.2, aes(color=color)) +
        scale_colour_manual(values= c(`q.val>=0.05`="black", `q.val<0.05`="#2196F3"))+
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth=0.5) + 
        xlab("effect size") + ylab("") + ggtitle(var[i]) + theme_minimal()
      if(all(is.na(sub$ci_low))) ll_plot[[i]] <- pp
      else  ll_plot[[i]] <- pp + geom_errorbar(aes(x = est, xmin = ci_low, xmax = ci_high, 
                                                   color=color), width = 0.2, linewidth=0.5) 
      
    }
    
    ggsave_adaptSize(paste(path, savingName,  "_res_all.png", sep=""), ll_plot)
  }
  
  
  ### Forest plot on PRS results
  if(nrow(sig_res_PRS)>0 & !("Xinputs" %in% colnames(sig_res_PRS))){
    cat(green("Create forest plot for significant PRS \n"))
    sig_res_PRS$color <- case_when(sig_res_PRS$q.value < 0.05 ~  "q.val<0.05",
                                   sig_res_PRS$q.value >= 0.05 ~ "q.val>=0.05")
    
    xlab <- unique(sig_res_PRS$contrast)
    
    if(is.null(interact_var)) cond2 <- FALSE else cond2 <- interact_var %in% colnames(sig_res_PRS)
    
    if(interaction_in_model & cond2){
      var <- unique(sig_res_PRS[, interact_var])
      ll_plot <- list()
      for(i in seq(var)){
        sub <- subset(sig_res_PRS, get(interact_var)==var[i])
        sub$feature <- factor(sub$feature, levels = rev(sub$feature))
        pp <- ggplot(sub, aes(x=estimate, y=feature)) +  geom_point(size=1.2, aes(color=color)) +
          scale_colour_manual(values= c(`q.val>=0.05`="black", `q.val<0.05`="#2196F3"))+
          geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth=0.5) +
          xlab(xlab) + ylab("") + ggtitle(paste(interact_var, ": ", var[i], sep="")) + theme_minimal()
        if(all(is.na(sub$lower.CL))) ll_plot[[i]] <- pp
        else  ll_plot[[i]] <- pp + geom_errorbar(aes(x = estimate, xmin = lower.CL, xmax = upper.CL,
                                                     color=color), width = 0.2, linewidth=0.5)
      }
      ggsave_adaptSize(paste(path, savingName,  "_res_PRS.png", sep=""), ll_plot)
      
    } else{
      sub <- sig_res_PRS
      sub$feature <- factor(sub$feature, levels = rev(sub$feature))
      pp <- ggplot(sub, aes(x=estimate, y=feature)) + geom_point(size=1.2, aes(color=color)) +
        scale_colour_manual(values= c(`q.val>=0.05`="black", `q.val<0.05`="#2196F3"))+
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth=0.5) +
        xlab(xlab) + ylab("") + theme_minimal()
      pp  <- pp + geom_errorbar(aes(x = estimate, xmin = lower.CL, xmax = upper.CL,
                                    color=color), width = 0.2, linewidth=0.5)
      ggsave(paste(path, savingName,  "_res_PRS.png", sep=""),  pp, width=6, height=6)
    }
  }
  
  
  ### var significant from res_anov
  ## Enough to plot var significant from res_all below
  # if(nrow(sig_res_anov)>0) {
  #   cat(green("plots var significant from res_anov\n"))
  #   pp_anov <- plots_varsignificant2(sig_res_anov, DF, savingName="res_anov", counts=counts,
  #                                    trans=trans_plot,
  #                                    path=paste(path, "Plots/", sep=""), ...)
  #   ggsave_adaptSize(paste(path, "Plots/sigMarker_anov.png", sep=""), pp_anov)
  # }
  
  
  ### var significant from res_all
  if(nrow(sig_res_all)>0) {
    cat(green("plots var significant from res_all\n"))
    pp_res_all <- plots_varsignificant2(sig_res_all, DF, savingName="res_All", counts=counts,
                                        trans=trans_plot,
                                        path=paste(path, "Plots/", sep=""), ...)
    if(length(pp_res_all) > 56) pp_res_all <- pp_res_all[1:56]
    ggsave_adaptSize(paste(path, "Plots/sigMarker_resAll.png", sep=""), pp_res_all)
    
  }
  
  
  ### var significant from res_PRS
  if(nrow(sig_res_PRS)>0) {
    cat(green("plots var significant from res_PRS\n"))
    # sig_res_PRS$names <- compar_var
    sig_res_PRS$variable <- compar_var
    pp_PRS <- plots_varsignificant2(sig_res_PRS, DF, savingName="res_PRS", counts=counts,
                                    trans=trans_plot,
                                    path=paste(path, "Plots/", sep=""),
                                    interact_var=interact_var,...)
    if(length(pp_PRS) > 56) pp_PRS <- pp_PRS[1:56]
    ggsave_adaptSize(paste(path, "Plots/sigMarker_PRS.png", sep=""), pp_PRS)
  }
}


library(qqplotr) # qqplot
library(ggplotify)
check_model <- function(mod, outcome, DF, var_check, print=T, save=F, path=NULL, 
                        verbose=F, resid_panel_pck="ggResidpanel", saving_name=NULL){
  
  if(any(class(mod) %in% c("cpglm", "lm"))) REmod <- F
  else REmod <- T
  
  if(is.null(saving_name)) saving_name <- outcome
  
  ##### Residuals
  # qqplot/ hist residual
  # marginal res vs X (check linearity of y with covariates)
  # conditional res vs pred (homoscedasticity)
  # conditional res vs y (outliers)
  # y vs pred
  # plot_redres(lmer_mod, type="pearson_mar", xvar="Alter")# marginal res vs X 
  # plot_redres(lmer_mod, type="pearson_cond") # conditional res vs pred
  # plot_redres(lmer_mod, type = "pearson_cond", xvar = "Species:Gemmiger_formicilis") # cond res vs y
  
  # output <- list()
  
  
  if (!file.exists(paste(wd0, "/", path, "model_check/", sep=""))){
    dir.create(paste(wd0, "/", path, "model_check/", sep=""), recursive =T)
  }
  
  #setwd(paste(wd0, "/", path, "model_check/", sep=""))
  #path2 <- paste(path, "model_check/", sep="")
  # path2 <- "model_check/"
  
  
  ######################c
  # Extract model info #
  ######################c
  # extract residuals
  # type response=raw
  if(any(class(mod) == "MixMod")){
    # "mean_subject": only on the fixed effects (marginal residual)
    # "subject_specific": on both the fixed and random effects (conditional residual)
    # "marginal": only on the fixed effects with marginalized β coefficients (marginal residual)
    resid <- data.frame(subject_specific=residuals(mod, type="subject_specific"),
                        mean_subject=residuals(mod, type="mean_subject"),
                        marginal=residuals(mod, type="marginal"))
    resid$Sample_ID <- rownames(resid)
  } 
  else{
    resid <- data.frame(pearson=residuals(mod, type="pearson")) #conditional residual
    resid$Sample_ID <- DF$Sample_ID
    # for "zinbmm", "zigmm", "glmmTMB", "lmer" in same order than DF but to check otherwise
  }
  DF <- merge(DF, resid, by="Sample_ID")
  
  # extract fitted values
  fitted <- fitted(mod)
  
  # extract random effect
  if(REmod){
    if("glmmTMB" %in% class(mod))  RE <- ranef(mod)[[1]]$Participant_ID
    else if(any(class(mod) %in% c("zinbmm", "zigmm", "MixMod"))) RE <-  as.data.frame(ranef(mod))
    else RE <- ranef(mod)[[1]]
  }
  
  
  ##########c
  # DHARMa #
  ##########c
  if(!(any(class(mod) %in% c("zinbmm", "zigmm", "cpglm")))){
    
    if(verbose) cat(green(" - run DHARMa: simulate Residuals\n"))
    simRes <- simulateResiduals(fittedModel = mod, plot = F) #DHARMa
    if(verbose) cat(green(" - run DHARMa: simulate conditional Residuals \n"))
    simRes_cond <- simulateResiduals(fittedModel = mod, plot = F, re.form = NULL)
    
    
    if(print) try(plot(simRes), silent=T)
    if(save){
      # file= paste(saving_name, "_DHARMa_main.png", sep="")
      file= paste(path, "model_check/", saving_name, "_DHARMa_main.png", sep="")
      png(file=file, width=7, height=5, units="in", res=300)
      try(plot(simRes), silent=T)
      dev.off()
      
      #ggsave(paste(path2, saving_name, "_DHARMa_main.png", sep=""), width=7, height=5)
      
    }
    
    if(print){
      par(mfrow = c(2, 3))
      for(i in seq(length(var_check))){
        plotResiduals(simRes, form = DF[ ,var_check[i]], xlab=var_check[i])
        mtext(var_check[i], 1, line=2)
      }
    }
    if(save){
      if(length(var_check)<=4){
        width=7 
        height=7
      } else {
        width=10 
        height=10
      }
      
      # file= paste(saving_name, "_DHARMa_var.png", sep="")
      file= paste(path, "model_check/", saving_name, "_DHARMa_var.png", sep="")
      png(file, width=width, height=height, units="in", res=300)
      par(mfrow = n2mfrow(length(var_check)))
      for(i in seq(length(var_check))){
        plotResiduals(simRes, form = DF[ ,var_check[i]], xlab=var_check[i])
        mtext(var_check[i], 1, line=2)
      }
      dev.off()
      
      #ggsave(paste(path2, saving_name, "_DHARMa_var.png", sep=""), width=7, height=5)
      
    }  
  } else {
    simRes <- NULL
    simRes_cond <- NULL
  }
  
  
  
  ############################c
  # Standard residuals plots #
  ############################c
  if(verbose) cat(green(" - Standard residuals plots \n"))
  
  
  if(!all(class(mod) %in% c("lm", "glm", "lme", "lmer", "lmerTest", "lmerModLmerTest", "glmer"))){
    resid_panel_pck = "none"
  }
  if(resid_panel_pck == "ggResidpanel"){
    if(print) resid_panel(mod, "all") # Residuals from ggResidpanel (residuals are divided by sd for some models = resid(model, type = "response")/summary(model)$sigma))
    if(save){
      ggsave(paste(path, "model_check/", outcome, "_resid_panel.png", sep=""),
             resid_panel(mod, "all"), width=12, height=7)
    }
  }else{
    for(i in seq(ncol(resid)-1)){
      resid_i <- colnames(resid)[i]
      qq_resid <- ggplot(data = DF, aes(sample = .data[[resid_i]])) + 
        stat_qq_band(bandType = "pointwise", distribution = "norm", fill = "#FBB4AE", alpha = 0.4) +
        stat_qq_line(distribution = "norm", colour = "#FBB4AE") + 
        qqplotr::stat_qq_point(distribution = "norm") + 
        xlab("Normal quantiles") + theme_bw() + ylab(resid_i) + ggtitle("Q-Q Plot")
      
      sd_resid <- sd(DF[, resid_i])
      hist_resid <- ggplot(data = DF, aes(x = .data[[resid_i]])) + 
        geom_histogram(aes(y = after_stat(density)), color="white", fill = "grey30", bins=30) + theme_bw()+ 
        stat_function(fun = dnorm, color = "blue", args = list(mean = 0, sd = sd_resid))+ ggtitle("Histogram")
      
      pred_resid <- ggplot(data = DF, mapping = aes(x = fitted, y = .data[[resid_i]])) + geom_point() + 
        geom_abline(slope = 0, intercept = 0, color = "blue") + 
        labs(x = "Predicted Values",  y = resid_i)+ theme_bw() + ggtitle("Residual Plot")
      
      y_resid <- ggplot(data = DF, mapping = aes(x = .data[[outcome]], y = .data[[resid_i]])) + geom_point() + 
        geom_abline(slope = 0, intercept = 0, color = "blue") + 
        labs(x = outcome,  y = resid_i)+ theme_bw()+ ggtitle("Residual vs response")
      
      y_pred <- ggplot(data = DF, mapping = aes(x = fitted, y = .data[[outcome]])) + geom_point() + 
        geom_abline(slope = 1, intercept = 0, color = "blue") + 
        labs(x = "Predicted Values",  y = outcome)+ theme_bw()+ ggtitle("Response vs Predicted")
      
      hist_species <- ggplot(data = DF, aes(x = .data[[outcome]])) + 
        geom_histogram(color="white", fill = "grey30", bins=30) + theme_bw()+ 
        labs(x = "pseudo_log10")+ ggtitle(paste("Distribution", outcome))+ 
        scale_x_continuous(trans=scales::pseudo_log_trans(base = 10))
      
      ll_plot <- list()
      ll_plot[[1]] <- pred_resid
      ll_plot[[2]] <- y_resid
      ll_plot[[3]] <- y_pred
      ll_plot[[4]] <- qq_resid
      ll_plot[[5]] <- hist_resid
      ll_plot[[6]] <- hist_species
      
      if(print) print(wrap_plots(ll_plot, ncol=3, nrow=2))
      if(save) ggsave(paste(path, "model_check/", saving_name, "_", resid_i, ".png", sep=""), 
                      wrap_plots(ll_plot, ncol=3, nrow=2), width=12, height=6)
    }
  }
  
  
  
  ########################c
  # Random effects check #
  ########################c
  
  if(REmod){
    if(verbose) cat(green(" - Random effects check \n"))
    
    RE$Participant_ID <- rownames(RE)
    RE <- dplyr::rename(RE, RE_int = `(Intercept)`)
    DF_RE <- merge(DF, RE, by="Participant_ID")
    if("lmerMod" %in% class(mod)){
      suppressWarnings(output$qq_ranef <- plot_ranef(mod)) # Random effect distribution from redres 
    }else{
      qq_ranef <- ggplot(data = RE, aes(sample = RE_int)) + 
        stat_qq_band(bandType = "pointwise", distribution = "norm", fill = "#FBB4AE", alpha = 0.4) +
        stat_qq_line(distribution = "norm", colour = "#FBB4AE") + 
        qqplotr::stat_qq_point(distribution = "norm") + 
        xlab("Normal quantiles") + theme_bw() + ylab("ranef")
    }
    ll_plot <- list()
    ll_plot[[1]] <- qq_ranef
    for(i in seq(length(var_check))){
      if(any(class(DF_RE[,var_check[i]])=="factor")){
        ll_plot[[i+1]] <- ggplot(DF_RE, aes(x=.data[[var_check[i]]], y=RE_int)) + geom_boxplot()+ 
          geom_text_repel(data=DF_RE[is_outlier(DF_RE$RE_int),], aes(label=Sample_ID))
      } 
      else ll_plot[[i+1]] <- ggplot(DF_RE, aes(x=.data[[var_check[i]]], y=RE_int)) + geom_point() + 
          geom_text_repel(data=DF_RE[is_outlier(DF_RE$RE_int),], aes(label=Sample_ID))
    }
    if(print){
      nb_plot <- length(ll_plot)
      ncol <- ceiling(sqrt(nb_plot))
      nrow <- ceiling(nb_plot/ncol)
      print(wrap_plots(ll_plot, ncol=ncol, nrow=nrow))
    }
    if(save) ggsave_adaptSize(paste(path, "model_check/", saving_name, "_ranef.png", sep=""), ll_plot, 4, 3)
  }
  
  
  ######################c
  # Influential points #
  ######################c
  if(all(class(mod) %in% c("lmerMod", "lme", "lmerTest", "lmerModLmerTest"))){
    
    if(verbose) cat(green(" - Influential points \n"))
    
    inf_obs <- hlm_influence(mod, data=DF)
    inf_obs <- cbind(inf_obs, Timepoint=DF[as.numeric(inf_obs$id), "Timepoint"])
    inf_obs$Part_Time <- paste(inf_obs$Participant_ID, inf_obs$Timepoint, sep="_")
    inf_part <- hlm_influence(mod, level="Participant_ID")
    
    cutoff_obs <- internal_cutoff_mod(inf_obs$cooksd, "cooks.distance")
    cutoff_part <- internal_cutoff_mod(inf_part$cooksd, "cooks.distance")
    if(nrow(inf_obs)>0) inf_obs_plot <- dotplot_diag(inf_obs$cooksd, name = "cooks.distance", cutoff = cutoff_obs, index=inf_obs$Part_Time) + ggtitle("Influential observations")
    if(nrow(inf_part)>0) inf_part_plot <- dotplot_diag(inf_part$cooksd, name = "cooks.distance", cutoff = cutoff_part, index=inf_part$Participant_ID) + ggtitle("Influential participant")
    
    inf_obs <- inf_obs[order(inf_obs$cooksd, decreasing=TRUE),]
    inf_obs <- inf_obs[inf_obs$cooksd>cutoff_obs,] 
    if(nrow(inf_obs)>10 ) inf_obs <- inf_obs[1:10,] 
    
    inf_part <- inf_part[order(inf_part$cooksd, decreasing=TRUE),]
    inf_part <- inf_part[inf_part$cooksd>cutoff_part,] 
    if(nrow(inf_part)>10 ) inf_part <- inf_part[1:10,] 
    
    if(print) wrap_plots(inf_obs_plot, inf_part_plot)
    if(save){
      if(nrow(inf_obs)>0 & nrow(inf_part)>0){
        ggsave(paste(path, "model_check/", saving_name, "_inf_points.png", sep=""),
               wrap_plots(inf_obs_plot, inf_part_plot), width=9, height=7)
      } else if(nrow(inf_obs)>0){
        ggsave(paste(path, "model_check/", saving_name, "_inf_points.png", sep=""),
               inf_obs_plot, width=4.5, height=7)
      }else if(nrow(inf_part)>0){
        ggsave(paste(path, "model_check/", saving_name, "_inf_points.png", sep=""),
               inf_part_plot, width=4.5, height=7)
      }
    }
  } else {
    inf_obs <- data.frame()
    inf_part <- data.frame()
  }
  
  
  
  
  ################c
  # Effect plots #
  ################c
  if(!(any(class(mod) %in% c("zinbmm", "zigmm", "MixMod", "cpglm")))){
    if(verbose) cat(green(" - Effect plots \n"))
    
    pred_effects <- predictorEffects(mod) # Predicted effects
    if(print) plot(pred_effects)
    if(save){
      png(paste(path, "model_check/", saving_name, "_pred_effect.png", sep=""),
          width=8, height=8, units="in", res=300)
      try(plot(pred_effects), T)
      dev.off()
    }
  }
  return(list(simRes=simRes, simRes_cond=simRes_cond, 
              inf_obs=inf_obs, inf_part=inf_part))
}




library(lmerTest)
library(ggResidpanel)
library(redres)
library(DHARMa)
library(effects)
library(HLMdiag)
result_mixedmod <- function(outcome, groupvars, DF, compar_var, plot=F, verbose=T, form=NULL, stopCode=F){
  
  output <- list()
  
  if(!is.null(form)) form <- paste(outcome, form, sep=" ~ ")
  if(is.null(form)) form <- paste(outcome, paste(groupvars, collapse=" + "), sep=" ~ ")
  
  # Try interaction with clusters
  # form <- paste(outcome, "B*group_info + P + (1|Participant_ID)", sep=" ~ ")
  
  if(verbose) print(paste("Model: ", form, sep=""))
  
  var_check <- groupvars[!(groupvars %in% grep("[|]", groupvars, value = TRUE))]
  DF_mod <<- DF[, c("Participant_ID", outcome, var_check)] # otherwise predictorEffects() does not want to work
  lmer_mod <- lmer(form, data=DF_mod, REML = TRUE)
  
  lm_summary <- coef(summary(lmer_mod))
  para <- as.data.frame(lm_summary)[-1, -c(3:4)]
  para$name <- rownames(lm_summary)[-1]
  para$outcomes <- outcome
  rownames(para) <- NULL
  para <- para[, c(5,4,1:3)]
  
  output$para <- para
  output$model <- lmer_mod
  output$residuals <- residuals(lmer_mod)
  output$fitted <- fitted(lmer_mod)
  output$simRes <- simulateResiduals(fittedModel = lmer_mod, plot = F) #DHARMa
  output$simRes_cond <- simulateResiduals(fittedModel = lmer_mod, plot = F, re.form = NULL)
  RE <- ranef(lmer_mod)[[1]]
  RE$Participant_ID <- rownames(RE)
  RE <- dplyr::rename(RE, RE_int = `(Intercept)`)
  output$RE <- RE
  output$ANOVA <- car::Anova(lmer_mod, type="II") # type="III"
  output$emm <- emmeans(lmer_mod, compar_var)
  output$PRS <- pairs(output$emm, infer=T)
  
  inf_obs <- hlm_influence(lmer_mod, data=DF)
  inf_obs <- cbind(inf_obs, Timepoint=DF[as.numeric(inf_obs$id), "Timepoint"])
  inf_obs$Part_Time <- paste(inf_obs$Participant_ID, inf_obs$Timepoint, sep="_")
  if(stopCode) inf_part <- tibble(cooksd=numeric(), Participant_ID=numeric())
  else  inf_part <- hlm_influence(lmer_mod, level="Participant_ID")
  
  cutoff_obs <- internal_cutoff_mod(inf_obs$cooksd, "cooks.distance")
  cutoff_part <- internal_cutoff_mod(inf_part$cooksd, "cooks.distance")
  
  if(plot){
    
    output$model <- lmer_mod 
    output$res_plot <- resid_panel(lmer_mod, "all") # Residual from ggResidpanel
    output$pred_effects <- predictorEffects(lmer_mod) # Predicted effects
    suppressWarnings(output$ranef <- plot_ranef(lmer_mod)) # Random effect distribution from redres 
    if(nrow(inf_obs)>0) output$inf_obs_plot <- dotplot_diag(inf_obs$cooksd, name = "cooks.distance", cutoff = cutoff_obs, index=inf_obs$Part_Time) + ggtitle("Influential observations")
    if(nrow(inf_part)>0) output$inf_part_plot <- dotplot_diag(inf_part$cooksd, name = "cooks.distance", cutoff = cutoff_part, index=inf_part$Participant_ID) + ggtitle("Influential participant")
  } 
  
  
  inf_obs <- inf_obs[order(inf_obs$cooksd, decreasing=TRUE),]
  inf_obs <- inf_obs[inf_obs$cooksd>cutoff_obs,] 
  if(nrow(inf_obs)>10 ) inf_obs <- inf_obs[1:10,] 
  output$inf_obs <- inf_obs
  
  
  inf_part <- inf_part[order(inf_part$cooksd, decreasing=TRUE),]
  inf_part <- inf_part[inf_part$cooksd>cutoff_part,] 
  if(nrow(inf_part)>10 ) inf_part <- inf_part[1:10,] 
  output$inf_part <- inf_part
  
  
  
  return(output)
}

# modified from HLMdiag:::internal_cutoff to be more stric
internal_cutoff_mod <- function(x, name){
  q3 <- quantile(x, p = 0.90) # initially p=0.75
  x.iqr <- IQR(x)
  cutoff <- q3 + 3 * x.iqr
  
  return(cutoff)
}

is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}





library(ggpubr) # ggtexttable()
library(rlang)
print_result_MM <- function(res, sig_threshold=0.05){
  t1 <- res$res_anov_all
  t2 <- res$emm_all
  t3 <- res$PRS_all
  
  # row = which(t1$`Pr(>Chisq)`<0.05)
  row = which(t1$q.value <sig_threshold)
  tg1 <- ggtexttable(t1, rows=NULL, theme = ttheme("minimal"))
  if(!is_empty(row)) tg1 <- table_cell_font(tg1, row=(row+1), column =c(1, 2, 6), face=  "bold") 
  tg1 <- tab_add_title(tg1, "Anova",face = "bold", size = 15)
  
  tg2 <- ggtexttable(t2, rows=NULL, theme = ttheme("minimal"))
  tg2 <- tab_add_title(tg2, "Marginal means",face = "bold", size = 15)
  
  row = which(t3$q.value <sig_threshold)
  tg3 <- ggtexttable(t3, rows=NULL, theme = ttheme("minimal"))
  if(!is_empty(row)) tg3 <- table_cell_font(tg3, row=(row+1), column =seq(ncol(t3)), face=  "bold")
  tg3 <- tab_add_title(tg3, "Pairwise differences",face = "bold", size = 15)
  
  return(list(tg1=tg1, tg2=tg2, tg3=tg3))
}




rearrange_Maaslin_res <- function(significant_results){
  
  significant_results$feature <- gsub("\\.",":", significant_results$feature)
  significant_results <- dplyr::rename(significant_results, names=metadata)
  significant_results <- dplyr::rename(significant_results, outcomes=feature)
  significant_results <- dplyr::rename(significant_results, q.value=qval)
  significant_results <- dplyr::rename(significant_results, p.value=pval)
  
  return(significant_results)
}


# sig_res already contains only significant results
# when imported from maaslin (significant_results.tsv) it correspond to all taxa with q.val < 0.25
save_plots_varsignificant <- function(sig_res, DF, savingName, tse, compar_var, second_facet=NULL, 
                                      line_diff=F, counts=F, path){
  
  if(counts) abund_values = "counts" else  abund_values = "relabundance"
  
  unit_height = 4
  if(!is.null(second_facet)){
    nb_elmnt <- length(unique(DF[, second_facet]))
    unit_height = nb_elmnt * unit_height /2
  } 
  
  # create one plot for each significant variable in "Plots_savingName/"
  pp <- plots_varsignificant(sig_res, DF, savingName, second_facet=second_facet, 
                             height=unit_height, line_diff=line_diff, path=path)
  # create one plot with all significant variables in path
  cat(green( paste("Create sigTaxa", second_facet, ".png \n", sep="")))
  ggsave_adaptSize(paste(path, "sigTaxa", second_facet,".png", sep=""), pp, unit_height=unit_height)
  
  sig_res_compar_var <- subset(sig_res,  names==compar_var)
  taxa <- sig_res_compar_var$outcomes
  
  
  # create 2 combined plots (by time and by treatment) in path for each compar_var significant variable 
  # Timepoint plot
  cat(green( paste("Create sigTaxaInter_byTime", second_facet, ".png \n", sep="")))
  pp_list <- list()
  for(i in seq(length(taxa))){
    # colData(tse)$abundance <- getAbundanceFeature(tse, feature_id = taxa[i], abund_values = abund_values)
    colData(tse)$abundance <- assay(tse, abund_values)[taxa[i],]
    pp_i <- plot_sample_measure_CO(tse, "abundance", taxa[i], trans = "psd_log", x="Timepoint", print_table=F, colorInter=T, counts=counts) 
    if(!is.null(second_facet)) pp_list[[i]] <- pp_i + facet_grid(paste(second_facet, "~ Gruppe", sep="")) 
    else pp_list[[i]] <- pp_i
  }
  ggsave_adaptSize(paste(path, "sigTaxaInter_byTime", second_facet, ".png", sep=""), pp_list, unit_width=5, unit_height=unit_height)
  
  # Treatment plot
  cat(green( paste("Create sigTaxaInter_byTreat", second_facet, ".png \n", sep="")))
  pp_list <- list()
  for(i in seq(length(taxa))){
    # colData(tse)$abundance <- getAbundanceFeature(tse, feature_id = taxa[i], abund_values = abund_values)
    colData(tse)$abundance <- assay(tse, abund_values)[taxa[i],]
    pp_i <- plot_sample_measure_CO(tse,  "abundance", taxa[i], trans = "psd_log", print_table=F, counts=counts)  
    if(!is.null(second_facet)) pp_list[[i]] <- pp_i + facet_grid(paste(second_facet, "~ Intervention", sep="")) 
    else pp_list[[i]] <- pp_i
  }
  ggsave_adaptSize(paste(path, "sigTaxaInter_byTreat", second_facet, ".png", sep=""), pp_list, unit_width=4.5, unit_height = unit_height)
}



plots_varsignificant <- function(sig_res, DF, savingName, second_facet=NULL, width=5, height=4, line_diff=F, path=NULL){
  
  if("outcomes" %in% colnames(sig_res))  sig_res <- dplyr::rename(sig_res, feature=outcomes)
  
  pp_list <- list()
  for(i in seq(nrow(sig_res))){
    variable <- sig_res[i,"names"]
    if(variable == "Alter.sc") variable <- "Alter"
    if(variable == "Timepoint:Intervention") variable <- "Timepoint"
    
    if(any(class(DF[, variable]) %in% "numeric")){
      pp <- ggplot(DF, aes(x=.data[[variable]], y=.data[[sig_res[i,"feature"]]], color=Inter_treat)) + 
        geom_point()  + scale_colour_manual(values= cols_Inter_treat) + theme_light() + 
        guides(color="none")
      if(!is.null(second_facet))  pp <- pp + facet_grid(rows= second_facet)
    } else{
      pp <-  plot_sample_measure_CO(DF, sig_res[i,"feature"], trans = NULL, line_diff = line_diff, 
                                    print_table = F, colorInter=T) 
      if(variable != "Treatment") pp <- pp + facet_grid(paste(second_facet, "~", variable, sep=""))
      if(variable == "Treatment") pp <- pp + facet_grid(paste(second_facet, "~ Intervention", sep=""))
      
    }
    
    title <- paste(sig_res[i,"feature"], "\n vs", variable)
    subtitle <- paste("q.value =",  round(sig_res[i,"q.value"], 4), " p.value =",  round(sig_res[i,"p.value"], 4))
    pp <- pp + labs(title=title, subtitle = subtitle)
    
    #if(!(variable %in% c(groupvars, "Alter"))) variable <- "Baseline"
    
    sig_res[i,"feature"] <- gsub(":", "\\.", sig_res[i,"feature"])
    
    if(!is.null(path)){
      cat(green(paste("Saving significant plots in ... ", path, "Plots_", savingName, "\n", sep="")))
      
      ggsave(paste(path, "Plots_", savingName,"/", second_facet, i, "_",  variable, "_", sig_res[i,"feature"], ".png", sep=""),
             pp, width=width, height=height)
    }
    
    pp_list[[i]] <- pp
  }
  
  return(pp_list)
}


plots_varsignificant2 <- function(sig_res, DF, savingName, width=5, height=4, counts=F, 
                                  trans="psd_log", path=NULL, interact_var=NULL, line_data=NULL,...){
  
  if("outcomes" %in% colnames(sig_res))  sig_res <- dplyr::rename(sig_res, feature=outcomes)
  # if(counts) trans <- "psd_log"
  # else trans <- ""
  
  pp_list <- list()
  interact_var0 <- interact_var
  for(i in seq(nrow(sig_res))){
    # variable <- sig_res[i,"names"] # "names" has been changed to "variable" in plot_results()
    variable <- sig_res[i,"variable"]
    if(!is.null(interact_var0)) if(interact_var0 == "Xinputs") interact_var <- sig_res[i,"Xinputs"] 
    
    pp <- plot_variable(sig_res[i,"feature"], variable, DF, trans=trans, counts=counts, interact_var=interact_var,...)
    
    if(is.null(interact_var)){
      
      if("Xinputs" %in% colnames(sig_res)){
        title <- paste(sig_res[i,"feature"], "\n vs ", variable, " (", sig_res[i, "Xinputs"], ")", sep="")
      } else  title <- paste(sig_res[i,"feature"], "\n vs", variable)
      
      if(any(class(DF[, variable]) %in% "factor")) sub_comment <- paste(" (level: ", sig_res[i,"label"], ")", sep="")
      else sub_comment <- NULL
      
    } else{
      if(interact_var0 == "Xinputs"){
        title <- paste(sig_res[i,"feature"], "in", variable, "\n vs", sig_res[i, interact_var0])     
      } else{
        title <- paste(sig_res[i,"feature"], "in", sig_res[i, interact_var0], "\n vs", variable)     
      }
      sub_comment <- paste(" (", sig_res[i,"contrast"], ")", sep="")
    }  
    
    subtitle <- paste("q.value =",  round(sig_res[i,"q.value"], 4), 
                      " p.value =",  round(sig_res[i,"p.value"], 4),
                      sub_comment, sep="")
    
    pp <- pp + labs(title=title, subtitle = subtitle) #, caption =caption )
    
    
    if(!is.null(path)){
      cat(green(paste("Saving significant plots in ... ", path, "Plots_", savingName, "\n", sep="")))
      
      
      if(is.null(interact_var)){
        if("Xinputs" %in% colnames(sig_res)) name_plot <- paste(i, "_",  variable, "_", sig_res[i,"feature"], "_", sig_res[i,"Xinputs"], sep="")
        else name_plot <- paste(i, "_",  variable, "_", sig_res[i,"feature"], sep="")
      } 
      else name_plot <- paste(i, "_",  variable, "_", interact_var, "_",sig_res[i,"feature"], sep="")
      if(nchar(name_plot) > 28) name_plot <- substr(name_plot, 1, 28)
      
      ggsave(paste(path, "Plots_", savingName,"/", name_plot, ".png", sep=""),
             pp, width=width, height=height)
    }
    
    pp_list[[i]] <- pp
  }
  
  return(pp_list)
}


plot_variable <- function(feature, variable, DF, trans="psd_log", counts=F, title=NULL, 
                          subtitle=NULL, interact_var=NULL, group=NULL,...){
  
  if(!("Inter_treat" %in% colnames(DF))) Inter_treat = NULL
  if(class(DF[, interact_var]) == "factor") group = interact_var else group=NULL
  
  if(any(class(DF[, variable]) %in% "numeric")){
    pp <- ggplot(DF, aes(x=.data[[variable]], y=.data[[feature]], color=Inter_treat)) + 
      geom_point()  + scale_colour_manual(values= cols_Inter_treat) + theme_light() + 
      guides(color="none")
    
    if(!is.null(trans)){
      if(trans == "psd_log"){
        pp <- pp + scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) 
        # compute stat AFTER axis transformation
        pp <- pp + ylab("psd_log10(y)")
      }
    }
  } else{
    
    if(variable == "B"){
      if(!is.null(interact_var) & is.null(group)){
        caption  <- "lines do not correspond to the fitted model"
        pp <-ggplot(DF, aes(x=.data[[interact_var]], y=.data[[feature]], color=Treatment)) + 
          geom_point()  + theme_minimal() + geom_line(aes(group=Participant_ID), color="gray") + geom_smooth(method=lm) + 
          labs(caption=caption)
      }else{
        DF_mod <- DF
        DF_mod$Participant_ID <- paste(DF_mod$Participant_ID, DF_mod$P )
        pp <- plot_sample_measure_CO(DF_mod, feature, x=variable, colorInter=F, 
                                     print_table=F, counts=counts, trans=trans, color_lines = "Intervention", ...) +
          scale_x_discrete(limits = c("Baseline", "After"))
        if(!is.null(group)){
          pp <- pp + facet_grid( ~  get(group))
        }
      }
    }
    if(variable == "P"){
      pp <- plot_sample_measure_CO(DF, feature, x=variable, colorInter=T, 
                                   print_table=F, line_diff=F, counts=counts, trans=trans) 
    }
    if(variable == "C1"){
      DF_mod <- subset(DF, Gruppe=="A")
      pp <- plot_sample_measure_CO(DF_mod, feature, x="Timepoint", colorInter=T, 
                                   print_table=F, counts=counts, trans=trans, ...)  
    }
    if(variable == "C2"){
      DF_mod <- subset(DF, Gruppe=="B")
      pp <- plot_sample_measure_CO(DF_mod, feature, x="Timepoint", colorInter=T, 
                                   print_table=F, counts=counts, trans=trans, ...)  
    }
    # Old code, use for analysis by group yoghurt oatmeal study
    # if(variable == "T"){
    #   if("group_info" %in% colnames(DF)){
    #     DF2 <- reshape2::dcast(DF,  Participant_ID + group_info ~ Inter_treat, value.var = feature)
    #     DF2b <- DF2[,c("Participant_ID", "group_info", "Before.Yogurt", "After.Yogurt")]
    #     DF2c <- DF2[,c("Participant_ID", "group_info", "Before.YogurtOatmeal", "After.YogurtOatmeal")]
    #   } else {
    #     DF2 <- reshape2::dcast(DF,  Participant_ID  ~ Inter_treat, value.var = feature)
    #     DF2b <- DF2[,c("Participant_ID", "Before.Yogurt", "After.Yogurt")]
    #     DF2c <- DF2[,c("Participant_ID", "Before.YogurtOatmeal", "After.YogurtOatmeal")]
    #   }
    #   DF2b$Diff_val <- DF2b$After.Yogurt - DF2b$Before.Yogurt
    #   DF2b$Diff_name <- "Yogurt-Baseline"
    #   DF2c$Diff_val <- DF2$After.YogurtOatmeal  - DF2$Before.YogurtOatmeal
    #   DF2c$Diff_name <- "YogurtOatmeal-Baseline"
    #   if("group_info" %in% colnames(DF)){
    #     var_toKeep <- c(c("Participant_ID", "group_info", "Diff_val", "Diff_name"))
    #   } else{
    #     var_toKeep <- c(c("Participant_ID", "Diff_val", "Diff_name"))
    #   }
    #   DF2 <- rbind(DF2b[, var_toKeep], DF2c[, var_toKeep])
    #   DF2$Diff_name <- factor(DF2$Diff_name)
    #   DF2 <- subset(DF2, !is.na(Diff_val))
    # 
    #   pp <- plot_sample_measure(DF2, "Diff_val", x="Diff_name", colorInter=F, counts=counts)+
    #     geom_hline(yintercept=0, color = "darkred", linetype="dashed", linewidth=1, ...)
    # }
    
    if(variable == "T"){
      inter1 <- as.character(unique(subset(DF, Gruppe=="A" & Timepoint=="T1")$Intervention))
      inter2 <- as.character(unique(subset(DF, Gruppe=="B" & Timepoint=="T1")$Intervention))
      treat0 <- as.character(unique(subset(DF, Timepoint=="T1")$Treatment))
      treat1 <- as.character(unique(subset(DF, Timepoint=="T2")$Treatment))
      treat0_inter1 <- paste(treat0, inter1, sep=".")
      treat0_inter2 <- paste(treat0, inter2, sep=".")
      treat1_inter1 <- paste(treat1, inter1, sep=".")
      treat1_inter2 <- paste(treat1, inter2, sep=".")
      
      var_sup <- NULL
      if(!is.null(group)) var_sup <- group
      if(!is.null(interact_var)) var_sup <- interact_var
      
      if(!is.null(var_sup)){
        formula <- as.formula(paste("Participant_ID +", var_sup, " ~ Inter_treat"))
        DF2 <- reshape2::dcast(DF, formula, value.var = feature)
        DF2b <- DF2[,c("Participant_ID", var_sup, treat0_inter1, treat1_inter1)]
        DF2c <- DF2[,c("Participant_ID", var_sup, treat0_inter2, treat1_inter2)]
        
      }else{
        DF2 <- reshape2::dcast(DF,  Participant_ID  ~ Inter_treat, value.var = feature)
        DF2b <- DF2[,c("Participant_ID", treat0_inter1, treat1_inter1)]
        DF2c <- DF2[,c("Participant_ID", treat0_inter2, treat1_inter2)]
      }
      
      DF2b$Diff_val <- DF2b[,treat1_inter1] - DF2b[,treat0_inter1]
      DF2b$Diff_name <- paste(inter1, treat0, sep="-")
      DF2c$Diff_val <- DF2c[,treat1_inter2] - DF2c[,treat0_inter2]
      DF2c$Diff_name <- paste(inter2, treat0, sep="-")
      
      if(!is.null(var_sup)){
        var_toKeep <- c(c("Participant_ID", var_sup, "Diff_val", "Diff_name"))
      } else{
        var_toKeep <- c(c("Participant_ID", "Diff_val", "Diff_name"))
      }
      
      DF2 <- rbind(DF2b[, var_toKeep], DF2c[, var_toKeep])
      DF2$Diff_name <- factor(DF2$Diff_name)
      DF2 <- subset(DF2, !is.na(Diff_val)) 
      
      if(!is.null(interact_var) & is.null(group)){
        caption  <- "lines do not correspond to the fitted model"
        pp <-ggplot(DF2, aes(x=.data[[interact_var]], y=Diff_val, color=Diff_name)) + 
          geom_point()  + theme_minimal() + geom_line(aes(group=Participant_ID), color="gray") + geom_smooth(method=lm) + 
          labs(caption=caption)
      }else {
        pp <- plot_sample_measure(DF2, "Diff_val", x="Diff_name", colorInter=F, counts=counts)+ 
          geom_hline(yintercept=0, color = "darkred", linetype="dashed", linewidth=1, ...)
        
        if(!is.null(group)){
          pp <- pp + facet_grid( ~  get(group))+
            scale_x_discrete(labels=function(x) abbreviate(x))
        }
      }
    }
    
    if(variable=="Treatment"){
      pp <-  plot_sample_measure_CO(DF, feature,x="Treatment", counts=counts, trans = trans, line_diff = T, 
                                    print_table = F, colorInter=T, ...) 
      
      if(!is.null(group)){
        pp <- pp + facet_grid( ~  get(group))
      }
    }
    
    if(variable=="Timepoint"){
      
      pp <-  plot_sample_measure_CO(DF, feature,x="Timepoint", counts=counts, trans = trans, line_diff = T, 
                                    print_table = F, colorInter=T, ...) 
      
      if(!is.null(group)){
        pp <- pp + facet_grid( ~  get(group))
      }
      
    }
    
    if(!(variable %in% c("B", "P", "C1", "C2", "T", "Treatment", "Timepoint"))){
      if(!is.null(interact_var) & is.null(group)){
        caption  <- "lines do not correspond to the fitted model"
        pp <-ggplot(DF, aes(x=.data[[interact_var]], y=.data[[feature]], color=.data[[variable]])) +
          geom_point()  + theme_minimal() + geom_smooth(method=lm) +
          labs(caption=caption)
      } else{
        pp <- plot_sample_measure_CO(DF, feature, x=variable, colorInter=T, 
                                     print_table=F, line_diff=F, counts=counts, trans=trans) 
        if(!is.null(group)){
          pp <- pp + facet_grid( ~  get(group))
        }
      }
    }
  }
  
  
  
  
  if(is.null(title)) title <- paste(feature, "\n vs", variable)
  pp <- pp + labs(title=title, subtitle = subtitle)
  
  return(pp)
}



plot_taxa_variable <- function(taxa, DF, trans="psd_log", counts=F){
  if(counts) title <- "Counts"
  else title <- 'Relative abundance'
  
  ll_plot <- list()
  for(var in c("Treatment", "T", "P", "C1", "C2", "Geschlecht", "Alter")){
    ll_plot[[var]] <- plot_variable(taxa, var, DF, trans=trans, counts=counts, title=var, subtitle = title)
  }
  ll_plot[[1]] <- ll_plot[[1]]+ facet_grid(. ~  Intervention)
  wrap_plots(ll_plot) + plot_layout(guides = "collect") #+ plot_annotation(title = title)
}


ggsave_adaptSize <- function(filepath, list_plots, unit_width=5, unit_height=4, guides_collect='collect'){
  
  nb_plot <- length(list_plots)
  ncol <- ceiling(sqrt(nb_plot))
  nrow <- ceiling(nb_plot/ncol)
  ggsave(filepath, 
         wrap_plots(list_plots) + plot_layout(guides = guides_collect), 
         width=ncol*unit_width, height = nrow*unit_height, limitsize = FALSE)
}




library(crayon)
extract_Rosner_res <- function(file, q.value_tresh=1, p.value_tresh=1, name_levelsig, group_name=NULL){
  res <- read.xlsx(file=file, sheetName="PRS_all")
  if("outcomes" %in% colnames(res)) res <- dplyr::rename(res, feature=outcomes)
  res <- subset(res, q.value < q.value_tresh)
  res <- subset(res, p.value < p.value_tresh)
  #if(nrow(res)>50) res <- res[1:50,]
  suppressWarnings(if(any(res$estimate)==0)cat(red("at least one estimate = 0 and will not be shown in the heatmap\n")))
  res$levelsig <- -log(res$q.value) * sign(res$estimate)
  
  if(!is.null(group_name)) res <- res[, c("feature", group_name, "levelsig", "estimate", "q.value")]
  else res <- res[, c("feature", "levelsig", "estimate", "q.value")]
  
  res <- dplyr::rename(res, !!name_levelsig:=levelsig)
  res <- dplyr::rename(res, !!paste(name_levelsig, "est", sep="."):=estimate)
  res <- dplyr::rename(res, !!paste(name_levelsig, "q.val", sep="."):=q.value)
  
  # res <- res[, c("feature", "levelsig", "estimate", "q.value")]
  # colnames(res) <- c("feature", name_levelsig, paste(name_levelsig, c("est", "q.val"), sep="."))
  return(res)
}





library("xlsx")
print_Rosner_model <- function(path_analysis, model, q.value_tresh=1, p.value_tresh=1,
                               inter1, inter2, additional_name=NULL, name=model, print=T,  
                               order=F, col_rot=0, dif_full="Y-YO", group_name=NULL, 
                               custom_order=NULL, remove_patterns = NULL, new_patterns=NULL, 
                               change_feature_name=NULL, toKeep=NULL){
  
  
  dif_inter1 <- paste(inter1, "-B", sep="")
  dif_inter2 <- paste(inter2, "-B", sep="")
  

  file_full <- paste(path_analysis , model, "_full", additional_name, "_results.xlsx", sep="")
  file_inter1 <- paste(path_analysis , model, "_", inter1, additional_name, "_results.xlsx", sep="")
  file_inter2 <- paste(path_analysis , model, "_", inter2, additional_name, "_results.xlsx", sep="")

  
  # Use results from pairwise comparison in PRS_all
  
  if(file.exists(file_full)){
    res <- extract_Rosner_res(file_full, q.value_tresh, p.value_tresh, dif_full, group_name)
    res_after_prt <- res[, colnames(res)!=dif_full]
    res_after <- res[, -tail(seq(ncol(res)),2)]
  }else{
    res_after <- data.frame(feature=NA, level_sig=NA)
    res_after <- dplyr::rename(res_after, !!dif_full:=level_sig)
    res_after_prt <- data.frame(feature=NA, est=NA, q.val=NA)
    res_after_prt <- dplyr::rename(res_after_prt, !!paste(dif_full, ".est", sep="") :=est)
    res_after_prt <- dplyr::rename(res_after_prt, !!paste(dif_full, ".q.val", sep=""):=q.val)
  }
  
  if(file.exists(file_inter1)){
    res <- extract_Rosner_res(file_inter1, q.value_tresh, p.value_tresh, dif_inter1, group_name)
    res_Int1_prt <- res[, colnames(res)!=dif_inter1]
    res_Int1 <- res[, -tail(seq(ncol(res)),2)]
  }else{
    res_Int1 <- data.frame(feature=NA, level_sig=NA)
    res_Int1 <- dplyr::rename(res_Int1, !!dif_inter1:=level_sig)
    res_Int1_prt <- data.frame(feature=NA, est=NA, q.val=NA)
    res_Int1_prt <- dplyr::rename(res_Int1_prt, !!paste(inter1, "-B.est", sep=""):=est)
    res_Int1_prt <- dplyr::rename(res_Int1_prt, !!paste(inter1, "-B.q.val", sep=""):=q.val)
  }
  
  if(file.exists(file_inter2)){
    res <- extract_Rosner_res(file_inter2, q.value_tresh, p.value_tresh, dif_inter2, group_name)
    res_Int2_prt <- res[, colnames(res)!=dif_inter2]
    res_Int2 <- res[, -tail(seq(ncol(res)),2)]
  }else{
    res_Int2 <- data.frame(feature=NA, level_sig=NA)
    res_Int2 <- dplyr::rename(res_Int2, !!dif_inter2:=level_sig)
    res_Int2_prt <- data.frame(feature=NA, est=NA, q.val=NA)
    res_Int2_prt <- dplyr::rename(res_Int2_prt, !!paste(inter2, "-B.est", sep=""):=est)
    res_Int2_prt <- dplyr::rename(res_Int2_prt, !!paste(inter2, "-B.q.val", sep=""):=q.val)
  }
  
  
  res <- merge(res_after, res_Int1, by=c("feature", group_name), all=T, sort=F)
  res <- merge(res, res_Int2, by=c("feature", group_name), all=T, sort=F)
  res <- subset(res, !is.na(feature))
  #res[is.na(res)] <- 0
  
  if(order) res <- res[order(res$feature),]
  if(!is.null(custom_order)) res <- arrange(res, match(feature, custom_order))
  if(!is.null(change_feature_name)) res[which(res$feature %in% names(change_feature_name)), "feature"] <- change_feature_name
  if(!is.null(remove_patterns)){
    nb_patterns <- length(remove_patterns)
    if(is.null(new_patterns)) new_patterns <- rep(" ", nb_patterns)
    for(i in seq(nb_patterns)){
      res$feature <- gsub(remove_patterns[i], new_patterns[i], res$feature )
    }
  }
  if(!is.null(toKeep)){
    res <- subset(res, feature %in% toKeep)
  }
  
  res_prt <- merge(res_after_prt, res_Int1_prt, by=c("feature", group_name), all=T, sort=F)
  res_prt <- merge(res_prt, res_Int2_prt, by=c("feature", group_name), all=T, sort=F)
  res_prt <- subset(res_prt, !is.na(feature))
  res_prt$model <- model
  if(order) res_prt <- res_prt[order(res_prt$feature),]
  if(!is.null(custom_order)) res_prt <- arrange(res_prt, match(feature, custom_order))
  if(!is.null(change_feature_name)) res_prt[which(res_prt$feature %in% names(change_feature_name)), "feature"] <- change_feature_name
  if(!is.null(remove_patterns)){
    nb_patterns <- length(remove_patterns)
    if(is.null(new_patterns)) new_patterns <- rep(" ", nb_patterns)
    for(i in seq(nb_patterns)){
      res_prt$feature <- gsub(remove_patterns[i], new_patterns[i], res_prt$feature )
    }
  }
  if(!is.null(toKeep)){
    res_prt <- subset(res_prt, feature %in% toKeep)
  }
  
  
  if(!is.null(group_name)){
    
    res <- split(res, res[,group_name])
    res <- lapply(res, function(x){
      rownames(x) <- x[,"feature"]
      x["feature"] <- NULL
      x[group_name] <- NULL
      ; x})
    
    res_prt <- split(res_prt, res_prt[,group_name])
    res_prt <- lapply(res_prt, function(x){
      rownames(x) <- x[,"feature"]
      x["feature"] <- NULL
      x[group_name] <- NULL
      ; x})
    
    HM <- list()
    for(i in seq(length(res))){
      name <- names(res[i])
      res_i <- res[[i]]
      res_prt_i <- res_prt[[i]]
      res.qval_i <- res_prt_i[, paste(c(dif_full, dif_inter1, dif_inter2), ".q.val", sep="")]
      res.qval_i[is.na(res.qval_i)] <- 1
      
      HM[[i]] <- list(create_heatmap(res_prt_i, res_i, res.qval_i, name, col_rot))
      names(HM[[i]]) <- name
    }
    HM <- unlist(HM)
  }else{
    rownames(res) <- res$feature
    res$feature <-NULL
    
    rownames(res_prt) <- res_prt$feature
    res_prt$feature <-NULL
    
    res.qval <- res_prt[, paste(c(dif_full, dif_inter1, dif_inter2), ".q.val", sep="")]
    res.qval[is.na(res.qval)] <- 1
    HM <-create_heatmap(res_prt, res, res.qval, name, col_rot)
  }
  
  
  # if(print) print(res_prt)
  if(print) return(list(HM=HM, res_prt=res_prt, res=res))
  else  return(HM)
}


create_heatmap <- function(res_prt, res, res.qval, name, col_rot){
  
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
  #return(list(HM=HM))
  return(HM)
}


