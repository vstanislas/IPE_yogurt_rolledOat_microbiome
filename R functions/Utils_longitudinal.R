
# the following variables are needed in obj:
# Participant_ID
# Treatment or Timepoint when using line_diff
library(ggnewscale)
plot_sample_measure <- function(DF, measure, subtitle="", x= "Treatment", colorInter=F, trans = NULL, 
                                lineMean=F, printMean=F, ind_filter=NULL, line_diff=T, 
                                linewidth_diff=1.4, ylim=NULL, counts=F, color_lines=NULL, ...){
  # scale_log = T replace by trans="log"
  guides_col <- "none"
  if(!("Inter_treat" %in% colnames(DF))) Inter_treat = NULL
  
  
  if(colorInter) cols <- cols_Inter_treat
  
  nb_ind <- length(unique(DF$Participant_ID))
  
  ## Calculate difference between time-points or treatments
  class_meas <- class(DF[,measure])
  if(line_diff & x %in% c("Treatment", "Timepoint") & class_meas!="factor"){
    DF <- computeDiff(DF, x, measure, ...)
  } 
  if(line_diff & (!(x %in% c("Treatment", "Timepoint")) | class_meas=="factor")){
    DF$Diff <- "None"
    DF$value <- 0
  } 
  if(is.null(color_lines)) color_lines <- "Diff"

  
  ## Filtering subjects
  if(class(ind_filter) == "call") ind_filter <-  subset(DF, eval(ind_filter))$Participant_ID
  if(!is.null(ind_filter)) {
    DF <- subset(DF, Participant_ID %in% ind_filter)
    nb_ind_filt <- length(unique(DF$Participant_ID))
    perc_ind <- round(nb_ind_filt/nb_ind*100, 1)
    subtitle <- paste(subtitle, " nb subjects: ", nb_ind_filt, " (" , perc_ind, "%)", sep="")
  }
  
  
  ## Plot
  pp <- ggplot(DF, aes(.data[[x]], .data[[measure]])) +
    stat_boxplot(geom = "errorbar",color="azure4", width=0.2) 
  
  if(line_diff){
    if(color_lines=="Diff"){
      line_sub <- DF
      col_lines <- c(Decrease="#700054", Increase="#02a0c8", Unchange="#ff8d3f")
    } else {
      line_sub <- DF
      col_lines <- c("green4", "gold2")
      guides_col <- NULL
    } 
    
    pp <- pp + geom_point(size=2, color="gray70") + 
      geom_line(data= line_sub, aes(group=Participant_ID, color=.data[[color_lines]]), linewidth=linewidth_diff, alpha=0.2) + 
      scale_colour_manual(values= col_lines)+ 
     guides(color=guides_col)
    
    
    if(colorInter){
      if(counts) pp <- pp +  new_scale_color() + 
          geom_violin(aes(fill=Inter_treat), color="gray5", linewidth=0.8, width= 0.6, alpha=0.5)
      else pp <- pp +  new_scale_color() + 
          geom_boxplot(aes(fill=Inter_treat), color="gray5", width= 0.3, alpha=0.5)
      
      pp <- pp + scale_fill_manual(values= cols)  
    }  
    else {
      if(counts) pp <- pp + geom_violin(color="gray5",fill="gray90", linewidth=0.8, width= 0.6, alpha=0.5)
      else pp <- pp + geom_boxplot(color="gray5",fill="gray90", width= 0.3, alpha=0.5)  
    }  
  } else {
    if(colorInter)  {
      pp <- pp + geom_jitter(aes(color=Inter_treat), width = 0.15, height=0, size=1.2, alpha=0.7)  +
        scale_colour_manual(values= cols) 
    }
    else  pp <- pp + geom_jitter(width = 0.2, height=0, size=1.2, alpha=0.7) 
    pp <- pp+ 
      geom_violin(color="gray60", linewidth=0.7, width= 0.6, alpha=0)   + 
      geom_boxplot(color="gray40", width= 0.3, alpha=0.5, outlier.shape = NA)  
  }
  
  if(lineMean){
    pp <- pp + 
      stat_summary(fun=mean, colour="darkred", geom="point", shape=4, size=1, alpha=0.7, stroke = 1.5)+
      stat_summary(group=1, fun=mean, colour="darkred", geom="line", linewidth=linewidth_diff,  alpha=0.7)
  } 
  
  if(printMean){
    pp <- pp + 
      stat_summary(fun=mean, geom="label", aes(label=round(..y.., 1)), size=4, label.size = 0, alpha=0.7,
                   col="gray30", vjust=0,  hjust=0.5, fontface = "bold",
                   position=position_dodge(width=1), show.legend=FALSE)
  } 
  pp <- pp  +
    labs(title=measure, subtitle=subtitle) + 
    xlab("") + theme_light() + 
    guides(color=guides_col, fill="none") +
    theme(strip.background = element_rect(fill="steelblue", linewidth=1.5, linetype="solid"),
          strip.text = element_text(face="bold"))
  
  
  if(is.null(trans)) trans = ""
  
  if(trans=="psd_log"){ 
    # consideration of zero/negative values
    # pseudo_log_trans(): smoothly transition to linear scale around 0. 
    # https://win-vector.com/2012/03/01/modeling-trick-the-signed-pseudo-logarithm/
    if(printMean | lineMean){
      pp <- pp + coord_trans(y = "pseudo_log", ylim=ylim) # compute stat before axis transformation
      pp <- pp + ylab("psd_ln(y)") # natural logarithm with base=exp(1), not possible to change base through coord_trans()
    } 
    else{
      pp <- pp + scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits= ylim) # compute stat AFTER axis transformation
      pp <- pp + ylab("psd_log10(y)")
    } 
  }   
  
  if(trans=="log"){
    # No consideration of zero/negative values
    if(printMean | lineMean) pp <- pp + coord_trans(y = "log10", ylim=ylim) 
    else pp <- pp + scale_y_continuous(trans="log10", limits= ylim) # compute stat AFTER axis transformation
    pp <- pp + ylab("log10(y)")
  }   
  
  
  if(trans=="sqrt"){
    if(printMean | lineMean) pp <- pp + coord_trans(y = "sqrt", ylim=ylim) # compute stat before axis transformation
    else pp <- pp + scale_y_continuous(trans="sqrt", limits= ylim) # compute stat AFTER axis transformation
    pp <- pp + ylab("sqrt(y)")
  }   
  
  
  return(pp)
}




plot_pair_measures <- function(obj, meas1, meas2, subject_name=T, facet_grid=T, scale_log=T){
  
  
  if(class(obj) == "TreeSummarizedExperiment"){
    df <- as.data.frame(merge(t(assays(obj)[["relabundance"]]), 
                              colData(obj), by="row.names"))
  }  else df <- obj
  
  
  p <- ggplot(df, aes(.data[[meas1]], .data[[meas2]])) +
    geom_point(aes(shape = Timepoint, color=Inter_treat), size=1.5, alpha=0.7, stroke = 2) +
    scale_color_manual(values=cols_Inter_treat) +
    scale_shape_manual(values = c("T1"= 1, "T2"= 16, "T3"= 2, "T4"= 17, 3)) + theme_light()  
  
  if(scale_log){
    p <- p + 
      scale_x_continuous(trans=scales::pseudo_log_trans(base = 10))+ 
      scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))
  }
  
  if(subject_name & facet_grid) p <- p  + guides(color = "none") 
  
  
  if(subject_name) p <- p  + new_scale_color() +
    scale_color_manual(values = cols_Participant_Id) + 
    scale_fill_manual(values = cols_Participant_Id) + 
    geom_text(aes(label=Participant_ID, colour=Participant_ID), check_overlap = T,  nudge_y = 0.002) + guides(color = "none")
  
  if(facet_grid) p <- p + guides(color = "none", shape="none")+ facet_grid(Treatment ~ Intervention) 
  
  return(p)
}



