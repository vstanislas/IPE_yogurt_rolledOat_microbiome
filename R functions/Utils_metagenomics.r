


# Function for low abondance OTUs < 0.5%
# Defining OTU filtering function
phyloseq_CleanUp<-function(PHYLOSEQDATA, thresh_abond=0.5, thresh_prev=0.5, norm=FALSE) {
  dim_OTU_0 <- dim(otu_table(PHYLOSEQDATA))
  
  # Delete OTU with no sample
  ps = prune_taxa(taxa_sums(PHYLOSEQDATA) > 0, PHYLOSEQDATA) 
  otu <- otu_table(ps)
  Filter1 <- dim_OTU_0 - dim(otu)
  
  # Transform the OTU table into relative abundance
  if(norm) otu.rel <- prop.table(otu, margin = 2)*100 else otu.rel <- otu
  
  # Abundance filtering: 
  # Delete OTU that correspond to less than thresh_abond % of the sample abundance over all samples. 
  keep <- apply(otu.rel, 1, max) >= thresh_abond
  clean_ps<-prune_taxa(keep, ps) 
  # Filter2 <- dim(otu) - dim(otu_table(clean_ps))
  Filter2 <- c(nrow(otu.rel) - sum(keep), 0)
  
  # Prevalence filtering:
  keep <- prevalence(otu.rel) >= thresh_prev
  clean_ps<-prune_taxa(names(keep[keep]), clean_ps) 
  Filter3 <- c(nrow(otu.rel) - sum(keep), 0)
  
  # Change NA in the taxonomy table as "Ambiguous_taxa"  
  tax_table(clean_ps)[is.na(tax_table(clean_ps))] <- 'Ambiguous_taxa'
  
  # Keep only OTU from Bacteria kingdom, which are not Mitochondria, chloroplast or from ambiguous taxa (NA)
  clean_taxatype=subset_taxa(ps, Kingdom == "Bacteria" & 
                               Family != "Mitochondria" & 
                               Class != "Chloroplast" &
                               Order != "Ambiguous_taxa")
  Filter4 <- dim(otu_table(ps)) - dim(otu_table(clean_taxatype))
  keep <- intersect(taxa_names(clean_ps), taxa_names(clean_taxatype))
  clean_subset <- prune_taxa(keep, ps) 
  
  dim_OTU_fin <- dim(otu_table(clean_subset))
  
  
  df_filt <- data.frame(dim_OTU_0, Filter1, Filter2, Filter3, Filter4, dim_OTU_fin)
  rownames(df_filt) <- c("OTU", "Samples")
  colnames(df_filt) <- c("dim_data", "del_noSample", paste("del_abond<", thresh_abond, "%", sep=""), 
                         paste("del_prev<", thresh_prev, "%", sep=""), "del_taxaType", "dim_dataClean")
  
  deleted_taxa = prune_taxa(!(taxa_names(PHYLOSEQDATA) %in% taxa_names(clean_subset)), PHYLOSEQDATA) 
  
  
  return(list(clean_subset =clean_subset, df_filt=df_filt, deleted_taxa=deleted_taxa))
}


# From ArchR package (https://rdrr.io/github/GreenleafLab/ArchR/src/R/ColorPalettes.R)
whitePurple = c("9"='#f7fcfd',"6"='#e0ecf4',"8"='#bfd3e6',"5"='#9ebcda',"2"='#8c96c6',"4"='#8c6bb1', "7"='#88419d',"3"='#810f7c',"1"='#4d004b')
solarExtra = c("5"='#3361A5', "7"='#248AF3', "1"='#14B3FF', "8"='#88CEEF', "9"='#C1D5DC', "4"='#EAD397', "3"='#FDB31A',"2"= '#E42A2A', "6"='#A31D1D')  #buencolors
# mycolors <- colorRampPalette(whitePurple)
mycolors <- colorRampPalette(c("#f7fcfd", '#9ebcda', '#8c6bb1', "#810f7c", "#F81894"))


plot_ordination <- function(tse, ordi_name, sd_threshold=NULL, rank_dominance="Family", 
                            genus_id=NULL, functional_id=NULL, data_type="Taxonomy", 
                            assay="relabundance", top_dom=T){
  
  if(data_type=="Taxonomy"){
   # tse <- addPerSampleDominantTaxa(tse, rank = rank_dominance, name="dominant_feature", assay_name=assay)
    tse <- addDominant(tse, group = rank_dominance, name="dominant_feature", assay.type=assay)
    
    if(!is.null(genus_id)){
      tse_gen <- agglomerateByRank(tse, rank="Genus")
      # colData(tse)$abundance  <- getAbundanceFeature(tse_gen, feature_id = genus_id, assay.type = assay)
      colData(tse)$abundance <- assay(tse_gen, assay)[genus_id,]
      title_plot <- genus_id
    }
  }
  
  if(data_type=="Functional"){
    
    if(top_dom){
      mat <- assay(tse, assay.type=assay) # was relabundance
      
      # find dominant functional feature
      idx <- apply(t(mat) == colMaxs(mat),1L, function(x) which(x)[1]) #in case of doublon take the first maximum id
      dom_path <- rownames(mat)[idx]
      names(dom_path) <- rep( names(idx), times = lengths(idx) )
      
      # find second dominant functional feature
      for(i in seq(length(idx))){
        mat[idx[i], names(idx[i])] = 0
      }
      idx <- apply(t(mat) == colMaxs(mat),1L, function(x) which(x)[1])
      dom_path_2 <- rownames(mat)[idx]
      names(dom_path_2) <- rep( names(idx), times = lengths(idx) )
      
      tse$dominant_feature <- dom_path
      tse$second_dominant_feature <- dom_path_2
    }
    
    if(!is.null(functional_id)){
      # colData(tse)$abundance  <- getAbundanceFeature(tse, feature_id = functional_id, assay.type = assay) # was CPM
      colData(tse)$abundance <- assay(tse, assay)[functional_id,]
      title_plot <- functional_id
    }
  }
  
  colData(tse)$Total_abundance  <- colSums(assays(tse)[[assay]]) # was first assay [[1]] considered
  
  df <- as.data.frame(cbind(reducedDim(tse, ordi_name),colData(tse)))
  
  if(!all(colnames(df)[1:2] == c("V1", "V2"))){
    if(grepl("PCA", ordi_name)) df <- dplyr::rename(df, V1 = PC1, V2 = PC2)
    if(grepl("TSNE", ordi_name)) df <- dplyr::rename(df, V1 = TSNE1, V2 = TSNE2) # needed after changing R version
    if(grepl("UMAP", ordi_name)) df <- dplyr::rename(df, V1 = UMAP1, V2 = UMAP2) # needed after changing R version
  }
  
  Participants_df <- subset(df, Timepoint=="T2")
  Participants_df$V2 <- Participants_df$V2 - (IQR(Participants_df$V2)/15)
  
  if(!is.null(sd_threshold)){
    Participants_toPrint <- df %>% group_by(Participant_ID) %>% dplyr::summarise(V1 = sd(V1), V2 = sd(V2))
    Participants_toPrint <- subset(Participants_toPrint, V1 >sd_threshold | V2>sd_threshold)$Participant_ID
    Participants_df_sub <- subset(Participants_df, Participant_ID %in% Participants_toPrint)
  }
  
  p_total_ab <-  ggplot(df, aes(x=V1, y=V2, shape = Timepoint, color=Total_abundance)) + 
    geom_point(size=2, alpha=0.7, stroke = 2 ) +  
    scale_shape_manual(values = c(1, 16, 2, 17, 3))+ 
    scale_color_gradientn(colours = mycolors(20)) +
    ggtitle(ordi_name) +
    theme_classic()
  
  p_gruppe <-  ggplot(df, aes(x=V1, y=V2, shape = Timepoint, color=Gruppe)) + 
    geom_point(size=2, alpha=0.7, stroke = 2 ) +  
    scale_shape_manual(values = c(1, 16, 2, 17, 3))+ 
    ggtitle(ordi_name) +
    stat_ellipse(aes(shape=NULL))+
    theme_classic()
  
  if(top_dom){
    p_topdom <-  ggplot(df, aes(x=V1, y=V2, shape = Timepoint, color=dominant_feature)) + 
      geom_point(size=2, alpha=0.7, stroke = 2 ) +  
      scale_shape_manual(values = c(1, 16, 2, 17, 3))+ 
      scale_color_manual(values=cols) + 
      ggtitle(ordi_name) +
      #stat_ellipse(aes(shape=NULL))+ #comment to remove circle
      theme_classic()
  }else p_topdom <- NULL
  
  if(data_type == "Functional" & top_dom){
    p_2ndtopdom <-  ggplot(df, aes(x=V1, y=V2, shape = Timepoint, color=second_dominant_feature)) + 
      geom_point(size=2, alpha=0.7, stroke = 2 ) +  
      scale_shape_manual(values = c(1, 16, 2, 17, 3))+ 
      scale_color_manual(values=cols) + 
      ggtitle(ordi_name) +
      stat_ellipse(aes(shape=NULL))+ #comment to remove circle
      theme_classic()
  } else p_2ndtopdom <- NULL
  
  if(!is.null(genus_id) | !is.null(functional_id)){
    p_1spec <-  ggplot(df, aes(x=V1, y=V2, shape = Timepoint, color=abundance)) +
      geom_point(size=2, alpha=0.7, stroke = 2 ) +
      scale_shape_manual(values = c(1, 16, 2, 17, 3))+
      scale_color_gradientn(colours = mycolors(20)) +
      ggtitle(ordi_name, title_plot) +
      # stat_ellipse(aes(shape=NULL))+  #comment to remove circle
      theme_classic()
  } else p_1spec <- NULL
  
  
  
  p_inter_treat <-  ggplot(df, aes(x=V1, y=V2, shape = Timepoint, color=Inter_treat)) + 
    geom_point(size=2, alpha=0.7, stroke = 2 ) +  
    scale_shape_manual(values = c(1, 16, 2, 17, 3))+ 
    scale_color_manual(values=cols_Inter_treat) + 
    ggtitle(ordi_name) +
    stat_ellipse(aes(shape=NULL))+
    theme_classic()
  
  
  
  p_part <- ggplot(df, aes(x=V1, y=V2)) +
    geom_point(aes(shape = Timepoint, color=Inter_treat), size=2, alpha=0.7, stroke = 2) +  
    scale_color_manual(values=cols_Inter_treat)+
    scale_shape_manual(values = c(1, 16, 2, 17, 3)) +
    new_scale_color() +
    # geom_polygon(aes(fill=Participant_ID), alpha=0.3, show.legend=F) +
    geom_text(data=Participants_df, aes(label=Participant_ID, colour=Participant_ID), check_overlap = F,  fontface = "bold")+
    scale_color_manual(values = cols_Participant_Id) + scale_fill_manual(values = cols_Participant_Id) +
    guides(color = "none") + 
    ggtitle(ordi_name) + theme_classic() 
  
  if(max(table(df$Participant_ID)) < 3) p_part <- p_part + geom_line(aes(group=Participant_ID, color=Participant_ID), show.legend=F)
  else p_part <-   p_part + geom_polygon(aes(fill=Participant_ID), alpha=0.3, show.legend=F)
  
  
  if(!is.null(sd_threshold)){
    p_hv_part <- ggplot(subset(df, Participant_ID %in% Participants_toPrint), aes(x=V1, y=V2)) +
      geom_point(aes(shape = Timepoint, color=Participant_ID), size=2, alpha=0.7, stroke = 1) +  
      scale_shape_manual(values = c(1, 16, 2, 17, 3)) +
      # geom_polygon(aes(fill=Participant_ID), alpha=0.3, show.legend=F)  +
      geom_text(data=Participants_df_sub, aes(label=Participant_ID, colour=Participant_ID), check_overlap = F,  fontface = "bold")+
      scale_color_manual(values = cols_Participant_Id) + scale_fill_manual(values = cols_Participant_Id)+
      guides(color = "none") + 
      ggtitle(ordi_name) + theme_classic() 
    
    if(max(table(df$Participant_ID)) < 3) p_hv_part <-  p_hv_part + geom_line(aes(group=Participant_ID, color=Participant_ID), show.legend=F)
    else p_hv_part <-   p_hv_part + geom_polygon(aes(fill=Participant_ID), alpha=0.3, show.legend=F)
    
    Participants_high_sd <- data.frame(Participant_ID = as.character(Participants_toPrint), ordi=ordi_name)
    
  } else {
    p_hv_part <- NULL
    Participants_high_sd <- NULL
  }
  
  return(list(p_gruppe=p_gruppe, p_topdom=p_topdom, p_2ndtopdom=p_2ndtopdom, 
              p_1spec=p_1spec, p_inter_treat=p_inter_treat, p_part=p_part, 
              p_hv_part=p_hv_part, p_total_ab=p_total_ab, 
              Participants_high_sd = Participants_high_sd))
}



