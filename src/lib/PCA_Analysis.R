PCA_Analysis <-function(Data,Norm_Meta_Data,Experiment,graph_path) {
  
  ifelse(!dir.exists(file.path(graph_path, paste(Experiment,"/",sep=""))), dir.create(file.path(graph_path, paste(Experiment,"/",sep=""))), FALSE)
  
  col_names<-colnames(Norm_Meta_Data[,-c(1,2)])
  # color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  color = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  for (var_name in col_names){
    colnames(Data)<-paste(t((Norm_Meta_Data))["Sample.ID",],
                          t((Norm_Meta_Data))["Batch",],t((Norm_Meta_Data))[var_name,],sep=' | ')
    res.pca<-PCA(Data, scale.unit = FALSE, ncp = 5, ind.sup = NULL,
                 quanti.sup = NULL, quali.sup = NULL, row.w = NULL,
                 col.w = NULL, graph = FALSE, axes = c(1,2))
    if (grepl('\\.', var_name)==TRUE){var_label<-gsub('\\.','_',var_name)}else{var_label<-var_name}
    
    graph_path_name<-paste(graph_path,"PCA_",Experiment,"_color_by_",var_label,".tiff",sep="")
    # tiff(file=graph_path_name, width = 1500, height =850, units = "px",res = 90)
    
    n<-length(unique(as.character(as.data.frame(Norm_Meta_Data)[[var_name]])))
    if(n>10){col=sample(color, n)}else{col=c( "#80000000", "#56B4E9","#00AFBB", "#E7B800", "#FC4E07","#5D6D7E","#33FFBB","#000000","#E333FF","#52FF33")}
    
    fviz_pca_var(res.pca,
                 habillage=as.data.frame(Norm_Meta_Data)[[var_name]],
                 # col.var = "cos2", # Color by contributions to the PC
                 # palette = c( "#80000000", "#56B4E9","#00AFBB", "#E7B800", "#FC4E07","#5D6D7E","#33FFBB","#000000","#E333FF","#52FF33"),
                 palette = col,
                 # gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE,     # Avoid text overlapping
                 geom = c("text","point"),
                 title=paste("PCA of Variables with Sample ID|Batch and ",var_name,sep="")
                
                 
    ) +  theme_gray(base_size=14)
    # dev.off()
    ggsave(graph_path_name,device='tiff',
           width = 17, height = 10, dpi = 350, units = "in",limitsize = FALSE)

  }
  
}