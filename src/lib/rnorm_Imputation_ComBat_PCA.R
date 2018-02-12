rnorm_Imputation_ComBat_PCA <-function(Missing_Data_Fill,Norm_Meta_Data,graph_path) {
  
  col_names<-colnames(Norm_Meta_Data[,-c(1,2)])
  
  for (var_name in col_names){
    colnames(ComBat_Results)<-paste(t((Norm_Meta_Data))["Sample.ID",],
                                       t((Norm_Meta_Data))["Batch",],t((Norm_Meta_Data))[var_name,],sep=' | ')
    res.pca<-PCA(Missing_Data_Fill, scale.unit = FALSE, ncp = 5, ind.sup = NULL,
                 quanti.sup = NULL, quali.sup = NULL, row.w = NULL,
                 col.w = NULL, graph = FALSE, axes = c(1,2))
    graph_path_name<-paste(graph_path,"ComBat_Clustering_color_by_",var_name,".tiff",sep="")
    tiff(file=graph_path_name, width = 1500, height =850, units = "px",res = 90)
    fviz_pca_var(res.pca,
                 habillage=as.data.frame(Norm_Meta_Data)$Family_group,
                 #col.var = "cos2", # Color by contributions to the PC
                 palette = c( "#80000000", "#56B4E9","#00AFBB", "#E7B800", "#FC4E07","#5D6D7E","#33FFBB","#000000","#E333FF","#52FF33"),
                 # gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE,     # Avoid text overlapping
                 geom = c("text","point"),
                 title=paste("PCA of Variables with Sample ID|Batch and ",var_name,sep="")
                 
    )
    dev.off()
    
    
    colnames(ComBat_Results)<-paste(t((Norm_Meta_Data))["Sample.ID",],t((Norm_Meta_Data))["Batch",],t((Norm_Meta_Data))["Family_group",],sep=' | ')
    
    res.pca<-PCA(ComBat_Results, scale.unit = FALSE, ncp = 5, ind.sup = NULL,
                 quanti.sup = NULL, quali.sup = NULL, row.w = NULL,
                 col.w = NULL, graph = TRUE, axes = c(1,2))
    
    fviz_pca_var(res.pca,
                 habillage=as.data.frame(Norm_Meta_Data)$Family_group,
                 # col.var = "cos2", # Color by contributions to the PC
                 # gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 palette = c( "#80000000", "#56B4E9","#00AFBB", "#E7B800", "#FC4E07","#5D6D7E","#33FFBB","#000000","#E333FF","#52FF33"),
                 repel = TRUE,     # Avoid text overlapping
                 geom = c("text","point"),
                 title="PCA of Variables with Sample ID|Batch and Family Group",
                 labelsize=4
    )
    
    
    
    
    
  }
}