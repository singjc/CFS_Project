Clustering_Replacing_NaN_With_rnorm_mean_SD <-function(Missing_Data_Fill,graph_path){
  
  
  Dist_Matrix<-dist(t(Missing_Data_Fill),method="euclidean")
  
  CFS_cluster = hclust(Dist_Matrix)
  #Transpose so each row is protein
  # CFS_cluster_id <- myplclust(CFS_cluster, labels = as.character(Norm_Meta_Data["Sample.ID",1:ncol(Norm_Meta_Data)]),lab.col=rainbow_hcl(8)[as.fumeric(as.character(Norm_Meta_Data["Batch",1:ncol(Norm_Meta_Data)]))], main="Dendrogram of Total Sample Clustering, color by batch")
  #Transpose so each column is protein
  graph_path_name<-paste(graph_path,"rnorm_Clustering_color_by_batch.tiff",sep="")
  tiff(file=graph_path_name, width = 1529, height = 876, units = "px",res = 100)
  CFS_cluster_id <- myplclust(CFS_cluster, 
                              labels = as.character(paste(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Sample.ID"],Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Batch"],sep="|")),
                              lab.col=rainbow_hcl(8)[as.fumeric(as.character(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Batch"]))], 
                              main="Dendrogram of Total Sample Clustering, color by batch")
  
  dev.off()
  
  graph_path_name<-paste(graph_path,"rnorm_Clustering_color_by_Family.tiff",sep="")
  tiff(file=graph_path_name, width = 1529, height = 876, units = "px",res = 100)
  CFS_cluster_id <- myplclust(CFS_cluster, 
                              labels = as.character(paste(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Sample.ID"],Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Family"],sep="|")),
                              lab.col=rainbow_hcl(8)[as.fumeric(as.character(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Family"]))], 
                              main="Dendrogram of Total Sample Clustering, color by Family")
  dev.off()
  
  graph_path_name<-paste(graph_path,"rnorm_Clustering_color_by_Gender.tiff",sep="")
  tiff(file=graph_path_name, width = 1529, height = 876, units = "px",res = 100)
  CFS_cluster_id <- myplclust(CFS_cluster, 
                              labels = as.character(paste(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Sample.ID"],Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Gender"],sep="|")),
                              lab.col=rainbow_hcl(8)[as.fumeric(as.character(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Gender"]))], 
                              main="Dendrogram of Total Sample Clustering, color by Gender")
  dev.off()
  
  graph_path_name<-paste(graph_path,"rnorm_Clustering_color_by_Age.tiff",sep="")
  tiff(file=graph_path_name, width = 1529, height = 876, units = "px",res = 100)
  CFS_cluster_id <- myplclust(CFS_cluster, 
                              labels = as.character(paste(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Sample.ID"],Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Age"],sep="|")),
                              lab.col=rainbow_hcl(8)[as.fumeric(as.character(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Age"]))],
                              main="Dendrogram of Total Sample Clustering, color by Age")
  dev.off()
  
  graph_path_name<-paste(graph_path,"rnorm_Clustering_color_by_CFS.tiff",sep="")
  tiff(file=graph_path_name, width = 1529, height = 876, units = "px",res = 100)
  CFS_cluster_id <- myplclust(CFS_cluster, 
                              labels = as.character(paste(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Sample.ID"],Norm_Meta_Data[1:nrow(Norm_Meta_Data),"CFS"],sep="|")),
                              lab.col=rainbow_hcl(8)[as.fumeric(as.character(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"CFS"]))], 
                              main="Dendrogram of Total Sample Clustering, color by CFS")
  dev.off()
  
  graph_path_name<-paste(graph_path,"rnorm_Clustering_color_by_Disease_Status.tiff",sep="")
  tiff(file=graph_path_name, width = 1529, height = 876, units = "px",res = 100)
  CFS_cluster_id <- myplclust(CFS_cluster, 
                              labels = as.character(paste(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Sample.ID"],Norm_Meta_Data[1:nrow(Norm_Meta_Data),"disease.status"],sep="|")),
                              lab.col=rainbow_hcl(8)[as.fumeric(as.character(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"disease.status"]))],
                              main="Dendrogram of Total Sample Clustering, color by Disease Status")
  dev.off()
  
  graph_path_name<-paste(graph_path,"rnorm_Clustering_color_by_EDS.tiff",sep="")
  tiff(file=graph_path_name, width = 1529, height = 876, units = "px",res = 100)
  CFS_cluster_id <- myplclust(CFS_cluster, 
                              labels = as.character(paste(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Sample.ID"],Norm_Meta_Data[1:nrow(Norm_Meta_Data),"EDS"],sep="|")),
                              lab.col=rainbow_hcl(8)[as.fumeric(as.character(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"EDS"]))],
                              main="Dendrogram of Total Sample Clustering, color by EDS")
  dev.off()
  
  graph_path_name<-paste(graph_path,"rnorm_Clustering_color_by_Age_Group.tiff",sep="")
  tiff(file=graph_path_name, width = 1529, height = 876, units = "px",res = 100)
  CFS_cluster_id <- myplclust(CFS_cluster, 
                              labels = as.character(paste(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Sample.ID"],Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Age_group"],sep="|")),
                              lab.col=rainbow_hcl(8)[as.fumeric(as.character(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Age_group"]))], 
                              main="Dendrogram of Total Sample Clustering, color by Age Group")
  dev.off()
  
  graph_path_name<-paste(graph_path,"rnorm_Clustering_color_by_Family_Group.tiff",sep="")
  tiff(file=graph_path_name, width = 1529, height = 876, units = "px",res = 100)
  CFS_cluster_id <- myplclust(CFS_cluster, 
                              labels = as.character(paste(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Sample.ID"],Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Family_group"],sep="|")),
                              lab.col=rainbow_hcl(8)[as.fumeric(as.character(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Family_group"]))], 
                              main="Dendrogram of Total Sample Clustering, color by Family Group")
  dev.off()
  
  graph_path_name<-paste(graph_path,"rnorm_Clustering_color_by_TMT_Label.tiff",sep="")
  tiff(file=graph_path_name, width = 1529, height = 876, units = "px",res = 100)
  CFS_cluster_id <- myplclust(CFS_cluster, 
                              labels = as.character(paste(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Sample.ID"],Norm_Meta_Data[1:nrow(Norm_Meta_Data),"TMT.Label"],sep="|")),
                              lab.col=rainbow_hcl(8)[as.fumeric(as.character(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"TMT.Label"]))], 
                              main="Dendrogram of Total Sample Clustering, color by TMT Label")
  dev.off()
}