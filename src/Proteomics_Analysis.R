list.of.packages <- c("rstudioapi", "dplyr","VennDiagram","colorspace","ggplot2","reshape","rlist","eulerr","venn",
                      "rafalib","RColorBrewer","biganalytics","pca3d","sinkr","factoextra")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


library(rstudioapi)
library(dplyr)
library(VennDiagram)
library(colorspace)
library(ggplot2)
library(reshape)
library(rlist)
library(eulerr)
library(venn)
library(rafalib)
library(RColorBrewer)
library(biganalytics)
# library(pca3d)


# library(devtools)
# install_github("marchtaylor/sinkr")
library(sinkr)

#Gets project path based on location of source script
project_path <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
lib_path<-paste(project_path,'/src/lib',sep="")
lapply(list.files(pattern = "[.]R$", full.names=TRUE,recursive = FALSE,path=lib_path), source)
#Gets path to where data is stored
data_path <- paste(project_path,'/Data/Proteomics_Dataset',sep="")
#Gets paths of all the Protein data text files with Gene Names Column
files <- list.files(data_path, pattern="^Gene_Name", full.names=T, recursive=FALSE)
#Initiate path for saving Results/Graphs
graph_path<-paste(project_path,'/Results/Graphs/',sep="")

# trace(Protein_Dataset_Info_Extraction, exit = function().last_env <<- parent.frame())
# get(Results[4],.last_env)

# ---- Extract Protein Information ----
Protein_Dataset_Info_Extraction(files)

#Get Protein abundance ratio datasets for each batch
Protein_AR_Datasets <- ls(all.names=TRUE,pattern="*\\d_AR$")
#Temporarily store very first Protein abundance ratio data
tmp <- get(Protein_AR_Datasets[1])
i=1
#Find common proteins between all datasets using semi_join. Semi_join returns values that are matching in the same row for both variables.
for (kk in length(Protein_AR_Datasets)){
  if (kk==length(Protein_AR_Datasets)){break}
  tmp <- semi_join(tmp, get(Protein_AR_Datasets[kk+i]), by = 'Unique.Sequence.ID')
}
Common_Proteins <- tmp
rm(tmp,kk,i)

#Plotting venn diagram of the number of proteins found in each batch and common protein between batches
# grid.newpage()
# venn_protein = draw.triple.venn(area1 = 4887, area2 = 5593, area3 = 6618, n12= 4109, n23 = 4646, n13 = 4018, n123 = 3782, category = c('Set 1', 'Set 2', 'Set 3'), fill =  c("#E495A5" ,"#ABB065" ,"#39BEB1" ))

# ---- Filtering For Common Proteins Goes Here ----
Filtering_Common_Proteins(Protein_AR_Datasets,Variable_Name="_AR_Common_Filtered")

# ---- Transposing date frames for easier maniputation ----
Transpose_Data(Protein_AR_Datasets,Variable_Name="_AR_Transpose")

# ---- Plotting distributions and histograms of abundance ratios for each batch ----
Protein_AR_Transpose_Datasets<-ls(all.names=TRUE,pattern="^Protein\\d_AR_Transpose$")
for (kk in 1:length(Protein_AR_Transpose_Datasets)){
  tmp<-get(paste("Protein",kk,"_AR_Transpose",sep=""))
  assign(paste("Melt_Protein",kk,"_AR_Transpose",sep=""),melt(tmp))
  tmp2<-get(paste("Melt_Protein",kk,"_AR_Transpose",sep=""))
  graph_path_name<-paste(graph_path,"Normal_Distribution_Plot_of_Protein_Dataset_",kk,".tiff",sep="")
  graph_name = paste("Normal Distribution Plot of Protein Dataset ",kk,sep="")
  tmp3<-ggplot(data = tmp2, aes(x = value, col=variable)) + geom_density(na.rm=T)  + theme(legend.position = "none")+ggtitle(graph_name)
  ggsave(graph_path_name,device='tiff')
  # graph_name = paste("Histogram of Protein Dataset ",kk,sep="")
  # ggplot(data=tmp2,aes(x=value))+geom_histogram(bins=1000,na.rm = TRUE) + labs(title=graph_name) + labs(x="Protein Abundance Ratio", y="Frequency")
  # ggave(graph_name,device='tiff')
  # assign(paste("Protein",kk,"_NormPlot",sep=""),tmp3)
}
rm(tmp,tmp2,tmp3)

#Combining all three datasets to plot overall distribution
Combined_Protein_Transpose<-bind_rows(lapply(as.list(Protein_AR_Transpose_Datasets),get))
graph_name = paste("Histogram of Combined Protein Dataset",sep="")
ggplot(data=melt(Combined_Protein_Transpose[,-c(1)]),aes(x=value))+geom_histogram(bins=1000,na.rm = TRUE) + labs(title=graph_name) + labs(x="Protein Abundance Ratio", y="Frequency")
ggsave(graph_name,device='tiff')

#Removing proteins with any NA(missing) values
Combined_Protein_Transpose <- Combined_Protein_Transpose[sapply(Combined_Protein_Transpose, function(x) !any(is.na(x)))]
# variable<- Combined_Protein_Transpose[,c(1)]
# value<-colnames(Combined_Protein_Transpose)[2:length(Combined_Protein_Transpose)]
#Deletes columns 1 and 2
Combined_Protein_Transpose <- Combined_Protein_Transpose[,-c(1:2, 4:10)]
#Plotting distribution of abundance ratios for combined data
Melt_Protein_Transpose <- melt(Combined_Protein_Transpose)
graph_path_name<-paste(graph_path,"Normal_Distribution_Plot_of_Combined_Protein_Dataset",".tiff",sep="")
graph_name = paste("Normal Distribution Plot of Combined Protein Dataset",sep="")
Proteins_NormPlot <- ggplot(data = Melt_Protein_Transpose, aes(x=value, col=variable)) + geom_density(na.rm=TRUE)  + theme(legend.position = "none")+ggtitle(graph_name)
ggsave(graph_path_name,device="tiff")

#Boxplot for proteins
graph_path_name<-paste(graph_path,"BoxPlot_of_Combined_Protein_Dataset_",kk,".tiff",sep="")
graph_name = paste("BoxPlot of Combined Protein Dataset",sep="")
Proteins_BoxPlot = ggplot(data=Melt_Protein_Transpose, aes(x=variable, y=value)) + geom_boxplot(aes(x= variable, y = value), na.rm=TRUE) +theme(legend.position = "none")+ggtitle(graph_name)
ggsave(graph_path_name,device='tiff')

graph_path_name<-paste(graph_path,"BoxPlot_function_of_Combined_Protein_Dataset_",kk,".tiff",sep="")
tiff(file=graph_path_name)
Proteins_using_Boxplot_Func<- boxplot(Combined_Protein_Transpose,main="Boxplot of Combined Protein Dataset using boxplot func",xlab="Protein",ylab="Normalized Abundance Ratio")
dev.off()

# --- Import Patient Info ----
#Gets paths of all the Protein data text files
Patient_files <- list.files(data_path, pattern="*Sample_List_Information.csv", full.names=T, recursive=FALSE)
Patient_Info_Extraction(Patient_files)

List_of_Transposed_Data<-ls(all.names=TRUE,pattern="^Protein\\d_AR_Transpose$")
List_of_Patient_Info<-ls(all.names=TRUE,pattern="^Patient_Sample_Set\\d$")

# ---- Appending Patient Info to Protein Data ----
Patient_Info_Append(List_of_Transposed_Data,List_of_Patient_Info)

Protein_Datasets<-ls(all.names=TRUE,pattern="^Protein_Dataset\\d$")
# #Combining all datasets to plot overall distribution
# Combined_Data_with_Patient_Info<-bind_rows(lapply(as.list(Protein_Datasets),get))
# Filtered_Combined_Data <- Combined_Data_with_Patient_Info[sapply(Combined_Data_with_Patient_Info, function(x) !any(is.na(x)))] 

# ---- Normalising data (median centring) ----
Data_Normalization(Protein_Datasets,graph_path,plot=F)
FileName<-paste(project_path,'/Results/OutputFiles/','NormDataFrameset_1_to_6.csv',sep="")
write.csv(NormDataFrame,file=FileName)


Combined_Protein_Transpose <- NormDataFrame[,-c(1:12)]
#Plotting distribution of abundance ratios for combined data
Melt_Protein_Transpose <- melt(Combined_Protein_Transpose)
graph_path_name<-paste(graph_path,"Normalized_Normal_Distribution_Plot_of_Combined_Protein_Dataset",".jpeg",sep="")
graph_name = paste("Normal Distribution Plot of Combined Protein Dataset",sep="")
Proteins_NormPlot <- ggplot(data = Melt_Protein_Transpose, aes(x=value, col=variable)) + geom_density(na.rm=TRUE)  + theme(legend.position = "none")+ggtitle(graph_name)
ggsave(graph_path_name,device="jpeg")


# ---- Hierarchial Clustering ----
graph_path<-paste(graph_path,'Dendrograms/',sep='')
tmp_O<-NormDataFrame
# tmp2<- tmp[sapply(tmp, function(x) !any(is.na(x)))] 
# #Combining patient id with different groups so that patient id can be seen in dendrogram
# tmp3<-apply(tmp2[,11:ncol(tmp2)],2,as.numeric)
tmp<-as.matrix(tmp_O)
tmp[is.na(tmp)] <- 0
tmp<-as.data.frame(tmp)
tmp1<-as.data.frame(tmp_df[,c(1:12)])
tmp2[]<-as.data.frame(lapply(tmp_df[,-c(1:12)],function(x) as.numeric(as.character(x))))
tmp<-cbind(tmp1,tmp2)

# tmp<-as.data.frame(as.numeric(tmp[,-c(1:12)]))


Dist_Matrix<-dist(tmp[,-c(1:12)],method="euclidean")
CFS_cluster = hclust(Dist_Matrix)
tmp$fam_id <- paste(NormDataFrame$Sample.ID, NormDataFrame$Family_group, sep="|")
tmp$CFS_id <- paste(NormDataFrame$Sample.ID, NormDataFrame$CFS, sep="|")
tmp$age_id <- paste(NormDataFrame$Sample.ID, NormDataFrame$Age, sep="|")
tmp$gender_id <- paste(NormDataFrame$Sample.ID, NormDataFrame$Gender, sep="|")
tmp$EDS_id <- paste(NormDataFrame$Sample.ID, NormDataFrame$EDS, sep="|")
tmp$ds_id <- paste(NormDataFrame$Sample.ID, NormDataFrame$disease.status, sep="|")
tmp$batch_id <- paste(NormDataFrame$Sample.ID, NormDataFrame$Batch, sep="|")
tmp$channel_id <- paste(NormDataFrame$Sample.ID, NormDataFrame$TMT.Label, sep="|")

#Color coding dendrogram by different groups

# CFS_cluster_id <- myplclust(CFS_cluster, labels = as.character(tmp$Sample.ID),lab.col=rainbow_hcl(8)[as.fumeric(as.character(tmp$Batch))], main="Dendrogram of Total Sample Clustering, color by batch")

graph_path_name<-paste(graph_path,"Clustering_color_by_Family_Group.tiff",sep="")
tiff(file=graph_path_name, width = 1529, height = 876, units = "px",res = 100)
CFS_cluster_fam = myplclust(CFS_cluster, labels = as.character(tmp$fam_id), lab.col= rainbow_hcl(5)[as.fumeric(as.character(tmp$Family_group))], main="Dendrogram of Total Sample Clustering, color by Family Group")
dev.off()

graph_path_name<-paste(graph_path,"Clustering_color_by_CFS.tiff",sep="")
tiff(file=graph_path_name, width = 1529, height = 876, units = "px",res = 100)
CFS_cluster_CFS = myplclust(CFS_cluster, labels = as.character(tmp$CFS_id), lab.col= rainbow_hcl(3)[as.fumeric(as.character(tmp$CFS))], main="Dendrogram of Total Sample Clustering, color by CFS")
dev.off()

graph_path_name<-paste(graph_path,"Clustering_color_by_Age_Group.tiff",sep="")
tiff(file=graph_path_name, width = 1529, height = 876, units = "px",res = 100)
CFS_cluster_age = myplclust(CFS_cluster, labels = as.character(tmp$age_id), lab.col= rainbow_hcl(9)[as.fumeric(as.character(tmp$Age_group))], main="Dendrogram of Total Sample Clustering, color by Age Group")
dev.off()

graph_path_name<-paste(graph_path,"Clustering_color_by_Gender.tiff",sep="")
tiff(file=graph_path_name, width = 1529, height = 876, units = "px",res = 100)
CFS_cluster_gender = myplclust(CFS_cluster, labels = as.character(tmp$gender_id), lab.col= rainbow_hcl(2)[as.fumeric(as.character(tmp$Gender))], main="Dendrogram of Total Sample Clustering, color by Gender")
dev.off()

graph_path_name<-paste(graph_path,"Clustering_color_by_EDS.tiff",sep="")
tiff(file=graph_path_name, width = 1529, height = 876, units = "px",res = 100)
eds_col = c("#86B875", "#E495A5", "#7DB0DD")
CFS_cluster_EDS = myplclust(CFS_cluster, labels = as.character(tmp$EDS_id), lab.col= eds_col[as.fumeric(as.character(tmp$EDS))], main="Dendrogram of Total Sample Clustering, color by EDS")
dev.off()

graph_path_name<-paste(graph_path,"Clustering_color_by_Disease_Status.tiff",sep="")
tiff(file=graph_path_name, width = 1529, height = 876, units = "px",res = 100)
CFS_cluster_ds = myplclust(CFS_cluster, labels = as.character(tmp$ds_id), lab.col= rainbow_hcl(8)[as.fumeric(as.character(tmp$disease.status))], main="Dendrogram of Total Sample Clustering, color by Disease Status")
dev.off()

graph_path_name<-paste(graph_path,"Clustering_color_by_batch.tiff",sep="")
tiff(file=graph_path_name, width = 1529, height = 876, units = "px",res = 100)
CFS_cluster_batch = myplclust(CFS_cluster, labels = as.character(tmp$batch_id), lab.col= rainbow_hcl(3)[as.fumeric(as.character(tmp$Batch))])
dev.off()

graph_path_name<-paste(graph_path,"Clustering_color_by_Channel.tiff",sep="")
tiff(file=graph_path_name, width = 1529, height = 876, units = "px",res = 100)
CFS_cluster_channel = myplclust(CFS_cluster, labels = as.character(tmp$channel_id), lab.col= rainbow_hcl(9)[as.fumeric(as.character(tmp$TMT.Label))], main="Dendrogram of Total Sample Clustering, color by TMT Label")
dev.off()

Data_WD<-NormDataFrame[grep("WD", NormDataFrame$Sample.ID), ]
tmp<-as.matrix(Data_WD)
tmp[is.na(tmp)] <- 0
tmp<-as.data.frame(tmp)
Dist_Matrix<-dist(tmp[,13:ncol(tmp)],method="euclidean")
CFS_cluster = hclust(Dist_Matrix)
CFS_cluster_id <- myplclust(CFS_cluster, labels = as.character(tmp$Sample.ID),lab.col=rainbow_hcl(8)[as.fumeric(as.character(tmp$Batch))], main="Dendrogram of WD Samples, color by batch")



# Hierarchial Clustering with DINEOF  -------------------------------------

# Method two
library(sinkr)
# graph_path<-paste(graph_path,'Dendrograms/',sep='')
Norm_Data_Matrix <- (as.matrix(NormDataFrame[,-c(1:12)]))
Norm_Meta_Data <- (as.matrix(NormDataFrame[,c(1:12)]))

Missing_Data_Fill <- dineof(Norm_Data_Matrix)

Dist_Matrix<-dist(Missing_Data_Fill$Xa,method="euclidean")
CFS_cluster = hclust(Dist_Matrix)
#Transpose so each row is protein
# CFS_cluster_id <- myplclust(CFS_cluster, labels = as.character(Norm_Meta_Data["Sample.ID",1:ncol(Norm_Meta_Data)]),lab.col=rainbow_hcl(8)[as.fumeric(as.character(Norm_Meta_Data["Batch",1:ncol(Norm_Meta_Data)]))], main="Dendrogram of Total Sample Clustering, color by batch")
#Transpose so each column is protein
graph_path_name<-paste(graph_path,"Dineof_Clustering_color_by_batch.tiff",sep="")
tiff(file=graph_path_name, width = 1529, height = 876, units = "px",res = 100)
CFS_cluster_id <- myplclust(CFS_cluster, 
                            labels = as.character(paste(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Sample.ID"],Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Batch"],sep="|")),
                            lab.col=rainbow_hcl(8)[as.fumeric(as.character(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Batch"]))], 
                            main="Dendrogram of Total Sample Clustering, color by batch")

dev.off()

graph_path_name<-paste(graph_path,"Dineof_Clustering_color_by_Family.tiff",sep="")
tiff(file=graph_path_name, width = 1529, height = 876, units = "px",res = 100)
CFS_cluster_id <- myplclust(CFS_cluster, 
                            labels = as.character(paste(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Sample.ID"],Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Family"],sep="|")),
                            lab.col=rainbow_hcl(8)[as.fumeric(as.character(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Family"]))], 
                            main="Dendrogram of Total Sample Clustering, color by Family")
dev.off()

graph_path_name<-paste(graph_path,"Dineof_Clustering_color_by_Gender.tiff",sep="")
tiff(file=graph_path_name, width = 1529, height = 876, units = "px",res = 100)
CFS_cluster_id <- myplclust(CFS_cluster, 
                            labels = as.character(paste(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Sample.ID"],Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Gender"],sep="|")),
                            lab.col=rainbow_hcl(8)[as.fumeric(as.character(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Gender"]))], 
                            main="Dendrogram of Total Sample Clustering, color by Gender")
dev.off()

graph_path_name<-paste(graph_path,"Dineof_Clustering_color_by_Age.tiff",sep="")
tiff(file=graph_path_name, width = 1529, height = 876, units = "px",res = 100)
CFS_cluster_id <- myplclust(CFS_cluster, 
                            labels = as.character(paste(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Sample.ID"],Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Age"],sep="|")),
                            lab.col=rainbow_hcl(8)[as.fumeric(as.character(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Age"]))],
                            main="Dendrogram of Total Sample Clustering, color by Age")
dev.off()

graph_path_name<-paste(graph_path,"Dineof_Clustering_color_by_CFS.tiff",sep="")
tiff(file=graph_path_name, width = 1529, height = 876, units = "px",res = 100)
CFS_cluster_id <- myplclust(CFS_cluster, 
                            labels = as.character(paste(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Sample.ID"],Norm_Meta_Data[1:nrow(Norm_Meta_Data),"CFS"],sep="|")),
                            lab.col=rainbow_hcl(8)[as.fumeric(as.character(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"CFS"]))], 
                            main="Dendrogram of Total Sample Clustering, color by CFS")
dev.off()

graph_path_name<-paste(graph_path,"Dineof_Clustering_color_by_Disease_Status.tiff",sep="")
tiff(file=graph_path_name, width = 1529, height = 876, units = "px",res = 100)
CFS_cluster_id <- myplclust(CFS_cluster, 
                            labels = as.character(paste(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Sample.ID"],Norm_Meta_Data[1:nrow(Norm_Meta_Data),"disease.status"],sep="|")),
                            lab.col=rainbow_hcl(8)[as.fumeric(as.character(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"disease.status"]))],
                            main="Dendrogram of Total Sample Clustering, color by Disease Status")
dev.off()

graph_path_name<-paste(graph_path,"Dineof_Clustering_color_by_EDS.tiff",sep="")
tiff(file=graph_path_name, width = 1529, height = 876, units = "px",res = 100)
CFS_cluster_id <- myplclust(CFS_cluster, 
                            labels = as.character(paste(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Sample.ID"],Norm_Meta_Data[1:nrow(Norm_Meta_Data),"EDS"],sep="|")),
                            lab.col=rainbow_hcl(8)[as.fumeric(as.character(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"EDS"]))],
                            main="Dendrogram of Total Sample Clustering, color by EDS")
dev.off()

graph_path_name<-paste(graph_path,"Dineof_Clustering_color_by_Age_Group.tiff",sep="")
tiff(file=graph_path_name, width = 1529, height = 876, units = "px",res = 100)
CFS_cluster_id <- myplclust(CFS_cluster, 
                            labels = as.character(paste(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Sample.ID"],Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Age_group"],sep="|")),
                            lab.col=rainbow_hcl(8)[as.fumeric(as.character(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Age_group"]))], 
                            main="Dendrogram of Total Sample Clustering, color by Age Group")
dev.off()

graph_path_name<-paste(graph_path,"Dineof_Clustering_color_by_Family_Group.tiff",sep="")
tiff(file=graph_path_name, width = 1529, height = 876, units = "px",res = 100)
CFS_cluster_id <- myplclust(CFS_cluster, 
                            labels = as.character(paste(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Sample.ID"],Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Family_group"],sep="|")),
                            lab.col=rainbow_hcl(8)[as.fumeric(as.character(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Family_group"]))], 
                            main="Dendrogram of Total Sample Clustering, color by Family Group")
dev.off()

graph_path_name<-paste(graph_path,"Dineof_Clustering_color_by_TMT_Label.tiff",sep="")
tiff(file=graph_path_name, width = 1529, height = 876, units = "px",res = 100)
CFS_cluster_id <- myplclust(CFS_cluster, 
                            labels = as.character(paste(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"Sample.ID"],Norm_Meta_Data[1:nrow(Norm_Meta_Data),"TMT.Label"],sep="|")),
                            lab.col=rainbow_hcl(8)[as.fumeric(as.character(Norm_Meta_Data[1:nrow(Norm_Meta_Data),"TMT.Label"]))], 
                            main="Dendrogram of Total Sample Clustering, color by TMT Label")
dev.off()






# Hierarchial Clustering with Removal of more than 60% missing values -------------------

# Method two
library(sinkr)
# graph_path<-paste(graph_path,'Dendrograms/',sep='')
Norm_Data_Matrix <- (as.matrix(NormDataFrame[,-c(1:12)]))
Norm_Meta_Data <- (as.matrix(NormDataFrame[,c(1:12)]))

#Removing Missing Data that have more than 60% missing values per protein
data<-t(Norm_Data_Matrix)
Total_Col<-ncol(data)
remove_Rows=matrix()
k=1
for (i in 1:nrow(data)){
  
  sum_na<-sum(is.na(data[i,]))
  if ((sum_na/Total_Col)>=0.6){
   print(i)
    remove_Rows[k]<-i
    k=k+1
  }
  
}

data<-data[-remove_Rows,]
row_means<-rowMeans(data,na.rm=T)
row_SD<-apply(data,1,sd,na.rm=T)
rnorm(data,row_means,row_SD)

data<-as.data.frame(data)
for (i in 1:nrow(data)){
  
  indices <- which(is.na(data[i,]))
  
  data[indices]<- rnorm(length(indices), row_means[i], row_SD[i])
  
}
Missing_Data_Fill<-as.matrix(data)
Clustering_Replacing_NaN_With_rnorm_mean_SD(Missing_Data_Fill,graph_path)
Clustering_ComBat(Missing_Data_Fill,graph_path)  






# ---- Data Seperartion ----

Norm_Data_Matrix <- as.matrix(NormDataFrame[,-c(1:12)])
Norm_Meta_Data <- as.matrix(NormDataFrame[,c(1:12)])

# ---- PCA ----
library(sinkr)
library(pca3d)
library(factoextra)
library(FactoMineR)
# library(doParallel)

# NormData<-NormDataFrame[,13:ncol(NormDataFrame)]
NormData<-(Norm_Data_Matrix)
NormData[is.na(Norm_Data_Matrix)] <- 0



test<-dineof(t(Norm_Data_Matrix))
tmp<-as.data.frame(test$Xa)
colnames(tmp)<-paste(t((Norm_Meta_Data))["Sample.ID",],t((Norm_Meta_Data))["Batch",],t((Norm_Meta_Data))["Gender",],sep=' | ')
# tmp<-NormData
# rownames(tmp)<-NormDataFrame$Set_TMT_Label


pca_CFS<-prcomp(tmp)

# biplot(pca_CFS,main="PCA biplot of CFS sample proteins")

# pca3d(pca_CFS,legend=T)


fviz_eig(pca_CFS)
fviz_pca_ind(pca_CFS,
             col.ind = "cos2", # Color by the quality of representation
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
              repel = TRUE,     # Avoid text overlapping
              select.ind = list(cos2 = 0.99),
             title="PCA of independants when cos2 >= 0.99"

)

fviz_pca_var(pca_CFS,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,     # Avoid text overlapping
             geom = c("text","point"),
             title="PCA of Variables with Sample ID|Batch and Gender"
             
)

fviz_pca_biplot(pca_CFS, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969",  # Individuals color
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                # repel = TRUE,     # Avoid text overlapping
                select.ind = list(cos2 = 0.99),
                geom.var = c("text","point"),
                geom.ind = c("text","point"),
                title="PCA of biplot when cos2 >= 0.99"
                
)

# tmp1<-NormDataFrame[,1]
# tmp<-(Norm_Data_Matrix)
# tmp<-dineof(tmp)
# rownames(tmp$Xa)<-tmp1
# tmp<-as.data.frame(tmp$Xa)

cor.mat <- round(cor(tmp),2)
# install.packages("corrplot")
library("corrplot")
{plot.new(); 
corrplot(cor.mat, type="upper", 
         tl.col="black", tl.srt=45,tl.cex=0.5,
         title="Correlations plot of each sample in order of hclust",
         outline=T,order="hclust")
dev.off()}

# install.packages("PerformanceAnalytics")
library("PerformanceAnalytics")
chart.Correlation(tmp, histogram=TRUE, pch=19)

res_pca<-PCA(t(tmp))



# ---- Univariate Modeling ----
Protein_fill<-dineof(Norm_Data_Matrix)
Protein_Dataframe<-cbind(Norm_Meta_Data,as.data.frame(Protein_fill$Xa))
Uni_Model_Results<-Univariate_Linear_Model(Protein_Dataframe)
Adjusted_PVals<-FDR_Adjusted_PVal(Uni_Model_Results)
Adjusted_PVals = cbind(as.data.frame(colnames(Protein_Dataframe)[-c(1:12)]), Adjusted_PVals)
colnames(Adjusted_PVals)[1] = 'Gene_Name'
P_Value_Thresholding(Adjusted_PVals)
NormDataFrame[,c(1:12)]


  
fit<-lm(formula = y ~ Batch+Family_group+disease.status+Gender+Age,data = NormDataFrame)


