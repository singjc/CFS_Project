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
library(colorspace)
library(biganalytics)
library(pca3d)

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
  graph_name = paste("Histogram of Protein Dataset ",kk,sep="")
  ggplot(data=tmp2,aes(x=value))+geom_histogram(bins=1000,na.rm = TRUE) + labs(title=graph_name) + labs(x="Protein Abundance Ratio", y="Frequency")
  ggave(graph_name,device='tiff')
  assign(paste("Protein",kk,"_NormPlot",sep=""),tmp3)
}
rm(tmp,tmp2,tmp3)

#Combining all three datasets to plot overall distribution
Combined_Protein_Transpose<-bind_rows(lapply(as.list(Protein_AR_Transpose_Datasets),get))
graph_name = paste("Histogram of Combined Protein Dataset",sep="")
ggplot(data=melt(Combined_Protein_Transpose[,-c(1)]),aes(x=value))+geom_histogram(bins=1000,na.rm = TRUE) + labs(title=graph_name) + labs(x="Protein Abundance Ratio", y="Frequency")
ggave(graph_name,device='tiff')

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


# ---- Hierarchial Clustering ----
tmp<-NormDataFrame
# tmp2<- tmp[sapply(tmp, function(x) !any(is.na(x)))] 
# #Combining patient id with different groups so that patient id can be seen in dendrogram
# tmp3<-apply(tmp2[,11:ncol(tmp2)],2,as.numeric)
tmp<-as.matrix(tmp)
tmp[is.na(tmp)] <- 0
tmp<-as.data.frame(tmp)

Dist_Matrix<-dist(tmp[,13:ncol(tmp)],method="euclidean")
CFS_cluster = hclust(Dist_Matrix)
tmp$fam_id <- paste(NormDataFrame$sample_id, NormDataFrame$Family_group, sep="|")
tmp$CFS_id <- paste(NormDataFrame$sample_id, NormDataFrame$CFS, sep="|")
tmp$age_id <- paste(NormDataFrame$sample_id, NormDataFrame$Age, sep="|")
tmp$gender_id <- paste(NormDataFrame$sample_id, NormDataFrame$Gender, sep="|")
tmp$EDS_id <- paste(NormDataFrame$sample_id, NormDataFrame$EDS, sep="|")
tmp$ds_id <- paste(NormDataFrame$sample_id, NormDataFrame$disease.status, sep="|")
tmp$batch_id <- paste(NormDataFrame$sample_id, NormDataFrame$Batch, sep="|")
tmp$channel_id <- paste(NormDataFrame$sample_id, NormDataFrame$TMT.Label, sep="|")

#Color coding dendrogram by different groups

CFS_cluster_id <- myplclust(CFS_cluster, labels = as.character(tmp$Sample.ID),lab.col=rainbow_hcl(8)[as.fumeric(as.character(tmp$Batch))], main="Dendrogram of Total Sample Clustering, color by batch")


CFS_cluster_fam = myplclust(CFS_cluster, labels = as.character(tmp$fam_id), lab.col= rainbow_hcl(5)[as.fumeric(as.character(tmp$Family_group))])
CFS_cluster_CFS = myplclust(CFS_cluster, labels = as.character(tmp$CFS_id), lab.col= rainbow_hcl(3)[as.fumeric(as.character(tmp$CFS))])
CFS_cluster_age = myplclust(CFS_cluster, labels = as.character(tmp$age_id), lab.col= rainbow_hcl(9)[as.fumeric(as.character(tmp$Age_group))])
CFS_cluster_gender = myplclust(CFS_cluster, labels = as.character(tmp$gender_id), lab.col= rainbow_hcl(2)[as.fumeric(as.character(tmp$Gender))])
eds_col = c("#86B875", "#E495A5", "#7DB0DD")
CFS_cluster_EDS = myplclust(CFS_cluster, labels = as.character(tmp$EDS_id), lab.col= eds_col[as.fumeric(as.character(tmp$EDS))])
CFS_cluster_ds = myplclust(CFS_cluster, labels = as.character(tmp$ds_id), lab.col= rainbow_hcl(8)[as.fumeric(as.character(tmp$disease.status))])
CFS_cluster_batch = myplclust(CFS_cluster, labels = as.character(tmp$batch_id), lab.col= rainbow_hcl(3)[as.fumeric(as.character(tmp$Batch))])
CFS_cluster_channel = myplclust(CFS_cluster, labels = as.character(tmp$channel_id), lab.col= rainbow_hcl(9)[as.fumeric(as.character(tmp$TMT.Label))])

Data_WD<-NormDataFrame[grep("WD", NormDataFrame$Sample.ID), ]
tmp<-as.matrix(Data_WD)
tmp[is.na(tmp)] <- 0
tmp<-as.data.frame(tmp)
Dist_Matrix<-dist(tmp[,13:ncol(tmp)],method="euclidean")
CFS_cluster = hclust(Dist_Matrix)
CFS_cluster_id <- myplclust(CFS_cluster, labels = as.character(tmp$Sample.ID),lab.col=rainbow_hcl(8)[as.fumeric(as.character(tmp$Batch))], main="Dendrogram of WD Samples, color by batch")

# ---- Data Seperartion ----

Norm_Data_Matrix <- as.matrix(NormDataFrame[,-c(1:12)])
Norm_Meta_Data <- as.matrix(NormDataFrame[,c(1:12)])

# ---- PCA ----
library(sinkr)
library(pca3d)
library(factoextra)
library(doParallel)

# NormData<-NormDataFrame[,13:ncol(NormDataFrame)]
NormData<-(Norm_Data_Matrix)
NormData[is.na(Norm_Data_Matrix)] <- 0



test<-dineof(NormData)
tmp<-as.data.frame(test$Xa)
tmp<-NormData
rownames(tmp)<-NormDataFrame$Set_TMT_Label


pca_CFS<-prcomp(t(tmp))

biplot(pca_CFS,main="PCA biplot of CFS sample proteins")

pca3d(pca_CFS,legend=T)


fviz_eig(pca_CFS)
fviz_pca_ind(pca_CFS,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             
)

fviz_pca_var(pca_CFS,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             
)

# fviz_pca_biplot(pca_CFS, repel = TRUE,
#                 col.var = "#2E9FDF", # Variables color
#                 col.ind = "#696969"  # Individuals color
# )

no_cores <- detectCores() - 1
# Create cluster with desired number of cores
cl <- makeCluster(3)
# Register cluster
registerDoParallel(cl)
# Find out how many cores are being used
getDoParWorkers()
# Stopping cluster
stopCluster(cl)

# ---- Linear Modeling ----

NormDataFrame$Batch<-factor(NormDataFrame$Batch)
NormDataFrame$Family_group<-factor(NormDataFrame$Family_group)
NormDataFrame$disease.status<-factor(NormDataFrame$disease.status)
NormDataFrame$Gender<-factor(NormDataFrame$Gender)
NormDataFrame$Age<-factor(NormDataFrame$Age)

y<-NormDataFrame[,-c(1:12)]

fit<-lm(formula =  y ~ disease.status,data=NormDataFrame)

Protein_AR_LM<-t(NormDataFrame[,-c(1:12)])
colnames(Protein_AR_LM)<-t(NormDataFrame[,c(8)])
# tmp<-bind_rows(Predictor_Variables,Protein_AR_LM)

fit<-lm(y~.,data=as.data.frame(Protein_AR_LM))
  
fit<-lm(formula = y ~ Batch+Family_group+disease.status+Gender+Age,data = NormDataFrame)


