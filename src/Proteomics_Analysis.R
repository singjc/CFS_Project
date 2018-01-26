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

# source('~/Users/justinsing/Documents/Hannes_Rost/CFS/CFS_Project/src/lib/Patient_Info_Extraction.R ')
# paste(Protein_AR_Datasets,'/src/lib')

#Gets project path based on location of source script
project_path <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
lib_path<-paste(project_path,'/src/lib',sep="")
lapply(list.files(pattern = "[.]R$", full.names=TRUE,recursive = FALSE,path=lib_path), source)
#Gets path to where data is stored
data_path <- paste(project_path,'/Data/Proteomics_Dataset',sep="")
#Gets paths of all the Protein data text files
files <- list.files(data_path, pattern="*Proteins.txt", full.names=T, recursive=FALSE)
#Initiate path for saving Results/Graphs
graph_path<-paste(project_path,'/Results/Graphs/',sep="")

# trace(Protein_Dataset_Info_Extraction, exit = function().last_env <<- parent.frame())
# get(Results[4],.last_env)
#Extract Protein Information
Protein_Dataset_Info_Extraction(files)

#Get Protein abundance ratio datasets for each batch
Protein_AR_Datasets <- ls(all.names=TRUE,pattern="*\\d_AR$")
#Temporarily store very first Protein abundance ratio data
tmp <- get(Protein_AR_Datasets[1])
i=1
#Find common proteins between all datasets using semi_join. Semi_join returns values that are matching in the same row for both variables.
for (kk in length(Protein_AR_Datasets)){
  if (kk==length(Protein_AR_Datasets)){break}
  tmp <- semi_join(tmp, get(Protein_AR_Datasets[kk+i]), by = 'Unique Sequence ID')
}
Common_Proteins <- tmp
rm(tmp,kk,i)

#Plotting venn diagram of the number of proteins found in each batch and common protein between batches
# grid.newpage()
# venn_protein = draw.triple.venn(area1 = 4887, area2 = 5593, area3 = 6618, n12= 4109, n23 = 4646, n13 = 4018, n123 = 3782, category = c('Set 1', 'Set 2', 'Set 3'), fill =  c("#E495A5" ,"#ABB065" ,"#39BEB1" ))

#Filtering For Common Proteins Goes Here
Filtering_Common_Proteins(Protein_AR_Datasets,Variable_Name="_AR_Common_Filtered")

#Transposing date frames for easier maniputation
Transpose_Data(Protein_AR_Datasets,Variable_Name="_AR_Transpose")

#Plotting distributions of abundance ratios for each batch
Protein_AR_Transpose_Datasets<-ls(all.names=TRUE,pattern="^Protein\\d_AR_Transpose$")
for (kk in 1:length(Protein_AR_Transpose_Datasets)){
  tmp<-get(paste("Protein",kk,"_AR_Transpose",sep=""))
  assign(paste("Melt_Protein",kk,"_AR_Transpose",sep=""),melt(tmp))
  tmp2<-get(paste("Melt_Protein",kk,"_AR_Transpose",sep=""))
  graph_path_name<-paste(graph_path,"Normal_Distribution_Plot_of_Protein_Dataset_",kk,".tiff",sep="")
  graph_name = paste("Normal Distribution Plot of Protein Dataset ",kk,sep="")
  tmp3<-ggplot(data = tmp2, aes(x = value, col=variable)) + geom_density(na.rm=T)  + theme(legend.position = "none")+ggtitle(graph_name)
  tmp3
  ggsave(graph_path_name,device='tiff')
  assign(paste("Protein",kk,"_NormPlot",sep=""),tmp3)
}
rm(tmp,tmp2,tmp3)

#Combining all three datasets to plot overall distribution
Combined_Protein_Transpose<-bind_rows(lapply(as.list(Protein_AR_Transpose_Datasets),get))

#Removing proteins with any NA(missing) values
Combined_Protein_Transpose <- Combined_Protein_Transpose[sapply(Combined_Protein_Transpose, function(x) !any(is.na(x)))]
# variable<- Combined_Protein_Transpose[,c(1)]
# value<-colnames(Combined_Protein_Transpose)[2:length(Combined_Protein_Transpose)]
#Deletes columns 1 and 2
Combined_Protein_Transpose <- Combined_Protein_Transpose[,-c(1:2, 4:10)]
#Plotting distribution of abundance ratios for combined data
Melt_Protein_Transpose <- melt(Combined_Protein_Transpose)
graph_path_name<-paste(graph_path,"Normal_Distribution_Plot_of_Combined_Protein_Dataset_",kk,".tiff",sep="")
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

#Import Patient Info
#Gets paths of all the Protein data text files
Patient_files <- list.files(data_path, pattern="*Sample_List_Information.csv", full.names=T, recursive=FALSE)
Patient_Info_Extraction(Patient_files)

List_of_Transposed_Data<-ls(all.names=TRUE,pattern="^Protein\\d_AR_Transpose$")
List_of_Patient_Info<-ls(all.names=TRUE,pattern="^Patient_Sample_Set\\d$")

# Appending Patient Info to Protein Data
for (dataset in 1:length(List_of_Transposed_Data)){
  tmp<-get(List_of_Transposed_Data[dataset])  
  #Appending Abundance Ratio string in column one, and appending rest of data
  tmp2<-cbind((row.names(tmp)),tmp,all=TRUE)
  colnames(tmp2)[1]<-'TMT.Label'
  #TMT.Label Cleaning up column to to get only the specific TMT label for each sample.
  tmp2$TMT.Label <- gsub(pattern='^Abundance\\sRatio\\s\\(log2\\)\\:\\s', '', tmp2$TMT.Label)
  tmp2$TMT.Label <- gsub(pattern='.....126.', '', tmp2$TMT.Label)
  tmp2$TMT.Label <- gsub(pattern="\\(", '', tmp2$TMT.Label)
  #Converts values in TMT.Label as factors to ensure merging with patient data works through easily
  tmp2$TMT.Label<-as.factor(tmp2$TMT.Label)
  rm(tmp)
  tmp<-get(List_of_Patient_Info[dataset])
  #Merging patient information with protein level dataset
  tmp2 <- merge(tmp, tmp2, by = 'TMT.Label',all=TRUE)
  #Adding batch number to dataset
  tmp2 <- cbind(as.data.frame(c(dataset)), tmp2)
  #Renaming batch number to 'set'
  colnames(tmp2)[1] = 'Batch'
  tmp2$Gender<-as.character(tmp2$Gender)
  assign(paste("Protein_Dataset",dataset,sep = ""),tmp2,envir=.GlobalEnv)
}
rm(tmp2)

Protein_Datasets<-ls(all.names=TRUE,pattern="^Protein_Dataset\\d$")
#Combining all three datasets to plot overall distribution
Combined_Data_with_Patient_Info<-bind_rows(lapply(as.list(Protein_Datasets),get))
Filtered_Combined_Data <- Combined_Data_with_Patient_Info[sapply(Combined_Data_with_Patient_Info, function(x) !any(is.na(x)))] 

#Normalising data (median centring)
median_proteins = aggregate(Filtered_Combined_Data[,-c(1:10)], by = list(Filtered_Combined_Data[,2]), FUN = median, na.rm = F)
Normalized_Protein_Data = Filtered_Combined_Data
Num_Unique_Batches<-length(unique(Filtered_Combined_Data$Batch))
Num_Unique_TMT_Labels<-length(unique(Filtered_Combined_Data$TMT.Label))
k<-Num_Unique_TMT_Labels
for (i in 1:(Num_Unique_Batches*Num_Unique_TMT_Labels)){

Normalized_Protein_Data[i:k,-c(1:10)] = Normalized_Protein_Data[i:k,-c(1:10)] - median_proteins[,-1]
k=k+1
}
rm(i,k)


# Hierarchial Clustering
tmp = Normalized_Protein_Data

#Combining patient id with different groups so that patient id can be seen in dendrogram
tmp2<-apply(tmp[,11:ncol(tmp)],2,as.numeric)
tmp3<- tmp2[sapply(tmp2, function(x) !any(is.na(x)))] 
CFS_cluster = hclust(dist(tmp3))
tmp$fam_id <- paste(Normalized_Protein_Data$sample_id, Normalized_Protein_Data$Family_group, sep="|")
tmp$CFS_id <- paste(Normalized_Protein_Data$sample_id, Normalized_Protein_Data$CFS, sep="|")
tmp$age_id <- paste(Normalized_Protein_Data$sample_id, Normalized_Protein_Data$Age, sep="|")
tmp$gender_id <- paste(Normalized_Protein_Data$sample_id, Normalized_Protein_Data$Gender, sep="|")
tmp$EDS_id <- paste(Normalized_Protein_Data$sample_id, Normalized_Protein_Data$EDS, sep="|")
tmp$ds_id <- paste(Normalized_Protein_Data$sample_id, Normalized_Protein_Data$disease.status, sep="|")
tmp$batch_id <- paste(Normalized_Protein_Data$sample_id, Normalized_Protein_Data$Batch, sep="|")
tmp$channel_id <- paste(Normalized_Protein_Data$sample_id, Normalized_Protein_Data$TMT.Label, sep="|")

#Color coding dendrogram by different groups
eds_col = c("#86B875", "#E495A5", "#7DB0DD")
CFS_cluster_id = myplclust(CFS_cluster, labels = as.character(tmp$sample_id))
CFS_cluster_fam = myplclust(CFS_cluster, labels = as.character(tmp$fam_id), lab.col= rainbow_hcl(5)[as.fumeric(as.character(tmp$Family_group))])
CFS_cluster_CFS = myplclust(CFS_cluster, labels = as.character(tmp$CFS_id), lab.col= rainbow_hcl(3)[as.fumeric(as.character(tmp$CFS))])
CFS_cluster_age = myplclust(CFS_cluster, labels = as.character(tmp$age_id), lab.col= rainbow_hcl(9)[as.fumeric(as.character(tmp$Age_group))])
CFS_cluster_gender = myplclust(CFS_cluster, labels = as.character(tmp$gender_id), lab.col= rainbow_hcl(2)[as.fumeric(as.character(tmp$Gender))])
CFS_cluster_EDS = myplclust(CFS_cluster, labels = as.character(tmp$EDS_id), lab.col= eds_col[as.fumeric(as.character(tmp$EDS))])
CFS_cluster_ds = myplclust(CFS_cluster, labels = as.character(tmp$ds_id), lab.col= rainbow_hcl(8)[as.fumeric(as.character(tmp$disease.status))])
CFS_cluster_set = myplclust(CFS_cluster, labels = as.character(tmp$set_id), lab.col= rainbow_hcl(3)[as.fumeric(as.character(tmp$set))])
CFS_cluster_channel = myplclust(CFS_cluster, labels = as.character(tmp$channel_id), lab.col= rainbow_hcl(9)[as.fumeric(as.character(tmp$TMT.Label))])






