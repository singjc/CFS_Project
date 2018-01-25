library(rstudioapi)
library(dplyr)
library(VennDiagram)
library(colorspace)
library(ggplot2)
library(reshape)
library(rlist)

source('~/Users/justinsing/Documents/Hannes_Rost/CFS/CFS_Project/src/lib/Patient_Info_Extraction.R ')

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
#Extract Protein Information
Dataset_Extraction<-Protein_Dataset_Info_Extraction(files)
# get(Results[4],.last_env)

#Get Protein abundance ratio datasets for each batch
Protein_AR_Datasets <- ls(all.names=TRUE,pattern="*\\d_AR$")
#Temporarily store very first Protein abundance ratio data
tmp <- get(Protein_AR_Datasets[1])
i=1
#Find common proteins between all 3 datasets using semi_join. Semi_join returns values that are matching in the same row for both variables.
for (kk in length(Protein_AR_Datasets)){
  if (kk==length(Protein_AR_Datasets)){break}
  tmp <- semi_join(tmp, get(Protein_AR_Datasets[kk+i]), by = 'Unique Sequence ID')
}
Common_Proteins <- tmp
rm(tmp,kk,i)

#Plotting venn diagram of the number of proteins found in each batch and common protein between batches
grid.newpage()
venn_protein = draw.triple.venn(area1 = 4887, area2 = 5593, area3 = 6618, n12= 4109, n23 = 4646, n13 = 4018, n123 = 3782, category = c('Set 1', 'Set 2', 'Set 3'), fill =  c("#E495A5" ,"#ABB065" ,"#39BEB1" ))

#Filtering protein abundances for only commonly found proteins across all batches
i=1
num_datasets = 1:length(Protein_AR_Datasets)
for(kk in length(Protein_AR_Datasets)){
  #Temporarily store very first Protein abundance ratio data
  tmp <- get(Protein_AR_Datasets[kk]) 
  if(kk==length(Protein_AR_Datasets)){break}
  #Get numbers that are to be compared against
  tmp2 <- num_datasets[num_datasets!=kk]
  for (jj in (tmp2)){
    tmp3 <- semi_join(tmp,get(get(Protein_AR_Datasets[jj]), by = 'Unique Sequence ID'))
    assign((paste("Protein",kk,"_AR",sep="")),tmp3)
    tmp<-get(Protein_AR_Datasets[kk])
  }
}
rm(tmp,tmp2,tmp3,kk,i,jj)

#Transposing data such that protein id and subject information can be found in individual columns, making it easier to commbine datasets later on
for (kk in 1:length(Protein_AR_Datasets)){
  tmp<-get(Protein_AR_Datasets[kk])
  #Obtain Sequence column information
  Seq_Col<-tmp[,c(1,2)]
  #Obtain Abundance Ration Sample IDS
  AR_ID<-colnames(tmp[-c(1,2)])
  assign((paste("Protein",kk,"_AR_Transpose",sep="")),as.data.frame(t(tmp[,-c(1,2)])))
  tmp2<-get(paste("Protein",kk,"_AR_Transpose",sep=""))
  #Convert dataframe to a matrix, to convert numbers to numeric data
  tmp2_matrix<-as.matrix(tmp2)
  row.names(tmp2_matrix)<-NULL
  tmp2_matrix<-apply(tmp2_matrix,2,as.numeric)
  tmp2<-as.data.frame(tmp2_matrix)
  row.names(tmp2)<-AR_ID
  
  # colnames(tmp2)<- Seq_Col[,1]
  # assign((paste("Protein",kk,"_AR_Transpose",sep="")),tmp2)
  # tmp3<-colnames(get(paste("Protein",kk,"_AR",sep="")))
  # tmp4<-cbind(as.data.frame(tmp3[-c(1:2)]),tmp2)
  # assign((paste("Protein",kk,"_AR_Transpose",sep="")),tmp4)
  # tmp5<-get(paste("Protein",kk,"_AR_Transpose",sep=""))
  
  colnames(tmp2)[1]<- 'Sample_ID'
  assign((paste("Protein",kk,"_AR_Transpose",sep="")),tmp2)
} 
rm(tmp,tmp2,tmp2_matrix,kk)

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
