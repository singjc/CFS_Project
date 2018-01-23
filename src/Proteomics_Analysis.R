library(rstudioapi)
library(dplyr)
library(VennDiagram)
library(colorspace)

#Gets project path based on location of source script
project_path <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
#Gets path to where data is stored
data_path <- paste(project_path,'/Data/Proteomics_Dataset',sep="")
#Gets paths of all the Protein data text files
files <- list.files(data_path, pattern="*Proteins.txt", full.names=T, recursive=FALSE)
i = 1
for (file in files ){
  #Extract only useful and relevant data from protein data file
  assign(paste("Protein",i,sep=""),(read.table(file))[,-c(1,2,8:49, 89:94)])
  #Store current protein data into a static temporary variable
  tmp <- get(paste("Protein",i,sep=""))
  #Rename column headers baser on names in the first row
  colnames(tmp) <-as.character(unlist(tmp[1, ])) 
  tmp <- tmp[-1, ]
  assign((paste("Protein",i,sep="")),tmp)
  #Extract only Abundance Ratios
  assign(paste("Protein",i,"_AR",sep=""),get(paste("Protein",i,sep=""))[c(1,5,16:24)])
  #Store current protein dataset into a static temporary variable
  tmp <- get(paste("Protein",i,"_AR",sep=""))
  #Convert Unique Sequence ID to character vectors
  tmp[["Unique Sequence ID"]] <- as.character(tmp[["Unique Sequence ID"]])
  #Re-assign converted tmp protein data to current Protein_AR variable
  assign((paste("Protein",i,"_AR",sep="")),tmp)
  #Store current protein dataset into a static temporary variable
  tmp <- get(paste("Protein",i,"_AR",sep=""))
  #Convert Sequence to character vectors
  tmp[["Sequence"]] <- as.character(tmp[["Sequence"]])
  #Re-assign converted tmp protein data to current Protein_AR variable
  assign((paste("Protein",i,"_AR",sep="")),tmp)
 i=i+1
 #Removing tmp variable from environment
 rm(tmp)
}
#Removing i variable from environment
rm(i)

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

for (kk in length(Protein_AR_Datasets)){
  if (kk==1){cat("One")}
  else if(kk==2){cat("Two")}
  else if(kk==3){cat("THREE")}
  
  tmp<-get(Protein_AR_Datasets[kk])
  Seq_Col<-tmp[,c(1,2)]
  assign((paste("Protein",kk,"_AR_Transpose",sep="")),as.data.frame(t(tmp[,-c(1,2)])))
  tmp2<-get(paste("Protein",kk,"_AR_Transpose",sep=""))
  colnames(tmp2)<- Seq_Col[,1]
  assign((paste("Protein",kk,"_AR_Transpose",sep="")),tmp2)
  tmp3<-colnames(get(paste("Protein",kk,"_AR",sep="")))
  tmp4<-cbind(as.data.frame(tmp3[-c(1:2)]),tmp2)
  assign((paste("Protein",kk,"_AR_Transpose",sep="")),tmp4)
  tmp5<-get(paste("Protein",kk,"_AR_Transpose",sep=""))
  colnames(tmp5)[1]<- 'Sample_ID'
  assign((paste("Protein",kk,"_AR_Transpose",sep="")),tmp5)
  
} 
rm(tmp,tmp2,tmp3,tmp4,tmp5)

