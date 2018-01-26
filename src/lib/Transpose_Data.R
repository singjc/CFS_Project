
Transpose_Data <-function(Protein_AR_Datasets,Variable_Name){
  #Transposing data such that protein id and subject information can be found in individual columns, making it easier to commbine datasets later on
  for (kk in 1:length(Protein_AR_Datasets)){
    tmp<-get(Protein_AR_Datasets[kk])
    #Obtain Sequence column information
    Seq_Col<-tmp[,c(1,2)]
    #Obtain Abundance Ration Sample IDS
    AR_ID<-colnames(tmp[-c(1,2)])
    assign((paste("Protein",kk,Variable_Name,sep="")),as.data.frame(t(tmp[,-c(1,2)])))
    tmp2<-get(paste("Protein",kk,Variable_Name,sep=""))
    #Convert dataframe to a matrix, to convert numbers to numeric data
    tmp2_matrix<-as.matrix(tmp2)
    row.names(tmp2_matrix)<-NULL
    tmp2_matrix<-apply(tmp2_matrix,2,as.numeric)
    tmp2<-as.data.frame(tmp2_matrix)
    row.names(tmp2)<-AR_ID
    colnames(tmp2)<- Seq_Col[,1]
    assign((paste("Protein",kk,Variable_Name,sep="")),tmp2,envir=.GlobalEnv)
  } 
  rm(tmp,tmp2,tmp2_matrix,kk)
  return(Seq_Col)
}