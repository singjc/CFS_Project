Filtering_Common_Proteins <-function(Protein_AR_Datasets,Variable_Name){
  #Filtering protein abundances for only commonly found proteins across all batches
  i=1
  num_datasets = 1:length(Protein_AR_Datasets)
  for(kk in 1:length(Protein_AR_Datasets)){
    #Temporarily store very first Protein abundance ratio data
    tmp <- get(Protein_AR_Datasets[kk]) 
    # if(kk==length(Protein_AR_Datasets)){break}
    #Get numbers that are to be compared against
    tmp2 <- num_datasets[num_datasets!=kk]
    for (jj in (tmp2)){
      tmp3 <- semi_join(tmp,(get(Protein_AR_Datasets[jj])), by = 'Unique Sequence ID')
      assign((paste("Protein",kk,Variable_Name,sep="")),tmp3,envir = .GlobalEnv)
      tmp<-get(Protein_AR_Datasets[kk])
    }
  }
  rm(tmp,tmp2,tmp3,kk,i,jj)
}