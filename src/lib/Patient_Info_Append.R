Patient_Info_Append<-function(List_of_Transposed_Data,List_of_Patient_Info){
  
  for (dataset in 1:length(List_of_Transposed_Data)){
    tmp<-get(List_of_Transposed_Data[dataset])  
    #Appending Abundance Ratio string in column one, and appending rest of data
    tmp2<-cbind((row.names(tmp)),tmp,all=TRUE)
    colnames(tmp2)[1]<-'TMT.Label'
    #TMT.Label Cleaning up column to to get only the specific TMT label for each sample.
    tmp2$TMT.Label <- gsub(pattern='^Abundance.Ratio..log2....', '', tmp2$TMT.Label)
    tmp2$TMT.Label <- gsub(pattern='.....126.', '', tmp2$TMT.Label)
    # tmp2$TMT.Label <- gsub(pattern="\\(", '', tmp2$TMT.Label) For original file
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
  
  
  
}