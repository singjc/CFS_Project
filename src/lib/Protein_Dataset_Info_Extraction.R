Protein_Dataset_Info_Extraction <- function(files) {
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
  
  Protein_Datasets <- ls(all.names=TRUE,pattern="^Protein\\d$")
  
  #Working on Returing variables from function
  return(sapply(Protein_Datasets,get))
}