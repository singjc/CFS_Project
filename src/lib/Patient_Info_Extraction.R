#This function is used to extract patient info form csv files.
Patient_Info_Extraction <- function(Patient_files) {
  i<-1
  for (patient_file in Patient_files){
    tmp<-read.csv(patient_file,header=TRUE,colClasses = c("character"))
    assign(paste('Patient_Sample_Set',i,sep=""),tmp,envir=.GlobalEnv)
    i=i+1
  }
  rm(i,tmp)
}


