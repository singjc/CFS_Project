FDR_Adjusted_PVal <- function(Results,variables) {
  # x$BH_Y = 
  #   p.adjust(Results$yes_raw.p, 
  #            method = "BH")
  # x$BH_U = 
  #   p.adjust(Results$u_raw.p, 
  #            method = "BH")
  # 
  adjusted_PVal<-data.frame(matrix(NA,nrow=nrow(Results),ncol=ncol(Results)))
  factor_levels<-colnames(Results)
  for (i in 1:length(factor_levels)){
    adjusted_PVal[,i]<-p.adjust(Results[,i],method="BH")
  }
  variables<-strsplit(colnames(Results),"_Raw_P-Value") 
  colnames(adjusted_PVal)<-paste(variables,"_Adjusted_P-value",sep="")
  return(adjusted_PVal)
}