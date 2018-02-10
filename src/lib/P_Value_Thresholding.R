P_Value_Thresholding<-function(Adjusted_PVals){
  
  factor_levels<-colnames(Adjusted_PVals)[colnames(Adjusted_PVals)!="Gene_Name"]
  for (i in 1:length(factor_levels)){
    tmp<-data.frame(matrix(NA,nrow=length(Adjusted_PVals[(Adjusted_PVals[,factor_levels[i]]<0.05),"Gene_Name"]),ncol=2))
    tmp[,1]<-as.data.frame(Adjusted_PVals[(Adjusted_PVals[,factor_levels[i]]<0.05),"Gene_Name"])
    tmp[,2]<-as.data.frame(Adjusted_PVals[(Adjusted_PVals[,factor_levels[i]]<0.05),factor_levels[i]])
    variable=strsplit(factor_levels[i],"_Adjusted_P-value")
    colnames(tmp)<-c("Gene_Name",paste(variable,"_Sig.PVal",sep=""))
    assign(paste(variable,"_PVal_Sig",sep=""),tmp,envir=.GlobalEnv)
  }

}