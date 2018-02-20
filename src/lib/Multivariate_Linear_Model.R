Multivariate_Linear_Model <- function(Protein_Dataframe) {
  # Protein_Dataframe<-NormDataFrame
  # Protein_Dataframe<-Missing_Data_Fill
  tmp = list()
  for(i in names(Protein_Dataframe)[-c(1:8)]){
    print(i)
    tmp[[i]] <- lm(Protein_Dataframe[,i] ~  CFS + Age_group + Family_group + Gender + Batch
                   , data=Protein_Dataframe, na.action=na.exclude)
  }
  # tmp1 = data.frame()
  # for (i in 1:length(tmp)) {
  #   tmp1[i,1] = summary(tmp[[i]])$coefficients[2,4]
  # }
  # for (i in 1:length(tmp)) {
  #   tmp1[i,2] = summary(tmp[[i]])$coefficients[3,4]
  # }
  # colnames(tmp1) = c('yes_raw.p', 'u_raw.p')
  # return(tmp1)
  
  tmp2=data.frame()
  for (i in 1:length(tmp)){
    factor_levels<-rownames(summary(tmp[[i]])$coefficients)
    for (k in 1:length(factor_levels)){
      tmp2[i,k]<-summary(tmp[[i]])$coefficients[k,4]
    }
  }
  colnames(tmp2)<-c(paste(factor_levels,"Raw_P-Value",sep="_"))
  return(tmp2)
}