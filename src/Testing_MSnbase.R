# source("https://bioconductor.org/biocLite.R")
# biocLite("MSnbase")
library("MSnbase")

exprs=t(data)
fdata=as.data.frame(t(rownames(data)))
# colnames(fdata)<-'Gene_Name'
pdata=as.data.frame(Norm_Meta_Data[,-c(1)])

MS_Set_Data<-MSnSet(exprs,fdata,pdata)

library("limma")

design <- model.matrix(~ CFS + Age_group + Family_group + Gender + Batch, data = Protein_Dataframe)
v <- voom(t(res.comp$completeObs), design, plot = TRUE)