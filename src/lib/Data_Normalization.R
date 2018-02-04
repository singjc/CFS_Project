Data_Normalization<-function(Protein_Datasets,graph_path,plot){
  # Combine all datasets into one data frame
  Combined_Data<-bind_rows(lapply(as.list(Protein_Datasets),get))
  tmp<-Combined_Data
  # Make a new label column to indentify individual sets and respective TMT labels
  Set_TMT_Label<-do.call(paste, c(tmp[,c(1,3)], sep = "set, "))
  data<-cbind(Set_TMT_Label,tmp)
  if (plot == T){
    data_reshaped <- melt(data,id.vars = "Set_TMT_Label",measure.vars = -c(1:12))
    graph_name<-"Boxplot of all runs, non-normalized"
    ggplot(data_reshaped, aes(x=Set_TMT_Label, y=value)) + geom_boxplot() + theme(axis.text.x  = element_text(angle=45, vjust=0.5))+ggtitle(graph_name)+xlab("Channel per Batch run")+theme(plot.title = element_text(hjust = 0.5))+ylab("Protein Abundance Ratio")
    graph_path_name<-paste(graph_path,graph_name,".tiff",sep="")
    ggsave(graph_path_name,device='tiff')
  }
  # Obtain mean for each individual batch/run/dataset
  Single_Batch_Mean<-as.data.frame(rowMeans(as.matrix(tmp[,-c(1:12)]),na.rm=TRUE))
  # Obtain the mean between each of the batches
  Mean_btwn_Batches<-colMeans(Single_Batch_Mean)
  # Get an adjusted mean for each bach-sample
  Adjusted_Means<-as.matrix(sweep(Single_Batch_Mean,2,Mean_btwn_Batches,"-"))
  dataMatrix<-as.matrix(tmp[,-c(1:12)])
  NormData<-sweep(dataMatrix,1,Adjusted_Means,'-')
  NormDataFrame<<-cbind(data[,c(1:12)],NormData)
  
  if (plot==T){
    NormPlotData<-melt(NormDataFrame,id.vars = "Set_TMT_Label",measure.vars = -c(1:12))
    graph_name<-"Boxplot of all runs, Normalized"
    ggplot(NormPlotData, aes(x=Set_TMT_Label, y=value)) + geom_boxplot() + theme(axis.text.x  = element_text(angle=45, vjust=0.5))+ggtitle(graph_name)+xlab("Channel per Batch run")+theme(plot.title = element_text(hjust = 0.5))+ylab("Protein Abundance Ratio")
    graph_path_name<-paste(graph_path,graph_name,".tiff",sep="")
    ggsave(graph_path_name,device='tiff')
  }
  
}