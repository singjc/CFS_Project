
i<-bind_rows(Protein_Dataset1,Protein_Dataset2,Protein_Dataset3,Protein_Dataset4,Protein_Dataset5,Protein_Dataset6)

boxplot(t(i[,-c(1:11)]),main="Boxplot of Non-Normalized Combined Protein Dataset",xlab="Protein",ylab="Abundance Ratio",names=do.call(paste, c(i[,1:2], sep = "set, ")) ,las=3)

dat<-i
new_col<-do.call(paste, c(i[,1:2], sep = "set, "))
dat<-cbind(new_col,dat)

data <- melt(dat,id.vars = "new_col",measure.vars = -c(1:12))
graph_name<-"Boxplot of all runs, non-normalized"
ggplot(data, aes(x=new_col, y=value)) + geom_boxplot() + theme(axis.text.x  = element_text(angle=45, vjust=0.5))+ggtitle(graph_name)+xlab("Channel per Batch run")+theme(plot.title = element_text(hjust = 0.5))+ylab("Protein Abundance Ratio")

med<- aggregate(data[,-c(1:2)], by = list(data[,1]), FUN = median, na.rm = TRUE)
eee<-dat[,-c(1:12)]-med[-1]

df_new <- sweep(dat[,-c(1:12)],1,med[-1],"-")

#Attemp1
# centering with 'scale()'
center_scale <- function(x) {
  scale(x, scale = FALSE)
}
NormDATA<-as.data.frame( center_scale(dat[,-c(1:12)]))
NormDataFrame<-cbind(dat[,c(1:12)],NormDATA)
NormPlotData<-melt(NormDataFrame,id.vars = "new_col",measure.vars = -c(1:12))
graph_name<-"Boxplot of all runs, Normalized"
ggplot(NormPlotData, aes(x=new_col, y=value)) + geom_boxplot() + theme(axis.text.x  = element_text(angle=45, vjust=0.5))+ggtitle(graph_name)+xlab("Channel per Batch run")+theme(plot.title = element_text(hjust = 0.5))+ylab("Protein Abundance Ratio")

# Attempt2
# center with 'apply()'
center_apply <- function(x) {
  apply(x, 2, function(y) y - mean(y))
}
NormDATA<-as.data.frame( center_apply(dat[,-c(1:12)]))
NormDataFrame<-cbind(dat[,c(1:12)],NormDATA)
NormPlotData<-melt(NormDataFrame,id.vars = "new_col",measure.vars = -c(1:12))
graph_name<-"Boxplot of all runs, Normalized"
ggplot(NormPlotData, aes(x=new_col, y=value)) + geom_boxplot() + theme(axis.text.x  = element_text(angle=45, vjust=0.5))+ggtitle(graph_name)+xlab("Channel per Batch run")+theme(plot.title = element_text(hjust = 0.5))+ylab("Protein Abundance Ratio")

# Attemp3
# center with 'colMeans()'
center_colmeans <- function(x) {
  xcenter = colMeans(x)
  x - rep(xcenter, rep.int(nrow(x), ncol(x)))
}

NormDATA<-as.data.frame( center_colmeans(dat[,-c(1:12)]))
NormDataFrame<-cbind(dat[,c(1:12)],NormDATA)
NormPlotData<-melt(NormDataFrame,id.vars = "new_col",measure.vars = -c(1:12))
graph_name<-"Boxplot of all runs, Normalized"
ggplot(NormPlotData, aes(x=new_col, y=value)) + geom_boxplot() + theme(axis.text.x  = element_text(angle=45, vjust=0.5))+ggtitle(graph_name)+xlab("Channel per Batch run")+theme(plot.title = element_text(hjust = 0.5))+ylab("Protein Abundance Ratio")

#Attempt 4
Single_Batch_Mean<-as.data.frame(rowMeans(as.matrix(dat[,-c(1:12)]),na.rm=TRUE))
col_Mean<-colMeans(Single_Batch_Mean)
new_mean<-as.matrix(sweep(Single_Batch_Mean,2,col_Mean,"-"))
dataTest<-as.matrix(dat[,-c(1:12)])
NormDATA<-apply(dataTest,MARGIN=2,function(y) y-new_mean)
NormData<-sweep(dataTest,1,new_mean,'-')
NormDataFrame<-cbind(dat[,c(1:12)],NormDATA)
NormPlotData<-melt(NormDataFrame,id.vars = "new_col",measure.vars = -c(1:12))
graph_name<-"Boxplot of all runs, Normalized"
ggplot(NormPlotData, aes(x=new_col, y=value)) + geom_boxplot() + theme(axis.text.x  = element_text(angle=45, vjust=0.5))+ggtitle(graph_name)+xlab("Channel per Batch run")+theme(plot.title = element_text(hjust = 0.5))+ylab("Protein Abundance Ratio")

