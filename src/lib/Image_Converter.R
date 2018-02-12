Image_Converter<-function(path){
  
  library("jpeg")
  library("tiff")
  
  # path<-"/Users/justinsing/Documents/Hannes_Rost/CFS/CFS_Project/Results/Graphs/Dendrograms/rnorm_Clustering/"
  files<-list.files(path,pattern="*.tiff$")
  
  for (file in files){
    path_file<-paste(path,file,sep='')
    print(path_file)
    img <- readTIFF(path_file, native=TRUE)
    name<-paste(path,strsplit(file,'.tiff'),".jpeg",sep='')
    writeJPEG(img, target = name, quality = 1)
  }
}