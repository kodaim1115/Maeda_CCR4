#Create csv file for non scaled data; Maeda CyTOF

library(data.table)

setwd("~/proj/CyTOF/")

dir <- "data/maeda/processed/"

fname <- c("OUH04",
           "OUH05",
           "OUH10",
           "OUH16")

for(i in 1:length(fname)){
  load(paste0(dir,fname[i],".RData"))
  write.table(savedat,file=paste0(dir,fname[i],"_no_scale.csv"),
              sep=",",row.names=FALSE,col.names=FALSE)
}
