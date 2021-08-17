#Process calibrated Maeda CCR4 CyTOF data

library(data.table)

setwd("~/proj/CyTOF/")

dir <- "data/maeda/processed/"

load(paste0(dir,"OUH04.RData"))
head(savedat)

row <- rownames(savedat)
col <- colnames(savedat)

fname <- c("OUH04_calib.csv",
           "OUH05_imputed.csv",
           "OUH10_calib.csv",
           "OUH16_calib.csv")

savename <- c("OUH04_calib.RData",
              "OUH05_calib.RData",
              "OUH10_calib.RData",
              "OUH16_calib.RData")


for(i in 1:length(fname)){
  savedat <- as.matrix(fread(paste0(dir,fname[i]),sep=","))
  rownames(savedat) <- row
  colnames(savedat) <- col
  
  save(savedat,file=paste0(dir,savename[i]))
}



