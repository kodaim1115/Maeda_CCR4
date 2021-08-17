#Maeda CyTOF preprocess; Scaling 0 to 1
setwd("~/proj/CyTOF/")
library(matrixStats)
library(data.table)

patients <- c("OUH04","OUH05","OUH16","OUH10")

Y <- list()
for(i in 1:length(patients)){
  load(paste0("data/maeda/processed/",patients[i],".RData"))
  Y[[i]] <- savedat
}
sapply(Y,nrow)

markers <- colnames(Y[[1]])
# [1] "EOMES"  "CD40"   "CD45RA" "CD8a"   "CD11c"  "CD14"   "pLck"   "OX40"  
# [9] "CD107a" "CCR4"   "CD3"    "PD1"    "CD86"   "pSTAT3" "Foxp3"  "Tbet"  
# [17] "CTLA4"  "CD80"   "CXCR3"  "CXCR5"  "LAG3"   "pNFkB"  "CCR7"   "CD40L" 
# [25] "CD25"   "CCR8"   "CD20"   "pS6"    "PDL1"   "CD4"


#Check distribution of each marker before scaling
breaks <- seq(0,10,by=0.1)
for(i in 1:length(markers)){
  setEPS()
  postscript(paste0("data/maeda/check_distribution_",markers[i],".eps"))
  par(mfrow=c(2,2))
  for(j in 1:length(Y)){
    hist(Y[[j]][,markers[i]],main="",xlim=c(0,10),breaks=breaks)
  }
  dev.off()
}

for(i in 1:length(Y)){
  data <- Y[[i]]
  rng <- colQuantiles(data,probs=c(0.01,0.99))
  savedat <- t((t(data)-rng[,1])/(rng[,2]-rng[,1]))
  savedat[savedat<0] <- 0
  savedat[savedat>1] <- 1
  save(savedat,file=paste0("data/maeda/processed/",patients[i],"_scale.RData"))
}
hist(savedat[,2])

for(i in 1:length(patients)){
  load(paste0("data/maeda/processed/",patients[i],"_scale.RData"))
  write.table(savedat,file=paste0("data/maeda/processed/",patients[i],".csv"),
              sep=",", row.names=FALSE, col.names=FALSE)
}

