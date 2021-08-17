#Analysis on supcomp result for Maeda CCR4 CyTOF; NO SCALE; CALIBRATED
#Figure
#LATEST VERSION

setwd("~/proj/CyTOF")
library(Rcpp)
library(gplots)
library(pheatmap)
library(tidyverse)
library(RColorBrewer)
library(umap)
sourceCpp("calc_SSE.cpp")

set.seed(123)

fname <- list.files("data/maeda/processed",full.names=TRUE)
fname <- fname[grep("no_scale_calib.RData",fname)]
fname

#min clusters
min <- 2
#max clusters
max <- 40
step <- max-min+1
L <- seq(min,max,length=step)
Y <- t_id <- list()
for(i in 1:length(fname)){
  load(fname[i])
  Y[[i]] <- savedat
  t_id[[i]] <- as.integer(rownames(savedat)) 
}
markers <- colnames(Y[[1]])
# markers <- gsub("....._","",markers)
# markers <- gsub("..._","",markers)
# markers <- gsub("CD4.CD4","CD4",markers)
markers

D <- length(Y)
wks <- sort(unique(t_id[[1]]))
T <- length(sort(unique(t_id[[1]])))
kmY <- do.call(rbind,Y)
nrow(kmY)

SSE <- c()
for(i in 1:step){
  load(paste0("data/supcomp/20200123_Maeda_CyTOF_",i+1,".RData"))
  mergeZ <- do.call(rbind,result$Z)
  mergeW <- do.call(rbind,result$W)
  SSE[i] <- calc_SSE_v2(kmY,mergeW,L[i])
}

setEPS()
postscript("data/fig/20200207_Maeda_CyTOF_SSE.eps")
par(pin=c(4.5,4.5))
plot(c(min:max),SSE,xlab="Number of clusters", ylab="SSE",
     cex.lab=3,cex.axis=2,type="b")
dev.off()

#cluster number
nclust <- 38
load(paste0("data/supcomp/20200123_Maeda_CyTOF_",nclust,".RData"))
pi <- result$pi
mu <- result$mu

rownames(mu) <- markers
colnames(mu) <- as.character(c(1:ncol(mu)))

breaks <- seq(-1,1,by=0.01)
color <- colorRampPalette(rev(brewer.pal(n=7 , name="RdYlBu")))(length(breaks))

setEPS()
postscript(paste0("data/fig/20200207_Maeda_CyTOF_heatmap_",nclust,"_clusters.eps"),height=10,width=12)
pheatmap(mu,color=color,breaks=breaks,scale="row",
         display_numbers=TRUE,
         legend_breaks=c(-1,1),legend_labels=c("Low","High"))
dev.off()


#pi
junk_pi <- do.call(rbind,pi)
setEPS()
postscript(paste0("data/fig/20200207_Maeda_CyTOF_pi_",nclust,"_clusters.eps"),height=15,width=15)
par(mfrow=c(5,8),mar=c(3,4,2,1)+0.1,oma=c(0,0,0,0), mgp=c(2.5,1,0),las=1)

for(i in 1:nclust){
  if(i==1){
    plot(0,0,xlim=c(1,T),ylim=c(0,1.2*max(junk_pi[,i])),main=paste0("cluster ",i),
         xlab="",ylab="",type="n",xaxt="n",cex.lab=2,cex.axis=1.3)
    legend("topright",col=1:D,lwd=1,legend=c("OUH04","OUH05","OUH10","OUH16"))
  } else{
    plot(0,0,xlim=c(1,T),ylim=c(0,1.2*max(junk_pi[,i])),main=paste0("cluster ",i),
         xlab="",ylab="",type="n",xaxt="n",cex.lab=2,cex.axis=1.3)
  }
  axis(1,at=c(1:T),labels=c("pre","post"),cex.axis=2)
  for(j in 1:D){
    lines(pi[[j]][,i],lwd=1,col=j)
  }
  
}
dev.off()


#remove abberant cluster
remove <- c(7,13,15,18,20,23,25,26,29,30,31,37)
sub_mu <- mu[,-remove]
colnames(sub_mu) <- as.character(c(1:ncol(sub_mu)))
setEPS()
postscript(paste0("data/fig/20200207_Maeda_CyTOF_heatmap_",nclust,"_sub_clusters.eps"),height=10,width=12)
pheatmap(sub_mu,color=color,breaks=breaks,scale="row",
         display_numbers=TRUE,
         legend_breaks=c(-1,1),legend_labels=c("Low","High"))
dev.off()

#pi with boxplot
sub_pi <- list()
for(i in 1:length(pi)){
  sub_pi[[i]] <- pi[[i]][,-remove] 
}
junk_sub_pi <- do.call(rbind,sub_pi)
  
setEPS()
postscript(paste0("data/fig/20200207_Maeda_CyTOF_pi2_",nclust,"_clusters.eps"),height=15,width=15)
par(mfrow=c(5,6),mar=c(3,4,2,1)+0.1,oma=c(0,0,0,0), mgp=c(2.5,1,0),las=1)
box <- list()
for(i in 1:ncol(sub_mu)){
  #if(i %in% remove)
  #  next
  box[[1]] <- junk_sub_pi[,i][c(1,3,5,7)]
  box[[2]] <- junk_sub_pi[,i][c(2,4,6,8)]
  boxplot(box,
          #xlim=c(1,T),
          ylim=c(0,1.2*max(junk_sub_pi[,i])),main=paste0("cluster ",i),
          xlab="",ylab="",type="n",xaxt="n",cex.main=2,cex.lab=1.5,cex.axis=1.3)
  #plot(0,0,xlim=c(1,T),ylim=c(0,1.2*max(junk_pi[,i])),main=paste0("cluster ",i),
  #      xlab="",ylab="",type="n",xaxt="n",cex.lab=2,cex.axis=1.3)
  axis(1,at=c(1:T),labels=c("pre","post"),cex.axis=1.5)
  for(j in 1:D){
    lines(sub_pi[[j]][,i],lwd=1,col=j)
  }
  ttest <- t.test(box[[1]], box[[2]], paired=TRUE)
  if(ttest$p.value<0.05)
    text(x=T, y=max(junk_sub_pi[,i]), labels="*", pos=3, cex=2)
  else if(ttest$p.value<0.01)
    text(x=T, y=max(junk_sub_pi[,i]), labels="**", pos=3, cex=2)
  else if(ttest$p.value<0.001)
    text(x=T, y=max(junk_sub_pi[,i]), labels="***", pos=3, cex=2)
}
dev.off()

#Ratio
pi_ratio <- pi
for(i in 1:length(pi)){
  tmp_ratio <- pi[[i]]
  ratio <- tmp_ratio[2,]/tmp_ratio[1,]
  pi_ratio[[i]] <- ratio
}
pi_ratio

pi_ratio_mat <- matrix(0,length(pi_ratio),nclust)
for(i in 1:length(pi_ratio)){
  pi_ratio_mat[i,] <- pi_ratio[[i]]
}


setEPS()
postscript("data/fig/20200207_Maeda_CyTOF_pi_bar.eps",height=15,width=15)
par(mfrow=c(5,6),mar=c(3,4,2,1)+0.1,oma=c(0,0,0,0), mgp=c(2.5,1,0),las=1)

for(i in 1:nclust){
  if(i==1){
    barplot(pi_ratio_mat[,i],col=1:D,main=paste0("cluster ",i))
    legend("topright",col=1:D,pch=15,legend=c("OUH04","OUH05","OUH10","OUH16"))
  } else{
    barplot(pi_ratio_mat[,i],col=1:D,main=paste0("cluster ",i))
  }
  axis(1,at=c(1:T),labels=c("pre","post"),cex.axis=1.5)
  abline(h=1,lty=2)
}
dev.off()


#Explore result

len <- nrow(Y[[1]])
col <- rep(1:4,each=len)
mergeZ <- do.call(rbind,result$Z)
colnames(mergeZ) <- markers
mergeW <- do.call(rbind,result$W)

fun <- function(x) which(x==1)
sel_id <- apply(mergeW,1,fun)
cell_id <- which(sel_id == 10)

sel_dat <- mergeZ[cell_id,]
plot(sel_dat[,"CD3"],sel_dat[,"CD20"],pch=20,cex=0.1,col=col[cell_id])

#umap
sub_id <- sub_kmY <- sub_W <- list()
for(i in 1:length(Y)){
  sub_id[[i]] <- c(sample(c(1:121378),50000,replace=FALSE), sample(c(121379:242756),50000,replace=FALSE))
  sub_kmY[[i]] <- Y[[i]][sub_id[[i]],]
  sub_W[[i]] <- result$W[[i]][sub_id[[i]],]
}
sub_kmY <- do.call(rbind,sub_kmY)
dim(sub_kmY)

cluster <- do.call(rbind,sub_W)
fun <- function(x) which(x==1)
cluster <- apply(cluster,1,fun)

umap = umap(sub_kmY, method='umap-learn')
#save(umap, file="data/fig/umap.RData")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

col = gg_color_hue(length(unique(cluster)))

#All cell
setEPS()
postscript("data/fig/20200207_Maeda_CyTOF_umap.eps",height=5,width=5)
plot(umap$layout, pch=20, cex=0.1, col=col[cluster], xlab="UMAP 1", ylab="UMAP 2")
dev.off()


#Remove abberant
remove_id <- which(cluster %in% remove)
sub_umap <- umap$layout[-remove_id,]
dim(sub_umap)  
sub_cluster <- cluster[-remove_id]
temp <- rep(0,length(sub_cluster))

#Make sub_cluster continuous
for(i in 1:length(unique(sub_cluster))){
  id <- which(sub_cluster==sort(unique(sub_cluster))[i])
  temp[id] <- i
}
sort(unique(temp))
sub_cluster <- temp

col <- gg_color_hue(length(unique(sub_cluster)))
#unique <- sort(unique(sub_cluster))
#col <- match(sub_cluster, unique)

setEPS()
postscript("data/fig/20200207_Maeda_CyTOF_sub_umap.eps",height=5,width=5)
plot(sub_umap, pch=20, cex=0.1, col=col[cluster], xlab="UMAP 1", ylab="UMAP 2")
legend("right", legend=as.character(c(1:ncol(sub_mu))), 
       col=col, fill=col, pch=1, cex=0.5)
dev.off()

#Pre post

pre_id <- which(rownames(sub_kmY)==0)
pre_id <- pre_id[-which(pre_id %in% remove_id)]
length(pre_id)

post_id <- which(rownames(sub_kmY)==1)
post_id <- post_id[-which(post_id %in% remove_id)]
length(post_id)

length(c(pre_id,post_id))
dim(sub_umap)

pre_umap <- umap$layout[pre_id,]
post_umap <- umap$layout[post_id,]

setEPS()
postscript("data/fig/20200207_Maeda_CyTOF_sub_umap.eps",height=5,width=5)
plot(pre_umap, pch=20, cex=0.1, col=col[cluster], xlab="UMAP 1", ylab="UMAP 2", main="Pre")
legend("right", legend=as.character(c(1:nclust)[-remove]), 
       col=col[sub_cluster], fill=col[-remove], pch=1, cex=0.5)
dev.off()

setEPS()
postscript("data/fig/20200207_Maeda_CyTOF_sub_umap.eps",height=5,width=5)
plot(post_umap, pch=20, cex=0.1, col=col[cluster], xlab="UMAP 1", ylab="UMAP 2", main="Pre")
legend("right", legend=as.character(c(1:nclust)[-remove]), 
       col=col[sub_cluster], fill=col[-remove], pch=1, cex=0.5)
dev.off()

