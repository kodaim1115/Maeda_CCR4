#For supcomp analysis; Maeda CCR4 CyTOF data; CALIBRATED
#LATEST VERSION
library(Rcpp)
library(RcppArmadillo)

sourceCpp("CYBERTRACK2.cpp")

file <- c("OUH04_calib.RData",
          "OUH05_calib.RData",
          "OUH10_calib.RData",
          "OUH16_calib.RData")


Y <- t_id <- list()
samplesize <- c()
for(i in 1:length(file)){
  load(file[i])
  Y[[i]] <- savedat
  t_id[[i]] <- as.integer(rownames(savedat)) 
  samplesize[i] <- length(t_id[[i]])
}
D <- length(Y)
T <- length(sort(unique(t_id[[1]])))
cumsumsize <- cumsum(samplesize)

L_vec <- c(2:40)
res <- NULL
for(k in 1:length(L_vec)){
  res <- rbind(res,cbind(L_vec[k]))
}
L <- as.numeric(res[ugeid,1])

kmY <- do.call(rbind,Y)

kminit <- function(y,L,seed = sample.int(.Machine$integer.max, 1)){
  set.seed(seed)
  kmres <- kmeans(y,L,iter.max=100,algorithm="Lloyd")
  list(mean=t(kmres$centers),
       var=simplify2array(lapply(split(as.data.frame(y),kmres$cluster),var)),
       cluster=kmres$cluster)
}
kmpar <- kminit(kmY,L,123)
12
K <- ncol(kmY)
num_iter <- 100
num_iter_refine <- 50
wis_iter <- 20
tau <- 1e-5
nu <- K+2
Lambda <- diag(K)
piini <- list()
for(d in 1:D) piini[[d]] <- matrix(1/L,T,L)
alphaini <- c(rep(1,T))

Wini <- list()
for(d in 1:D){
  Wini[[d]] <- matrix(0,nrow(Y[[d]]),L)
  Wini[[d]][,1] <- 1 
}
muini <- kmpar$mean
Sigmaini <- kmpar$var
for(l in 1:L) Sigmaini[,,l] <- Sigmaini[,,l]+(1e-5*diag(K))

gzi=0.01
P <- 1

result <- LONGINUS(Y,L,D,P,Wini,piini,alphaini,muini,Sigmaini,tau,nu,gzi,Lambda,num_iter,num_iter_refine,wis_iter,t_id)
save(result,file=paste0("/archive/data/rwdc/kodaim/20191231_Maeda_CyTOF_",L,".RData"))

