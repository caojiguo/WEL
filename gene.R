
#smoothing and computing rho
dyncor.ir <- function(Y1,Y2,time)
{ #Y1= X1+err1; Y2= X2+err2; time=tt; bdw <- h[l]
  n <- dim(Y1)[1]
  m <- length(time)
  Y1_pred <- matrix(NA, nrow=n, ncol=m)
  Y2_pred <- matrix(NA, nrow=n, ncol=m)
  for(i in 1:n)
  { #i = 1
    # cat("i = ", i, "\n")
    index <- which(!is.na(Y1[i,]))
    u <- time[index]
    cvBwSel <-  gam(as.numeric(Y1[i,index]) ~ s(u, bs="cr"), method="REML")
    Y1_pred[i,] <-  predict(cvBwSel, newdata= data.frame(u = time))
    cvBwSel <- gam(as.numeric(Y2[i,index]) ~ s(u, bs="cr"), method="REML")
    Y2_pred[i,] <- predict(cvBwSel, newdata= data.frame(u = time))
  }
  
 
  for(j in 1:m)
  {
    
    Y1_pred[,j] <- Y1_pred[,j]-mean(Y1_pred[,j])
    Y2_pred[,j] <- Y2_pred[,j]-mean(Y2_pred[,j])
  }
  
  for(j in 1:m)
  {
    Y1_pred[,j] <- Y1_pred[,j] - apply(Y1_pred,1,trapz,x=time)
    Y2_pred[,j] <- Y2_pred[,j] - apply(Y2_pred,1,trapz,x=time)
  }
  
  M1 <- Y1_pred^2
  M2 <- Y2_pred^2
  
  for(j in 1:m)
  {
    Y1_pred[,j] <- Y1_pred[,j]/sqrt(apply(M1,1,trapz,x=time))
    Y2_pred[,j] <- Y2_pred[,j]/sqrt(apply(M2,1,trapz,x=time))
  }
  
  rho <- apply(Y1_pred*Y2_pred,1,trapz,x=time)
  return (mean(rho))
  
}

#===========================================
#data analysis
genedat <- read.table("gene58rep34.txt",head=T)
gene9 <- as.matrix(genedat[genedat$gene=="F8",])
gene45 <- as.matrix(genedat[genedat$gene=="F44",])

library(fda)
#obs <- c(0,2,4,6,8,18,24,32,48,72)/72 
windows()
pdf("gene9prof.pdf")
fda::matplot(x=c(0,2,4,6,8,18,24,32,48,72), y=t(genedat[genedat$gene=="F8",-1]), type="l", 
             col="black", ylab="Gene expression", xlab="Time")
matpoints(x=c(0,2,4,6,8,18,24,32,48,72), y=t(genedat[genedat$gene=="F8",-1]), pch=4)
dev.off()

pdf("gene45prof.pdf")
fda::matplot(x=c(0,2,4,6,8,18,24,32,48,72), y=t(genedat[genedat$gene=="F44",-1]), type="l", 
             col="black", ylab="Gene expression", xlab="Time")
matpoints(x=c(0,2,4,6,8,18,24,32,48,72), y=t(genedat[genedat$gene=="F44",-1]), pch=4)
dev.off()

source("CI.R")
source("scel.R")
tau=10^(-4);

library(locpol)
library(caTools)


obs <- c(0,2,4,6,8,18,24,32,48,72)/72 
tt <- 0:72/72

#define local linear smoothing


m <- length(tt)
n <- nrow(gene9)
gene9_pred <- gene45_pred <- matrix(0,nrow=n,ncol=m)



for(i in 1:n)
{ #i = 1
  cvBwSel <- gam(as.numeric(gene9[i,-1]) ~ s(obs, bs="cr"), method="REML")
  gene9_pred[i,] <- predict(cvBwSel, newdata= data.frame(obs = tt))
  cvBwSel <- gam(as.numeric(gene45[i,-1]) ~ s(obs, bs="cr"), method="REML")
  gene45_pred[i,] <- predict(cvBwSel, newdata= data.frame(obs = tt))
}


Sd.gene9 <- gene9_pred
Sd.gene45 <- gene45_pred

for(j in 1:m)
{
  Sd.gene9[,j] <- gene9_pred[,j]-mean(gene9_pred[,j])
  Sd.gene45[,j] <- gene45_pred[,j]-mean(gene45_pred[,j])
  
}

for(j in 1:m)
{
  Sd.gene9[,j] <- Sd.gene9[,j] - apply(Sd.gene9,1,trapz,x=tt)
  Sd.gene45[,j] <- Sd.gene45[,j] - apply(Sd.gene45,1,trapz,x=tt)
}

M1 <- Sd.gene9^2
M2 <- Sd.gene45^2


for(j in 1:m)
{
  Sd.gene9[,j] <- Sd.gene9[,j]/sqrt(apply(M1,1,trapz,x=tt))
  Sd.gene45[,j] <- Sd.gene45[,j]/sqrt(apply(M2,1,trapz,x=tt))
}



rhog9g45 <- apply(Sd.gene9*Sd.gene45,1,trapz,x=tt)

#point estimate
mean(rhog9g45)
#0.36

#confidence interval based on empirical likelihood
CIrhog9g45 <- CIlen(rhog9g45,bs_value=qchisq(0.95,1))$CI
#0.16 0.52


#bootstrap analysis

for(i in 1:n)
{ #i = 1
  cvBwSel <- gam(as.numeric(gene9[i,-1]) ~ s(obs, bs="cr"), method="REML")
  gene9_pred[i,] <- predict(cvBwSel, newdata= data.frame(obs = tt))
  cvBwSel <- gam(as.numeric(gene45[i,-1]) ~ s(obs, bs="cr"), method="REML")
  gene45_pred[i,] <- predict(cvBwSel, newdata= data.frame(obs = tt))
}

g9 <- g45 <- matrix(NA, nrow=n, ncol=m)
g9[, match(obs, tt)] <- as.numeric(gene9[,-1])
g45[, match(obs, tt)] <- as.numeric(gene45[,-1])


resid_g9 <- g9 - gene9_pred
resid_g45 <- g45 - gene45_pred


B <- 500
rhog9g45.boot<-numeric(B)
set.seed(123)
for(b in 1:B)
{#b = 1
  cat("b = ", b, "\n")
  index1 <- sample(1:n, size=n, replace=T)
  X1 <- gene9_pred[index1,]
  X2 <- gene45_pred[index1,]
  
  index2 <- sample(1:n, size=n, replace=T)
  err1 <- resid_g9[index2,]
  err2 <- resid_g45[index2,]
  
  rhog9g45.boot[b] <- dyncor.ir(Y1= X1+err1, Y2= X2+err2, time=tt)
}

quantile(rhog9g45.boot, probs=c(0.025,0.975))
#0.10 0.51




