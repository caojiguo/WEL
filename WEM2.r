
#=======================================
#empirical likelihood
source("scel.R")
source("CI.R")
tau=10^(-4);

#======================================================
tt <- seq(0,1,length=100)
m <- length(tt)
sigma_11 <- diag(c(1,1/2,1/3))
sigma_12 <- diag(c(1/3,1/4,1/6))
sigma_22 <- diag(c(1/2,1/3,1/4))
Sigma <- cbind(rbind(sigma_11,sigma_12),rbind(sigma_12,sigma_22))
eta_0 <- rep(1,m)
eta_1 <- 2*sqrt(3)*(tt-1/2)
eta_2 <- 6*sqrt(5)*(tt-1/2)^2 - 1/2*sqrt(5)
h <- seq(0.01,0.41,by=0.025)
L <- length(h)


#====================
library(mvtnorm)
library(KernSmooth)
library(caTools)
n <- 100  #repeat for n = 50
f_1 <- f_2 <- matrix(NA,nrow=n,ncol=m)
Y_1 <- Y_2 <- matrix(NA,nrow=n,ncol=m)
Y1_pred <- Y2_pred <- Y_1
f1_hat <- f2_hat <- matrix(NA,nrow=n,ncol=m)
rho.target <- 0.5
K <- 100
logelr0 <- numeric(K)
count <- numeric(L)
CI.len <- numeric(L)
len <- numeric(K)
set.seed(20170615)

#==============================================
#scenario 2: regular design

Sys.time()
for(l in 1:L)
{ #l = 9
  for(k in 1:K)
  { #k = 1
    eps <- rmvnorm(n, sigma = Sigma)
    e <- matrix(rnorm(n*m*2,mean=0,sd=1/2),nrow=n*2)
    
    for(j in 1:m)
    {
      f_1[,j] <- 1 + eps[,1] + eta_1[j]*eps[,2] + eta_2[j]*eps[,3]
      f_2[,j] <- 1 + eps[,4] + eta_1[j]*eps[,5] + eta_2[j]*eps[,6]
      Y_1[,j] <- f_1[,j]  + e[1:n,j]
      Y_2[,j] <- f_2[,j]  + e[-(1:n),j]
    }
    
    for(i in 1:n)
{ #i = 1
      fit1 <- locpoly(x=tt, y=Y_1[i,],bandwidth=h[l],gridsize=m)
      Y1_pred[i,] <- fit1$y
      fit2 <- locpoly(x=tt, y=Y_2[i,],bandwidth=h[l],gridsize=m)
      Y2_pred[i,] <- fit2$y
}
#      windows()
#      plot(x=tt, y=Y_1[1,],type="l",col="red")
#       lines(x=tt, y=Y1_pred[1,], col="blue")
    
    
    for(j in 1:m)
    {
      f1_hat[,j] <- Y1_pred[,j]-mean(Y1_pred[,j])
      f2_hat[,j] <- Y2_pred[,j]-mean(Y2_pred[,j])
    }
    
    for(j in 1:m)
    {
      f1_hat[,j] <- f1_hat[,j] - apply(f1_hat,1,trapz,x=tt)
      f2_hat[,j] <- f2_hat[,j] - apply(f2_hat,1,trapz,x=tt)
    }
    
    M1 <- f1_hat^2
    M2 <- f2_hat^2
    
    for(j in 1:m)
    {
      f1_hat[,j] <- f1_hat[,j]/sqrt(apply(M1,1,trapz,x=tt))
      f2_hat[,j] <- f2_hat[,j]/sqrt(apply(M2,1,trapz,x=tt))
    }
    
    rho <- apply(f1_hat*f2_hat,1,trapz,x=tt)
    logelr0[k]=-2*emplik(rho,rho.target)$logelr
     if(logelr0[k] < qchisq(0.95,df=1))
      count[l] <- count[l] + 1
     len[k] <- CIlen(rho,bs_value=qchisq(0.95,1))$len
  }
  CI.len[l] <- mean(len)
}
Sys.time()

#in addition, we've kept logelr0 of K = 100 simulation runs when h = 0.21
#for qq plot
