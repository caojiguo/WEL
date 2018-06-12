source("CIwt.R")
tau=10^(-4);


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


library(mvtnorm)
library(KernSmooth)
library(caTools)
n <- 50  #repeat for n = 100
num <- rep(0, n)
f_1 <- f_2 <- matrix(NA,nrow=n,ncol=m)
Y_1 <- Y_2 <- matrix(NA,nrow=n,ncol=m)
Y1_pred <- Y2_pred <- Y_1
f1_hat <- f2_hat <- matrix(NA,nrow=n,ncol=m)
rho.target <- 0.5
K <- 100

#scenario 1: irregular design

interv <- 25:100 #number of obs. in sparse trajectories

#logelr0 <- 0
#logelr0 <- numeric(K)
count <- numeric(L)
CI.len <- numeric(L)
len <- numeric(K)
set.seed(20170615)



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
    
    # index <- sample(1:n, n*pi)
    for (i in 1:n)
    { #i = 49
      num[i] <- sample(interv, size=1)
      ni <- sample(1:m, num[i])
      fit1 <- locpoly(x=tt[ni], y=Y_1[i,ni],bandwidth=h[l],gridsize=m, range.x=c(0,1))
      Y1_pred[i,] <- fit1$y
      fit2 <- locpoly(x=tt[ni], y=Y_2[i,ni],bandwidth=h[l],gridsize=m, range.x=c(0,1))
      Y2_pred[i,] <- fit2$y
      
    }
    
    
    
    
    for(j in 1:m)
    {
      Y1_pred[,j][is.na(Y1_pred[,j])] <- mean(Y1_pred[,j], na.rm=T)
      f1_hat[,j] <- Y1_pred[,j]-mean(Y1_pred[,j])
      Y2_pred[,j][is.na(Y2_pred[,j])] <- mean(Y2_pred[,j], na.rm=T)
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
    weig <- sum(num)/(n*n*num)
    logelr0 = -2*emplik.wt(u=cbind(1,rho),ds=rep(1/n,n), qs=weig, tx=c(1,rho.target))
    if(logelr0 < qchisq(0.95,df=1))
      count[l] <- count[l] + 1
    len[k] <- CIlen.wt(x=rho,ds=rep(1/n,n),qs=weig,bs_value=qchisq(0.95,1))$leng
  }
  CI.len[l] <- mean(len)
}

#in addition, we've kept logelr0 of K = 100 simulation runs when h = 0.21
# for qq plot