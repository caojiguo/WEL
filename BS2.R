#=========================================
#second scenario: regular design


#function used to estimate rho under regular design
dyncor <- function(Y1,Y2,time,bdw)
{ #Y1= X1+err1; Y2= X2+err2; time=tt; bdw <- h[l]
   n <- dim(Y1)[1]
   m <- length(time)
   Y1_pred <- Y1
   Y2_pred <- Y2
   for(i in 1:n)
  { #i = 1
    fit1 <- locpoly(x=time, y=Y1[i,],bandwidth=bdw, gridsize=m)
    Y1_pred[i,] <- fit1$y
    fit2 <- locpoly(x=time, y=Y2[i,],bandwidth=bdw, gridsize=m)
    Y2_pred[i,] <- fit2$y
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

#====================
library(mvtnorm)
#library(locpol)
library(KernSmooth)
library(caTools)
#library(emplik)

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

n <- 50
K <- 100
f_1 <- f_2 <- matrix(NA,nrow=n,ncol=m)
Y_1 <- Y_2 <- matrix(NA,nrow=n,ncol=m)
Y1_pred <- Y2_pred <- Y_1
f1_hat <- f2_hat <- matrix(NA,nrow=n,ncol=m)
rho.target <- 0.5
#rho.est <- matrix(0, nrow=K, ncol=L)

#logelr0 <- matrix(0, nrow=K, ncol=L)
count <- numeric(L)
bootCI.len <- numeric(L)
len <- numeric(K)
#resid1 <- resid2 <- matrix(NA,nrow=n,ncol=m)
B <- 500
rho.boot <- numeric(B)


#scenario 1: dense and regular design
set.seed(20170615)
Sys.time()
for(l in 1:L)
{  #l=9
  
  for(k in 1:K)
  { #k = 1
    eps <- rmvnorm(n, sigma = Sigma)
    e <- matrix(rnorm(n*m*2,mean=0,sd=1/2),nrow=n*2)
    
    for(j in 1:m)
    {
      Y_1[,j] <- 1 + eps[,1] + eta_1[j]*eps[,2] + eta_2[j]*eps[,3] +  e[1:n,j]
      Y_2[,j] <- 1 + eps[,4] + eta_1[j]*eps[,5] + eta_2[j]*eps[,6] + e[-(1:n),j]
    }
   
    for(i in 1:n)
    { #i = 1; 
      fit1 <- locpoly(x=tt, y=Y_1[i,],bandwidth=h[l], gridsize=m)
      Y1_pred[i,] <- fit1$y
      fit2 <- locpoly(x=tt, y=Y_2[i,],bandwidth=h[l], gridsize=m)
      Y2_pred[i,] <- fit2$y
    }
    resid1 <- Y_1 - Y1_pred
    resid2 <- Y_2 - Y2_pred
    
    

    for(b in 1:B)
    {#b = 1
      cat("b = ", b, "\n")
      index1 <- sample(1:n, size=n, replace=T)
      X1 <- Y1_pred[index1,]
      X2 <- Y2_pred[index1,]
      index2 <- sample(1:n, size=n, replace=T)
      err1 <- resid1[index2,]
      err2 <- resid2[index2,]
      rho.boot[b] <- dyncor(Y1= X1+err1, Y2= X2+err2, time=tt, bdw <- h[l])
    }
    
   if(rho.target>quantile(rho.boot,probs=0.025) & rho.target<quantile(rho.boot,probs=0.975))
      count[l] <- count[l] + 1
    len[k] <- quantile(rho.boot,probs=0.975) - quantile(rho.boot,probs=0.025)
    
  }
   bootCI.len[l] <- mean(len)
}

#repeat the above procedure for n=100 