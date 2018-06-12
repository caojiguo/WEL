#==============================================
#first scenario: irregular design


#function used to estimate rho under irregular design
dyncor2 <- function(Y1,Y2,time,bdw)
{ #Y1= X1+err1; Y2= X2+err2; time=tt; bdw <- h[l]
  n <- dim(Y1)[1]
  m <- length(time)
  Y1_pred <- matrix(NA, nrow=n, ncol=m)
  Y2_pred <- matrix(NA, nrow=n, ncol=m)
  for(i in 1:n)
  { #i = 2
    # cat("i = ", i, "\n")
    index <- which(!is.na(Y1[i,]))
    fit1 <- locpoly(x=time[index], y=Y1[i,index],bandwidth=bdw, gridsize=m, range.x=c(0,1))
    Y1_pred[i,] <- fit1$y
    fit2 <- locpoly(x=time[index], y=Y2[i,index],bandwidth=bdw, gridsize=m, range.x=c(0,1))
    Y2_pred[i,] <- fit2$y
  }
  
  #   plot(x=time[index], y=Y1[i,index], type="l", col="red")
  #   lines(x=tt, y=Y1_pred[i,], col="blue")
  
  for(j in 1:m)
  {
    if(anyNA(Y1_pred[,j]))
    {
      Y1_pred[,j][is.na(Y1_pred[,j])] <- mean(Y1_pred[,j], na.rm=T)
      Y1_pred[,j] <- Y1_pred[,j]-mean(Y1_pred[,j])
      Y2_pred[,j][is.na(Y2_pred[,j])] <- mean(Y2_pred[,j], na.rm=T)
      Y2_pred[,j] <- Y2_pred[,j]-mean(Y2_pred[,j])
    }
    
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

n <- 100
K <- 100
f_1 <- f_2 <- matrix(NA,nrow=n,ncol=m)
Y_1 <- Y_2 <- matrix(NA,nrow=n,ncol=m)
Y1_pred <- Y2_pred <- Y_1
f1_hat <- f2_hat <- matrix(NA,nrow=n,ncol=m)
rho.target <- 0.5
#rho.est <- matrix(0, nrow=K, ncol=L)

interv <- 25:100 #number of obs. in sparse trajectories


count <- numeric(L)
bootCI.len <- numeric(L)
len <- numeric(K)
B <- 500
rho.boot <- numeric(B)



set.seed(20170615)
#Sys.time()
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
    
    
    
    for (i in 1:n)
    { 
      num <- sample(interv, size=1)
      ni <- sample(1:m, num)
      Y_1[i, -ni] <- NA
      Y_2[i, -ni] <- NA
      fit1 <- locpoly(x=tt[ni], y=Y_1[i,ni],bandwidth=h[l],gridsize=m, range.x=c(0,1))
      Y1_pred[i,] <- fit1$y
      fit2 <- locpoly(x=tt[ni], y=Y_2[i,ni],bandwidth=h[l],gridsize=m, range.x=c(0,1))
      Y2_pred[i,] <- fit2$y
      
    }
    
    
    
    for(j in 1:m)
    {
      Y1_pred[,j][is.na(Y1_pred[,j])] <- mean(Y1_pred[,j], na.rm=T)
      #f1_hat[,j] <- Y1_pred[,j]-mean(Y1_pred[,j])
      Y2_pred[,j][is.na(Y2_pred[,j])] <- mean(Y2_pred[,j], na.rm=T)
      #f2_hat[,j] <- Y2_pred[,j]-mean(Y2_pred[,j])
    }
    
    resid1 <- Y_1 - Y1_pred
    resid2 <- Y_2 - Y2_pred
    
    for(b in 1:B)
    {
      #b = 1
      #cat("b = ", b, "\n")
      index1 <- sample(1:n, size=n, replace=T)
      X1 <- Y1_pred[index1,]
      X2 <- Y2_pred[index1,]
      index2 <- sample(1:n, size=n, replace=T)
      err1 <- resid1[index2,]
      err2 <- resid2[index2,]
      rho.boot[b] <- dyncor2(Y1= X1+err1, Y2= X2+err2, time=tt, bdw <- h[l])
    }
    
    Sys.time()
    if(rho.target>quantile(rho.boot,probs=0.025) & rho.target<quantile(rho.boot,probs=0.975))
      count[l] <- count[l] + 1
    len[k] <- quantile(rho.boot,probs=0.975) - quantile(rho.boot,probs=0.025)
    
  }
  bootCI.len[l] <- mean(len)
}

#repeat the above procedure for n=50, 
