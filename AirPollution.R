pm25 <- read.csv("pm25_2000.csv",header=T)
no2 <- read.csv("no2_2000.csv",header=T)

head(pm25)

dim(pm25)



library(locpol)
library(caTools)

source("CIwt.R")
tau=10^(-4);


tt <- 1:366/366
n <- nrow(pm25)
m <- ncol(pm25)
num <- numeric(n)
pm25_pred <- pm25
no2_pred <- no2

St.pm25 <- pm25
St.no2 <- no2

deg <- 1
kernel <- EpaK
est <- function(bw, dat, x) return(locPolSmootherC(dat$x,dat$y, x, bw, deg,
                                                   kernel)$beta0)

for(i in 1:n)
{
   # i = 1
  index1 <- which(!is.na(pm25[i,]))
  cvBwSel <- regCVBwSelC(tt[index1], pm25[i,index1], deg, kernel, interval = c(0.01, 0.25))
  pm25_pred[i,] <-  est(cvBwSel, data.frame(x=tt[index1], y=as.numeric(pm25[i,index1])), tt)
 
  index2 <- which(!is.na(no2[i,]))
  cvBwSel <- regCVBwSelC(tt[index2], no2[i,index2], deg, kernel, interval = c(0.01, 0.25))
  no2_pred[i,] <- est(cvBwSel, data.frame(x=tt[index2], y=as.numeric(no2[i,index2])), tt)

  num[i] <- sqrt(length(index1)*length(index2))
}

summary(num)


for(j in 1:m)
{
  pm25_pred[,j][is.na(pm25_pred[,j])] <- mean(pm25_pred[,j], na.rm=T)
  St.pm25[,j] <- pm25_pred[,j]-mean(pm25_pred[,j])
  no2_pred[,j][is.na(no2_pred[,j])] <- mean(no2_pred[,j], na.rm=T)
  St.no2[,j] <- no2_pred[,j]-mean(no2_pred[,j])
}

anyNA(St.pm25); anyNA(St.no2)



for(j in 1:m)
{
  St.pm25[,j] <- St.pm25[,j] - apply(St.pm25,1,trapz,x=tt)
  St.no2[,j] <- St.no2[,j] - apply(St.no2,1,trapz,x=tt)
  
}

M1 <- St.pm25^2
M2 <- St.no2^2


for(j in 1:m)
{
  St.pm25[,j] <- St.pm25[,j]/sqrt(apply(M1,1,trapz,x=tt))
  St.no2[,j] <- St.no2[,j]/sqrt(apply(M2,1,trapz,x=tt))
 
}

# set.seed(2017720)
# randum <- sample(1:n, size=1)
# #windows()
# pdf("standardized plot of two airpollutants.pdf")
# par(lwd=2, font.lab=2, cex.lab=1.2, mar=c(5,4.6,4,2)+0.1)
# plot(x=tt*366, y=St.pm25[randum,], type="l", lty=1, xlab="Day", ylab="Value")
# lines(x=tt*366, y=St.no2[randum,], lty=2)
# legend("topleft", legend=c("PM2.5", "NO2"), lty=1:2,lwd=2)
# dev.off()


rho <- apply(St.pm25*St.no2,1,trapz,x=tt)
mean(rho) #point est
#0.411, 
weig <- sum(num)/(n*n*num)

CI <- CIlen.wt(x=rho,ds=rep(1/n,n),qs=weig,bs_value=qchisq(0.95,1))$ci
#0.3455525 0.4681721




#bootstrap analysis
est <- function(bw, dat, x) return(locPolSmootherC(dat$x,dat$y, x, bw, deg,
                                                   kernel)$beta0)
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
    # fit1 <- locpoly(x=time[index], y=Y1[i,index],bandwidth=bdw, gridsize=m, range.x=c(0,1))
    # Y1_pred[i,] <- fit1$y
    cvBwSel <- regCVBwSelC(time[index], Y1[i,index], 1, EpaK, interval = c(0.01, 0.25))
    #est(cvBwSel, data.frame(x=tt, y=C5[i,]), tt)
    Y1_pred[i,] <- est(cvBwSel, data.frame(x=time[index], y=as.numeric(Y1[i,index])), time)
    # fit2 <- locpoly(x=time[index], y=Y2[i,index],bandwidth=bdw, gridsize=m, range.x=c(0,1))
    # Y2_pred[i,] <- fit2$y
    
    index2 <- which(!is.na(Y2[i,]))
    cvBwSel <- regCVBwSelC(time[index2], Y2[i,index2], 1, EpaK, interval = c(0.01, 0.25))
    Y2_pred[i,] <- est(cvBwSel, data.frame(x=time[index2], y=as.numeric(Y2[i,index2])), time)
  }
  
  #   plot(x=time[index], y=Y1[i,index], type="l", col="red")
  #   lines(x=tt, y=Y1_pred[i,], col="blue")
  
  # for(j in 1:m)
  # {
  #   if(anyNA(Y1_pred[,j]))
  #   {
  #     Y1_pred[,j][is.na(Y1_pred[,j])] <- mean(Y1_pred[,j], na.rm=T)
  #     Y1_pred[,j] <- Y1_pred[,j]-mean(Y1_pred[,j])
  #     Y2_pred[,j][is.na(Y2_pred[,j])] <- mean(Y2_pred[,j], na.rm=T)
  #     Y2_pred[,j] <- Y2_pred[,j]-mean(Y2_pred[,j])
  #   }
  #   
  # }
  
  for(j in 1:m)
  {
    Y1_pred[,j][is.na(Y1_pred[,j])] <- mean(Y1_pred[,j], na.rm=T)
    Y1_pred[,j] <- Y1_pred[,j]-mean(Y1_pred[,j])
    Y2_pred[,j][is.na(Y2_pred[,j])] <- mean(Y2_pred[,j], na.rm=T)
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


resid1 <- pm25 - pm25_pred
resid2 <- no2 - no2_pred


B <- 500
rho.boot  <- numeric(B)
set.seed(123)

for(b in 1:B)
{
  #b = 1
  cat("b = ", b, "\n")
  index1 <- sample(1:n, size=n, replace=T)
  X1 <- pm25_pred[index1,]
  X2 <- no2_pred[index1,]
  index2 <- sample(1:n, size=n, replace=T)
  err1 <- resid1[index2,]
  err2 <- resid2[index2,]
  rho.boot[b] <- dyncor.ir(Y1= X1+err1, Y2= X2+err2, time=tt)
}

quantile(rho.boot, probs=c(0.025,0.975))
#0.2598284 0.4748692 