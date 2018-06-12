
#smoothing and computing rho
dyncor.boot <- function(Y1,Y2,time)
{ #Y1= X1+err1; Y2= X2+err2; time=tt; bdw <- h[l]
  n <- dim(Y1)[1]
  m <- length(time)
  Y1_pred <- Y1
  Y2_pred <- Y2
  for(i in 1:n)
  { #i = 1
    #fit1 <- locpoly(x=time, y=Y1[i,],bandwidth=bdw, gridsize=m)
    cvBwSel <- regCVBwSelC(time, Y1[i,], 1, EpaK, interval = c(0.01, 0.25))
    #est(cvBwSel, data.frame(x=tt, y=C5[i,]), tt)
    Y1_pred[i,] <- est(cvBwSel, data.frame(x=time, y=Y1[i,]), time)
    cvBwSel <- regCVBwSelC(time, Y2[i,], 1, EpaK, interval = c(0.01, 0.25))
    Y2_pred[i,] <- est(cvBwSel, data.frame(x=time, y=Y2[i,]), time)
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
load("C5C3P8.RData", verbose=T)
windows()
plot(x=0:255, y=apply(C5, 2, mean), type="l", col="red", xlab="", 
      ylab="Value", ylim=range(apply(P8, 2, mean)))
lines(x=0:255, y=apply(C3, 2, mean), type="l", col="blue")
lines(x=0:255, y=apply(P8, 2, mean), type="l", col="green")


source("CI.R")
source("scel.R")
tau=10^(-4);

library(locpol)
library(caTools)



tt <- 0:255/255

#define local linear smoothing
deg <- 1
kernel <- EpaK
est <- function(bw, dat, x) return(locPolSmootherC(dat$x,dat$y, x, bw, deg,
                                                   kernel)$beta0)

m <- length(tt)
n <- nrow(C5)

C5_pred <- C5
C3_pred <- C3
P8_pred <- P8

for(i in 1:n)
{ #i = 1
  cvBwSel <- regCVBwSelC(tt, C5[i,], deg, kernel, interval = c(0.01, 0.25))
  C5_pred[i,] <- est(cvBwSel, data.frame(x=tt, y=C5[i,]), tt)
  cvBwSel <- regCVBwSelC(tt, C3[i,], deg, kernel, interval = c(0.01, 0.25))
  C3_pred[i,] <- est(cvBwSel, data.frame(x=tt, y=C3[i,]), tt)
  cvBwSel <- regCVBwSelC(tt, P8[i,], deg, kernel, interval = c(0.01, 0.25))
  P8_pred[i,] <- est(cvBwSel, data.frame(x=tt, y=P8[i,]), tt)
}



# windows()
# plot(x=tt, y=C5[1,],type="l",col="red")
# lines(x=tt, y=C5_pred[1,], col="blue", lty=2)
# lines(x=tt, y=FC5[1,],type="l",col="blue")
# lines(x=tt, y=FC5_pred[1,], col="blue", lty=2)
# lines(x=tt, y=P8[1,],type="l",col="green")
# lines(x=tt, y=P8_pred[1,], col="green", lty=2)

Sd.C5 <- C5
Sd.C3 <- C3
Sd.P8 <- P8



for(j in 1:m)
{
  Sd.C5[,j] <- C5_pred[,j]-mean(C5_pred[,j])
  Sd.C3[,j] <- C3_pred[,j]-mean(C3_pred[,j])
  Sd.P8[,j] <- P8_pred[,j]-mean(P8_pred[,j])
}

for(j in 1:m)
{
  Sd.C5[,j] <- Sd.C5[,j] - apply(Sd.C5,1,trapz,x=tt)
  Sd.C3[,j] <- Sd.C3[,j] - apply(Sd.C3,1,trapz,x=tt)
  Sd.P8[,j] <- Sd.P8[,j] - apply(Sd.P8,1,trapz,x=tt)
}

M1 <- Sd.C5^2
M2 <- Sd.C3^2
M3 <- Sd.P8^2

for(j in 1:m)
{
  Sd.C5[,j] <- Sd.C5[,j]/sqrt(apply(M1,1,trapz,x=tt))
  Sd.C3[,j] <- Sd.C3[,j]/sqrt(apply(M2,1,trapz,x=tt))
  Sd.P8[,j] <- Sd.P8[,j]/sqrt(apply(M3,1,trapz,x=tt))
}



rhoC5C3 <- apply(Sd.C5*Sd.C3,1,trapz,x=tt)
rhoC5P8 <- apply(Sd.C5*Sd.P8,1,trapz,x=tt)
rhoC3P8 <- apply(Sd.C3*Sd.P8,1,trapz,x=tt)

#point estimate
mean(rhoC5C3); mean(rhoC5P8); mean(rhoC3P8)
#0.826, 0.439, 0.369

#confidence interval based on empirical likelihood
CIC5C3 <- CIlen(rhoC5C3,bs_value=qchisq(0.95,1))$CI
#0.772 0.863
CIC5P8 <- CIlen(rhoC5P8,bs_value=qchisq(0.95,1))$CI
#0.357 0.511
CIC3P8 <- CIlen(rhoC3P8,bs_value=qchisq(0.95,1))$CI
#0.289 0.444


#bootstrap analysis

for(i in 1:n)
{ #i = 1; 
  cvBwSel <- regCVBwSelC(tt, C5[i,], deg, kernel, interval = c(0.01, 0.25))
  C5_pred[i,] <- est(cvBwSel, data.frame(x=tt, y=C5[i,]), tt)
  cvBwSel <- regCVBwSelC(tt, C3[i,], deg, kernel, interval = c(0.01, 0.25))
  C3_pred[i,] <- est(cvBwSel, data.frame(x=tt, y=C3[i,]), tt)
  cvBwSel <- regCVBwSelC(tt, P8[i,], deg, kernel, interval = c(0.01, 0.25))
  P8_pred[i,] <- est(cvBwSel, data.frame(x=tt, y=P8[i,]), tt)
}
residC5 <- C5 - C5_pred
residC3 <- C3 - C3_pred
residP8 <- P8 - P8_pred
  

B <- 500
rhoC5C3.boot <- rhoC5P8.boot <- rhoC3P8.boot <- numeric(B)
set.seed(123)
for(b in 1:B)
{#b = 1
  cat("b = ", b, "\n")
  index1 <- sample(1:n, size=n, replace=T)
  X1 <- C5_pred[index1,]
  X2 <- C3_pred[index1,]
  X3 <- P8_pred[index1,]
  index2 <- sample(1:n, size=n, replace=T)
  err1 <- residC5[index2,]
  err2 <- residC3[index2,]
  err3 <- residP8[index2,]
  rhoC5C3.boot[b] <- dyncor.boot(Y1= X1+err1, Y2= X2+err2, time=tt)
  rhoC5P8.boot[b] <- dyncor.boot(Y1= X1+err1, Y2= X3+err3, time=tt)
  rhoC3P8.boot[b] <- dyncor.boot(Y1= X2+err2, Y2= X3+err3, time=tt)
}

quantile(rhoC5C3.boot, probs=c(0.025,0.975))
#0.7702790 0.8584177 
quantile(rhoC5P8.boot, probs=c(0.025,0.975))
#0.3542459 0.5123570 
quantile(rhoC3P8.boot, probs=c(0.025,0.975))
#0.2804984 0.4413344 



##======================
#data analysis for air pollution data


pm25 <- read.csv("pm25_2000.csv",header=T)
no2 <- read.csv("no2_2000.csv",header=T)

source("CIwt.R")



tt.air <- 1:366/366
n.air <- nrow(pm25)
m.air <- ncol(pm25)
num.air <- numeric(n.air)
pm25_pred <- pm25
no2_pred <- no2

St.pm25 <- pm25
St.no2 <- no2



for(i in 1:n.air)
{
  # i = 1
  index1 <- which(!is.na(pm25[i,]))
  cvBwSel <- regCVBwSelC(tt.air[index1], pm25[i,index1], deg, kernel, interval = c(0.01, 0.25))
  pm25_pred[i,] <-  est(cvBwSel, data.frame(x=tt.air[index1], y=as.numeric(pm25[i,index1])), tt.air)
  
  index2 <- which(!is.na(no2[i,]))
  cvBwSel <- regCVBwSelC(tt.air[index2], no2[i,index2], deg, kernel, interval = c(0.01, 0.25))
  no2_pred[i,] <- est(cvBwSel, data.frame(x=tt.air[index2], y=as.numeric(no2[i,index2])), tt.air)
  
  num.air[i] <- sqrt(length(index1)*length(index2))
}

summary(num.air)

j = 2
any(is.infinite(pm25_pred[,j]))

vec <- apply(pm25_pred, 2, function(x) any(is.infinite(x)))
any(vec)

for(j in 1:m.air)
{
  pm25_pred[,j][is.na(pm25_pred[,j])] <- mean(pm25_pred[,j], na.rm=T)
  St.pm25[,j] <- pm25_pred[,j]-mean(pm25_pred[,j])
  no2_pred[,j][is.na(no2_pred[,j])] <- mean(no2_pred[,j], na.rm=T)
  St.no2[,j] <- no2_pred[,j]-mean(no2_pred[,j])
}

anyNA(St.pm25); anyNA(St.no2)



for(j in 1:m.air)
{
  St.pm25[,j] <- St.pm25[,j] - apply(St.pm25,1,trapz,x=tt.air)
  St.no2[,j] <- St.no2[,j] - apply(St.no2,1,trapz,x=tt.air)
  
}

M1 <- St.pm25^2
M2 <- St.no2^2


for(j in 1:m.air)
{
  St.pm25[,j] <- St.pm25[,j]/sqrt(apply(M1,1,trapz,x=tt.air))
  St.no2[,j] <- St.no2[,j]/sqrt(apply(M2,1,trapz,x=tt.air))
  
}

rho <- apply(St.pm25*St.no2,1,trapz,x=tt.air)
mean(rho) #point est
#0.411, 
weig <- sum(num.air)/(n*n*num.air)

CI <- CIlen.wt(x=rho,ds=rep(1/n.air,n.air),qs=weig,bs_value=qchisq(0.95,1))$ci
# 0.3758870 0.4431477





set.seed(2017820)
randum1 <- sample(1:n, size=1)
randum2 <- sample(1:n.air, size=1)


pdf("standardized plot.pdf", height=7, width=9)
par(mfrow=c(1,2), font.lab=2, cex.lab=1.2, lwd=1,mai=c(1.02,1.01,0.82,0.42))
plot(x=tt.air*366, y=St.pm25[randum2,], type="l", lty=1, xlab="Day", ylab="Value")
lines(x=tt.air*366, y=St.no2[randum2,], lty=2)
legend("bottomleft", legend=c("PM2.5", "NO2"), lty=1:2,lwd=1)
plot(x=tt, y=Sd.C5[randum1,], type="l", lty=1, xlab="Time", ylab="Recordings")
lines(x=tt, y=Sd.C3[randum1,], lty=2)
lines(x=tt, y=Sd.P8[randum1,], lty=3)
legend("bottomleft", legend=c("C5", "C3", "P8"), lty=1:3,lwd=1)
dev.off()
