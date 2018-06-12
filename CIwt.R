
##########################################
# R Function for WKL Lagrange Multiplier #
#                                        #
# Input: u=(x_1,x_2,...,x_n)             #
#        ds=(d_1,d_2,...,d_n)            #
#        (Design Wights: d_i=1/pi_i)     #
#        ds=(1,1,...1) for iid data      #
#        qs=(q_1,q_2,...,q_n)            #
#        (q-weights)                     #
#        tx: benchmark totals for x      #
# Output: likelihood ratio               #
#                                        #
# inspired by Changbao Wu  
##########################################
emplik.wt=function(u,ds,qs,tx)
{ #tx=1, sum q^{-1} = n, d_i = 1/n, x_1 = (1, x_1), list by row
  n=length(ds)
  M=0*tx
  dif=1
  tol=1e-8
  k=0
  while(dif>tol & k<=100){
    D1=t(u)%*%((ds/(1+qs*u%*%M))*rep(1,n))-tx
    DD=-t(u)%*%(c((qs*ds/(1+qs*(u%*%M))^2))*u)
    aa=abs(det(DD))
    if(aa<=0.00000001) k=101
    if(aa>0.00000001){
      D2=solve(DD,D1,tol=1e-40)
      dif=max(abs(D2))
      rule=1
      while(rule>0){
        rule=0
        if(min(1+qs*(u%*%(M-D2)))<=0) rule=rule+1
        if(rule>0) D2=D2/2
      }
    }
    M=M-D2
    k=k+1
  }
  if(k>=100) M=0*tx+10000
  #return(as.vector(M))
  
  ws=ds/(1+qs*(u%*%M))
  WKL1=sum((ds/qs)*(log(ws/ds)-ws/ds+1))
  ws2=rep(1/n,n)
  WKL2=sum((ds/qs)*(log(ws2/ds)-ws2/ds+1))
  return(WKL1-WKL2)
}







##########################################
# R Function for WKL Lagrange Multiplier #
#                                        #
# Input: x=(x_1,x_2,...,x_n)             #
#        ds=(d_1,d_2,...,d_n)            #
#        (Design Wights: d_i=1/pi_i)     #
#        ds=(1,1,...1) for iid data      #
#        qs=(q_1,q_2,...,q_n)            #
#        (q-weights)                     #

# Output: length of CI                   #
#                                        #
# (Modify code provided by Changbao wu)  #
##########################################


#specify tau first

CIlen.wt=function(x,ds,qs,bs_value)
{
  #tx=1, sum q^{-1} = n, d_i = 1/n, x_1 = (1, x_1), list by row
  barx=mean(x)
  b=barx
  mu=b+1;
  
  #u=cbind(1,dat),ds=rep(1/n,n),qs=rep(1,n),tx=c(1,mu)
  logelr=-2*emplik.wt(u=cbind(1,x),ds,qs,tx=c(1,mu))
  while(logelr<=bs_value)
    
 {mu=mu+0.5;
  
  logelr=-2*emplik.wt(u=cbind(1,x),ds,qs,tx=c(1,mu))}
  
  c=mu
  
  while(c-b>tau)
    
  {int=(c+b)/2
  
  logelr=-2*emplik.wt(u=cbind(1,x),ds,qs,tx=c(1,int))
  
  if (logelr<=bs_value) b=int
  
  else c=int
  
  }
  
  
  d1=(c+b)/2;
  
  b=barx
  mu=b-1;
  logelr=-2*emplik.wt(u=cbind(1,x),ds,qs,tx=c(1,mu))
  
  while(logelr<=bs_value)
    
  {mu=mu-0.5;
  
  logelr=-2*emplik.wt(u=cbind(1,x),ds,qs,tx=c(1,mu))}
  
  a=mu
  
  while(b-a>tau)
    
  {int=(a+b)/2
  
  logelr=-2*emplik.wt(u=cbind(1,x),ds,qs,tx=c(1,int));
  
  if(logelr<=bs_value) b=int
  else a=int }
  
  d2=(a+b)/2
  
  list(leng=d1-d2, ci=c(d2, d1))
}