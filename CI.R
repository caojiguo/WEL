CIlen=function(x,bs_value)
{
#x <- Z; bs_value=sprho
barx=mean(x)
b=barx
mu=b+1;
logelr=-2*emplik(x,mu)$logelr

while(logelr<=bs_value)
  
 {mu=mu+0.5;

 logelr=-2*emplik(x,mu)$logelr}

c=mu

while(c-b>tau)
  
  {int=(c+b)/2

   logelr=-2*emplik(x,int)$logelr

   if (logelr<=bs_value) b=int

   else c=int
  
  }


d1=(c+b)/2;

b=barx
mu=b-1;
logelr=-2*emplik(x,mu)$logelr

while(logelr<=bs_value)
  
 {mu=mu-0.5;

  logelr=-2*emplik(x,mu)$logelr}

a=mu

while(b-a>tau)
  
  {int=(a+b)/2

logelr=-2*emplik(x,int)$logelr;

if(logelr<=bs_value) b=int
else a=int }

d2=(a+b)/2

return(list(len=d1-d2, CI=c(d2,d1)))
}


