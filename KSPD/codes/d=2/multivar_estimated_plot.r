library(rgl)
set.seed(20)

x1=c(.2,.22,.21,.23,.24,.2,.2,.65,.66,.11,.111,.112,.114,.115,.2,.2,.2,.2,.2,.2,.2,.2,.2,.2,.2,.2,.2)
u<- seq(0.01,1,length.out = 200)
X=c(x1,u)
X=matrix(X)

v<- seq(10,11,length.out=200)
v2=c(5,8,9,13,12,14,13,7,9,10,10.1,10.1,10.05,10.03,4,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4,4,4,4,4)
X2=c(v2,v)
X2=matrix(X2)

noise <- rnorm(length(u),0,10)

y1=c(-800000000000,-95000000.1,-92400000,-87900000,-5400000.8,-650000.2,-91000000,-830000.9,-100.88,-2000000000000,-2000000000000.05,-2000000000000.06,-2000000000000.09,-2000000000000.02,-10000000000,-10000000000001,-100000002,-1000000003,-1000000004,-1000000005,-800000006,-1000000000007,-1000000000007,-1000000000007,-1000000000007,-1000000000007,-1000000000007)
y2=-100*(8*exp(v)*10*exp(10*u))+as.matrix(noise)
Fun=-100*(8*exp(X2)*10*exp(10*X))
Y=c(y1,y2)

Nadarya_Watson_method=Robust_Method=Final=basic=G=matrix(0, nrow = length(X), ncol = 3)
K=P=matrix(0, nrow = length(X), ncol = length(X))
M=T=D=Z=W=basicMethod=matrix(0, nrow = length(X), ncol = 1)

G=cbind(X,X2,Y)

sigma=c(89558784586,97896,500,6000,789,12345,658967,58586,11)
h=c(seq(0.01,0.1,length.out=10),seq(0.2,5,length.out=10))

Robust_dist=Nadaraya_Watson_dist=matrix(0, nrow = length(sigma), ncol = length(h))
h=matrix(h)
for(k in c(1:length(sigma))){
  for(j in c(1:length(X))){
    K[,j]=exp(-rowSums((sweep(G,2,G[j,],FUN="-"))^2)/(sigma[k])^2)
  }
  
  for(n in c(1:length(X))){
    for(i in c(1:length(X))){
      T[i]=exp(-sum((G[i,]-G[n,])^2)/(sigma[k])^2)
      D[i]= sqrt(1+1-2*T[i])
      if(D[i]==0)
      {Z[i]=0
      }else{
        Z[i]=1/D[i]}
    }
    P=1+K-as.vector(T)
    P=sweep(P,2,as.vector(T),FUN="-")
    W[n]=1-(1/(length(X)))*sqrt(t(Z)%*%P%*%Z)
  }
  
  for(j in c(1:length(h))){
    for(i in c(1:length(X))){
      temp=(sqrt((X[,1]-X[i,1])^2+(X2[,1]-X2[i,1])^2))/h[j,1]
      Ker=1/(sqrt(2*pi))*exp(-(temp^2)/2)
      num1=t(as.matrix(Ker))%*%Y
      denom1=sum(Ker)
      
      basicMethod[i,1]=num1/denom1
      
      num=t(as.matrix(Ker*W))%*%Y
      denom=(sum(Ker*W))
      M[i,1]=num/denom
    }
    
    Robust_dist[k,j]=(sum((Fun-M)^2))/length(Fun)
    Nadaraya_Watson_dist[k,j]=(sum((Fun-basicMethod)^2))/length(Fun)
  }
}
ind=which(Robust_dist == min(Robust_dist), arr.ind = TRUE)
index=which(Nadaraya_Watson_dist == min(Nadaraya_Watson_dist), arr.ind = TRUE)

sig=c(sigma[ind[1,1]],sigma[index[1,1]])
H=matrix(0, nrow = 1, ncol = 1)
mfrow3d(nr=1,nc=2,sharedMouse = TRUE)

for(k in c(1:length(sig))){
  for(j in c(1:length(X))){
    K[,j]=exp(-rowSums((sweep(G,2,G[j,],FUN="-"))^2)/(sig[k])^2)
  }
  
  for(n in c(1:length(X))){
    for(i in c(1:length(X))){
      T[i]=exp(-sum((G[i,]-G[n,])^2)/(sig[k])^2)
      D[i]= sqrt(1+1-2*T[i])
      if(D[i]==0)
      {Z[i]=0
      }else{
        Z[i]=1/D[i]}
    }
    P=1+K-as.vector(T)
    P=sweep(P,2,as.vector(T),FUN="-")
    W[n]=1-(1/(length(X)))*sqrt(t(Z)%*%P%*%Z)
  }
  if(k==1)
    H[1]=h[ind[1,2]]
  else
    H[1]=h[index[1,2]]
  for(j in c(1:length(H))){
    for(i in c(1:length(X))){
      temp=(sqrt((X[,1]-X[i,1])^2+(X2[,1]-X2[i,1])^2))/H[j,1]
      Ker=1/(sqrt(2*pi))*exp(-(temp^2)/2)
      num1=t(as.matrix(Ker))%*%Y
      denom1=sum(Ker)
      
      basicMethod[i,1]=num1/denom1
      
      num=t(as.matrix(Ker*W))%*%Y
      denom=(sum(Ker*W))
      M[i,1]=num/denom
    }
    if (k==1)
     Robust_Method= Final=cbind(X,X2,M)
    else
  Nadarya_Watson_method=  basic=cbind(X,X2,basicMethod)
  
    
    plot3d(Robust_Method)
    plot3d(Nadarya_Watson_method)
    rglwidget
    
  }
}
  

cat("Robust-distance-min = ",min(Robust_dist),"\n")


cat("h_min =",h[ind[1,2]],"\n")
cat("sigma_min =",sigma[ind[1,1]],"\n")

cat("Nadarya-Watson-distance-min = ",min(Nadaraya_Watson_dist),"\n")

cat("h_min =",h[index[1,2]],"\n")
cat("sigma_min =",sigma[index[1,1]],"\n")




