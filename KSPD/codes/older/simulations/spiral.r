library(MASS)
library(rgl)
set.seed(20)

x1=c(0.2,0.22,0.21,0.23,0.24,0.2,0.2,0.65,0.66,0.1,0.111,0.112,0.114,0.115,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.22,0.21,0.23,0.24,0.2,0.2,0.65,0.66,0.1,0.111,0.112,0.114,0.115,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2)

t=seq(9,20,length.out = 400)
u<-100* sin(t)
v<- 100*cos(t)
v2=c(0,0,0,0,0,0,0,0,0,0,0.1,0.1,0.05,0.03,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1,0.1,0.05,0.03,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0,0,0,0,0)
X=c(x1,u)
X=matrix(X)
noise <- rnorm(length(u),0,0.3)
M=basicMethod=matrix(0, nrow = length(X), ncol = 1)
y1=c(300,300,300,300,300,300,300,300,300,300,300.1,300.1,300.05,300.03,300,300.1,300.2,300.3,300.4,300.5,300.6,300.7,300,300,300,300,300,300,300,300,300,300,300,300,300,300,300,300.1,300.1,300.05,300.03,300,300.1,300.2,300.3,300.4,300.5,300.6,300.7,300,300,300,300,300)
Fun=30*t
y2=Fun+as.matrix(noise)

Y=c(y1,y2)
X2=c(v2,v)
X2=matrix(X2)

Final=basic=G=matrix(0, nrow = length(X), ncol = 3)
for(i in c(1:length(X))){
  for(j in c(1:3)){
    if(j==1){
      G[i,j]=X[i]}
    if(j==2)
      G[i,j]=X2[i]
    else if(j==3){
      G[i,j]=Y[i]}
  }}
sigma=80000
K=matrix(0, nrow = length(X), ncol = length(X))
W=matrix(0, nrow = length(X), ncol = 1)

for(j in c(1:length(X))){
  K[,j]=exp(-rowSums((sweep(G,2,G[j,],FUN="-"))^2)/sigma^2)
}

T=matrix(0, nrow = length(X), ncol = 1)
D=matrix(0, nrow = length(X), ncol = 1)
Z=matrix(0, nrow = length(X), ncol = 1)
for(n in c(1:length(X))){
  for(i in c(1:length(X))){
    T[i]=exp(-sum((G[i,]-G[n,])^2)/sigma^2)
    
    D[i]= sqrt(1+1-2*T[i])
    if(D[i]==0)
    {
      Z[i]=0
    }else{
      Z[i]=1/D[i]
    }
  }
  P=matrix(0, nrow = length(X), ncol = length(X))
  P=1+K-as.vector(T)
  
  for(j in c(1:length(X))){
    P[,j]=P[,j]-T[j]
  }
  
  W[n]=1-1/(length(X))*sqrt(t(Z)%*%P%*%Z)}
h=seq(.08,.08,1)
h=matrix(h)
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
  
  for(i in c(1:length(X))){
    for(j in c(1:3)){
      if(j==1){
        basic[i,j]=X[i]}
      if(j==2)
        basic[i,j]=X2[i]
      else if(j==3){
        basic[i,j]=basicMethod[i]}
    }}
  
  
  for(i in c(1:length(X))){
    for(j in c(1:3)){
      if(j==1){
        Final[i,j]=X[i]}
      if(j==2)
        Final[i,j]=X2[i]
      else if(j==3){
        Final[i,j]=M[i]}
    }}
  Robust_dist=(sum((Fun-M[c((length(x1)+1):length(X))])^2))/length(Fun)
  Nadaraya_Watson_dist=(sum((Fun-basicMethod[c((length(x1)+1):length(X))])^2))/length(Fun)
  cat("Robust_distance ", Robust_dist,"\n")
  cat("Nadarya_Watson_distance",Nadaraya_Watson_dist,"\n")
  
  mfrow3d(nr=1,nc=2,sharedMouse = TRUE)
  plot3d(Final)
  plot3d(basic)
  rglwidget()
  
}



