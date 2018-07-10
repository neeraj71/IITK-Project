library(MASS)
set.seed(20)

x1=c(.1,.12,.11,.13,.14,.7,.7,.65,.66,.11,.111,.112,.114,.115)
u<- seq(0.01,1,length.out = 200)
v<- seq(10,11,length.out=200)
v2=c(5,8,9,13,12,14,13,7,9,10,10.1,10.1,10.05,10.03)
X=c(x1,u)
X=matrix(X)
noise <- rnorm(length(u),0,0.3)
M=basicMethod=matrix(0, nrow = length(X), ncol = 1)
y1=c(-80,-95.1,-924.-879,-54.8,-65.2,-91,-83.9,-98.9,-100.88,20,20.05,20.06,20.09,20.02)
y2=1000*(5000*(u)^3-3*u^4+3*u^5-u^6+v^6+1000*sin(v)+100*sqrt(v))+as.matrix(noise)
Y=c(y1,y2)
X2=c(v2,v)
X2=matrix(X2)
#X=matrix(u)
#X2=matrix(v)
#Y=matrix(y2)

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
sigma=8000000000000
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
h=seq(0.03,0.03,1)
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
  #pairs(G)
  for(i in c(1:length(X))){
    for(j in c(1:3)){
      if(j==1){
        basic[i,j]=X[i]}
      if(j==2)
        basic[i,j]=X2[i]
      else if(j==3){
        basic[i,j]=basicMethod[i]}
    }}
 # pairs(basic)
  
  for(i in c(1:length(X))){
    for(j in c(1:3)){
      if(j==1){
        Final[i,j]=X[i]}
      if(j==2)
        Final[i,j]=X2[i]
      else if(j==3){
        Final[i,j]=M[i]}
    }}
library(rgl)
  library(plot3D)
  mfrow3d(nr=1,nc=2,sharedMouse = TRUE)
  plot3d(Final)
  plot3d(basic)
  rglwidget
  #surf3D (X, X2, M, colvar = M, phi = 40, theta = 40, col = NULL, NAcol = "white", breaks = NULL, border = NA, facets = TRUE, colkey = NULL, panel.first = NULL, clim = NULL, clab = NULL, bty = "n", lighting = FALSE, shade = NA, ltheta = -135, lphi = 0, inttype = 1, add = FALSE, plot = TRUE)
  
  #pairs(Final)
  
}



