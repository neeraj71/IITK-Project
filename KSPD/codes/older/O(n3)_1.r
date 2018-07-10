library(MASS)
set.seed(20)

x1=c(.1,.12,.11,.13,.14,.7,.7,.71,.65,.67,.675,.66,.4,.41,.42,.43,.44,.45,.46)
u<- seq(0.01,1,length.out = 200)
X=c(x1,u)
X=matrix(X)
noise <- rnorm(length(u),0,0.3)
Y=M=basicMethod=matrix(0, nrow = length(X), ncol = 1)
y1=c(5,5.1,4.9,4.8,5.2,-1,-.9,-.92,-.93,-.93,-.9,-.88,-8001,-8001,-8000,-8003,-8004,-8005,-8006)
y2=400*(100*u^3+50*sin(u)-200)+as.matrix(noise)
Y=c(y1,y2)

G=matrix(0, nrow = length(X), ncol = 2)
for(i in c(1:length(X))){
  for(j in c(1:2)){
    if(j==1){
      G[i,j]=X[i]}else{
        G[i,j]=Y[i]}
  }}
sigma=8000000
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
h=seq(0.04,0.04,1)
h=matrix(h)
for(j in c(1:length(h))){
  for(i in c(1:length(X))){
    temp=(X[,1]-X[i,1])/h[j,1]
    Ker=1/(sqrt(2*pi))*exp(-(temp^2)/2)
    num1=t(as.matrix(Ker))%*%Y
    denom1=sum(Ker)
    
    basicMethod[i,1]=num1/denom1
    
    num=t(as.matrix(Ker*W))%*%Y
    denom=(sum(Ker*W))
    M[i,1]=num/denom
  }
  
  par(mfrow=c(2,1))
  plot(X,Y,ylim=range(c(Y,M)))
  par(new=TRUE)
  plot(X,M, ylim=range(c(Y,M)), axes = FALSE, xlab = "", ylab = "")
  
  plot(X,Y,ylim=range(c(Y,basicMethod)))
  par(new=TRUE)
  plot(X,basicMethod, ylim=range(c(Y,basicMethod)), axes = FALSE, xlab = "", ylab = "")
  
}



