set.seed(20)

x1=c(.1,.12,.11,.13,.14,.7,.7,.65,.66)
u<- seq(0,1,length.out = 200)
X=c(x1,u)
X=matrix(X)
noise <- rnorm(length(u),0,0.3)
Y=M=basicMethod=matrix(0, nrow = length(X), ncol = 1)
y1=c(5,5.1,4.9,4.8,5.2,-1,-.9,-.9,-.88)
y2=300*(u^3-3*u^4+3*u^5-u^6)+as.matrix(noise)
Y=c(y1,y2)

G=matrix(0, nrow = length(X), ncol = 2)
for(i in c(1:length(X))){
  for(j in c(1:2)){
    if(j==1){
      G[i,j]=X[i]}else{
        G[i,j]=Y[i]}
  }}
sigma=8
K=matrix(0, nrow = length(X), ncol = length(X))
W=matrix(0, nrow = length(X), ncol = 1)

for(j in c(1:length(X))){
  K[,j]=exp(-rowSums((sweep(G,2,G[j,],FUN="-"))^2)/sigma^2)
}

T=matrix(0, nrow = length(X), ncol = 1)
D=matrix(0, nrow = length(X), ncol = 1)
Z=matrix(0, nrow = length(X), ncol = 1)
for(n in c(1:length(X))){
  
  T[,1]=exp(-rowSums((sweep(G,2,G[n,],FUN="-"))^2)/sigma^2)
  D[,1]=sqrt(1+1-2*T)
  Z=1/D
  Z[n]=0
   P=matrix(0, nrow = length(X), ncol = length(X))
  P=1+K-as.vector(T)
  
  P=sweep(P,2,as.vector(T),FUN="-")
  
  W[n]=1-1/(length(X))*sqrt(t(Z)%*%P%*%Z)
  }
h=seq(0.03,0.03,1)
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



