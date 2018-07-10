set.seed(20)

x1=c(.1,.12,.11,.13,.14,.7,.7,.65,.66)
u<- seq(0,1,length.out = 200)
X=c(x1,u)
X=matrix(X)
noise <- rnorm(length(u),0,0.3)
Y=matrix(0, nrow = length(X), ncol = 1)
M=basicMethod=matrix(0, nrow = length(X)-1, ncol = length(X))
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
W=matrix(0, nrow = length(X)-1, ncol = 1)
g=matrix(0, nrow = length(X)-1, ncol = 2)
Ki=matrix(0, nrow = length(X)-1, ncol = length(X)-1)
for(j in c(1:length(X))){
  K[,j]=exp(-rowSums((sweep(G,2,G[j,],FUN="-"))^2)/sigma^2)
}
n=1
T=matrix(0, nrow = length(X)-1, ncol = 1)
D=matrix(0, nrow = length(X)-1, ncol = 1)
Z=matrix(0, nrow = length(X)-1, ncol = 1)
for(n in c(1:length(X))){
for(m in c(1:(length(X)-1))){
  if(n==1)
   { g=G[c((n+1):length(X)),]
   Ki=K[c((n+1):(length(X))),c((n+1):(length(X)))]}
  if(n>1 && n<length(X))
   { g=G[c(1:(n-1),(n+1):length(X)),]
   Ki=K[(c((1:(n-1)),((n+1):length(X)))),(c((1:(n-1)),((n+1):length(X))))]}
  if(n==length(X))
   {g=G[c(1:(n-1)),]
   Ki=K[c(1:(n-1)),c(1:(n-1))]}
  #if(n!=length(X))
  T[,1]=exp(-rowSums((sweep(g,2,g[m,],FUN="-"))^2)/sigma^2)
  #else if(n==length(x))
   # T[,1]=exp(-rowSums((sweep(g,2,g[n-1,],FUN="-"))^2)/sigma^2)
  D[,1]=sqrt(1+1-2*T)
  Z=1/D
 # if(n!=length(X))
  Z[m]=0
  P=matrix(0, nrow = length(X)-1, ncol = length(X)-1)
  P=1+Ki-as.vector(T)
  
  P=sweep(P,2,as.vector(T),FUN="-")
  #if(n!=length(X))
  W[m]=1-1/(length(X))*sqrt(t(Z)%*%P%*%Z)
}
h=seq(0.03,0.03,1)
h=matrix(h)
for(j in c(1:length(h))){
  for(i in c(1:(length(X)-1))){
    temp=(g[,1]-g[i,1])/h[j,1]
    Ker=1/(sqrt(2*pi))*exp(-(temp^2)/2)
    num1=t(as.matrix(Ker))%*%(g[,2])
    denom1=sum(Ker)
    
    basicMethod[i,n]=num1/denom1
    
    num=t(as.matrix(Ker*W))%*%g[,2]
    denom=(sum(Ker*W))
    M[i,n]=num/denom
  }
  }
}


