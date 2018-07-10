set.seed(1)

u<- seq(0,1,length.out = 200)
x1=c(.1,.12,.11,.13,.14,.7,.7,.65,.66)
X=c(x1,u)
X=matrix(X)
Fun=300*(X^3-3*X^4+3*X^5-X^6) 

W=T=D=Z=M=basicMethod=matrix(0, nrow = length(X), ncol = 1)
K=P=matrix(0, nrow = length(X), ncol = length(X))
G=matrix(0, nrow = length(X), ncol = 2)

y1=c(5,5.1,4.9,4.8,5.2,-1,-.9,-.9,-.88)
noise <- rnorm(length(u),0,0.3)
y2=300*(u^3-3*u^4+3*u^5-u^6) +as.matrix(noise)
Y=c(y1,y2)
Y=matrix(Y)

for(i in c(1:length(X))){
  for(j in c(1:2)){
    if(j==1){
      G[i,j]=X[i]}
    else{
      G[i,j]=Y[i]}
  }
}

#sigma=8
#sigma=c(.8,8,80,800,4000,10000,20000,30000,40000,50000,60000)
sigma=2^{-4:7}
h=c(seq(0.01,0.1,length.out=10),seq(0.2,5,length.out=10))

Robust_dist=Nadaraya_Watson_dist=matrix(0, nrow = length(sigma), ncol = length(h))
h=matrix(h)

for(k in c(1:length(sigma))){
for(j in c(1:length(X))){
  K[,j]=exp(-rowSums((sweep(G,2,G[j,],FUN="-"))^2)/(sigma[k])^2)
}

for(n in c(1:length(X))){
  
  T[,1]=exp(-rowSums((sweep(G,2,G[n,],FUN="-"))^2)/(sigma[k])^2)
  D[,1]=sqrt(1+1-2*T)
  Z=1/D
  Z[n]=0
 
  P=1+K-as.vector(T)
  
  P=sweep(P,2,as.vector(T),FUN="-")
  
  W[n]=1-1/(length(X))*sqrt(t(Z)%*%P%*%Z)
}

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
  
  
  Robust_dist[k,j]=(sum((Fun-M)^2))/length(Fun)
  Nadaraya_Watson_dist[k,j]=(sum((Fun-basicMethod)^2))/length(Fun)

  
}
}
cat("Robust-distance-min= ",min(Robust_dist),"\n")

ind=which(Robust_dist == min(Robust_dist), arr.ind = TRUE)
cat("h_min=",h[ind[1,2]],"\n")
cat("sigma_min=",sigma[ind[1,1]],"\n")

cat("Nadarya-Watson-distance-min= ",min(Nadaraya_Watson_dist),"\n")
index=which(Nadaraya_Watson_dist == min(Nadaraya_Watson_dist), arr.ind = TRUE)
cat("h_min=",h[index[1,2]],"\n")
cat("sigma_min=",sigma[index[1,1]],"\n")

