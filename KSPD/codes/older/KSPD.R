#library(mvtnorm)
library(MASS)

set.seed(1)
n=50
cont=-50

X <- rnorm(n,0,1)
y <- 1 + 2*X + 4*X^2
#y <- 5 + 2*(class[,1]) + 4*(class[,1])^2
#+ 60*class[,1]
noise <- rnorm(n, mean=0, sd=0.5)
noisy.y <- y + noise
#plot(class[,1],y)

G=cbind(noisy.y,X) #,class[,2]
#pairs(G)
#G[n,]=rep(1,dim(G)[2])
G[n,]=c(cont,0)
plot(G[,2],G[,1])

sigma=5
K=matrix(0, nrow = n, ncol = n)
for(i in c(1:n)){
  for(j in c(1:n)){
    K[i,j]=exp(-sum((G[i,]-G[j,])^2)/sigma^2)
  }
}

W=T=D=Z=matrix(0, nrow = n, ncol = 1)
for(N in c(1:n)){

for(i in c(1:n)){
  T[i]=exp(-sum((G[i,]-G[N,])^2)/sigma^2)
  D[i]=sqrt(1+1-2*T[i])
  if(D[i]==0){Z[i]=0}
  else{Z[i]=1/D[i]}
}

P=matrix(0, nrow = n, ncol = n)
for(i in c(1:n)){
  for(j in c(1:n)){
    P[i,j]=1+K[i,j]-T[i]-T[j]
  }
}

W[N]=1-1/(n)*sqrt(t(Z)%*%P%*%Z)
}
W1=W-min(W)
W1=W1/sum(W1)
#W[n]=0
#print(W)

#G1=cbind(noisy.y,X,X^2)
#G1[n,]=c(cont,0,0)
G1=cbind(noisy.y,X^2)
G1[n,]=c(cont,0,0)
G2=data.frame(G1)
fit1=lm(noisy.y ~ ., data=G2)
print(fit1$coefficients)

fit2=lm(noisy.y ~ ., data=G2, weights=W1)
print(fit2$coefficients)

fit21 <- lm(noisy.y ~ poly(X,2), weights=W1)
print(fit21$coefficients)

fit3=rlm(noisy.y ~ ., data=G2)
print(fit3$coefficients)
