library(rgl)
set.seed(20)

library(readr)
 realrr <- read_csv("C:/Users/neeraj yadav/Downloads/realrr.csv")
 attach(realrr)
 size=dim(realrr)
 len=size[1]
#dim_data including y
#feature matrix containing x1,x2..,xn
 mat=matrix(0, nrow = size[1], ncol = size[2])
dim_data=size[2]
mat=cbind(X111,X222,X333,X444,X555);

Final=basic=G=matrix(0, nrow = size[1], ncol = dim_data)
K=P=matrix(0, nrow = len, ncol =len )
Y=M=T=D=Z=W=basicMethod=matrix(0, nrow = len, ncol = 1)
Y=X666


G=cbind(X111,X222,X333,X444,X555,X666)

sigma=2^{-2:18}
h=c(seq(0.01,0.1,length.out=10),seq(0.2,5,length.out=10))

Robust_dist=Nadaraya_Watson_dist=matrix(0, nrow = length(sigma), ncol = length(h))
h=matrix(h)
for(k in c(1:length(sigma))){
  for(j in c(1:len)){
    K[,j]=exp(-rowSums((sweep(G,2,G[j,],FUN="-"))^2)/(sigma[k])^2)
  }
  
  for(n in c(1:size[1])){
    for(i in c(1:size[1])){
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
    for(i in c(1:length(X111))){
    # remove->  temp=(sqrt((X[,1]-X[i,1])^2+(X2[,1]-X2[i,1])^2))/h[j,1]
      temp=rowSums((sweep(mat,2,mat[i,],FUN="-"))^2)/h[j,1]
      Ker=1/(sqrt(2*pi))*exp(-(temp)/2)
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
cat("Robust-distance-min = ",min(Robust_dist),"\n")

ind=which(Robust_dist == min(Robust_dist), arr.ind = TRUE)
cat("h_min =",h[ind[1,2]],"\n")
cat("sigma_min =",sigma[ind[1,1]],"\n")

cat("Nadarya-Watson-distance-min = ",min(Nadaraya_Watson_dist),"\n")
index=which(Nadaraya_Watson_dist == min(Nadaraya_Watson_dist), arr.ind = TRUE)
cat("h_min =",h[index[1,2]],"\n")
cat("sigma_min =",sigma[index[1,1]],"\n")





