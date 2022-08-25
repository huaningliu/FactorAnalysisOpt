
#This code implements the frequentist estimation procedure in Section 3.1

library(quantreg)

R <- 5
p <- 8
VAR <- 1

N <- 100
P <- 3  # T in my case

TAU <- 0.05

U <- matrix(runif(N*P,0,1),nrow=N,ncol=P)

LAM <- matrix(runif(P*R,-1,1),nrow=P,ncol=R)
FAC <- matrix(runif(N*R,0,2),nrow=N,ncol=R)

TFL <- FL <- matrix(0,nrow=N,ncol=P)

for(i in 1:N){
  for(j in 1:P){
    
    B <- LAM[j,]+0.1*U[i,j]
    if(U[i,j]<=0.2){FL[i,j] <- FAC[i,1:3]%*%B[1:3]}
    if(0.2<=U[i,j] && U[i,j]<=0.8){FL[i,j] <- FAC[i,1:4]%*%B[1:4]}
    if(0.8<=U[i,j]){FL[i,j] <- FAC[i,1:5]%*%B[1:5]}
    
    B <- LAM[j,]+0.1*TAU
    if(U[i,j]<=0.2){TFL[i,j] <- FAC[i,1:3]%*%B[1:3]}
    if(0.2<=U[i,j] && U[i,j]<=0.8){TFL[i,j] <- FAC[i,1:4]%*%B[1:4]}
    if(0.8<=U[i,j]){TFL[i,j] <- FAC[i,1:5]%*%B[1:5]}
    
  }
}


AX <- matrix(runif(p*P*N,0,1),nrow=P*N)

for(i in 1:N){
  for(j in 1:P){
    a <- FAC[i,1]^2; b <- LAM[j,1]+0.01*U[i,j]
    AX[(N*(j-1)+1):(N*j),1] <- AX[(N*(j-1)+1):(N*j),1]+0.02*a+0.02*b
    a <- FAC[i,2]^2; b <- LAM[j,2]+0.01*U[i,j]
    AX[(N*(j-1)+1):(N*j),3] <- AX[(N*(j-1)+1):(N*j),3]-0.01*a+0.02*b
    a <- FAC[i,3]^2; b <- LAM[j,3]+0.01*U[i,j]
    AX[(N*(j-1)+1):(N*j),5] <- AX[(N*(j-1)+1):(N*j),5]+0.02*a+0.02*b
  }
}

TXB <- XB <- matrix(0,nrow=N,ncol=P)

for(i in 1:N){
  for(j in 1:P){
    
    X <- AX[(N*(j-1)+1):(N*j),]
    B <- c(-1,1,-1,1,-1,1,-1,-1)+j/P+0.1*U[i,j]
    if(U[i,j]<=0.2){XB[i,j] <- X[i,1:p]%*%B[1:p]}
    if(0.2<=U[i,j] && U[i,j]<=0.8){XB[i,j] <- X[i,1:p]%*%B[1:p]}
    if(0.8<=U[i,j]){XB[i,j] <- X[i,1:p]%*%B[1:p]}
    
    X <- AX[(N*(j-1)+1):(N*j),]
    B <- c(-1,1,-1,1,-1,1,-1,-1)+j/P+0.1*TAU
    if(U[i,j]<=0.2){TXB[i,j] <- X[i,1:p]%*%B[1:p]}
    if(0.2<=U[i,j] && U[i,j]<=0.8){TXB[i,j] <- X[i,1:p]%*%B[1:p]}
    if(0.8<=U[i,j]){TXB[i,j] <- X[i,1:p]%*%B[1:p]}
    
  }
}

TB <- matrix(0,p,P)
for(j in 1:P){TB[,j] <- c(-1,1,-1,1,-1,1,-1,-1)+j/P+0.1*TAU}


TXBFL <- TXB+TFL+qnorm(TAU,0,VAR)

AY <- XB+FL+qnorm(U,0,VAR)

if(TAU==0.05){R <- 3}
if(TAU==0.50){R <- 4}
if(TAU==0.95){R <- 5}