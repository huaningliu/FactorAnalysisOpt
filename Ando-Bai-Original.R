#===================================#
#Estimation 
#===================================#

XB <- matrix(0,nrow=N,ncol=P)
FL <- matrix(0,nrow=N,ncol=P)
B <- matrix(0,nrow=p+1,ncol=P)

for(j in 1:P){
  y <- AY[,j]; X <- AX[(N*(j-1)+1):(N*j),]
  fit <- conquer(X = X, Y = y, tau=tau)$coeff # rq(y~X,tau=TAU)
  B[,j] <- fit# (fit$coefficients)
  XB[,j] <- cbind(1,X)%*%B[,j]
}

Z <- AY-XB
VEC <- eigen(Z%*%t(Z))$vectors; 
F <- sqrt(N)*(VEC)[,1:R];
L <- t(solve(t(F)%*%F)%*%t(F)%*%Z);
FL <- F%*%t(L)

B.old <- B
FL.old <- FL

for(ITE in 1:100){
  
  for(i in 1:N){
    y <- AY[i,]; X <- L
    fit <- rq(y-XB[i,]~X+0,tau=TAU)
    # print((fit$coefficients)[1:R])
    F[i,] <- (fit$coefficients)[1:R]
  }
  
  QRF <- qr(F);
  F <- sqrt(N)*qr.Q(QRF);
  QRL <- qr(qr.R(QRF)%*%t(L));
  F <- sqrt(N)*qr.Q(QRF)%*%qr.Q(QRL);
  
  for(j in 1:P){
    y <- AY[,j]; X <- cbind(AX[(N*(j-1)+1):(N*j),],F)
    fit <- rq(y~X,tau=TAU) #conquer(X = X, Y = y, tau=tau)# a$coeff # rq(y~X,tau=TAU)
    B[,j] <- (fit$coeff)[1:(p+1)]
    L[j,] <- (fit$coeff)[(p+2):(p+R+1)]
    XB[,j] <- cbind(1,AX[(N*(j-1)+1):(N*j),])%*%B[,j]
    FL[,j] <- F%*%L[j,]
  }
  
  Up <- sum((B.old-B)^2)/(length(B))+sum((FL.old-FL)^2)/(length(FL))
  
  B.old <- B
  FL.old <- FL
  
  print(Up)
  
  if(Up<=0.001){break}
  
}

standard = matrix(0, nrow = N, ncol = T)
for (t in 1:T) {
  for (i in 1:N) {
    x_it = AX[(N*(t-1)+1+i-1),]
    x_it = c(x_it, 1)
    b_i = B[,i]
    f_t = F[t,]
    lamb_i = L[i,]
    standard[i,t] = x_it %*% b_i + f_t %*% lamb_i
  }
}


correct_B = B