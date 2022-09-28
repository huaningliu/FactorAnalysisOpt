# read data
asia <- read.csv("C:/Users/liuhn/Desktop/test-data/asia.csv")
euro <- read.csv("C:/Users/liuhn/Desktop/test-data/euro.csv")
globe <- read.csv("C:/Users/liuhn/Desktop/test-data/globe.csv")
namer <- read.csv("C:/Users/liuhn/Desktop/test-data/namer.csv")
japan <- read.csv("C:/Users/liuhn/Desktop/test-data/japan.csv")
exus <- read.csv("C:/Users/liuhn/Desktop/test-data/exus.csv")

asia = data.matrix(asia)[,c('Mkt.RF','SMB','HML','RMW','CMA', 'RF')]
euro = data.matrix(euro)[,c('Mkt.RF','SMB','HML','RMW','CMA', 'RF')]
globe = data.matrix(globe)[,c('Mkt.RF','SMB','HML','RMW','CMA', 'RF')]
namer = data.matrix(namer)[,c('Mkt.RF','SMB','HML','RMW','CMA', 'RF')]
japan = data.matrix(japan)[,c('Mkt.RF','SMB','HML','RMW','CMA', 'RF')]
exus = data.matrix(exus)[,c('Mkt.RF','SMB','HML','RMW','CMA', 'RF')]

colnames(asia) = c(1,2,3,4,5,6)
colnames(euro) = c(1,2,3,4,5,6)
colnames(globe) = c(1,2,3,4,5,6)
colnames(namer) = c(1,2,3,4,5,6)
colnames(japan) = c(1,2,3,4,5,6)
colnames(exus) = c(1,2,3,4,5,6)

# period to examine
period = 2
tau = 0.05
r_tau = r = R = 5

if (period == 1) {
  start = 4306
  end = 4914
} else if (period == 2) {
  start = 5002
  end = 5871
} else if (period == 3) {
  start = 5872
  end = 6457
} else if (period == 4) {
  start = 4306
  end = 6457
} else {
  start = 1
  end = 1000
}


asia = asia[(start:end),]
euro = euro[start:end,]
globe = globe[start:end,]
namer = namer[start:end,]
japan = japan[start:end,]

# order: Global, North America, Europe, Asia Pacific

N = 5
T = length(asia[,1])
p = k = 5

TAU <- 0.05

AY = matrix(0, nrow = N, ncol = T)
AY[1,] = globe[,1]
AY[2,] = namer[,1]
AY[3,] = euro[,1]
AY[4,] = asia[,1]
AY[5,] = japan[,1]
# AY[6,] = exus[,6]

y = AY
X = array(0, c(N, T, p))
X[1,,] = globe[,(2:6)]
X[2,,] = namer[,(2:6)]
X[3,,] = euro[,(2:6)]
X[4,,] = asia[,(2:6)]
X[5,,] = japan[,(2:6)]
# X[6,,] = exus[,(1:5)]

AX = matrix(0, nrow = N*T, ncol = p)
for (i in 1:T) {
  AX[(N*(i-1)+1):(N*i),] = X[,i,]
}

#===================================#
#Estimation 
#===================================#

P = T

XB <- matrix(0,nrow=N,ncol=P)
FL <- matrix(0,nrow=N,ncol=P)
B <- matrix(0,nrow=p+1,ncol=P)

for(j in 1:P){
  y <- AY[,j]; X <- AX[(N*(j-1)+1):(N*j),]
  fit <- rq(y~X,tau=TAU)
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