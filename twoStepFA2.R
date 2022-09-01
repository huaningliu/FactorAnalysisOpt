library(rstiefel)
library(ggplot2)
library(quantreg)
library(conquer)

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

asia = asia[1:500,]
euro = euro[1:500,]
globe = globe[1:500,]
namer = namer[1:500,]
japan = japan[1:500,]

# order: Global, North America, Europe, Asia Pacific

AY = matrix(0, nrow = N, ncol = T)
AY[1,] = globe[,6]
AY[2,] = namer[,6]
AY[3,] = euro[,6]
AY[4,] = asia[,6]
AY[5,] = japan[,6]
# AY[6,] = exus[,6]

y = AY
X = array(0, c(N, T, p))
X[1,,] = globe[,(1:5)]
X[2,,] = namer[,(1:5)]
X[3,,] = euro[,(1:5)]
X[4,,] = asia[,(1:5)]
X[5,,] = japan[,(1:5)]
# X[6,,] = exus[,(1:5)]

N = 5
T = length(asia[,1])
p = k = 5
r_tau = r = 5

delta = 0.001
tau = 0.05

# quantile residual
for (i in 1:N) {
  y_i = y[i,]
  X_i = X[i,,]
  qr_fit = conquer(X_i, y_i, tau=tau)
  qr_nres = qr_fit$res * (qr_fit$res <= 0)
  y[i,] = qr_nres/tau + (y_i - qr_fit$res)
}

F = (T^0.5)*rustiefel(T, r_tau)
LAMUP = matrix(1, nrow = r_tau, ncol = r_tau)
LAMUP[upper.tri(LAMUP)] <- 0 
LAMDOWN = matrix(1, nrow = N-r_tau, ncol = r_tau)
LAM = rbind(LAMUP, LAMDOWN)
# LAM = matrix(1, nrow = min(r_tau, N), ncol = r_tau)
# LAM[upper.tri(LAM)] <- 0 
B = matrix(1, nrow = N, ncol = p+1)

beta = rep(0, k)
M_F = diag(T) - F %*% t(F) / T

maxIter = 10
conv_count = 0   # converge
prev_loss = Inf  # converge

# X: N*T*p
# y: N*T

for (i in 1:maxIter) {
  # projection matrix
  M_F = diag(T) - F %*% t(F) / T
  # update beta(F)
  mult_mat = matrix(0, nrow = k, ncol = k)
  cum_mat = rep(0, k)
  for (i in 1:N) {
    # print("1")
    mult_mat = mult_mat + t(X[i,,]) %*% M_F %*% X[i,,]  # TBC
    cum_mat = cum_mat + t(X[i,,]) %*% M_F %*% y[i,]     # TBC
  }
  beta = solve(mult_mat) %*% cum_mat
  # calc W as error w/o latent factor
  W = matrix(0, nrow = T, ncol = N)
  for (i in 1:N) {
    # print("2")
    W[,i] = y[i,] - X[i,,] %*% beta   # TBC
  }
  # decomp = eigen(W %*% t(W))
  F = svd(W %*% t(W))$u[,1:r]*sqrt(T)# decomp$vectors[,1:r]
  # update \Lambda
  lambda = (t(W) %*% F) / T
  loss = 0
  # calc loss
  for (i in 1:N) {
    # print("3")
    fac = y[i,] - X[i,,] %*% beta - F %*% lambda[i,]  # TBC
    loss = loss + t(fac) %*% fac
  }
  # another loss = tr(t(W) %*% M_F %*% W)
  # check convergence
  if (abs(loss - prev_loss) < 0.001) {
    conv_count = conv_count + 1
  } else {
    conv_count = 0
  }
  if (conv_count >= 3) {
    print('Converged!')
    print(loss)
    break
  }
  prev_loss = loss
  print(loss)
}


ESF = F