library(rstiefel)
library(ggplot2)
library(quantreg)

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

N = 5
T = length(asia[,1])
p = 5
r_tau = 5

delta = 0.001
tau = 0.05

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

# initiate
F = (T^0.5)*rustiefel(T, r_tau)
LAMUP = matrix(1, nrow = r_tau, ncol = r_tau)
LAMUP[upper.tri(LAMUP)] <- 0 
LAMDOWN = matrix(1, nrow = N-r_tau, ncol = r_tau)
LAM = rbind(LAMUP, LAMDOWN)
# LAM = matrix(1, nrow = min(r_tau, N), ncol = r_tau)
# LAM[upper.tri(LAM)] <- 0 
B = matrix(1, nrow = N, ncol = p+1)

# quantile-based algo
K = 100
record = c()
for (i in 1:K) {
  # concat x and f together for b and lambda
  old_B = B
  old_LAM = LAM
  old_F = F
  # regress on f
  for (j in 1:T) {
    x_t = X[,j,]
    y_t = y[,j]
    # b_i = B[j,]
    cur_y = y_t - rowSums(cbind(x_t, rep(1, length(x_t[,1]))) * B)
    # sec_fit = conquer(X = LAM, Y = cur_y, tau=tau)#$coeff
    # x_0 = rep(0.1, r_tau)
    # print("haha")
    sec_fit = rq(cur_y~0+LAM,tau=tau)$coefficients
    # print("nono")
    # f_t = sec_fit$coeff
    # print(f_t)
    F[j,] = sec_fit
  }
  # y_i hat error calc #
  for (i in 1:N) { # N
    y_i = y[i,]
    x_i = X[i,,]
    X_i = cbind(x_i, F)
    print("Nothing happend")
    first_fit = rq(y_i~X_i,tau=tau)$coefficients
    # first_fit_2 = conquer(X = X_i, Y = y_i, tau=tau)$coeff
    l = length(first_fit)
    b_i = c(first_fit[2:(2+p-1)], first_fit[1]) # c(first_fit[(2+r_tau):l], first_fit[1])
    lambda_i = first_fit[(2+p):l] # first_fit[2:(2+r_tau-1)]
    B[i,] = b_i
    LAM[i,] = lambda_i
  }
  
  # QR decomp for F
  QR_F = qr(F)
  Q_F = qr.Q(QR_F)
  R_F = qr.R(QR_F)
  
  # QR decomp for RLam
  QR = qr(R_F %*% t(LAM))
  Q_lambda = qr.Q(QR)
  R_lambda = qr.R(QR)
  
  # update F and LAM
  F = sqrt(T)*(Q_F %*% Q_lambda)
  LAM = t(R_lambda)
  
  # calculate converge criteria
  criteria = sum((B - old_B)^2) / (length(B)) + sum(((F %*% t(LAM)) - (old_F %*% t(old_LAM)))^2) / (N*T)
  # print(LAM)
  record = c(record, criteria)
  if (criteria < delta) {
    print(criteria)
    print("the algorithm is done!")
    break
  }
  print(criteria)
}

quantileF = F