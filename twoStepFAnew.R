library(rstiefel)
library(ggplot2)
library(quantreg)
library(conquer)
library(adaHuber)

test_period = 4 # could be 1,2,3,4
# period to examine
tau = 0.05
r_tau = r = 6 # number of latent factor

# read data, periods
factor_data_p1 <- read.csv("C:/Users/liuhn/Desktop/test-data/factor_data_p1.csv")
portfolios_p1 <- read.csv("C:/Users/liuhn/Desktop/test-data/portfolio_data_p1.csv")

factor_data_p2 <- read.csv("C:/Users/liuhn/Desktop/test-data/factor_data_p2.csv")
portfolios_p2 <- read.csv("C:/Users/liuhn/Desktop/test-data/portfolio_data_p2.csv")

factor_data_p3 <- read.csv("C:/Users/liuhn/Desktop/test-data/factor_data_p3.csv")
portfolios_p3 <- read.csv("C:/Users/liuhn/Desktop/test-data/portfolio_data_p3.csv")

factor_data_p4 <- read.csv("C:/Users/liuhn/Desktop/test-data/factor_data_p4.csv")
portfolios_p4 <- read.csv("C:/Users/liuhn/Desktop/test-data/portfolio_data_p4.csv")

if (test_period == 1) {
  factor_data = factor_data_p1
  portfolios = portfolios_p1
} else if (test_period == 2) {
  factor_data = factor_data_p2
  portfolios = portfolios_p2
} else if (test_period == 3) {
  factor_data = factor_data_p3
  portfolios = portfolios_p3
} else if (test_period == 4) {
  factor_data = factor_data_p4
  portfolios = portfolios_p4
} else {
  stop("Period number is not valid")
  quit()
}

# --- read data end, start process data

colnames(factor_data) = c(1,2,3,4,5)
colnames(portfolios) = c(1,2,3,4,5,6)
factor_data = data.matrix(factor_data)[,c(1,2,3,4,5)]
portfolios = data.matrix(portfolios)

# --- processing ends, setting constants

# order: 

N = 6
T = dim(portfolios)[1]
p = k = 5   # number of common factor

AY = t(portfolios)

y = AY
X = array(0, c(N, T, p))
for (i in 1:N) {
  X[i,,] = factor_data
}

delta = 0.001

# quantile FA
# initiate
F = (T^0.5)*rustiefel(T, r_tau)
if(N >= r_tau) {
  LAMUP = matrix(1, nrow = r_tau, ncol = r_tau)
  LAMUP[upper.tri(LAMUP)] <- 0 
  LAMDOWN = matrix(1, nrow = N-r_tau, ncol = r_tau)
  LAM = rbind(LAMUP, LAMDOWN)#matrix(runif(N*r_tau,-1,1),nrow=N,ncol=r_tau)# rbind(LAMUP, LAMDOWN)
} else {
  # LAMUP = matrix(1, nrow = N, ncol = N)
  # LAMUP[upper.tri(LAMUP)] <- 0 
  # LAMDOWN = matrix(0, nrow = N, ncol = r_tau-N)
  LAM = cbind(LAMUP, LAMDOWN)#matrix(runif(N*r_tau,-1,1),nrow=N,ncol=r_tau)#cbind(LAMUP, LAMDOWN)
}
# LAM = matrix(1, nrow = min(r_tau, N), ncol = r_tau)
# LAM[upper.tri(LAM)] <- 0 
B = matrix(0, nrow = N, ncol = p+1)
# F <- matrix(runif(N*r_tau,0,2),nrow=N,ncol=r_tau)

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
    cur_y = y_t - B%*%c(x_t[1,],1)
    # sec_fit = conquer(X = LAM, Y = cur_y, tau=tau)#$coeff
    # x_0 = rep(0.1, r_tau)
    sec_fit = rq(cur_y~0+LAM,tau=tau)$coefficients
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

qF_step1 = F
L_step1 = LAM
B_step1 = B

# calculate residual
for (i in 1:N) {
  for (t in 1:T) {
    y[i,t] = y[i,t] - c(X[i,t,],1) %*% B_step1[i,] - qF_step1[t,] %*% L_step1[i,]
  }
}

# evaluate ES
for (i in 1:N) {
  res = y[i,]
  x = X[i,,]
  fitting = adaHuber.reg(x, res, method='adaptive')
  y[i,] = cbind(1, x) %*% fitting$coef
}

F = (T^0.5)*rustiefel(T, r_tau)
if(N >= r_tau) {
  LAMUP = matrix(1, nrow = r_tau, ncol = r_tau)
  LAMUP[upper.tri(LAMUP)] <- 0 
  LAMDOWN = matrix(1, nrow = N-r_tau, ncol = r_tau)
  LAM = rbind(LAMUP, LAMDOWN)
} else {
  LAMUP = matrix(1, nrow = N, ncol = N)
  LAMUP[upper.tri(LAMUP)] <- 0 
  LAMDOWN = matrix(0, nrow = N, ncol = r_tau-N)
  LAM = cbind(LAMUP, LAMDOWN)
}
B = matrix(1, nrow = N, ncol = p+1)

beta = rep(0, k)
M_F = diag(T) - F %*% t(F) / T

maxIter = 400
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
  if (conv_count >= 2) {
    print('Converged!')
    print(loss)
    break
  }
  prev_loss = loss
  print(loss)
}


EF_step2 = F