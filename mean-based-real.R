library(rstiefel)
library(ggplot2)
library(quantreg)
library(conquer)

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

factor_data = data.matrix(factor_data)[,c(1,2,3,4,5)]
portfolios = data.matrix(portfolios)

N = 6
T = dim(portfolios)[1]
p = k = 5   # number of common factor

# read data
AY = t(portfolios)

y = AY
X = array(0, c(N, T, p))
for (i in 1:N) {
  X[i,,] = factor_data
}

delta = 0.001
tau = 0.05

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
  if (conv_count >= 3) {
    print('Converged!')
    print(loss)
    break
  }
  prev_loss = loss
  print(loss)
}


baseF = F