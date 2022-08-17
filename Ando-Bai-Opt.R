library(rstiefel)
library(conquer)

# constants
N = 100     # number of units
T = 100    # terminal time
tau = 0.05 # or 0.95 in the paper
r_tau = 3  # number of common factors, different by Table[1] in Ando-Bai work
p = 8      # dim of x as data
delta = 0.001  # stopping criterion

# initiate
F = (T^0.5)*rustiefel(T, r_tau)
LAMUP = matrix(1, nrow = r_tau, ncol = r_tau) # matrix(runif(r_tau^2, -1, 1), nrow = r_tau, ncol = r_tau)
LAMUP[upper.tri(LAMUP)] <- 0 
LAMDOWN = matrix(1, nrow = N-r_tau, ncol = r_tau) # matrix(runif((N-r_tau)*r_tau, -1, 1), nrow = N-r_tau, ncol = r_tau)
LAM = rbind(LAMUP, LAMDOWN)
B = matrix(0, nrow = N, ncol = p+1)

# data
y = AY
X = array(c(AX), c(N, T, p))

optim_f <- function(f, Y, X) {
  mu = Y - X %*% f
  return(sum(mu*(tau - 1*(mu < 0))))
}  ## rq replace optim

# start the algorithm
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
    b_i = B[j,]
    cur_y = y_t - rowSums(cbind(x_t, rep(1, length(x_t[,1]))) * B)
    # sec_fit = conquer(X = LAM, Y = cur_y, tau=tau)#$coeff
    # x_0 = rep(0.1, r_tau)
    sec_fit = rq(cur_y~0+LAM,tau=tau)$coefficients# optim(x_0, optim_f, method = "BFGS", Y = cur_y, X = LAM)
    # f_t = sec_fit$coeff
    # print(f_t)
    F[j,] = sec_fit
  }
  # y_i hat error calc #
  for (i in 1:N) { # N
    y_i = y[i,]
    x_i = X[i,,]
    X_i = cbind(x_i, F)
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
# plot(record, type='l', main = 'Criteria Plot', xlab = 'Iteration', ylab = 'Loss')

# mean-based factor model, least squared

# PCA <-> lm
# ES <- QR-LS
library(verification)
# out calc
out = matrix(0, nrow = N, ncol = T)
for (t in 1:T) {
  for (i in 1:N) {
    x_it = X[i,t,]
    x_it = c(x_it, 1)
    b_i = B[i,]
    f_t = F[t,]
    lamb_i = LAM[i,]
    out[i,t] = x_it %*% b_i + f_t %*% lamb_i
  }
}