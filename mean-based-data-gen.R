library(pracma)

maxiter = 1
Case = 7
BETA=array(0, c(maxiter, 5,Case))
SIGMA=matrix(0, nrow = maxiter, ncol = Case)
NNN = matrix(0, nrow = maxiter, ncol = Case)
grandmean = rep(0, maxiter)
p = 5
TT = c(10,20,50,100,100,100,100)
NN = c(100,100,100,100,10,20,50)
BETA2 = array(0, c(maxiter, 10,Case))

for (kk in 1:Case) {
  T = TT[kk]
  N = NN[kk]
  oneT = rep(1, T)
  oneN = rep(1, N)
  M_T = diag(T) - oneT %*% t(oneT)/T
  M_N = diag(N) - oneN %*% t(oneN)/T
  
  for (niter in 1:maxiter) {
    r = 2
    p = 5
    k = 5
    
    # DGP
    lambda = matrix(rnorm(N*r), N, r)
    F = matrix(rnorm(T*r), T, r)
    e = 2*matrix(rnorm(T*N), T, N)
    XX = array(0, c(T,N,p))
    
    X1 = matrix(rnorm(T*N), T, N) + F %*% t(lambda) + matrix(1, nrow = T, ncol = r) %*% t(lambda) + F %*% matrix(1, nrow = r, ncol = N) + 1
    X2 = matrix(rnorm(T*N), T, N) + F %*% t(lambda) + matrix(1, nrow = T, ncol = r) %*% t(lambda) + F %*% matrix(1, nrow = r, ncol = N) + 1
    X3 = matrix(0, nrow = T, ncol = N)
    X4 = X3
    X5 = X3
    X3[1:T,1:N] = 1
    aa = rowSums(lambda) + rnorm(N)
    bb = rowSums(F) + rnorm(T)
    
    for (t in 1:T) {
      X4[t,] = aa
    }
    for (i in 1:N) {
      X5[,i] = bb
    }
    truebeta = c(1,3,5,2,4)
    XX[,,1] = X1
    XX[,,2] = X2
    XX[,,3] = X3
    XX[,,4] = X4
    XX[,,5] = X5
    X = XX[,,1:p]
    Y = F %*% t(lambda) + e
    for (k in 1:p) {
      Y = Y + X[,,k] * truebeta[k]
    }
    # DGP ends
    y = Y
    X = X
    
    # algorithm starts
    # F = (T^0.5)*rustiefel(T, r)
    # LAMUP = matrix(1, nrow = r, ncol = r)
    # LAMUP[upper.tri(LAMUP)] <- 0
    # LAMUP[lower.tri(LAMUP)] <- 0
    # LAMDOWN = matrix(0, nrow = N-r, ncol = r)
    # LAM = lambda #rbind(LAMUP, LAMDOWN)
    beta = rep(0, k)
    M_F = diag(T) - F %*% t(F) / T
    
    maxIter = 100
    conv_count = 0
    prev_loss = Inf
    
    for (i in 1:maxIter) {
      # projection matrix
      M_F = diag(T) - F %*% t(F) / T
      # update beta(F)
      mult_mat = matrix(0, nrow = k, ncol = k)
      cum_mat = rep(0, k)
      for (i in 1:N) {
        mult_mat = mult_mat + t(X[,i,]) %*% M_F %*% X[,i,]  # TBC
        cum_mat = cum_mat + t(X[,i,]) %*% M_F %*% y[,i]     # TBC
      }
      beta = solve(mult_mat) %*% cum_mat
      # calc W as error w/o latent factor
      W = matrix(0, nrow = T, ncol = N)
      for (i in 1:N) {
        W[,i] = y[,i] - X[,i,] %*% beta   # TBC
      }
      # decomp = eigen(W %*% t(W))
      F = svd(W %*% t(W))$u[,1:r]*sqrt(T)# decomp$vectors[,1:r]
      # update \Lambda
      lambda = (t(W) %*% F) / T
      loss = 0
      for (i in 1:N) {
        fac = y[,i] - X[,i,] %*% beta - F %*% lambda[i,]  # TBC
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
      # print(loss)
    }
    BETA[niter,,kk] = beta
    
  }
}