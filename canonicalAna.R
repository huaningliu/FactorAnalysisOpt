library(CCA)
library(CCP)

# we have: quantileF, ESF

# dim(quantileF) == dim(ESF)

# colnames(qF_step1) = c("q1", "q2", "q3", "q4", "q5")
# colnames(EF_step2) = c("e1", "e2", "e3", "e4", "e5")
# colnames(baseF) = c("b1", "b2", "b3", "b4", "b5")

# comp = CCA::matcor(baseF, ESF)

# baseF, qF_step1, EF_step2

cc1 <- cc(qF_step1, EF_step2)
rho <- cc1$cor
n <- dim(EF_step2)[1]
p <- length(qF_step1[1,])
q <- length(EF_step2[1,])

CCP::p.asym(rho, n, p, q, tstat = "Wilks")

cc1 <- cc(baseF, EF_step2)
rho <- cc1$cor
n <- dim(EF_step2)[1]
p <- length(baseF[1,])
q <- length(EF_step2[1,])

CCP::p.asym(rho, n, p, q, tstat = "Wilks")