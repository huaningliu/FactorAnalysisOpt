library(CCA)
library(CCP)

# we have: quantileF, ESF

dim(quantileF) == dim(ESF)

colnames(quantileF) = c("q1", "q2", "q3", "q4", "q5")
colnames(ESF) = c("e1", "e2", "e3", "e4", "e5")


comp = CCA::matcor(quantileF, ESF)

cc1 <- cc(quantileF, ESF)
rho <- cc1$cor
n <- dim(quantileF)[1]
p <- 5# length(quantileF)
q <- 5# length(ESF)

CCP::p.asym(rho, n, p, q, tstat = "Wilks")