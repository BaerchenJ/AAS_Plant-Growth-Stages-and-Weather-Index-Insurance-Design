##############################################
# File: Train-W.r 
# Author: [Jing Zou] 
# Created: [2023.06.10] 
# Description: [estimate a whole-cycle Generalized Additive Model using PSANOVA method] 
##############################################

WCPSA <- function(y, x1, x2, nseg){
  
  n = length(y)
  
  ## definitions of the following three parameters please refer to 
  ## Lee, D. J., Durbán, M., & Eilers, P. (2013). Efficient two-dimensional smoothing with P-spline ANOVA mixed models and nested bases. Computational Statistics & Data Analysis, 61, 22-37.
  
  bdeg = 3  
  
  div = 3
  
  pord = 2
  
  ########################################################### first phase for x1 and x2
  MM1 <- MM.basis(x1, min(x1) - 0.01, max(x1) + 0.01, nseg, bdeg, pord)  # bases for x1
  X1 <- MM1$X
  G1 <- MM1$Z
  d1 <- MM1$d
  B1 <- MM1$B
  
  MM2 <- MM.basis(x2, min(x2) - 0.01, max(x2) + 0.01, nseg, bdeg, pord)  # bases for x2
  X2 <- MM2$X
  G2 <- MM2$Z
  d2 <- MM2$d
  B2 <- MM2$B
  
  MM1n <- MM.basis(x1, min(x1) - 0.01, max(x1) + 0.01, nseg/div, bdeg, pord)  # Nested bases for x1
  G1n <- MM1n$Z
  d1n <- MM1n$d
  B1n <- MM1n$B
  
  MM2n <- MM.basis(x2, min(x2) - 0.01, max(x2) + 0.01, nseg/div, bdeg, pord)  # Nested bases for x2
  G2n <- MM2n$Z
  d2n <- MM2n$d
  B2n <- MM2n$B
  
  c1 <- ncol(B1)
  c2 <- ncol(B2)
  
  c1n <- ncol(B1n)
  c2n <- ncol(B2n)
  
  one1. <- matrix(X1[, 1], ncol = 1)
  one2. <- matrix(X2[, 1], ncol = 1)
  
  x1. <- matrix(X1[, 2], ncol = 1)
  x2. <- matrix(X2[, 2], ncol = 1)
  
  X <- rowtens(X2, X1)  # -> Fixed effects
  d12 <- c(rep(1, c2n - pord) %x% d1n + d2n %x% rep(1, c1n - pord))
  Delta1 <- diag(1/sqrt(d1))
  Delta2 <- diag(1/sqrt(d2))
  Delta12 <- diag(1/sqrt(d12))
  
  # random effects matrices
  Z1 <- G1 %*% Delta1  # smooth random comp. fx1
  Z2 <- G2 %*% Delta2  # smooth random comp. fx2
  Z1x2 <- rowtens(x2., G1n)  # linear:smooth comp. x2:fx1
  Z2x1 <- rowtens(G2n, x1.)  # linear:smooth comp. fx2:x1
  Z12 <- rowtens(G2n, G1n) %*% Delta12  # smooth interaction  fx1:fx2
  
  # Random effects matrix
  Z <- cbind(Z1, Z2, Z1x2, Z2x1, Z12)
  
  # Compute cross-products
  XtX. <- crossprod(X)
  XtZ. <- crossprod(X, Z)
  ZtZ. <- crossprod(Z)
  Xty. <- crossprod(X, y)
  Zty. <- crossprod(Z, y)
  yty <- sum(y^2)
  
  V <- construct.block2(XtX., XtZ., ZtZ.)
  u <- c(Xty., Zty.)
  
  np <- c(ncol(X), ncol(Z1), ncol(Z2), ncol(Z1x2), ncol(Z2x1), ncol(Z12))
  
  nblocks <- 6  # number of components (1 fixed component + 5 random effects)
  
  idx <- rep(1:nblocks, np)
  end1 <- proc.time()[3]
  nn <- sum(np)
  
  # Initialize algorithm
  la <- c(0, c(rep(1, nblocks)))
  
  thr <- 1e-06
  
  #library(MASS)
  library(svd)
  # Run Schall's EM algorithm
  
  for (it in 1:1000) {
    
    # Build penalty matrix
    P <- diag(la[idx])
    pv <-  V + P
    Vu <- cbind(V, u)
    
    # solve equations
    # if (is.singular(pv) == TRUE){
    s <- propack.svd(pv)
    bb <- s$v %*% diag(1/s$d) %*% t(s$u)
    Hb <-  bb %*% Vu
    #Hb <- ginv(V + P) %*% Vu
    HH <- Hb[, 1:nn]  
    b <- Hb[, nn + 1]
    
    
    # Compute effective dimensions and variances
    ed <- tapply(diag(HH), idx, sum)
    ssb <- tapply(b^2, idx, sum)
    tau <- ssb/ed + thr
    ssr <- yty - crossprod(b, 2 * u - V %*% b)
    sig2 <- (ssr/(n - sum(ed)))
    
    # New lambdas and convergence check
    lanew <- c(0, rep(sig2, nblocks - 1))/tau
    dla <- mean(abs(log10(la[2:nblocks]) - log10(lanew[2:nblocks])))
    la <- lanew
    if (dla < 1e-06) 
      break
  } 
  M <- cbind(X,Z)
  Fit <- M %*% b
  Fit <- as.numeric(Fit)
  residuals = y - Fit
  res <- list(b = b, Fit = Fit, M = M, idx = idx, dla = dla, it = it)
  res
}
