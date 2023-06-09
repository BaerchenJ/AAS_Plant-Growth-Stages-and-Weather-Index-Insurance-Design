
PDPSA <- function(y, x1, x2, x3, x4, x5, x6, x7, x8, nseg){
  
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
  
  XS1 <- rowtens(X2, X1)  # -> Fixed effects
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
  ZS1 <- cbind(Z1, Z2, Z1x2, Z2x1, Z12)
  
  #######################################################  second phase for x3 and x4
  MM3 <- MM.basis(x3, min(x3) - 0.01, max(x3) + 0.01, nseg, bdeg, pord)  # bases for x3
  X3 <- MM3$X
  G3 <- MM3$Z
  d3 <- MM3$d
  B3 <- MM3$B
  
  MM4 <- MM.basis(x4, min(x4) - 0.01, max(x4) + 0.01, nseg, bdeg, pord)  # bases for x4
  X4 <- MM4$X
  G4 <- MM4$Z
  d4 <- MM4$d
  B4 <- MM4$B
  
  MM3n <- MM.basis(x3, min(x3) - 0.01, max(x3) + 0.01, nseg/div, bdeg, pord)  # Nested bases for x3
  G3n <- MM3n$Z
  d3n <- MM3n$d
  B3n <- MM3n$B
  
  MM4n <- MM.basis(x4, min(x4) - 0.01, max(x4) + 0.01, nseg/div, bdeg, pord)  # Nested bases for x4
  G4n <- MM4n$Z
  d4n <- MM4n$d
  B4n <- MM4n$B
  
  c3 <- ncol(B3)
  c4 <- ncol(B4)
  
  c3n <- ncol(B3n)
  c4n <- ncol(B4n)
  
  one3. <- matrix(X3[, 1], ncol = 1)
  one4. <- matrix(X4[, 1], ncol = 1)
  
  x3. <- matrix(X3[, 2], ncol = 1)
  x4. <- matrix(X4[, 2], ncol = 1)
  
  XS2 <- rowtens(X4, X3)  # -> Fixed effects
  d34 <- c(rep(1, c4n - pord) %x% d3n + d4n %x% rep(1, c3n - pord))
  Delta3 <- diag(1/sqrt(d3))
  Delta4 <- diag(1/sqrt(d4))
  Delta34 <- diag(1/sqrt(d34))
  
  # random effects matrices
  Z3 <- G3 %*% Delta3  # smooth random comp. fx3
  Z4 <- G4 %*% Delta4  # smooth random comp. fx4
  Z3x4 <- rowtens(x4., G3n)  # linear:smooth comp. x4:fx3
  Z4x3 <- rowtens(G4n, x3.)  # linear:smooth comp. fx4:x3
  Z34 <- rowtens(G4n, G3n) %*% Delta34  # smooth interaction  fx3:fx4
  
  # Random effects matrix
  ZS2 <- cbind(Z3, Z4, Z3x4, Z4x3, Z34)
  
  ##################################################### third phase for x5 and x6
  MM5 <- MM.basis(x5, min(x5) - 0.01, max(x5) + 0.01, nseg, bdeg, pord)  # bases for x5
  X5 <- MM5$X
  G5 <- MM5$Z
  d5 <- MM5$d
  B5 <- MM5$B
  
  MM6 <- MM.basis(x6, min(x6) - 0.01, max(x6) + 0.01, nseg, bdeg, pord)  # bases for x6
  X6 <- MM6$X
  G6 <- MM6$Z
  d6 <- MM6$d
  B6 <- MM6$B
  
  MM5n <- MM.basis(x5, min(x5) - 0.01, max(x5) + 0.01, nseg/div, bdeg, pord)  # Nested bases for x5
  G5n <- MM5n$Z
  d5n <- MM5n$d
  B5n <- MM5n$B
  
  MM6n <- MM.basis(x6, min(x6) - 0.01, max(x6) + 0.01, nseg/div, bdeg, pord)  # Nested bases for x6
  G6n <- MM6n$Z
  d6n <- MM6n$d
  B6n <- MM6n$B
  
  c5 <- ncol(B5)
  c6 <- ncol(B6)
  
  c5n <- ncol(B5n)
  c6n <- ncol(B6n)
  
  one5. <- matrix(X5[, 1], ncol = 1)
  one6. <- matrix(X6[, 1], ncol = 1)
  
  x5. <- matrix(X5[, 2], ncol = 1)
  x6. <- matrix(X6[, 2], ncol = 1)
  
  XS3 <- rowtens(X6, X5)  # -> Fixed effects
  d56 <- c(rep(1, c6n - pord) %x% d5n + d6n %x% rep(1, c5n - pord))
  Delta5 <- diag(1/sqrt(d5))
  Delta6 <- diag(1/sqrt(d6))
  Delta56 <- diag(1/sqrt(d56))
  
  # random effects matrices
  Z5 <- G5 %*% Delta5  # smooth random comp. fx5
  Z6 <- G6 %*% Delta6  # smooth random comp. fx6
  Z5x6 <- rowtens(x6., G5n)  # linear:smooth comp. x6:fx5
  Z6x5 <- rowtens(G6n, x5.)  # linear:smooth comp. fx6:x5
  Z56 <- rowtens(G6n, G5n) %*% Delta56  # smooth interaction  fx5:fx6
  
  # Random effects matrix
  ZS3 <- cbind(Z5, Z6, Z5x6, Z6x5, Z56)
  
  ##################################################### fourth phase for x7 and x8
  MM7 <- MM.basis(x7, min(x7) - 0.01, max(x7) + 0.01, nseg, bdeg, pord)  # bases for x7
  X7 <- MM7$X
  G7 <- MM7$Z
  d7 <- MM7$d
  B7 <- MM7$B
  
  MM8 <- MM.basis(x8, min(x8) - 0.01, max(x8) + 0.01, nseg, bdeg, pord)  # bases for x8
  X8 <- MM8$X
  G8 <- MM8$Z
  d8 <- MM8$d
  B8 <- MM8$B
  
  MM7n <- MM.basis(x7, min(x7) - 0.01, max(x7) + 0.01, nseg/div, bdeg, pord)  # Nested bases for x7
  G7n <- MM7n$Z
  d7n <- MM7n$d
  B7n <- MM7n$B
  
  MM8n <- MM.basis(x8, min(x8) - 0.01, max(x8) + 0.01, nseg/div, bdeg, pord)  # Nested bases for x8
  G8n <- MM8n$Z
  d8n <- MM8n$d
  B8n <- MM8n$B
  
  c7 <- ncol(B7)
  c8 <- ncol(B8)
  
  c7n <- ncol(B7n)
  c8n <- ncol(B8n)
  
  one7. <- matrix(X7[, 1], ncol = 1)
  one8. <- matrix(X8[, 1], ncol = 1)
  
  x7. <- matrix(X7[, 2], ncol = 1)
  x8. <- matrix(X8[, 2], ncol = 1)
  
  XS4 <- rowtens(X8, X7)  # -> Fixed effects
  d78 <- c(rep(1, c8n - pord) %x% d7n + d8n %x% rep(1, c7n - pord))
  Delta7 <- diag(1/sqrt(d7))
  Delta8 <- diag(1/sqrt(d8))
  Delta78 <- diag(1/sqrt(d78))
  
  # random effects matrices
  Z7 <- G7 %*% Delta7  # smooth random comp. fx7
  Z8 <- G8 %*% Delta8  # smooth random comp. fx8
  Z7x8 <- rowtens(x8., G7n)  # linear:smooth comp. x8:fx7       
  Z8x7 <- rowtens(G8n, x7.)  # linear:smooth comp. fx8:x7
  Z78 <- rowtens(G8n, G7n) %*% Delta78  # smooth interaction  fx7:fx8
  
  # Random effects matrix
  ZS4 <- cbind(Z7, Z8, Z7x8, Z8x7, Z78)
  
  
  ## combine four phases together
  X <- cbind(XS1, XS2, XS3, XS4)
  Z <- cbind(ZS1, ZS2, ZS3, ZS4)
  
  # Compute cross-products
  XtX. <- crossprod(X)
  XtZ. <- crossprod(X, Z)
  ZtZ. <- crossprod(Z)
  Xty. <- crossprod(X, y)
  Zty. <- crossprod(Z, y)
  yty <- sum(y^2)
  
  V <- construct.block2(XtX., XtZ., ZtZ.)
  u <- c(Xty., Zty.)
  
  np <- c(ncol(XS1), ncol(XS2), ncol(XS3), ncol(XS4), ncol(Z1), ncol(Z2), ncol(Z1x2), ncol(Z2x1), ncol(Z12), ncol(Z3), ncol(Z4), ncol(Z3x4), ncol(Z4x3), ncol(Z34), ncol(Z5), ncol(Z6), ncol(Z5x6), ncol(Z6x5), ncol(Z56), ncol(Z7), ncol(Z8), ncol(Z7x8), ncol(Z8x7), ncol(Z78))
  nblocks <- 24
  
  idx <- rep(1:nblocks, np)
  end1 <- proc.time()[3]
  nn <- sum(np)
  
  # Initialize algorithm
  la <- c(0, 0, 0, 0, c(rep(1, nblocks - 4 )))
  
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
    lanew <- c(0, 0, 0, 0, rep(sig2, nblocks - 4))/tau
    dla <- mean(abs(log10(la[5:nblocks]) - log10(lanew[5:nblocks])))
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