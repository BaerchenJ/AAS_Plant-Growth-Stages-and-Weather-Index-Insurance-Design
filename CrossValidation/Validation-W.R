##################################################################### Test model

NEWXZ <- function(x1, x2, nseg) {
  bdeg = 3
  div = 3
  pord = 2
  # Build Mixed Model Bases
  
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
  
  #####################
  X <- rowtens(X2, X1)  # -> Fixed effects
  d3 <- c(rep(1, c2n - pord) %x% d1n + d2n %x% rep(1, c1n - pord))
  Delta1 <- diag(1/sqrt(d1))
  Delta2 <- diag(1/sqrt(d2))
  Delta3 <- diag(1/sqrt(d3))
  
  # random effects matrices
  Z1 <- G1 %*% Delta1  # smooth random comp. fx1
  Z2 <- G2 %*% Delta2  # smooth random comp. fx2
  Z1x2 <- rowtens(x2., G1n)  # linear:smooth comp. x2:fx1
  Z2x1 <- rowtens(G2n, x1.)  # linear:smooth comp. fx2:x1
  Z12 <- rowtens(G2n, G1n) %*% Delta3  # smooth interaction  fx1:fx2
  #####################
  
  # Random effects matrix
  Z <- cbind(Z1, Z2, Z1x2, Z2x1, Z12)
  
  M <- cbind(X, Z)
  
  M
}

