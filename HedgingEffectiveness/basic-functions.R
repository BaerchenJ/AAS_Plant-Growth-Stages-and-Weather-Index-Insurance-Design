library(splines)
library(fields)

# Compute B-spline base matrix
bspline <- function(X., XL., XR., NDX., BDEG.){
  dx <- (XR. - XL.)/NDX.
  knots <- seq(XL. - BDEG. * dx, XR. + BDEG. * dx, by = dx)
  B <- spline.des(knots, X., BDEG. + 1, 0 * X.)$design
  B
}

# row-tensor product (or Box product of two matrices)
rowtens <- function(X, B){
  one.1 <- matrix(1, 1, ncol(X))
  one.2 <- matrix(1, 1, ncol(B))
  kronecker(X, one.2) * kronecker(one.1, B)
}

# Mixed Model Basis
MM.basis <- function(x, xl, xr, ndx, bdeg, pord){
  
  B <- bspline(x, xl, xr, ndx, bdeg)
  m <- ncol(B)
  
  D <- diff(diag(m), differences = pord)
  
  P.svd <- svd(t(D) %*% D)
  U <- (P.svd$u)[, 1:(m - pord)]
  d <- (P.svd$d)[1:(m - pord)]
  
  Z <- B %*% U
  
  X <- NULL
  for (i in 0:(pord - 1)){
    X <- cbind(X, x^i)
  }
  list(X = X, Z = Z, d = d, B = B, m = m, D = D)
}

# Construct 2 x 2 block symmetric matrices:
construct.block2 <- function(A1, A2, A4){
  block <- rbind(cbind(A1, A2), cbind(t(A2), A4))
  return(block)
}

# function to convert fitted values back to original scale
back <- function(x, y){
  z = x * (max(y)-min(y)) + min(y)
  return(z)
}
