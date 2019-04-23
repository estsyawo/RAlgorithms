# Conjugate gradient
#===========================================================>
# Author: Clara-Christina Gerstner
# Purpose: To implement the conjugate gradient method for
# solving a system of linear equations.
#===========================================================>

# Algorithm

Conjgrad_fx <- function(A,b,x = rep(0, ncol(A))){
  r = b - A %*% x
  p = r
  k = 0
  dpr = t(r) %*% r
  
  dev=10
  tol=1e-10
  maxk = 1000
  
  while(dev > tol & k <= maxk){
    Ap = A %*% p
    alpha = t(r) %*% r / (t(p) %*% Ap)
    x = x + as.numeric(alpha) * p
    r = r - as.numeric(alpha) * Ap
    dev = max(abs(r))
    beta = t(r) %*% r/ dpr
    p = r + as.numeric(beta) * p
    k=k+1
  }
  return(x)
}

#===========================================================>

# Example 

A <- matrix(c(7,3,5,7), 2)
b <- c(1,2)

Conjgrad_fx(A,b)
