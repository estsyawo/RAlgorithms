# Cholesky decomposition 
#===========================================================>
# Author: Clara-Christina Gerstner
# Purpose: To take the product of a lower triangular matrix 
# and its conjugate transpose.
#===========================================================>

# Algorithm

Chol_fx <- function(A){
  L <-matrix(0:0,ncol(A),ncol(A)) 
  L[1,1] = sqrt(A[1,1])
  for (k in 2:ncol(A)){
    L[k,1:(k-1)] = solve(L[1:(k-1),1:(k-1)],A[1:(k-1),k])
    L[k,k] = sqrt(A[k,k]-t(L[k,1:(k-1)])%*%(L[k,1:(k-1)]))
  }
  return(L)
}

#===========================================================>

# Example

A <- matrix(round(runif(9)*10), 3)

Chol_solution <- t(chol(A))

Chol_fx(A)

all.equal(Chol_fx(A),Chol_solution)