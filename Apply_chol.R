# Appy Chol decomposition
#===========================================================>
# Author: Clara-Christina Gerstner
# Purpose: Solving system of equations using Cholesky decomposition
#===========================================================>

# Algorithm

forward_fx <- function(L,b){
  x <- rep(NA, length(b))
  x[1] <- b[1]/L[1,1]
for (i in 2: length(b)){
  x[i] <- b[i]
    for (j in 1:(i-1)){
    x[i] <- x[i] - L[i,j]*x[j]
    }
  x[i] <- x[i]/L[i,i]
}
  return(x)
}

backward_fx <- function(U,b){
  n <- length(b)
  x <- rep(NA, n)
  x[n] <- b[n]/U[n,n]
  for (i in (n-1):1){
    x[i] <- b[i]
    for (j in (i+1):n){
      x[i] <- x[i] - U[i,j]*x[j]
    }
    x[i] <- x[i]/U[i,i]
  }
  return(x)
}

#===========================================================>

# Example 1

L <- matrix(c(2,-1,3,0,3,0.5,0,0,-1),nrow = 3)
U <- matrix(c(3,0,0,1,3,0,-1,-1,1),nrow = 3)
b1 <- c(1,2,3)
b2 <- c(3,2,1)

forward_fx(L,b1)
backward_fx(U,b2)

# Example 2

AL <- matrix(c(7,3,0,7), nrow = 2)
AU <- matrix(c(7,0,5,7), nrow = 2)
b1 <- c(1,2)

forward_fx(AL,b1)
backward_fx(AU,b1)

# Compare to Conjgrad (with rounding error)

A <- matrix(c(7,3,5,7), 2)
b <- c(1,2)

Conjgrad_fx(A,b)
