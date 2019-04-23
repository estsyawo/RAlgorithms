# QR decomposition 
#===========================================================>
# Author: Clara-Christina Gerstner
# Purpose: To decompose matrix into an orthogonal matrix Q 
# and an upper triangular matrix R
#===========================================================>

# Step-by-step example

require(Matrix)

# Example 1
A <- matrix(c(12,6,-4,-51,167,24,4,-68,-41),nrow = 3)
A

# Example 2
A <- matrix(c(2,2,1,-2,1,2,18,0,0),nrow = 3)
A

#===========================================================>

# Run one Q at a time

H <- list()
m <- nrow(A)

# get Q1

k=1
a1 <- A[k:m,k]
e1 <- as.matrix(c(1, rep(0, length(a1)-1)))
v <- sign(a1[1]) * sqrt(sum(a1^2)) * e1 + a1
Q1 <- diag(length(a1)) - 2*as.vector((v)%*%t(v)) / (t(v) %*% (v))

H[[k]] <- Q1

A <- Q1%*%A
round(A)

# get Q2

k=2
a1 <- A[k:m,k]
e1 <- as.matrix(c(1, rep(0, length(a1)-1)))
v <- sign(a1[1]) * sqrt(sum(a1^2)) * e1 + a1
Q2 <- diag(length(a1)) - 2*as.vector((v)%*%t(v)) / (t(v) %*% (v))

Q2 <- rbind(c(0,0),Q2)
Q2 <- cbind(c(1,0,0),Q2)

H[[k]] <- Q2

A <- Q2%*%A
round(A)

# get Q3

k=3
a1 <- A[k:m,k]
e1 <- as.matrix(c(1, rep(0, length(a1)-1)))
v <- sign(a1[1]) * sqrt(sum(a1^2)) * e1 + a1
Q3 <- diag(length(a1)) - 2*as.vector((v)%*%t(v)) / (t(v) %*% (v))

Q3 <- rbind(c(0),c(0),Q3)
Q3 <- cbind(c(1,0,0),c(0,1,0),Q3)

H[[k]] <- Q3

A <- Q3%*%A
round(A)

# get Results

Q <- Q1 %*% Q2 %*% Q3 

R <- A

round(R)

#===========================================================>

# Algorithm

qr_fx <- function(A) {
  A2 <- as.matrix(A)
  
  m <- nrow(A)
  H <- list() 
  
  for (k in 1:m) {
    a <- A2[k:m,k]
    e <- as.matrix(c(1, rep(0, length(a)-1)))
    v <- sign(a[1]) * sqrt(sum(a^2)) * e + a
    
    # Compute Q matrix
    Q <- diag(length(a)) - 2 * as.vector(v %*% t(v)) / (t(v) %*% v)
    if (k > 1) {
      Q <- bdiag(diag(k-1), Q)
    }
    
    # Store result
    H[[k]] <- Q
    
    A2 <- Q %*% A2
  }
  
  Q1 <- Reduce("%*%", H) 
  Result <- list('Q'=round(Q1),'R'=round(A2))
  return(Result)
}

qr_fx(A)
