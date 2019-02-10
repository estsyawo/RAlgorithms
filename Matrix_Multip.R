# Matrix multiplication 
#===========================================================>
# Author: Clara-Christina Gerstner
# Purpose: Matrix product manipulations without explicit 
# transpose.
#===========================================================>

### XtX ###

# Algorithm: input: x; output: x'x 

fxtx = function(X){
  XX <-matrix(0:0,ncol(X),ncol(X)) 
  for(i in 1:ncol(X)){
    for(j in 1:i){
      for (k in 1:nrow(X)){
        XX[i,j]<-XX[i,j]+X[k,i]*X[k,j] 
      }
      if(i!=j){
        XX[j,i]=XX[i,j]} 
    }
  }
  return(XX)
}

#===========================================================>

# Example

X <- matrix(runif(9), 3); Y = matrix(rnorm(9), 3)

XtX<- t(X) %*% X

all.equal(fxtx(X),XtX)

#===========================================================>

### XtY ###

# Algorithm:  input: x,y  output: x'y

fxty <- function(X,Y){
  XY <-matrix(0:0,ncol(X),ncol(Y)) 
  for(i in 1:ncol(X)){
  for(j in 1:ncol(Y)){
    for (k in 1:nrow(Y)){
      XY[i,j]<-XY[i,j]+X[k,i]*Y[k,j] 
    }
  }
  }
  return(XY)
}

#===========================================================>

# Example

XtY <- t(X) %*% Y

all.equal(fxty(X, Y),XtY)

#===========================================================>

### XYt ###

# Algorithm:  input: x,y  output: xy'

fxyt <- function(X,Y){
XY2 <-matrix(0:0,ncol(X),ncol(Y)) 
for(i in 1:ncol(X)){
  for(j in 1:nrow(X)){
    for (k in 1:nrow(Y)){
      XY2[i,j]<-XY2[i,j]+X[i,k]*Y[j,k] 
    }
  }
}
return(XY2)
}

#===========================================================>

# Example

XYt <- X %*% t(Y)

all.equal(fxyt(X, Y),XYt)