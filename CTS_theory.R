# Classic True Score Theory

#===========================================================>
# Author: Clara-Christina Gerstner
# Purpose: Summary of functions relevant to classic test theory.
#===========================================================>

source("Basic_functions.R")

# Estimated True Score 

Est_True <- function(x,n){
  mean_x <- mean_fx(x,n)
  r <- pear_fx(x,x,n)
  T <- rep(0,length(x)) 
  for (i in 1:length(x)){
  T[i] <- r*(x[i]-mean_x) + mean_x
  }
  return(T)
}

# Example

Est_True(X, length(X))

#===========================================================>

# Deviate

dev_fx <- function(x, y, n){
  mean_x <- mean_fx(x, n)
  mean_y <- mean_fx(y, n)
  sum_fx(((x-(mean_x))*(y-(mean_y))),n)
}

# Example
dev_fx(X, X, length(X))

#===========================================================>

# Pearson Product Moment

pear_fx <- function(x,y,n){
mean_x <- mean_fx(x, n)
mean_y <- mean_fx(y, n)
  dev_fx(X, Y, n) / 
  (sqrt(dev_fx(X, X, length(X)) * dev_fx(Y, Y, length(X))))
}

# Example

X <- c(runif(10,1,10))
Y <- c(runif(10,1,10))

pear_fx(X,Y, length(X))

all.equal(pear_fx(X,Y, length(X)),cor(X,Y, method = "pearson"))

#===========================================================>

# Covariance

cov_fx <- function(x,y,n){
  (dev_fx(x,y,n)) / (n-1)
}

# Example

cov_fx(X,Y,length(X))

all.equal(cov_fx(X,Y,length(X)),cov(X,Y, method = "pearson"))

#===========================================================>

# Standard Error of Measurement

sem_fx <- function(x,n){
  sem <- SD_fx(x,n)*sqrt(1-pear_fx(x,x,n))
    return(sem)
}

# Example

sem_fx(X, length(X))

#===========================================================>

# Confidence limits for SEM

sem_ci <- function(x,z,n){
  S <-matrix(0:0,ncol=3,nrow=length(x)) 
  sem_xi <- abs(sem_fx(x,n))
  for (i in 1:length(x)){
    xi <- x[i]
    ul <- xi + z*sem_xi 
    ll <- xi - z*sem_xi
    S[i,] <- c(ll,xi,ul)
  }
  colnames(S) <- c("Lower limit", "X","Upper limit")
  return(S)
}

# Example

sem_ci(X,1.96,length(X))

#===========================================================>

# Standard Error of Estimate

see_fx <- function(x,n){
  see <- pear_fx(x,x,n) * sem_fx(x,n)
  return(see)
}

# Example

see_fx(X, length(X))

#===========================================================>

# Confidence limits for SEE

see_ci <- function(x,z,n){
  S <-matrix(0:0,ncol=3,nrow=length(x)) 
  tx <- Est_True(x,n)
  see_tx <- abs(see_fx(x,n))
  for (i in 1:length(x)){
    txi <- tx[i]
    ul <- txi + z*see_tx 
    ll <- txi - z*see_tx
    S[i,] <- c(ll,txi,ul)
  }
  colnames(S) <- c("Lower limit", "Est.True","Upper limit")
  return(S)
}

# Example

see_ci(X,1.96,length(X))

#===========================================================>

# Biserial correlation

bis_fx <- function(b){
  n <- nrow(b)
  mean_y1 <- mean_fx(B[,2][B[,1]==1],length(B[,2][B[,1]==1]))
  mean_y <- mean_fx(b[,2],n)
  sd_y <- SD_fx(b[,2],n)
  px <- mean_fx(b[,1],n)
  fz <- dnorm(qnorm(1-px))
  ((mean_y1 - mean_y)/sd_y)*(px/fz)
}

# Example
biserial(Y,X)
bis_fx(B)

#===========================================================>

# Point Biserial correlation

pbis_fx <- function(b){
  n <- nrow(b)
  mean_y1 <- mean_fx(B[,2][B[,1]==1],length(B[,2][B[,1]==1]))
  mean_y <- mean_fx(b[,2],n)
  sd_y <- SD_fx(b[,2],n)
  px <- mean_fx(b[,1],n)
  ((mean_y1 - mean_y)/sd_y)*sqrt(px/(1-px))
}

# Example
X <- c(0,0,0,0,0,1,1,1,1,1,1,1,1,1,1)   # Binary
Y <- c(1,4,0,5,7,4,5,7,9,8,6,9,7,6,5)   # Continuous
B <- as.matrix(cbind(X,Y))

pbis_fx(B)

#===========================================================>

# Phi correlation

phi_fx <- function(Con){
  a <- Con[1,1]
  b <- Con[1,2]
  c <- Con[2,1]
  d <- Con[2,2]
  n <- a+b+c+d
  pc <- d/n
  px <- (b+d)/n
  py <- (c+d)/n
  (pc - (px*py))/sqrt(px*(1-px)*py*(1-py))
}

# Example Allen & Yen p.37
Con_table <- matrix(c(8,12,4,16),2,2)  

all.equal(phi(Con_table),round(phi_fx(Con_table),2))

#===========================================================>

# Coefficient alpha

alpha_fx <- function(x,n){
  (n/(n-1)) * (1-(sum_fx((Var_fx(x,n)),n) / Var_fx((apply(x,1,sum_fx,n=n)),n) ))
}

# Example

M <- matrix(runif(100,1,10),10,10)

library(psych)
alpha_fx(M,nrow(M))
alpha_M <- alpha(M)
rawalpha_M <- alpha_M[[1]]$raw_alpha

all.equal(alpha_fx(M,nrow(M)), rawalpha_M)

#===========================================================>

# KR-20

kr20_fx <- function(x,n){
  p <- mean_fx(x,n)
  (n/(n-1)) * (1-(sum_fx((p*(1-p)),n) / Var_fx((apply(x,1,sum_fx,n=n)),n) ))
}

# Example

M <- matrix(round(runif(100,0,1)),10,10)
kr20_fx(M,nrow(M))

#===========================================================>