# Classic True Score Theory

#===========================================================>
# Author: Clara-Christina Gerstner
# Purpose: Summary of functions relevant to classic test theory.
#===========================================================>

source("Basic_functions.R")
library(psych)

# Estimated True Score 

Est_True <- function(x){
  r <- pear_fx(x,x)
  T <- rep(0,length(x)) 
  for (i in 1:length(x)){
  T[i] <- r*(x[i]-mean_fx(x)) + mean_fx(x)
  }
  return(T)
}

# Example
Est_True(X)

#===========================================================>

# Deviate

dev_fx <- function(x, y){
  sum_fx(((x-(mean_fx(x)))*(y-(mean_fx(y)))))
}

# Example
dev_fx(X, X)

#===========================================================>

# Pearson Product Moment

pear_fx <- function(x,y){
  dev_fx(x, y) / 
  (sqrt(dev_fx(x, x) * dev_fx(y, y)))
}

# Example 1
X <- c(0,0,0,0,0,1,1,1,1,1,1,1,1,1,1)   # Binary
Y <- c(1,4,0,5,7,4,5,7,9,8,6,9,7,6,5)   # Continuous

pear_fx(X,Y)

all.equal(pear_fx(X,Y),cor(X,Y, method = "pearson"))

# Example 2: Pearson correlation table

A <- TEST2[,1:5]

cor(A,A, method = "pearson")

#===========================================================>

# Covariance

cov_fx <- function(x,y){
  (dev_fx(x,y)) / (length(x)-1)
}

# Example
cov_fx(X,Y)

all.equal(cov_fx(X,Y),cov(X,Y, method = "pearson"))

#===========================================================>

# Standard Error of Measurement

sem_fx <- function(x,x){
  sem <- SD_fx(x)*sqrt(1-pear_fx(x,x))
    return(sem)
}

# Example
sem_fx(X)

#===========================================================>

# Confidence limits for SEM

sem_ci <- function(x,z){
  S <-matrix(0:0,ncol=3,nrow=length(x)) 
  sem_xi <- abs(sem_fx(x))
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
sem_ci(X,1.96)

#===========================================================>

# Standard Error of Estimate

see_fx <- function(x){
  see <- pear_fx(x,x) * sem_fx(x)
  return(see)
}

# Example
see_fx(X)

#===========================================================>

# Confidence limits for SEE

see_ci <- function(x,z){
  S <-matrix(0:0,ncol=3,nrow=length(x)) 
  tx <- Est_True(x)
  see_tx <- abs(see_fx(x))
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
see_ci(X,1.96)

#===========================================================>

# Biserial correlation

bis_fx <- function(b){
  n <- nrow(b)
  mean_y1 <- mean_fx(B[,2][B[,1]==1])
  mean_y <- mean_fx(b[,2])
  sd_y <- SD_fx(b[,2])
  px <- mean_fx(b[,1])
  fz <- dnorm(qnorm(1-px))
  ((mean_y1 - mean_y)/sd_y)*(px/fz)
}

# Example
X <- c(0,0,0,0,0,1,1,1,1,1,1,1,1,1,1)   # Binary
Y <- c(1,4,0,5,7,4,5,7,9,8,6,9,7,6,5)   # Continuous
B <- as.matrix(cbind(X,Y))

bis_fx(B)

#===========================================================>

# Point Biserial correlation

pbis_fx <- function(b){
  n <- nrow(b)
  mean_y1 <- mean_fx(B[,2][B[,1]==1])
  mean_y <- mean_fx(b[,2])
  sd_y <- SD_fx(b[,2])
  px <- mean_fx(b[,1])
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

alpha_fx <- function(x){
  n <- ncol(x)
  (n/(n-1)) * (1-(sum_fx((Var_fx(x))) / Var_fx((apply(x,1,sum_fx))) ))
}

# Example
TEST2 <- read.csv("TEST2.csv", header = TRUE)

alpha_fx(TEST2)
alpha_M <- cronbach(TEST2)
rawalpha_M <- alpha_M[[1]]

all.equal(alpha_fx(TEST2), rawalpha_M)

#===========================================================>

# KR-20

kr20_fx <- function(x){
  n <- nrow(x)
  p <- mean_fx(x)
  (n/(n-1)) * (1-(sum_fx((p*(1-p))) / Var_fx((apply(x,1,sum_fx))) ))
}

# Example
kr20_fx(TEST2)

#===========================================================>

# Correction for attenuation

corr_fx <- function(x,y){
  pear_fx(x,y) / sqrt((pear_fx(x,x)*pear_fx(y,y)))
}

# Example
X <- c(runif(10,1,10))
Y <- c(runif(10,1,10))

corr_fx(X,Y)
pear_fx(X,Y)

#===========================================================>

# Spearman-Brown Formula

spearman_fx <- function(x,y,a){
  n <- length(x)
  f <- (n+a)/n
  (f*pear_fx(x,y)) / (1+(f-1)*pear_fx(x,y))
}

# Example
X <- c(runif(10,1,10))
Y <- c(runif(10,1,10))

spearman_fx(X,Y, 5) # add 5 items to test
pear_fx(X,Y)

#===========================================================>

# Cohen's kappa

cohen_fx <- function(Con){
  a <- Con[1,1]
  b <- Con[1,2]
  c <- Con[2,1]
  d <- Con[2,2]
  n <- a+b+c+d
  p0 <- (a+d)/n
  pc <- ((((a+c)*(a+b))/n)+(((b+d)*(c+d))/n))/n
  ((p0 - pc)/(1 - pc))
}

# Example Allen & Yen p.37
Con_table <- matrix(c(22,4,2,11),2,2)  

cohen_fx(Con_table)
