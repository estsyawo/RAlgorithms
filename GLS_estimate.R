# GLS estimate 
#===========================================================>
# Author: Clara-Christina Gerstner
# Purpose: Use lower level functions to estimate parameters
# using GLS.
#===========================================================>

source("Matrix_Multip.R")
source("Conj_grad.R")
source("Chol_decom.R")

# Load data
data(longley)

# Get OLS estimates
lm_estimate <- lm(Employed~GNP+Population, data = longley)
summary(lm_estimate)
# OLS estimates do not take into account error correlation

# Using nlme package to get GLS estimates and error correlation
library(nlme)
gls_solution <- gls(Employed~GNP+Population,
                    correlation=corAR1(form=~Year),data=longley)
gls_summary <- summary(gls_solution)
(gls_estimate <- as.matrix(gls_solution[[4]]))
(st.errors <- sqrt(diag(gls_solution$varBeta)))

# correlated Errors by year (Phi in gls_solution output)
error <- 0.6441692

#===========================================================>

# Obtain GLS estimates using the closed form solution

# Covariance matrix V of error
V <- diag(length(longley$Employed))
V <- error^abs(row(V)-col(V))

# calculate generalized least square estimates
# beta = (X'V-1X)-1 X'V-1y
X <- model.matrix(lm_estimate)
Y <- as.matrix(longley$Employed)
V.inv <- solve(V)

(gls_estimate1 <- solve(t(X) %*% V.inv %*% X)%*% t(X)%*% V.inv %*% Y)
all.equal(gls_estimate, gls_estimate1)

#===========================================================>

# Alternative: obtain GLS estimate using Cholesky decomposition 
# and the conjugate gradient

K <- Chol_fx(V)
K.inv <- solve(t(K))

all.equal(solve(V), (K.inv)%*%t(K.inv))

# Mutliply X and Y by lower triangle of the cholesky decomp
B <- (K.inv)%*%X
Z <- (K.inv)%*%Y

# Obtain B'B and B'Z
A <- fxtx(B)
b <- fxty(B,Z)

# Solve set of linear equations
(gls_estimate2 <- Conjgrad_fx(A,b))
# results will be slightly different due to numerical rounding errors