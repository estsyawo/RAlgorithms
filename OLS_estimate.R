# OLS estimate 
#===========================================================>
# Author: Clara-Christina Gerstner
# Purpose: Use lower level functions to estimate parameters
# using OLS.
#===========================================================>

source("Matrix_Multip.R")
source("Conj_grad.R")
source("Chol_decom.R")

# Load data
X1 <- as.matrix(women$height)
X0 <- cbind(1,X1)                
Y1 <- as.matrix(women$weight)

# Obtain X'X and X'Y
A2 <- fxtx(X0)
b2 <- fxty(X0,Y1)

# Solve system of linear equations
Conjgrad_fx(A2,b2)

# Compare results to lm function
lm_solution <- summary(lm(Y1~X1))

lm_solution <- as.matrix(lm_solution$coefficients[1:2])

all.equal(Conjgrad_fx(A2,b2), lm_solution)