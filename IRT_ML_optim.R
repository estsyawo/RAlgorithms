# IRT Rasch Model Maximum likelihood
#===========================================================>
# Author: Clara-Christina Gerstner
# Purpose: Obtain theta value and probability for theta in a
# one-factor Rasch model using maximum likelihood with optim.
#===========================================================>

#-------------------------------------------------------->
# STEP I: Load sample data 
# Test 2 Data: 25 items, N= 660

TEST2 <- read.csv("TEST2.csv", header = TRUE)

#-------------------------------------------------------->
# STEP II: Run functions

# Maximum likelihood function
# Input: Test data, idP = person ID
# Output: Theta
Rasch_ML <- function(theta,TESTdat, idP){
  x <- TESTdat[idP,]
  t <- sum(x)
  b <- apply(TESTdat,2,mean)
  -(t*theta - sum(x*b) - sum(log(1+exp(theta-b))))
}

# Optim function
theta_obj <- optimize(Rasch_ML, c(-10,10), TESTdat=TEST2, id=5)

# Probability function
# Input: Testdata, idP = person id, idI = item id
# Output: probability for theta
P_theta <- function(theta,TESTdat, idP, idI){
  x <- TESTdat[idP,]
  t <- sum(x)
  b <- apply(TESTdat,2,mean)
  exp(theta-b[idI])/ (1+ exp(theta-b[idI]))
}

#-------------------------------------------------------->
# STEP III: Test function

(theta_out <- theta_obj[[1]])
P_theta(theta_obj$minimum, TEST2, 1, 4)
