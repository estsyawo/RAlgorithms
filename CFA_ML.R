# CFA Model Maximum likelihood
#===========================================================>
# Author: Clara-Christina Gerstner
# Purpose: Obtain lambda values for confirmatory factor model
# using maximum likelihood with optim. Apply approach to TIMSS 
# data set motivation variables.
#===========================================================>

#-------------------------------------------------------->
# STEP I: Create Covariance matrix from TIMSS

Data <- cbind(Dat$BSBS23A, Dat$BSBS23B, Dat$BSBS23C,Dat$BSBS23D, Dat$BSBS23E, Dat$BSBS23F, Dat$BSBS23G, Dat$BSBS23H, 
              Dat$BSBS24A, Dat$BSBS24B, Dat$BSBS24C, Dat$BSBS24D, Dat$BSBS24E, Dat$BSBS24F, Dat$BSBS24G, Dat$BSBS24H, Dat$BSBS24I,
              Dat$BSBS21A, Dat$BSBS21B, Dat$BSBS21C, Dat$BSBS21D, Dat$BSBS21E, Dat$BSBS21F, Dat$BSBS21G, Dat$BSBS21H, Dat$BSBS21I)

# Deal with missing data
rowi <- c(0); k=0
for (i in 1:nrow(Data)){
  if (any(is.na(Data[i,]))){
    k=k+1
    rowi[k]=i
  }
}

# create covariance matrix S
S <- cov(Data[-rowi,])

# Alternative: load S matrix
# 9 variables, 3 constructs with 3 variables each (Intrinsic motivation, utility motivation, self concept)
S <- matrix(c(0.5814956, 0.4246274, 0.4235168 ,0.2382251, 0.2186985, 0.2024063, 0.3624870, 0.2508839, 0.3714332,
            0.4246274 ,0.7198141 ,0.5449004, 0.2773675 ,0.2737452, 0.2376515 ,0.4361281 ,0.3233333 ,0.4587497,
            0.4235168, 0.5449004 ,0.8090091 ,0.2935289 ,0.3119527, 0.2564495 ,0.4329468 ,0.3101380 ,0.4578905,
            0.2382251 ,0.2773675 ,0.2935289 ,0.7162658 ,0.5209469, 0.4270696 ,0.3606683 ,0.3077300 ,0.3762689,
            0.2186985 ,0.2737452 ,0.3119527 ,0.5209469 ,0.8719185, 0.4819982 ,0.3330789 ,0.2650756 ,0.3481156,
            0.2024063 ,0.2376515 ,0.2564495 ,0.4270696 ,0.4819982, 0.7812140 ,0.2874521 ,0.2491076 ,0.3126798,
            0.3624870 ,0.4361281 ,0.4329468 ,0.3606683 ,0.3330789, 0.2874521 ,0.7954722 ,0.4868525 ,0.6780372,
            0.2508839 ,0.3233333 ,0.3101380 ,0.3077300 ,0.2650756, 0.2491076 ,0.4868525 ,0.6677226 ,0.5296753,
            0.3714332 ,0.4587497, 0.4578905 ,0.3762689, 0.3481156, 0.3126798 ,0.6780372 ,0.5296753, 0.8555101),ncol=9, nrow=9)

#-------------------------------------------------------->
# STEP II: Load functions

# Estimate Sigma (theta) and plot output
# Input: theta
# Output: Sigma
CFA_Sigma <- function(theta){
  lam <- theta[1:9]
  phi <- theta[10:12]
  delta <- theta[13:21]
  lammat <- matrix(0, ncol=3, nrow=9)
  lammat[1:3,1] <- lam[1:3]
  lammat[4:6,2] <- lam[4:6]
  lammat[7:9,3] <- lam[7:9]
  phimat <- diag(3)
  phimat[2,1] <- phi[1]; phimat[1,2] = phimat[2,1]
  phimat[3,1] <- phi[2]; phimat[1,3] = phimat[3,1]
  phimat[3,2] <- phi[3]; phimat[2,3] = phimat[3,2]
  deltamat <- diag(delta)
  Sig <- lammat%*%phimat%*%t(lammat)+deltamat  
  return(Sig)
}

# Output function
# Input: theta
# Output: lambda matrix, phi matrix, delta matrix
CFA_output <- function(theta){
  lam <- theta[1:9]
  phi <- theta[10:12]
  delta <- theta[13:21]
  lammat <- matrix(0, ncol=3, nrow=9)
  lammat[1:3,1] <- lam[1:3]
  lammat[4:6,2] <- lam[4:6]
  lammat[7:9,3] <- lam[7:9]
  phimat <- diag(3)
  phimat[2,1] <- phi[1]; phimat[1,2] = phimat[2,1]
  phimat[3,1] <- phi[2]; phimat[1,3] = phimat[3,1]
  phimat[3,2] <- phi[3]; phimat[2,3] = phimat[3,2]
  deltamat <- diag(delta)
  return(list(lammat=lammat, phimat=phimat, deltamat=deltamat))
}

# Maximum likelihood function
# Input: theta and S
CFA_ML <- function(theta, S){
  Sig <- CFA_Sigma(theta)
  log(det(Sig)) + sum(diag(S%*%solve(Sig)))
}

#-------------------------------------------------------->
# STEP III: Run functions

# Try different Starting values
th1 <- c(0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.6,0.6,0.6,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25)
th2 <- c(0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.1,0.1,0.1,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5)
th3 <- runif(21)

# Run ML and print output
ML_output <- optim(par=th1, fn=CFA_ML, S=S)

CFA_output(ML_output$par)

#-------------------------------------------------------->
# STEP IV: get loglike values

# Calculate CFA log Function 1
CFA_loglik <- function(theta, S){
  Sig <- CFA_Sigma(theta)
  log(det(Sig)) + sum(diag(S%*%solve(Sig))) - log(det(S)) - nrow(S) 
}
CFA_loglik(th1, S)

# Calculate CFA log Function 2
CFA_loglik2 <- function(theta, S, N){
  Sig <- CFA_Sigma(theta)
  (((N-1)/2) * (log(det(Sig)) + sum(diag(S%*%solve(Sig))))) 
}
2*CFA_loglik2(th1, S, 8133) # number of observations

# Finding: th1 provides the best log likelihood values

#-------------------------------------------------------->
# Compare to LISREL output:

# includes all three outpu matrices
th_LISREL <- c(0.5799,0.7344, 0.7373,0.7032, 0.7425, 0.6242, 0.7963, 0.6167, 0.8520,
               0.5421, 0.7372, 0.6036,0.2452,0.1805 ,0.2654 ,0.2218 , 0.3206, 0.3916, 
               0.1613,0.2875, 0.1295)

CFA_loglik(th_LISREL, S) 
# smaller than optim result

2*CFA_loglik2(th_LISREL, S, 8133)
# smaller than optim result

# Conclusion: LISREL provides more accurate estimates, 
# but our function gives reasonable results