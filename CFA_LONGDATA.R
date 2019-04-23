# CFA Model Longitudional data
#===========================================================>
# Author: Clara-Christina Gerstner
# Purpose: CFA model for longitudional data with constraint factor loadings.
#===========================================================>

#-------------------------------------------------------->
# STEP I: Set correlation matrix for S
source("long_data.R")
S <-  longFacCorMatFull

#-------------------------------------------------------->
# STEP II: Load functions

# Estimate Sigma (theta) and plot output
# Input: theta 
# Output: Sigma

CFA_Sigma_L <- function(theta, S, C){
  NS <- ncol(S)
  NC <- length(C)
  t <- NS/NC
  lam <- theta[1:t]
  delta <- theta[(t+1):(t+NS)]
  phi <- theta[(t+NS+1):(length(theta))]
  lammat <- matrix(0, ncol=NC, nrow=NS)
  for (k in 1:t){
    m <- k
    for (i in 1:(NS/t)){
      lammat[m,i] <- lam[k]
      m <- m+3
    }
  }
  phimat <- diag(NC)
  cnt=0
  for (j in 1:(NC-1)){
    for (i in (j+1):NC){
      cnt=cnt+1
      phimat[i,j] <- phi[cnt]; phimat[j,i] = phimat[i,j]
    }
  }
  deltamat <- diag(delta)
  Sig <- lammat%*%phimat%*%t(lammat)+deltamat  
  return(Sig)
}

# Calculate CFA log Function 
CFA_loglik <- function(theta, S, N){
  Sig <- CFA_Sigma_L(theta,S,C)
  (((N-1)/2) * (log(abs(det(Sig))) + sum(diag(S%*%solve(Sig)))))  
}
2*CFA_loglik(ML_output$par, S, 126)

# Estimation method I: Maximum likelihood
# Input: starting values theta, S
# Output: theta optimized
F_ML <- function(theta, S){
  Sig <- CFA_Sigma_L(theta,S,C)
  log(det(Sig)) + sum(diag(S%*%solve(Sig)))
}

# Estimation method II: Unweighted results
# Input: starting values theta, S
# Output: theta optimized
F_ULS <- function(theta, S){
  Sig <- CFA_Sigma_L(theta,S,C)
  Z <- (S-Sig)
  0.5*sum(diag(Z%*%Z))
}

# Output results
CFA_output_L <- function(theta,N){
  NS <- ncol(S)
  NC <- length(C)
  t <- NS/NC
  lam <- theta[1:t]
  delta <- theta[(t+1):(t+NS)]
  phi <- theta[(t+NS+1):(length(theta))]
  lammat <- matrix(0, ncol=NC, nrow=NS)
  for (k in 1:t){
    m <- k
    for (i in 1:(NS/t)){
      lammat[m,i] <- lam[k]
      m <- m+3
    }
  }
  phimat <- diag(NC)
  cnt=0
  for (j in 1:(NC-1)){
    for (i in (j+1):NC){
      cnt=cnt+1
      phimat[i,j] <- phi[cnt]; phimat[j,i] = phimat[i,j]
    }
  }
  deltamat <- diag(delta)
  loglik <- 2*CFA_loglik(theta, S, N)
  return(list(theta=theta,lammat=lammat, phimat=phimat, deltamat=deltamat, loglik=loglik))
}

#-------------------------------------------------------->
# STEP III: Get estimates for factor loading, factor correlations, residual correlations

# Different factor structures
C <- c(3,3,3,3) # 3 factors measured over 4 occasions

# Generate starting values
# add residual correlations (in delta)
# constrain lambda correlations
theta_fun <- function(S,C){
  th <- c(rep(0.8, (ncol(S)/length(C))), rep(0.5, (ncol(S))), rep(0.4,((length(C)^2-length(C))/2)))
}

# try different staring values
theta <- theta_fun(S,C)

# Run ML and print output
ML_output <- optim(par=theta, fn=F_ULS, S=S,control = list(maxit=50000))

CFA_output <- CFA_output_L(ML_output$par,126)
CFA_output

# Compare model and Saturated model
Sat_model <- (126-1)*(log(det(S)) + ncol(S))

q=abs(Sat_model-CFA_output$loglik)
df=78-length(theta)      # number of free parameters
tail_prob <- 1-(pchisq(q, df, ncp = 0, lower.tail = TRUE, log.p = FALSE))
tail_prob 