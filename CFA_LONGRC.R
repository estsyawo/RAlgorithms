# CFA Model Longitudional data
#===========================================================>
# Author: Clara-Christina Gerstner
# Purpose: CFA model for longitudional data with constraint factor loadings
# and added residual correlations.
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

CFA_Sigma_L2 <- function(theta, S, C){
  NS <- ncol(S)
  NC <- length(C)
  t <- NS/NC
  lam <- theta[1:t]
  delta <- theta[(t+1):(t+NS+(NS*t/2))]
  phi <- theta[((t+NS+(NS*t/2))+1):(length(theta))]
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
  deltamat <- diag(delta[1:ncol(S)])
  cnt=12
  for (j in 1:(NS-t)){
    vec <- j+3*c(1:t)
    v = vec[vec<=NS]
    #(j+t):t:NS
    for (i in v){
      cnt=cnt+1
      deltamat[i,j] <- delta[cnt]; deltamat[j,i] = deltamat[i,j]
      i 
    }
  }
  Sig <- lammat%*%phimat%*%t(lammat)+deltamat  
  return(Sig)
}

# Calculate CFA log Function 
CFA_loglik2 <- function(theta, S, N){
  Sig <- CFA_Sigma_L2(theta,S,C)
  ((N-1)/2) * (log(abs(det(Sig))) + sum(diag(S%*%solve(Sig))))
}

# Estimation method I: Maximum likelihood
# Input: starting values theta, S
# Output: theta optimized
F_ML <- function(theta, S){
  Sig <- CFA_Sigma_L2(theta,S,C)
  log(det(Sig)) + sum(diag(S%*%solve(Sig)))
}

# Estimation method II: Unweighted results
# Input: starting values theta, S
# Output: theta optimized
F_ULS <- function(theta, S){
  S <- as.matrix(S)
  Sig <- CFA_Sigma_L2(theta,S,C)
  Z <- (S-Sig)
  0.5*sum(diag(Z%*%Z))
}

# Estimation method III: GLS 
# Input: starting values theta and S
# Output: theta optimized
F_GLS <- function(theta, S){
  Sig <- CFA_Sigma_L2(theta,S,C)
  W <- solve(Sig)   # or S
  Z <- (S-Sig)%*%W  # or diag(1,ncol(S)) for ULS
  0.5*sum(diag(Z%*%Z))
}

# Output results
CFA_output_L2 <- function(theta,N){
  NS <- ncol(S)
  NC <- length(C)
  t <- NS/NC
  lam <- theta[1:t]
  delta <- theta[(t+1):(t+NS+(NS*t/2))]
  phi <- theta[((t+NS+(NS*t/2))+1):(length(theta))]
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
  deltamat <- diag(delta[1:ncol(S)])
  cnt=12
  for (j in 1:(NS-t)){
    vec <- j+3*c(1:t)
    v = vec[vec<=NS]
    #(j+t):t:NS
    for (i in v){
      cnt=cnt+1
      deltamat[i,j] <- delta[cnt]; deltamat[j,i] = deltamat[i,j]
      i 
    }
  }
  loglik <- 2*CFA_loglik2(ML_output$par, S, N)
  return(list(lammat=lammat, phimat=phimat, deltamat=deltamat, loglik=abs(loglik)))
}

#-------------------------------------------------------->
# STEP III: Get estimates for factor loading, factor correlations, residual correlations

# Different factor structures
C <- c(3,3,3,3)  # 3 factors measured over 4 occasions

# Generate starting values
# add residual correlations (in delta)
# constrain lambda correlations
theta_fun <- function(S,C){
  th <- c(rep(0.8, (ncol(S)/length(C))), rep(0.1, (ncol(S))), rep(0.0,((ncol(S)*C[1])/2)), rep(0.4,((length(C)^2-length(C))/2)))
}

theta <- theta_fun(S,C)

# Run ML and print output
ML_output <- optim(par=theta, fn=F_ULS, S=S, method = "CG", control = list(maxit=50000))   # Compare ML and ULS

ML_output <- nlm(p=theta, f=F_ML, S=S)

CFA_output <- CFA_output_L2(ML_output$par,126)
CFA_output
ML_output$iterations

# Compare RC model and Saturated model
Sat_model <- (126-1)*(log(det(S)) + ncol(S))

q=abs(Sat_model-CFA_output$loglik)
df=78-39      # number of free parameters
tail_prob <- 1-(pchisq(q, df, ncp = 0, lower.tail = TRUE, log.p = FALSE))
tail_prob 

# Conclusion: 3-factor model is not significantly different than Null model