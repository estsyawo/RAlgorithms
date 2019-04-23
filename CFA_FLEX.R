# CFA Model Different estimation methods
#===========================================================>
# Author: Clara-Christina Gerstner
# Purpose: Build a flexible version of confirmatory factor analysis.
#===========================================================>

#-------------------------------------------------------->
# STEP I: Set correlation matrix for S
source("cfa_data.R")
S

#-------------------------------------------------------->
# STEP II: Load functions

# Estimate Sigma (theta) and plot output
# Input: theta 
# Output: Sigma

CFA_Sigma_F <- function(theta, S, C){
  NS <- ncol(S)
  NC <- length(C)
  lam <- theta[1:NS]
  delta <- theta[(NS+1):(2*NS)]
  phi <- theta[(2*NS+1):((2*NS)+((NC^2-NC)/2))]
  lammat <- matrix(0, ncol=NC, nrow=NS)
  Ci <- C[1]
  lammat[1:Ci,1] <- lam[1:Ci]
  for (i in 1:(NC-1)){
    lammat[(Ci+1):(Ci+C[i+1]),(i+1)] <- lam[(Ci+1):(Ci+C[i+1])]
    Ci <- (Ci+C[i+1])
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
  S <- as.matrix(S)
  Sig <- CFA_Sigma_F(theta,S,C)
  (((N-1)/2) * (log(det(Sig)) + sum(diag(S%*%solve(Sig)))))  
}

# Estimation method I: Maximum likelihood
# Input: starting values theta, S
# Output: theta optimized
F_ML <- function(theta, S){
  Sig <- CFA_Sigma_F(theta,S,C)
  log(det(Sig)) + sum(diag(S%*%solve(Sig)))
}

# Estimation method II: Unweighted results
# Input: starting values theta, S
# Output: theta optimized
F_ULS <- function(theta, S){
  S <- as.matrix(S)
  Sig <- CFA_Sigma_F(theta,S,C)
  Z <- (S-Sig)
  0.5*sum(diag(Z%*%Z))
}

# Estimation method III: GLS 
# Input: starting values theta and S
# Output: theta optimized
F_GLS <- function(theta, S){
  Sig <- CFA_Sigma_F(theta,S,C)
  W <- solve(Sig)   # or S
  Z <- (S-Sig)%*%W  # or diag(1,ncol(S)) for ULS
  0.5*sum(diag(Z%*%Z))
}

# Output results
CFA_output_F <- function(theta, N){
  NS <- ncol(S)
  NC <- length(C)
  lam <- theta[1:NS]
  delta <- theta[(NS+1):(2*NS)]
  phi <- theta[(2*NS+1):((2*NS)+((NC^2-NC)/2))]
  lammat <- matrix(0, ncol=NC, nrow=NS)
  Ci <- C[1]
  lammat[1:Ci,1] <- lam[1:Ci]
  for (i in 1:(NC-1)){
    lammat[(Ci+1):(Ci+C[i+1]),(i+1)] <- lam[(Ci+1):(Ci+C[i+1])]
    Ci <- (Ci+C[i+1])
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
  return(list(lammat=lammat, phimat=phimat, deltamat=deltamat, loglik=loglik))
}

#-------------------------------------------------------->
# STEP III: Run functions

# Different factor structures
C <- c(5,5,5) # 3 factors, 5 variables each, correct model
# C <- c(2,3,4,6) # 4 factors, different number of variables

# Generate starting values
theta_fun <- function(S,C){
  th <- c(rep(0.8, ncol(S)), rep(0.5, ncol(S)), rep(0.1,((length(C)^2-length(C))/2)))
}

# Run ML and print output
ML_output <- optim(par=theta_fun(S,C), fn=F_ULS, S=S, method = "BFGS", control = list(maxit=50000)) # compare ML, ULS, GLS

# 3 factor model, N=150
CFA_output <- CFA_output_F(ML_output$par,150)
CFA_output

# 4 factor model, N=150
CFA_output2 <- CFA_output_F(ML_output$par, 150)
CFA_output2

#-------------------------------------------------------->
# STEP IV: Compare models

# Compare 3-factor and 4-factor model
q=CFA_output2$loglik-CFA_output$loglik
df= 36-33      
tail_prob <- 1-(pchisq(q, df, ncp = 0, lower.tail = TRUE, log.p = FALSE))
tail_prob 

# Conclusion: 3-factor model is clearly a better fit than 4-factor model

# Compare 3-factor and Saturated model
Sat_model <- (N-1)*(log(det(S)) + ncol(S))

q=CFA_output$loglik-Sat_model
df=120-33      # number of free parameters
tail_prob <- 1-(pchisq(q, df, ncp = 0, lower.tail = TRUE, log.p = FALSE))
tail_prob 

# Conclusion: 3-factor model is significantly different than Saturated model
# LISREL output tailprob is above 0.05

#-------------------------------------------------------->
# STEP V: Fit statistics

N <- 150
p <- ncol(S)  # number of variables
t <- length(theta_fun(S,C))  # number of estimated parameters
df <- ((p*(p+1))/2) - t

# AIC
(AIC_estimate <- 2*t + CFA_output$loglik)

# BIC
(BIC_estimate <- t*log10(N) + CFA_output$loglik)

# NNFI = TLI
Chib <- CFA_output$loglik
Chim <- CFA_output$loglik - Sat_model
dfb <- 105
dfm <- df
Delm <- max((Chim-dfm),0)
Delb <- max((Chib-dfb),0)

(TLI <- ((Chib/dfb)-(Chim/dfm))/((Chib/dfb)-1)) # Bollen p.273

# RMSEA
(RMSEA_estimate <- sqrt(max((((CFA_output$loglik/df) - 1)/(N - 1)) , 0)))

(RMSEA_estimate2 <- sqrt(Delm/(dfm*(N-1))))

# CFI
(CFI <- 1 - (Delm/Delb))

# SRMR
Sig <- CFA_Sigma_F(ML_output$par,S,C)
lobs <-  S[!lower.tri(S)]
limp <-  Sig[!lower.tri(Sig)]

(SRMR_estimate <- sqrt(mean((limp - lobs)^2)))

#-------------------------------------------------------->
# Compare Output lavaan
CFA.model <- '  Conflict =~ NA*c1 + c2 + c3 + c4 + c5 
Ambivalence =~ NA*a1 + a2 + a3 + a4 + a5
Maintenance   =~ NA*m1 + m2 + m3 + m4 + m5 
Conflict ~~ 1*Conflict
Ambivalence ~~ 1*Ambivalence
Maintenance ~~ 1*Maintenance'

library(lavaan)
S = data.frame(S)
row.names(S) = c("c1","c2","c3","c4","c5","a1","a2","a3","a4","a5","m1","m2","m3","m4","m5")
names(S) = c("c1","c2","c3","c4","c5","a1","a2","a3","a4","a5","m1","m2","m3","m4","m5")
S <- as.matrix(S)

fit <- cfa(CFA.model, sample.cov = as.matrix(S), sample.nobs = 150)
summary(fit, fit.measures=TRUE)
fitMeasures(fit)
