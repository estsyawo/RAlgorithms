# CFA Model Different estimation methods
#===========================================================>
# Author: Clara-Christina Gerstner
# Purpose: Obtain estimates for confirmatory factor model
# using maximum likelihood with optim. Apply approach to bk4 data set.
#===========================================================>

#-------------------------------------------------------->
# STEP I: Set correlation matrix for S

# 15 variables, 3 constructs with 5 variables each (conflict, ambivalence, maintenance)
S <- matrix(c(1.00000, 0.58532, 0.60861, 0.63821, 0.60726 , 0.19634,  0.20470,  0.16537,  0.15152,  0.12835,  0.06796, 0.19126,  0.18235,  0.06326, 0.12748,
               0.58532, 1.00000, 0.58533, 0.57862, 0.59925,  0.30879,  0.28136 , 0.25502,  0.28721,  0.24919,  0.15468, 0.21003,  0.17531,  0.15608, 0.15254,
               0.60861, 0.58533, 1.00000, 0.60880, 0.58689 , 0.29356 , 0.33663,  0.33472,  0.35463,  0.34669,  0.11121, 0.13194,  0.06865,  0.12054, 0.13924,
               0.63821, 0.57862, 0.60880, 1.00000, 0.59694,  0.28820,  0.30869 , 0.28084,  0.33568,  0.31658,  0.09564, 0.21647,  0.19861,  0.12256, 0.14209,
               0.60726, 0.59925, 0.58689, 0.59694, 1.00000 , 0.28436 , 0.26227,  0.24877,  0.34172,  0.35052,  0.07250, 0.16643,  0.10479,  0.11099, 0.07444,
               0.19634, 0.30879, 0.29356, 0.28820, 0.28436,  1.00000,  0.62700 , 0.66242,  0.58551,  0.67884, -0.08351, 0.07671,  0.04059, -0.04901, 0.03420,
               0.20470, 0.28136, 0.33663, 0.30869, 0.26227 , 0.62700 , 1.00000,  0.66926,  0.65772,  0.60823, -0.12091, 0.01050, -0.08253, -0.16360, 0.03734,
               0.16537, 0.25502, 0.33472, 0.28084, 0.24877,  0.66242,  0.66926 , 1.00000,  0.66722,  0.63258, -0.08820, 0.07863,  0.05481, -0.08883, 0.14695,
               0.15152, 0.28721, 0.35463, 0.33568, 0.34172 , 0.58551 , 0.65772,  0.66722,  1.00000,  0.62151, -0.02940, 0.03248,  0.03006, -0.04960, 0.08688,
               0.12835, 0.24919, 0.34669, 0.31658, 0.35052,  0.67884,  0.60823 , 0.63258,  0.62151,  1.00000, -0.10307, 0.06131, -0.08248, -0.06565, 0.05354,
               0.06796, 0.15468, 0.11121, 0.09564, 0.07250 ,-0.08351 ,-0.12091, -0.08820, -0.02940, -0.10307,  1.00000, 0.43092,  0.50106,  0.46820, 0.45281,
               0.19126, 0.21003, 0.13194, 0.21647, 0.16643,  0.07671,  0.01050 , 0.07863,  0.03248,  0.06131,  0.43092, 1.00000,  0.42783,  0.46365, 0.50013,
               0.18235, 0.17531, 0.06865, 0.19861, 0.10479 , 0.04059 ,-0.08253,  0.05481,  0.03006, -0.08248,  0.50106, 0.42783,  1.00000,  0.56989, 0.56869,
               0.06326, 0.15608, 0.12054, 0.12256, 0.11099, -0.04901, -0.16360 ,-0.08883, -0.04960, -0.06565,  0.46820, 0.46365,  0.56989,  1.00000, 0.60948,
               0.12748, 0.15254, 0.13924, 0.14209, 0.07444 , 0.03420 , 0.03734,  0.14695,  0.08688,  0.05354,  0.45281, 0.50013,  0.56869,  0.60948, 1.00000),ncol=15, nrow=15)

#-------------------------------------------------------->
# STEP II: Load functions

# Estimate Sigma (theta) and plot output
# Input: theta (33 values)
# Output: Sigma
CFA_Sigma <- function(theta){
  lam <- theta[1:15]
  phi <- theta[16:18]
  delta <- theta[19:33]
  lammat <- matrix(0, ncol=3, nrow=15)
  lammat[1:5,1] <- lam[1:5]
  lammat[6:10,2] <- lam[6:10]
  lammat[11:15,3] <- lam[11:15]
  phimat <- diag(3)
  phimat[2,1] <- phi[1]; phimat[1,2] = phimat[2,1]
  phimat[3,1] <- phi[2]; phimat[1,3] = phimat[3,1]
  phimat[3,2] <- phi[3]; phimat[2,3] = phimat[3,2]
  deltamat <- diag(delta)
  Sig <- lammat%*%phimat%*%t(lammat)+deltamat  
  return(Sig)
}
CFA_Sigma(th2)

# Output function
# Input: theta
# Output: lambda matrix, phi matrix, delta matrix
CFA_output <- function(theta){
  lam <- theta[1:15]
  phi <- theta[16:18]
  delta <- theta[19:33]
  lammat <- matrix(0, ncol=3, nrow=15)
  lammat[1:5,1] <- lam[1:5]
  lammat[6:10,2] <- lam[6:10]
  lammat[11:15,3] <- lam[11:15]
  phimat <- diag(3)
  phimat[2,1] <- phi[1]; phimat[1,2] = phimat[2,1]
  phimat[3,1] <- phi[2]; phimat[1,3] = phimat[3,1]
  phimat[3,2] <- phi[3]; phimat[2,3] = phimat[3,2]
  deltamat <- diag(delta)
  return(list(lammat=lammat, phimat=phimat, deltamat=deltamat))
}

# Maximum likelihood function
# Input: theta and S
F_ML <- function(theta, S){
  Sig <- CFA_Sigma(theta)
  log(det(Sig)) + sum(diag(S%*%solve(Sig)))
}

# FGLS function
# Input: theta and S
F_GLS <- function(theta, S){
  S <- as.matrix(S)
  Sig <- CFA_Sigma(theta)
  W <- solve(Sig)   # or S
  Z <- (S-Sig)%*%W  # or diag(1,ncol(S)) for ULS
  0.5*sum(diag(Z%*%Z))
}
F_GLS(th2,S)

# FGLS2 function
# Input: theta and S
F_GLS2 <- function(theta, S){
  S <- as.matrix(S)
  Sig <- CFA_Sigma(theta)
  Z <- (S-Sig)  # or diag(1,ncol(S)) for ULS
  W <- solve(Z%*%Z)
  0.5*sum(diag(Z%*%W%*%Z))
}
F_GLS2(th2,S)

# FULS function
# Input: theta and S
F_ULS <- function(theta, S){
  S <- as.matrix(S)
  Sig <- CFA_Sigma(theta)
  Z <- (S-Sig)
  0.5*sum(diag(Z%*%Z))
}
F_ULS(th2,S)

#-------------------------------------------------------->
# STEP III: Run functions

# Try different Starting values
th1 <- c(0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.6,0.6,0.6,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25, 0.25,0.25,0.25,0.25,0.25,0.25)
th2 <- c(0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.01,0.01,0.01,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5)
th3 <- runif(33)

# Run ML and print output
ML_output <- optim(par=th2, fn=F_ML, S=S)

CFA_output(ML_output$par)

# Run GLS and print output
GLS_output <- optim(par=th2, fn=F_GLS2, S=S)

CFA_output(GLS_output$par)

# Run GLS and print output
ULS_output <- optim(par=th2, fn=F_ULS, S=S)

CFA_output(ULS_output$par)

#-------------------------------------------------------->
# STEP IV: get log likelihood values

# Calculate CFA log Function 1
CFA_loglik <- function(theta, S){
  S <- as.matrix(S)
  Sig <- CFA_Sigma(theta)
  log(det(Sig)) + sum(diag(S%*%solve(Sig))) - log(det(S)) - nrow(S) 
}
CFA_loglik(th2, S)

# Calculate CFA log Function 2
CFA_loglik2 <- function(theta, S, N){
  S <- as.matrix(S)
  Sig <- CFA_Sigma(theta)
  (((N-1)/2) * (log(det(Sig)) + sum(diag(S%*%solve(Sig)))))  
}
L0 <- 2*CFA_loglik2(th2, S, 150) # number of observations

# Finding: th1 provides the best log likelihood values

#-------------------------------------------------------->
# Compare to LISREL output:

# includes all three outpu matrices
th_LISREL <- c(.7769,.7561,.7773,.7907,.7705,.7959,.8024,.8245,.7907,.7905,.6322,.6249,.7409,.7691,.7690,.4443,.2351,-.015,.3965,.4282,.3958,.3748,.4060,.3666,.3562,.3202,.3746,.3752,.6003,.6095,.4511,.4085,.4086)

CFA_loglik(th_LISREL, S) 

L1 <- 2*CFA_loglik2(th_LISREL, S, 150)

#-------------------------------------------------------->
# Fit statistics

df=((15*16)/2) - length(theta)
t <- length(th2)
N <- 150
p <- 15

# AIC
(AIC_estimate0 <- 2*t + L0)

(AIC_estimate1 <- 2*t + L1)

# BIC
(BIC_estimate0 <- t*log10(N) + L0)

(BIC_estimate1 <- t*log10(N) + L1)

# Chi square test
Null_model <- (N-1)*(log(det(S)) + p)

Diff_model0 <- L0 - Null_model

Diff_model1 <- L1 - Null_model

pchisq(Diff_model1, df=87, lower.tail=FALSE)
# not rejected at 0.05

pchisq(Diff_model0, df=df, lower.tail=FALSE)
# rejected at 0.05

# NNFI = TLI
Chib <- L1 # 1184.008 # Baseline model see Lavaan Output
Chim <- Diff_model1
df2 <- 105

TLI <- ((Chib/df2)-(Chim/df))/((Chib/df2)-1) # Bollen p.273
TLI

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
