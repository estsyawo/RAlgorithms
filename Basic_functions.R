# Basic functions

#===========================================================>
# Author: Clara-Christina Gerstner
# Purpose: To provide a set of basic functions to be called
# into higher level functions.
#===========================================================>

### Sum ###

# Function

sum_fx <- function(x){
  n <- length(x)
  sum_x = 0
  for (i in 1:n){
  sum_x <- sum_x + x[i]
  }
  sum_x
}

# Compare result to R function
X <- c(runif(10,1,10))
all.equal(sum_fx(X), sum(X))

#===========================================================>

### Mean ###

# Function

mean_fx <- function(x){
  if (ncol(x) == 1 || is.null(ncol(x))){
  n <- length(x)
  S <- sum_fx(x)/n
  }
  else{
    S <-matrix(0:0,ncol=1,nrow=nrow(x)) 
    for (i in 1:ncol(x)){
      n <- nrow(x)
      S[i] <- (sum_fx(x[,i]))/n
    }
    return(S)
  }
}

# Compare result to R function
all.equal(mean_fx(X), mean(X))

# Compare result to R function for Matrix
M <- matrix(runif(100,1,10),10,10)
all.equal(as.matrix(apply(M,2,mean)),mean_fx(M))

#===========================================================>

### Variance ###

# Function

Var_fx <- function(x) {
  if (ncol(x) == 1 || is.null(ncol(x))){
  n <- length(x)
  S <- sum_fx(((x-(mean_fx(x)))^2))/(n-1)
  }
  else{
    S <-matrix(0:0,ncol=1,nrow=nrow(x)) 
  for (i in 1:ncol(x)){
  n <- nrow(x)
  S[i] <- sum_fx(((x[,i]-(mean_fx(x[,i])))^2))/(n-1)
  }
  return(S)
  }
}

# Compare result to R function for Vector
all.equal(Var_fx(X), var(X))

# Compare result to R function for Matrix
M <- matrix(runif(100,1,10),10,10)
all.equal(as.matrix(apply(M,2,var)),Var_fx(M))

#===========================================================>

### Standard deviation ###

# Function

SD_fx <- function(x) {
  sqrt(Var_fx(x))
}

# Compare result to R function
all.equal(SD_fx(X), sd(X))

#===========================================================>

### Minimum ###

# Function

min_fx <- function(x) {
  n <- length(x)
  if (n > 1) {  
    min_x = x[1]
    for (i in 2:n) {
      if (min_x > x[i]) min_x = x[i]
    }
    return(min_x)
  }
  else{
    return(x[1])
  }
}

# Compare result to R function
all.equal(min_fx(X), min(X))

#===========================================================>

### Maximum ###

# Function

max_fx <- function(x) {
  n <- length(x)
  if (n > 1) {  
    max_x = x[1]
    for (i in 2:n) {
    if (max_x < x[i]) max_x = x[i]
  }
  return(max_x)
  }
  else{
    return(x[1])
  }
}

# Compare result to R function
all.equal(max_fx(X), max(X))

#===========================================================>

### Range ###

# Function

range_fx <- function(x) {
  return(c(min_fx(x), max_fx(x)))
}

# Compare result to R function
all.equal(range_fx(X), range(X))

#===========================================================>

### Median ###

median_fx <- function(l) {
  l <- sort(l); nl = length(l)
  if (nl%%2 != 0){
    return(l[(nl/2)+1])}
  else {
    return(0.5*(l[nl/2+1] + l[nl/2]))}
}

# Compare result to R function
all.equal(median_fx(X), median(X))
