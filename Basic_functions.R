# Basix functions

#===========================================================>
# Author: Clara-Christina Gerstner
# Purpose: To provide a set of basic functions to be called
# into higher level functions.
#===========================================================>

### Sum ###

# Function

sum_fx <- function(x, n){
  sum_x = 0
  for (i in 1:n){
  sum_x <- sum_x + x[i]
  }
  sum_x
}

# Compare result to R function
X <- c(runif(10,1,10))
all.equal(sum_fx(X,length(X)), sum(X))

#===========================================================>

### Mean ###

# Function

mean_fx <- function(x, n){
  sum_fx(x, n)/n
}

# Compare result to R function
all.equal(mean_fx(X,length(X)), mean(X))

#===========================================================>

### Variance ###

# Function

Var_fx <- function(x, n) {
  sum_fx(((x-(mean_fx(x, n)))^2),n)/(n-1)
}

# Compare result to R function
all.equal(Var_fx(X,length(X)), var(X))

#===========================================================>

### Standard deviation ###

# Function

SD_fx <- function(x,n) {
  sqrt(Var_fx(x, n))
}

# Compare result to R function
all.equal(SD_fx(X,length(X)), sd(X))

#===========================================================>

### Minimum ###

# Function

min_fx <- function(x,n) {
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
all.equal(min_fx(X,length(X)), min(X))

#===========================================================>

### Maximum ###

# Function

max_fx <- function(x,n) {
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
all.equal(max_fx(X,length(X)), max(X))

#===========================================================>

### Range ###

# Function

range_fx <- function(x,n) {
  return(c(min_fx(x,n), max_fx(x,n)))
}

# Compare result to R function
all.equal(range_fx(X,length(X)), range(X))

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