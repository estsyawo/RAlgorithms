# Test item correction
#===========================================================>
# Author: Clara-Christina Gerstner
# Purpose: Delete test items that do not contribute to the
# variance of a scale by evaluation alpha if deleted coefficients
# and item total correlations (r).
#===========================================================>

# Load sample data 
# Test 2 Data: 25 items, N= 660

TEST2 <- read.csv("TEST2.csv", header = TRUE)

#===========================================================>

# Compute item total correlation (r) and cronbach's alpha

library(multilevel)
Itemtot <- item.total(scale(TEST2))
Alph <- as.numeric(cronbach(scale(TEST2)))[1]

#===========================================================>

# Algorithm to drop items that reduce performance of scale
# Rule: each item needs to meet r=[0.2, 0.8]
# Algorithm drops item with lowest r one at a time

Drop_items_fx <- function(Test){
  Itemtot <- item.total(scale(Test))
  zg = 1
  while(zg==1){
  k= which.min(Itemtot[,2])
  if (Itemtot[k,2] < 0.2){
    print("Item dropped:")
    print(as.character(Itemtot[k,1]))
    Itemtot=Itemtot[-k,]
    Test <- Test[,-k]
    print(round(as.numeric(cronbach(scale(Test)))[1],5))
    }else{
    zg=-1
    }
  }
  zg = 1
  while(zg==1){
  k= which.max(Itemtot[,2])
  if (Itemtot[k,2] > 0.8){
    print("Item dropped:")
    print(as.character(Itemtot[k,1]))
    Itemtot=Itemtot[-k,]
    Test <- Test[,-k]
    print(round(as.numeric(cronbach(scale(Test)))[1],5))
    }else{
    zg=-1
    }
  }
  return(Itemtot)
}

#===========================================================>

# Example: Run algorithm for Test 2 data

Drop_items_fx(TEST2)
