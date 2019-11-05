### This file contains minimal code to generate the first line of results shown in Table 2 in 
### "Penalized and Constrained Optimization: An Application to High-Dimensional Website Advertising"

### The complete code for reproducibility can be found in "PaC_Reproducibility_Revision.R"
### However, this is intended to show reproducibility in more narrow context for ease of explanation


## Packages used in data and correlation matrix generation for this code
library(MASS)
library(MBESS)

## Packages used in the main functions of the "PAC Functions.R" source code
library(lars)
library(limSolve)
library(quadprog)

## Source "PAC_Functions_Revised.R"
source("PaC_Functions_Revised.R")


## Set number of training/testing sets required for code
## For complete reproducibility of the results in the paper, iter=100
## This can be adjusted in the code for faster runtimes if desired

iter = 100

## For reproducibility, the seed vector used to generate data
set.seed(1234)
seed = sample(1:10000, iter, replace = F)

###################################################################################
################### Creating PAC Coefficient Paths ################################
###################           (Table 2)            ################################
###################################################################################

## Set list of runs for first row of Table 2

fits=vector("list", iter)

## Test set size: 10,000 observations
## For all runs in the paper and this file:
## s=5 (the number of non-zero random uniform components) and
## sigma=1 (standard normal distribution)

n.test = 10000
s = 5
sigma = 1

## First setting in table

n = 100
p = 50
m = 5


## Run complete iterations of PAC fits with no correlation setting
for (i in 1:iter) {
  
  fits[[i]] = compare.to.lasso(n = n, p = p, m = m, cov.mat = NULL, sigma = sigma, trace = F, seed = seed[i])
  
}

###################################################################################
###################    Calculating Error Rates     ################################
###################       and Relaxed Fits         ################################
###################           (Table 2)            ################################
###################################################################################


### Section 5.1: Comparison to Existing Lasso Methods
### Evaluating Results, PAC vs. Lasso (both relaxed and unrelaxed)

## Create objects to store SSE values
## Note there are four error rates calculated (in order):
## PAC
## Lasso
## Relaxed PAC
## Relaxed Lasso

error.constrained.valid = rep(0,iter)
error.lars.valid = rep(0,iter)
error.constrained.valid.relaxed = rep(0,iter)
error.lars.valid.relaxed = rep(0,iter)


## Create lists for runs of each Lasso model
lars.coefs = vector("list",iter)

## For validation purposes, record index and lambda value of lowest error per run per method
valid.place = rep(0,iter)
valid.lam = rep(0,iter)
valid.place.lars = rep(0,iter)
valid.lam.lars = rep(0,iter)

valid.place.relaxed = rep(0,iter)
valid.lam.relaxed = rep(0,iter)
valid.place.lars.relaxed = rep(0,iter)
valid.lam.lars.relaxed = rep(0,iter)

#LASSO runs, saving coefficients
for (i in 1:iter) {
  
  
  lars.trial = lars(fits[[i]]$data$x, fits[[i]]$data$y, normalize = T, intercept = F)
  lars.coefs[[i]] = predict(lars.trial, s = fits[[i]]$constrained.fit$lambda, type = "coefficients", mode = "lambda")$coef
  
}

## Set relaxed coefficients to current coefficients (update later)

lars.coefs.relaxed = lars.coefs
fits.relaxed = fits

## Create the required number of validation/testing data sets
valid.set.x = vector("list",iter)
valid.set.y = matrix(0,n,iter)
test.set.x = vector("list",iter)
test.set.y = matrix(0,n.test,iter)

## Generate and store the required number of validation/testing data sets
for (i in 1:iter) {
  
  set.seed(seed[i]+1)
  valid.set.x[[i]] = matrix(rnorm(n*p),n,p)
  valid.set.y[,i] = as.vector(valid.set.x[[i]]%*%fits[[i]]$data$beta+sigma*rnorm(n))
  test.set.x[[i]] = matrix(rnorm(n.test*p),n.test,p)
  test.set.y[,i] = as.vector(test.set.x[[i]]%*%fits[[i]]$data$beta)
  
}

## Run relaxed PAC and relaxed Lasso
## Note the relaxed PAC requires a deconstructed formulation
## The X matrix (X.til) and Y vector (Y.til) are calculated by reformulating the PAC
## As shown in the PAC paper (e.g. Generalized Lasso and Appendix A)

for (i in 1:iter) {
  
  b = fits[[i]]$data$b
  y.data = fits[[i]]$data$y
  
  
  n.lam = length(fits[[i]]$constrained.fit$lambda)
  error = rep(0,n.lam)
  error.lars = rep(0,n.lam)
  error.relaxed = rep(0,n.lam)
  error.lars.relaxed = rep(0,n.lam)
  
  ## Calculate relaxed versions of each model at every lambda value
  for (k in 1:n.lam) {
    
    
    beta2.index = fits[[i]]$constrained.fit$b2index[,k]
    x2.data = fits[[i]]$data$x[,beta2.index]
    A2 = fits[[i]]$data$C.full[,beta2.index]
    Y.til = y.data-x2.data%*%solve(A2)%*%b
    
    total = c(1:p)
    
    beta1.index.adj = NULL
    
    ## Identify non-zero vector entries
    for (j in 1:p) {
      if (fits[[i]]$constrained.fit$coefs[,k][j]==0) {
        beta1.index.adj = c(beta1.index.adj,j)
      }
    }
    
    beta1.index = total[-c(beta2.index,beta1.index.adj)]
    
    if(length(beta1.index)==0) {
      
      fits.relaxed[[i]]$constrained.fit$coefs[,k]==fits[[i]]$constrained.fit$coefs[,k]
      
    } else {
      
      x1.star = fits[[i]]$data$x[,beta1.index]
      A1.star = fits[[i]]$data$C.full[,beta1.index]
      X.til = x1.star-x2.data%*%solve(A2)%*%A1.star
      
      ## Calculate relaxed PAC using lars function
      if (n<length(beta1.index)) {
        
        lars.trial = lars(X.til, Y.til, normalize = T, intercept = F)
        relax.fit=predict(lars.trial, s = fits[[i]]$constrained.fit$lambda, type = "coefficients", mode = "lambda")$coef
        fits.relaxed[[i]]$constrained.fit$coefs[,k][beta1.index] = relax.fit[k,]
        beta2 = solve(A2)%*%(b-A1.star%*%relax.fit[k,])
        
      } else {
        relax.fit = lm(Y.til~X.til-1)
        if(length(beta1.index)==1) {
          beta2 = solve(A2)%*%(b-A1.star*relax.fit$coef)
        } else {
          
          beta2 = solve(A2)%*%(b-A1.star%*%relax.fit$coef)
        }
        fits.relaxed[[i]]$constrained.fit$coefs[,k][beta1.index] = relax.fit$coef
      }
      
      fits.relaxed[[i]]$constrained.fit$coefs[,k][beta2.index] = beta2
      
    }
    
    relax.index = NULL
    
    ## Identify non-zero coefficient values for Lasso run
    for (j in 1:p) {
      
      if (lars.coefs[[i]][k,][j] != 0) {
        
        relax.index = c(relax.index,j)
        
      }
      
    }
    
    ## Calculate relaxed Lasso fit
    if(!is.null(relax.index)) {
      
      x1.lars = fits[[i]]$data$x[,relax.index]
      relax.fit = lm(y.data~x1.lars-1)
      lars.coefs.relaxed[[i]][k,][relax.index] = relax.fit$coef
      
    }
    
    ## Calculate SSE using pred.error() function for each of four methods at this lambda
    
    error[k] = pred.error(valid.set.x[[i]], valid.set.y[,i], fits[[i]]$constrained.fit$coefs[,k])
    error.lars[k] = pred.error(valid.set.x[[i]], valid.set.y[,i], lars.coefs[[i]][k,])
    error.relaxed[k] = pred.error(valid.set.x[[i]], valid.set.y[,i], fits.relaxed[[i]]$constrained.fit$coefs[,k])
    error.lars.relaxed[k] = pred.error(valid.set.x[[i]], valid.set.y[,i], lars.coefs.relaxed[[i]][k,])
    
  }
  
  ## Identify validation minimum SSE and the lambda at which it occurs for all four methods
  valid.place[i] = which.min(error)
  valid.lam[i] = fits[[i]]$constrained.fit$lambda[valid.place[i]]
  valid.place.lars[i] = which.min(error.lars)
  valid.lam.lars[i] = fits[[i]]$constrained.fit$lambda[valid.place.lars[i]]
  valid.place.relaxed[i] = which.min(error.relaxed)
  valid.lam.relaxed[i] = fits[[i]]$constrained.fit$lambda[valid.place.relaxed[i]]
  valid.place.lars.relaxed[i] = which.min(error.lars.relaxed)
  valid.lam.lars.relaxed[i] = fits[[i]]$constrained.fit$lambda[valid.place.lars.relaxed[i]]
  
  ## Calculate RMSE using minimum validation values applied to testing data
  error.constrained.valid[[i]] = sqrt(pred.error(test.set.x[[i]], test.set.y[,i], fits[[i]]$constrained.fit$coefs[,valid.place[i]])/10000)
  error.lars.valid[[i]] = sqrt(pred.error(test.set.x[[i]], test.set.y[,i], lars.coefs[[i]][valid.place.lars[i],])/10000)
  error.constrained.valid.relaxed[[i]] = sqrt(pred.error(test.set.x[[i]], test.set.y[,i], fits.relaxed[[i]]$constrained.fit$coefs[,valid.place.relaxed[i]])/10000)
  error.lars.valid.relaxed[[i]] = sqrt(pred.error(test.set.x[[i]], test.set.y[,i], lars.coefs.relaxed[[i]][valid.place.lars.relaxed[i],])/10000)
  
}


###################################################################################
###################    Calculating Error Rates     ################################
###################           (Table 2)            ################################
###################################################################################

## Table 2 Results

# Avg RMSE and SE, Setting 1, No Correlation
## (First row of table)

mean(error.lars.valid)
sd(error.lars.valid)/sqrt(iter)
mean(error.constrained.valid)
sd(error.constrained.valid)/sqrt(iter)
mean(error.lars.valid.relaxed)
sd(error.lars.valid.relaxed)/sqrt(iter)
mean(error.constrained.valid.relaxed)
sd(error.constrained.valid.relaxed)/sqrt(iter)








