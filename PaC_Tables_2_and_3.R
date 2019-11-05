### This file contains the complete code used to generate the results shown in
### "Penalized and Constrained Optimization: An Application to High-Dimensional Website Advertising"
### In most cases, the full code has been included to allow for step-by-step evaluation
### However, comments and explanations for the code are only included in the first example of a given block of code
### The other blocks run similarly

##  For a complete walkthrough of a single run of code in this file, please refer to
##  PaC_Simple_Table2.R, which provides a more complete, easier to follow example of this complete code

## Packages used in data and correlation matrix generation for this code
library(MASS)
library(MBESS)

## Packages used in the main functions of the "PAC Functions.R" source code
library(lars)
library(limSolve)
library(quadprog)

## Source "PAC Functions.R"
source("PaC_Functions_Revised.R")


## Set number of training/testing sets required for code
## For complete reproducibility of the results in the paper, iter=100
## This can be adjusted in the code for faster runtimes if desired

iter=100

## For reproducibility, the seed vector used to generate data
set.seed(1234)
seed=sample(1:10000, iter, replace=F)

###################################################################################
################### Creating PAC Coefficient Paths ################################
###################           (Table 2)            ################################
###################################################################################

## Set list of runs for each setting in Table 2
## Note fits.4 corresponds to Setting 1, fits.5 to Setting 2, and fits.7 to Setting 3
## Each list of fits has a corresponding correlation structure fit
fits.4=vector("list", iter)
fits.4.corr.dying=vector("list", iter)
fits.5=vector("list", iter)
fits.5.corr.dying=vector("list", iter)
fits.7=vector("list", iter)
fits.7.corr.dying=vector("list", iter)

### Section 5.1: PAC Comparison to Existing Lasso Methods
### Running PAC Fits

## Test set size: 10,000 observations
## For all runs in the paper and this file:
## s=5 (the number of non-zero random uniform components) and
## sigma=1 (standard normal distribution)

n.test=10000
s=5
sigma=1

## Setting 1

n=100
p=50
m=5

## Generate correlation matrix (dying-off correlation)
corr=matrix(0,p,p)
for (i in 1:p) {
  for (j in 1:p) {
    if (i==j) {
      corr[i,j]=1 }
    else {corr[i,j]=0.5^abs(i-j)}
  }
}

## Correlation matrix to covariance matrix
cov.mat2=cor2cov(corr,rep(1,p))

## Run complete iterations of PAC fits for correlation setting
for (i in 1:iter) {
  
  fits.4.corr.dying[[i]]=compare.to.lasso(n=n,p=p,m=m,cov.mat=cov.mat2,sigma=sigma,trace=F,seed=seed[i])
  
}

## Run complete iterations of PAC fits for no correlation setting
for (i in 1:iter) {
  
  fits.4[[i]]=compare.to.lasso(n=n,p=p,m=m,cov.mat=NULL,sigma=sigma,trace=F,seed=seed[i])
  
}


## Setting 2

n=50
p=500
m=10


#Generate correlation matrix
corr=matrix(0,p,p)
for (i in 1:p) {
  for (j in 1:p) {
    if (i==j) {
      corr[i,j]=1 }
    else {corr[i,j]=0.5^abs(i-j)}
  }
}

cov.mat2=cor2cov(corr,rep(1,p))

for (i in 1:iter) {
  
  fits.5.corr.dying[[i]]=compare.to.lasso(n=n,p=p,m=m,cov.mat=cov.mat2,sigma=sigma,trace=F,seed=seed[i])
  
}

for (i in 1:iter) {
  
  fits.5[[i]]=compare.to.lasso(n=n,p=p,m=m,cov.mat=NULL,sigma=sigma,trace=F,seed=seed[i])
  
}

## Setting 3

n=50
p=100
m=60

corr=matrix(0,p,p)
for (i in 1:p) {
  for (j in 1:p) {
    if (i==j) {
      corr[i,j]=1 }
    else {corr[i,j]=0.5^abs(i-j)}
  }
}

cov.mat2=cor2cov(corr,rep(1,p))

for (i in 1:iter) {
  
  fits.7.corr.dying[[i]]=compare.to.lasso(n=n,p=p,m=m,cov.mat=cov.mat2,sigma=sigma,trace=F,seed=seed[i])
  
}

for (i in 1:iter) {
  
  fits.7[[i]]=compare.to.lasso(n=n,p=p,m=m,cov.mat=NULL,sigma=sigma,trace=F,seed=seed[i])
  
}


###################################################################################
################### Creating PAC Coefficient Paths ################################
###################           (Table 3)            ################################
###################################################################################


### Section 5.2: Violations of Constraints
### Running PAC Fits

##Setting 3

n=50
p=100
m=60

## Four runs of this setting, one for each value of a
fits.test3.25=vector("list",iter)
fits.test3.5=vector("list",iter)
fits.test3.75=vector("list",iter)
fits.test31=vector("list",iter)

## PAC Runs with each error vector
for (i in 1:iter) {

set.seed(seed[i])
error.vec=runif(m,0,.25)
fits.test3.25[[i]]=compare.to.lasso(err=error.vec,n=n,p=p,m=m,sigma=1,seed=seed[i],backwards=F)
error.vec=runif(m,0,.5)
fits.test3.5[[i]]=compare.to.lasso(err=error.vec,n=n,p=p,m=m,sigma=1,seed=seed[i],backwards=F)
error.vec=runif(m,0,.75)
fits.test3.75[[i]]=compare.to.lasso(err=error.vec,n=n,p=p,m=m,sigma=1,seed=seed[i],backwards=F)
error.vec=runif(m,0,1.0)
fits.test31[[i]]=compare.to.lasso(err=error.vec,n=n,p=p,m=m,sigma=1,seed=seed[i],backwards=F)


}


## Setting 2

n=50
p=500
m=10


fits.test2.25=vector("list",iter)
fits.test2.5=vector("list",iter)
fits.test2.75=vector("list",iter)
fits.test21=vector("list",iter)

for (i in 1:iter) {

set.seed(seed[i])
error.vec=runif(m,0,.25)
fits.test2.25[[i]]=compare.to.lasso(err=error.vec,n=n,p=p,m=m,sigma=1,seed=seed[i],backwards=F)
error.vec=runif(m,0,.5)
fits.test2.5[[i]]=compare.to.lasso(err=error.vec,n=n,p=p,m=m,sigma=1,seed=seed[i],backwards=F)
error.vec=runif(m,0,.75)
fits.test2.75[[i]]=compare.to.lasso(err=error.vec,n=n,p=p,m=m,sigma=1,seed=seed[i],backwards=F)
error.vec=runif(m,0,1.0)
fits.test21[[i]]=compare.to.lasso(err=error.vec,n=n,p=p,m=m,sigma=1,seed=seed[i],backwards=F)


}

## Setting 1

n=100
p=50
m=5


fits.test1.25=vector("list",iter)
fits.test1.5=vector("list",iter)
fits.test1.75=vector("list",iter)
fits.test11=vector("list",iter)

for (i in 1:100) {

set.seed(seed[i])
error.vec=runif(m,0,.25)
fits.test1.25[[i]]=compare.to.lasso(err=error.vec,n=n,p=p,m=m,sigma=1,seed=seed[i],backwards=F)
error.vec=runif(m,0,.5)
fits.test1.5[[i]]=compare.to.lasso(err=error.vec,n=n,p=p,m=m,sigma=1,seed=seed[i],backwards=F)
error.vec=runif(m,0,.75)
fits.test1.75[[i]]=compare.to.lasso(err=error.vec,n=n,p=p,m=m,sigma=1,seed=seed[i],backwards=F)
error.vec=runif(m,0,1.0)
fits.test11[[i]]=compare.to.lasso(err=error.vec,n=n,p=p,m=m,sigma=1,seed=seed[i],backwards=F)


}

###################################################################################
###################    Calculating Error Rates     ################################
###################       and Relaxed Fits         ################################
###################           (Table 2)            ################################
###################################################################################


### Section 5.1: Comparison to Existing Lasso Methods
### Evaluating Results, PAC vs. Lasso (both relaxed and unrelaxed)

## Create objects to store SSE values at each setting/correlation structure
## Note there are four error rates calculated per group (in order):
## PAC
## Lasso
## Relaxed PAC
## Relaxed Lasso

error.constrained.valid.5.corr.dying=rep(0,iter)
error.lars.valid.5.corr.dying=rep(0,iter)
error.constrained.valid.5.corr.dying.relaxed=rep(0,iter)
error.lars.valid.5.corr.dying.relaxed=rep(0,iter)

error.constrained.valid.5=rep(0,iter)
error.lars.valid.5=rep(0,iter)
error.constrained.valid.5.relaxed=rep(0,iter)
error.lars.valid.5.relaxed=rep(0,iter)

error.constrained.valid.4=rep(0,iter)
error.lars.valid.4=rep(0,iter)
error.constrained.valid.4.relaxed=rep(0,iter)
error.lars.valid.4.relaxed=rep(0,iter)

error.constrained.valid.4.corr.dying=rep(0,iter)
error.lars.valid.4.corr.dying=rep(0,iter)
error.constrained.valid.4.corr.dying.relaxed=rep(0,iter)
error.lars.valid.4.corr.dying.relaxed=rep(0,iter)

error.constrained.valid.7.corr.dying=rep(0,iter)
error.lars.valid.7.corr.dying=rep(0,iter)
error.constrained.valid.7.corr.dying.relaxed=rep(0,iter)
error.lars.valid.7.corr.dying.relaxed=rep(0,iter)

error.constrained.valid.7=rep(0,iter)
error.lars.valid.7=rep(0,iter)
error.constrained.valid.7.relaxed=rep(0,iter)
error.lars.valid.7.relaxed=rep(0,iter)

## Create lists for runs of each Lasso model
lars.coefs.5.corr.dying=vector("list",iter)
lars.coefs.5=vector("list",iter)
lars.coefs.4.corr.dying=vector("list",iter)
lars.coefs.4=vector("list",iter)
lars.coefs.7.corr.dying=vector("list",iter)
lars.coefs.7=vector("list",iter)


## For validation purposes, record index and lambda value of lowest error per run per method
valid.place=rep(0,iter)
valid.lam=rep(0,iter)
valid.place.lars=rep(0,iter)
valid.lam.lars=rep(0,iter)

valid.place.relaxed=rep(0,iter)
valid.lam.relaxed=rep(0,iter)
valid.place.lars.relaxed=rep(0,iter)
valid.lam.lars.relaxed=rep(0,iter)

## Setting 1, Correlation


n=100
p=50
m=5

#LASSO runs, saving coefficients
for (i in 1:iter) {


lars.trial=lars(fits.4.corr.dying[[i]]$data$x,fits.4.corr.dying[[i]]$data$y,normalize=T,intercept=F)
lars.coefs.4.corr.dying[[i]]=predict(lars.trial,s=fits.4.corr.dying[[i]]$constrained.fit$lambda,type="coefficients",mode="lambda")$coef

}

## Set relaxed coefficients to current coefficients (update later)

lars.coefs.4.corr.dying.relaxed=lars.coefs.4.corr.dying
fits.4.corr.dying.relaxed=fits.4.corr.dying

## Create the required number of validation/testing data sets
valid.set.x=vector("list",iter)
valid.set.y=matrix(0,n,iter)
test.set.x=vector("list",iter)
test.set.y=matrix(0,n.test,iter)

## Generate and store the required number of validation/testing data sets
for (i in 1:iter) {

 set.seed(seed[i]+1)
 valid.set.x[[i]]=matrix(rnorm(n*p),n,p)
 valid.set.y[,i]=as.vector(valid.set.x[[i]]%*%fits.4.corr.dying[[i]]$data$beta+sigma*rnorm(n))
 test.set.x[[i]]=matrix(rnorm(n.test*p),n.test,p)
 test.set.y[,i]=as.vector(test.set.x[[i]]%*%fits.4.corr.dying[[i]]$data$beta)
 
}

## Run relaxed PAC and relaxed Lasso
## Note the relaxed PAC requires a deconstructed formulation
## The X matrix (X.til) and Y vector (Y.til) are calculated by reformulating the PAC
## As shown in the PAC paper (e.g. Generalized Lasso and Appendix A)

for (i in 1:iter) {

b=fits.4.corr.dying[[i]]$data$b
y.data=fits.4.corr.dying[[i]]$data$y


n.lam=length(fits.4.corr.dying[[i]]$constrained.fit$lambda)
error=rep(0,n.lam)
error.lars=rep(0,n.lam)
error.relaxed=rep(0,n.lam)
error.lars.relaxed=rep(0,n.lam)

## Calculate relaxed versions of each model at every lambda value
for (k in 1:n.lam) {


beta2.index=fits.4.corr.dying[[i]]$constrained.fit$b2index[,k]
x2.data=fits.4.corr.dying[[i]]$data$x[,beta2.index]
A2=fits.4.corr.dying[[i]]$data$C.full[,beta2.index]
Y.til=y.data-x2.data%*%solve(A2)%*%b

total=c(1:p)

beta1.index.adj=NULL

## Identify non-zero vector entries
for (j in 1:p) {
   if (fits.4.corr.dying[[i]]$constrained.fit$coefs[,k][j]==0) {
   beta1.index.adj=c(beta1.index.adj,j)
   }
  }

beta1.index=total[-c(beta2.index,beta1.index.adj)]

if(length(beta1.index)==0) {

	fits.4.corr.dying.relaxed[[i]]$constrained.fit$coefs[,k]==fits.4.corr.dying[[i]]$constrained.fit$coefs[,k]

} else {

	x1.star=fits.4.corr.dying[[i]]$data$x[,beta1.index]
	A1.star=fits.4.corr.dying[[i]]$data$C.full[,beta1.index]
	X.til=x1.star-x2.data%*%solve(A2)%*%A1.star

	## Calculate relaxed PAC using lars function
	if (n<length(beta1.index)) {

   		lars.trial=lars(X.til,Y.til,normalize=T,intercept=F)
   		relax.fit=predict(lars.trial,s=fits.4.corr.dying[[i]]$constrained.fit$lambda,type="coefficients",mode="lambda")$coef
   		fits.4.corr.dying.relaxed[[i]]$constrained.fit$coefs[,k][beta1.index]=relax.fit[k,]
   		beta2=solve(A2)%*%(b-A1.star%*%relax.fit[k,])

	} else {
		relax.fit=lm(Y.til~X.til-1)
		if(length(beta1.index)==1) {
			beta2=solve(A2)%*%(b-A1.star*relax.fit$coef)
		} else {

			beta2=solve(A2)%*%(b-A1.star%*%relax.fit$coef)
		}
		fits.4.corr.dying.relaxed[[i]]$constrained.fit$coefs[,k][beta1.index]=relax.fit$coef
	}

	fits.4.corr.dying.relaxed[[i]]$constrained.fit$coefs[,k][beta2.index]=beta2

}

relax.index=NULL

## Identify non-zero coefficient values for Lasso run
for (j in 1:p) {

   if (lars.coefs.4.corr.dying[[i]][k,][j]!=0) {

	relax.index=c(relax.index,j)

   }

}

## Calculate relaxed Lasso fit
if(!is.null(relax.index)) {

x1.lars=fits.4.corr.dying[[i]]$data$x[,relax.index]
relax.fit=lm(y.data~x1.lars-1)
lars.coefs.4.corr.dying.relaxed[[i]][k,][relax.index]=relax.fit$coef

}

## Calculate SSE using pred.error() function for each of four methods at this lambda

error[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],fits.4.corr.dying[[i]]$constrained.fit$coefs[,k])
error.lars[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],lars.coefs.4.corr.dying[[i]][k,])
error.relaxed[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],fits.4.corr.dying.relaxed[[i]]$constrained.fit$coefs[,k])
error.lars.relaxed[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],lars.coefs.4.corr.dying.relaxed[[i]][k,])

}

## Identify validation minimum SSE and the lambda at which it occurs for all four methods
valid.place[i]=which.min(error)
valid.lam[i]=fits.4.corr.dying[[i]]$constrained.fit$lambda[valid.place[i]]
valid.place.lars[i]=which.min(error.lars)
valid.lam.lars[i]=fits.4.corr.dying[[i]]$constrained.fit$lambda[valid.place.lars[i]]
valid.place.relaxed[i]=which.min(error.relaxed)
valid.lam.relaxed[i]=fits.4.corr.dying[[i]]$constrained.fit$lambda[valid.place.relaxed[i]]
valid.place.lars.relaxed[i]=which.min(error.lars.relaxed)
valid.lam.lars.relaxed[i]=fits.4.corr.dying[[i]]$constrained.fit$lambda[valid.place.lars.relaxed[i]]

## Calculate RMSE using minimum validation values applied to testing data
error.constrained.valid.4.corr.dying[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],fits.4.corr.dying[[i]]$constrained.fit$coefs[,valid.place[i]])/10000)
error.lars.valid.4.corr.dying[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],lars.coefs.4.corr.dying[[i]][valid.place.lars[i],])/10000)
error.constrained.valid.4.corr.dying.relaxed[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],fits.4.corr.dying.relaxed[[i]]$constrained.fit$coefs[,valid.place.relaxed[i]])/10000)
error.lars.valid.4.corr.dying.relaxed[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],lars.coefs.4.corr.dying.relaxed[[i]][valid.place.lars.relaxed[i],])/10000)

}


## Setting 1, No Correlation

#LASSO
for (i in 1:iter) {

lars.trial=lars(fits.4[[i]]$data$x,fits.4[[i]]$data$y,normalize=T,intercept=F)
lars.coefs.4[[i]]=predict(lars.trial,s=fits.4[[i]]$constrained.fit$lambda,type="coefficients",mode="lambda")$coef

}

lars.coefs.4.relaxed=lars.coefs.4
fits.4.relaxed=fits.4


for (i in 1:iter) {

 set.seed(seed[i]+1)
 valid.set.x[[i]]=matrix(rnorm(n*p),n,p)
 valid.set.y[,i]=as.vector(valid.set.x[[i]]%*%fits.4[[i]]$data$beta+sigma*rnorm(n))
 test.set.x[[i]]=matrix(rnorm(n.test*p),n.test,p)
 test.set.y[,i]=as.vector(test.set.x[[i]]%*%fits.4[[i]]$data$beta)
 
}

for (i in 1:iter) {

b=fits.4[[i]]$data$b
y.data=fits.4[[i]]$data$y


n.lam=length(fits.4[[i]]$constrained.fit$lambda)
error=rep(0,n.lam)
error.lars=rep(0,n.lam)
error.relaxed=rep(0,n.lam)
error.lars.relaxed=rep(0,n.lam)


for (k in 1:n.lam) {


beta2.index=fits.4[[i]]$constrained.fit$b2index[,k]
x2.data=fits.4[[i]]$data$x[,beta2.index]
A2=fits.4[[i]]$data$C.full[,beta2.index]
Y.til=y.data-x2.data%*%solve(A2)%*%b

total=c(1:p)

beta1.index.adj=NULL

for (j in 1:p) {
   if (fits.4[[i]]$constrained.fit$coefs[,k][j]==0) {
   beta1.index.adj=c(beta1.index.adj,j)
   }
  }

beta1.index=total[-c(beta2.index,beta1.index.adj)]

if(length(beta1.index)==0) {

	fits.4.relaxed[[i]]$constrained.fit$coefs[,k]==fits.4[[i]]$constrained.fit$coefs[,k]

} else {

	x1.star=fits.4[[i]]$data$x[,beta1.index]
	A1.star=fits.4[[i]]$data$C.full[,beta1.index]
	X.til=x1.star-x2.data%*%solve(A2)%*%A1.star

	if (n<length(beta1.index)) {

   		lars.trial=lars(X.til,Y.til,normalize=T,intercept=F)
   		relax.fit=predict(lars.trial,s=fits.4[[i]]$constrained.fit$lambda,type="coefficients",mode="lambda")$coef
   		fits.4.relaxed[[i]]$constrained.fit$coefs[,k][beta1.index]=relax.fit[k,]
   		beta2=solve(A2)%*%(b-A1.star%*%relax.fit[k,])

	} else {
		relax.fit=lm(Y.til~X.til-1)
		if(length(beta1.index)==1) {
			beta2=solve(A2)%*%(b-A1.star*relax.fit$coef)
		} else {

			beta2=solve(A2)%*%(b-A1.star%*%relax.fit$coef)
		}
		fits.4.relaxed[[i]]$constrained.fit$coefs[,k][beta1.index]=relax.fit$coef
	}

	fits.4.relaxed[[i]]$constrained.fit$coefs[,k][beta2.index]=beta2

}

relax.index=NULL

for (j in 1:p) {

   if (lars.coefs.4[[i]][k,][j]!=0) {

	relax.index=c(relax.index,j)

   }

}

if(!is.null(relax.index)) {

x1.lars=fits.4[[i]]$data$x[,relax.index]
relax.fit=lm(y.data~x1.lars-1)
lars.coefs.4.relaxed[[i]][k,][relax.index]=relax.fit$coef

}

error[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],fits.4[[i]]$constrained.fit$coefs[,k])
error.lars[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],lars.coefs.4[[i]][k,])
error.relaxed[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],fits.4.relaxed[[i]]$constrained.fit$coefs[,k])
error.lars.relaxed[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],lars.coefs.4.relaxed[[i]][k,])

}

valid.place[i]=which.min(error)
valid.lam[i]=fits.4[[i]]$constrained.fit$lambda[valid.place[i]]
valid.place.lars[i]=which.min(error.lars)
valid.lam.lars[i]=fits.4[[i]]$constrained.fit$lambda[valid.place.lars[i]]
valid.place.relaxed[i]=which.min(error.relaxed)
valid.lam.relaxed[i]=fits.4[[i]]$constrained.fit$lambda[valid.place.relaxed[i]]
valid.place.lars.relaxed[i]=which.min(error.lars.relaxed)
valid.lam.lars.relaxed[i]=fits.4[[i]]$constrained.fit$lambda[valid.place.lars.relaxed[i]]


error.constrained.valid.4[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],fits.4[[i]]$constrained.fit$coefs[,valid.place[i]])/10000)
error.lars.valid.4[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],lars.coefs.4[[i]][valid.place.lars[i],])/10000)
error.constrained.valid.4.relaxed[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],fits.4.relaxed[[i]]$constrained.fit$coefs[,valid.place.relaxed[i]])/10000)
error.lars.valid.4.relaxed[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],lars.coefs.4.relaxed[[i]][valid.place.lars.relaxed[i],])/10000)

}



## Setting 2, Correlation


n=50
p=500
m=10
n.test=10000
s=5
sigma=1

#LASSO
for (i in 1:iter) {


lars.trial=lars(fits.5.corr.dying[[i]]$data$x,fits.5.corr.dying[[i]]$data$y,normalize=T,intercept=F)
lars.coefs.5.corr.dying[[i]]=predict(lars.trial,s=fits.5.corr.dying[[i]]$constrained.fit$lambda,type="coefficients",mode="lambda")$coef

}

lars.coefs.5.corr.dying.relaxed=lars.coefs.5.corr.dying
fits.5.corr.dying.relaxed=fits.5.corr.dying

valid.set.x=vector("list",iter)
valid.set.y=matrix(0,n,iter)
test.set.x=vector("list",iter)
test.set.y=matrix(0,n.test,iter)

for (i in 1:iter) {

 set.seed(seed[i]+1)
 valid.set.x[[i]]=matrix(rnorm(n*p),n,p)
 valid.set.y[,i]=as.vector(valid.set.x[[i]]%*%fits.5.corr.dying[[i]]$data$beta+sigma*rnorm(n))
 test.set.x[[i]]=matrix(rnorm(n.test*p),n.test,p)
 test.set.y[,i]=as.vector(test.set.x[[i]]%*%fits.5.corr.dying[[i]]$data$beta)
 
}

for (i in 1:iter) {

b=fits.5.corr.dying[[i]]$data$b
y.data=fits.5.corr.dying[[i]]$data$y


n.lam=length(fits.5.corr.dying[[i]]$constrained.fit$lambda)
error=rep(0,n.lam)
error.lars=rep(0,n.lam)
error.relaxed=rep(0,n.lam)
error.lars.relaxed=rep(0,n.lam)


for (k in 1:n.lam) {


beta2.index=fits.5.corr.dying[[i]]$constrained.fit$b2index[,k]
x2.data=fits.5.corr.dying[[i]]$data$x[,beta2.index]
A2=fits.5.corr.dying[[i]]$data$C.full[,beta2.index]
Y.til=y.data-x2.data%*%solve(A2)%*%b

total=c(1:p)

beta1.index.adj=NULL

for (j in 1:p) {
   if (fits.5.corr.dying[[i]]$constrained.fit$coefs[,k][j]==0) {
   beta1.index.adj=c(beta1.index.adj,j)
   }
  }

beta1.index=total[-c(beta2.index,beta1.index.adj)]

if(length(beta1.index)==0) {

	fits.5.corr.dying.relaxed[[i]]$constrained.fit$coefs[,k]==fits.5.corr.dying[[i]]$constrained.fit$coefs[,k]

} else {

	x1.star=fits.5.corr.dying[[i]]$data$x[,beta1.index]
	A1.star=fits.5.corr.dying[[i]]$data$C.full[,beta1.index]
	X.til=x1.star-x2.data%*%solve(A2)%*%A1.star

	if (n<length(beta1.index)) {

   		lars.trial=lars(X.til,Y.til,normalize=T,intercept=F)
   		relax.fit=predict(lars.trial,s=fits.5.corr.dying[[i]]$constrained.fit$lambda,type="coefficients",mode="lambda")$coef
   		fits.5.corr.dying.relaxed[[i]]$constrained.fit$coefs[,k][beta1.index]=relax.fit[k,]
   		beta2=solve(A2)%*%(b-A1.star%*%relax.fit[k,])

	} else {
		relax.fit=lm(Y.til~X.til-1)
		if(length(beta1.index)==1) {
			beta2=solve(A2)%*%(b-A1.star*relax.fit$coef)
		} else {

			beta2=solve(A2)%*%(b-A1.star%*%relax.fit$coef)
		}
		fits.5.corr.dying.relaxed[[i]]$constrained.fit$coefs[,k][beta1.index]=relax.fit$coef
	}

	fits.5.corr.dying.relaxed[[i]]$constrained.fit$coefs[,k][beta2.index]=beta2

}

relax.index=NULL

for (j in 1:p) {

   if (lars.coefs.5.corr.dying[[i]][k,][j]!=0) {

	relax.index=c(relax.index,j)

   }

}

if(!is.null(relax.index)) {

x1.lars=fits.5.corr.dying[[i]]$data$x[,relax.index]
relax.fit=lm(y.data~x1.lars-1)
lars.coefs.5.corr.dying.relaxed[[i]][k,][relax.index]=relax.fit$coef

}

error[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],fits.5.corr.dying[[i]]$constrained.fit$coefs[,k])
error.lars[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],lars.coefs.5.corr.dying[[i]][k,])
error.relaxed[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],fits.5.corr.dying.relaxed[[i]]$constrained.fit$coefs[,k])
error.lars.relaxed[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],lars.coefs.5.corr.dying.relaxed[[i]][k,])

}

valid.place[i]=which.min(error)
valid.lam[i]=fits.5.corr.dying[[i]]$constrained.fit$lambda[valid.place[i]]
valid.place.lars[i]=which.min(error.lars)
valid.lam.lars[i]=fits.5.corr.dying[[i]]$constrained.fit$lambda[valid.place.lars[i]]
valid.place.relaxed[i]=which.min(error.relaxed)
valid.lam.relaxed[i]=fits.5.corr.dying[[i]]$constrained.fit$lambda[valid.place.relaxed[i]]
valid.place.lars.relaxed[i]=which.min(error.lars.relaxed)
valid.lam.lars.relaxed[i]=fits.5.corr.dying[[i]]$constrained.fit$lambda[valid.place.lars.relaxed[i]]


error.constrained.valid.5.corr.dying[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],fits.5.corr.dying[[i]]$constrained.fit$coefs[,valid.place[i]])/10000)
error.lars.valid.5.corr.dying[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],lars.coefs.5.corr.dying[[i]][valid.place.lars[i],])/10000)
error.constrained.valid.5.corr.dying.relaxed[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],fits.5.corr.dying.relaxed[[i]]$constrained.fit$coefs[,valid.place.relaxed[i]])/10000)
error.lars.valid.5.corr.dying.relaxed[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],lars.coefs.5.corr.dying.relaxed[[i]][valid.place.lars.relaxed[i],])/10000)

}


## Setting 2, No Correlation

#LASSO
for (i in 1:iter) {

lars.trial=lars(fits.5[[i]]$data$x,fits.5[[i]]$data$y,normalize=T,intercept=F)
lars.coefs.5[[i]]=predict(lars.trial,s=fits.5[[i]]$constrained.fit$lambda,type="coefficients",mode="lambda")$coef

}

lars.coefs.5.relaxed=lars.coefs.5
fits.5.relaxed=fits.5


for (i in 1:iter) {

 set.seed(seed[i]+1)
 valid.set.x[[i]]=matrix(rnorm(n*p),n,p)
 valid.set.y[,i]=as.vector(valid.set.x[[i]]%*%fits.5[[i]]$data$beta+sigma*rnorm(n))
 test.set.x[[i]]=matrix(rnorm(n.test*p),n.test,p)
 test.set.y[,i]=as.vector(test.set.x[[i]]%*%fits.5[[i]]$data$beta)
 
}

for (i in 1:iter) {

b=fits.5[[i]]$data$b
y.data=fits.5[[i]]$data$y


n.lam=length(fits.5[[i]]$constrained.fit$lambda)
error=rep(0,n.lam)
error.lars=rep(0,n.lam)
error.relaxed=rep(0,n.lam)
error.lars.relaxed=rep(0,n.lam)


for (k in 1:n.lam) {


beta2.index=fits.5[[i]]$constrained.fit$b2index[,k]
x2.data=fits.5[[i]]$data$x[,beta2.index]
A2=fits.5[[i]]$data$C.full[,beta2.index]
Y.til=y.data-x2.data%*%solve(A2)%*%b

total=c(1:p)

beta1.index.adj=NULL

for (j in 1:p) {
   if (fits.5[[i]]$constrained.fit$coefs[,k][j]==0) {
   beta1.index.adj=c(beta1.index.adj,j)
   }
  }

beta1.index=total[-c(beta2.index,beta1.index.adj)]

if(length(beta1.index)==0) {

	fits.5.relaxed[[i]]$constrained.fit$coefs[,k]==fits.5[[i]]$constrained.fit$coefs[,k]

} else {

	x1.star=fits.5[[i]]$data$x[,beta1.index]
	A1.star=fits.5[[i]]$data$C.full[,beta1.index]
	X.til=x1.star-x2.data%*%solve(A2)%*%A1.star

	if (n<length(beta1.index)) {

   		lars.trial=lars(X.til,Y.til,normalize=T,intercept=F)
   		relax.fit=predict(lars.trial,s=fits.5[[i]]$constrained.fit$lambda,type="coefficients",mode="lambda")$coef
   		fits.5.relaxed[[i]]$constrained.fit$coefs[,k][beta1.index]=relax.fit[k,]
   		beta2=solve(A2)%*%(b-A1.star%*%relax.fit[k,])

	} else {
		relax.fit=lm(Y.til~X.til-1)
		if(length(beta1.index)==1) {
			beta2=solve(A2)%*%(b-A1.star*relax.fit$coef)
		} else {

			beta2=solve(A2)%*%(b-A1.star%*%relax.fit$coef)
		}
		fits.5.relaxed[[i]]$constrained.fit$coefs[,k][beta1.index]=relax.fit$coef
	}

	fits.5.relaxed[[i]]$constrained.fit$coefs[,k][beta2.index]=beta2

}

relax.index=NULL

for (j in 1:p) {

   if (lars.coefs.5[[i]][k,][j]!=0) {

	relax.index=c(relax.index,j)

   }

}

if(!is.null(relax.index)) {

x1.lars=fits.5[[i]]$data$x[,relax.index]
relax.fit=lm(y.data~x1.lars-1)
lars.coefs.5.relaxed[[i]][k,][relax.index]=relax.fit$coef

}

error[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],fits.5[[i]]$constrained.fit$coefs[,k])
error.lars[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],lars.coefs.5[[i]][k,])
error.relaxed[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],fits.5.relaxed[[i]]$constrained.fit$coefs[,k])
error.lars.relaxed[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],lars.coefs.5.relaxed[[i]][k,])

}

valid.place[i]=which.min(error)
valid.lam[i]=fits.5[[i]]$constrained.fit$lambda[valid.place[i]]
valid.place.lars[i]=which.min(error.lars)
valid.lam.lars[i]=fits.5[[i]]$constrained.fit$lambda[valid.place.lars[i]]
valid.place.relaxed[i]=which.min(error.relaxed)
valid.lam.relaxed[i]=fits.5[[i]]$constrained.fit$lambda[valid.place.relaxed[i]]
valid.place.lars.relaxed[i]=which.min(error.lars.relaxed)
valid.lam.lars.relaxed[i]=fits.5[[i]]$constrained.fit$lambda[valid.place.lars.relaxed[i]]


error.constrained.valid.5[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],fits.5[[i]]$constrained.fit$coefs[,valid.place[i]])/10000)
error.lars.valid.5[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],lars.coefs.5[[i]][valid.place.lars[i],])/10000)
error.constrained.valid.5.relaxed[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],fits.5.relaxed[[i]]$constrained.fit$coefs[,valid.place.relaxed[i]])/10000)
error.lars.valid.5.relaxed[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],lars.coefs.5.relaxed[[i]][valid.place.lars.relaxed[i],])/10000)

}




## Setting 3, Correlation


n=50
p=100
m=60

#LASSO
for (i in 1:iter) {


lars.trial=lars(fits.7.corr.dying[[i]]$data$x,fits.7.corr.dying[[i]]$data$y,normalize=T,intercept=F)
lars.coefs.7.corr.dying[[i]]=predict(lars.trial,s=fits.7.corr.dying[[i]]$constrained.fit$lambda,type="coefficients",mode="lambda")$coef

}

lars.coefs.7.corr.dying.relaxed=lars.coefs.7.corr.dying
fits.7.corr.dying.relaxed=fits.7.corr.dying

valid.set.x=vector("list",100)
valid.set.y=matrix(0,n,100)
test.set.x=vector("list",100)
test.set.y=matrix(0,n.test,100)

for (i in 1:iter) {

 set.seed(seed[i]+1)
 valid.set.x[[i]]=matrix(rnorm(n*p),n,p)
 valid.set.y[,i]=as.vector(valid.set.x[[i]]%*%fits.7.corr.dying[[i]]$data$beta+sigma*rnorm(n))
 test.set.x[[i]]=matrix(rnorm(n.test*p),n.test,p)
 test.set.y[,i]=as.vector(test.set.x[[i]]%*%fits.7.corr.dying[[i]]$data$beta)
 
}

for (i in 1:iter) {

b=fits.7.corr.dying[[i]]$data$b
y.data=fits.7.corr.dying[[i]]$data$y


n.lam=length(fits.7.corr.dying[[i]]$constrained.fit$lambda)
error=rep(0,n.lam)
error.lars=rep(0,n.lam)
error.relaxed=rep(0,n.lam)
error.lars.relaxed=rep(0,n.lam)


for (k in 1:n.lam) {


beta2.index=fits.7.corr.dying[[i]]$constrained.fit$b2index[,k]
x2.data=fits.7.corr.dying[[i]]$data$x[,beta2.index]
A2=fits.7.corr.dying[[i]]$data$C.full[,beta2.index]
Y.til=y.data-x2.data%*%solve(A2)%*%b

total=c(1:p)

beta1.index.adj=NULL

for (j in 1:p) {
   if (fits.7.corr.dying[[i]]$constrained.fit$coefs[,k][j]==0) {
   beta1.index.adj=c(beta1.index.adj,j)
   }
  }

beta1.index=total[-c(beta2.index,beta1.index.adj)]

if(length(beta1.index)==0) {

	fits.7.corr.dying.relaxed[[i]]$constrained.fit$coefs[,k]==fits.7.corr.dying[[i]]$constrained.fit$coefs[,k]

} else {

	x1.star=fits.7.corr.dying[[i]]$data$x[,beta1.index]
	A1.star=fits.7.corr.dying[[i]]$data$C.full[,beta1.index]
	X.til=x1.star-x2.data%*%solve(A2)%*%A1.star

	if (n<length(beta1.index)) {

   		lars.trial=lars(X.til,Y.til,normalize=T,intercept=F)
   		relax.fit=predict(lars.trial,s=fits.7.corr.dying[[i]]$constrained.fit$lambda,type="coefficients",mode="lambda")$coef
   		fits.7.corr.dying.relaxed[[i]]$constrained.fit$coefs[,k][beta1.index]=relax.fit[k,]
   		beta2=solve(A2)%*%(b-A1.star%*%relax.fit[k,])

	} else {
		relax.fit=lm(Y.til~X.til-1)
		if(length(beta1.index)==1) {
			beta2=solve(A2)%*%(b-A1.star*relax.fit$coef)
		} else {

			beta2=solve(A2)%*%(b-A1.star%*%relax.fit$coef)
		}
		fits.7.corr.dying.relaxed[[i]]$constrained.fit$coefs[,k][beta1.index]=relax.fit$coef
	}

	fits.7.corr.dying.relaxed[[i]]$constrained.fit$coefs[,k][beta2.index]=beta2

}

relax.index=NULL

for (j in 1:p) {

   if (lars.coefs.7.corr.dying[[i]][k,][j]!=0) {

	relax.index=c(relax.index,j)

   }

}

if(!is.null(relax.index)) {

x1.lars=fits.7.corr.dying[[i]]$data$x[,relax.index]
relax.fit=lm(y.data~x1.lars-1)
lars.coefs.7.corr.dying.relaxed[[i]][k,][relax.index]=relax.fit$coef

}

error[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],fits.7.corr.dying[[i]]$constrained.fit$coefs[,k])
error.lars[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],lars.coefs.7.corr.dying[[i]][k,])
error.relaxed[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],fits.7.corr.dying.relaxed[[i]]$constrained.fit$coefs[,k])
error.lars.relaxed[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],lars.coefs.7.corr.dying.relaxed[[i]][k,])

}

valid.place[i]=which.min(error)
valid.lam[i]=fits.7.corr.dying[[i]]$constrained.fit$lambda[valid.place[i]]
valid.place.lars[i]=which.min(error.lars)
valid.lam.lars[i]=fits.7.corr.dying[[i]]$constrained.fit$lambda[valid.place.lars[i]]
valid.place.relaxed[i]=which.min(error.relaxed)
valid.lam.relaxed[i]=fits.7.corr.dying[[i]]$constrained.fit$lambda[valid.place.relaxed[i]]
valid.place.lars.relaxed[i]=which.min(error.lars.relaxed)
valid.lam.lars.relaxed[i]=fits.7.corr.dying[[i]]$constrained.fit$lambda[valid.place.lars.relaxed[i]]


error.constrained.valid.7.corr.dying[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],fits.7.corr.dying[[i]]$constrained.fit$coefs[,valid.place[i]])/10000)
error.lars.valid.7.corr.dying[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],lars.coefs.7.corr.dying[[i]][valid.place.lars[i],])/10000)
error.constrained.valid.7.corr.dying.relaxed[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],fits.7.corr.dying.relaxed[[i]]$constrained.fit$coefs[,valid.place.relaxed[i]])/10000)
error.lars.valid.7.corr.dying.relaxed[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],lars.coefs.7.corr.dying.relaxed[[i]][valid.place.lars.relaxed[i],])/10000)

}


## Setting 3, No Correlation

#LASSO
for (i in 1:iter) {

lars.trial=lars(fits.7[[i]]$data$x,fits.7[[i]]$data$y,normalize=T,intercept=F)
lars.coefs.7[[i]]=predict(lars.trial,s=fits.7[[i]]$constrained.fit$lambda,type="coefficients",mode="lambda")$coef

}

lars.coefs.7.relaxed=lars.coefs.7
fits.7.relaxed=fits.7


for (i in 1:iter) {

 set.seed(seed[i]+1)
 valid.set.x[[i]]=matrix(rnorm(n*p),n,p)
 valid.set.y[,i]=as.vector(valid.set.x[[i]]%*%fits.7[[i]]$data$beta+sigma*rnorm(n))
 test.set.x[[i]]=matrix(rnorm(n.test*p),n.test,p)
 test.set.y[,i]=as.vector(test.set.x[[i]]%*%fits.7[[i]]$data$beta)
 
}

for (i in 1:iter) {

b=fits.7[[i]]$data$b
y.data=fits.7[[i]]$data$y


n.lam=length(fits.7[[i]]$constrained.fit$lambda)
error=rep(0,n.lam)
error.lars=rep(0,n.lam)
error.relaxed=rep(0,n.lam)
error.lars.relaxed=rep(0,n.lam)


for (k in 1:n.lam) {


beta2.index=fits.7[[i]]$constrained.fit$b2index[,k]
x2.data=fits.7[[i]]$data$x[,beta2.index]
A2=fits.7[[i]]$data$C.full[,beta2.index]
Y.til=y.data-x2.data%*%solve(A2)%*%b

total=c(1:p)

beta1.index.adj=NULL

for (j in 1:p) {
   if (fits.7[[i]]$constrained.fit$coefs[,k][j]==0) {
   beta1.index.adj=c(beta1.index.adj,j)
   }
  }

beta1.index=total[-c(beta2.index,beta1.index.adj)]

if(length(beta1.index)==0) {

	fits.7.relaxed[[i]]$constrained.fit$coefs[,k]==fits.7[[i]]$constrained.fit$coefs[,k]

} else {

	x1.star=fits.7[[i]]$data$x[,beta1.index]
	A1.star=fits.7[[i]]$data$C.full[,beta1.index]
	X.til=x1.star-x2.data%*%solve(A2)%*%A1.star

	if (n<length(beta1.index)) {

   		lars.trial=lars(X.til,Y.til,normalize=T,intercept=F)
   		relax.fit=predict(lars.trial,s=fits.7[[i]]$constrained.fit$lambda,type="coefficients",mode="lambda")$coef
   		fits.7.relaxed[[i]]$constrained.fit$coefs[,k][beta1.index]=relax.fit[k,]
   		beta2=solve(A2)%*%(b-A1.star%*%relax.fit[k,])

	} else {
		relax.fit=lm(Y.til~X.til-1)
		if(length(beta1.index)==1) {
			beta2=solve(A2)%*%(b-A1.star*relax.fit$coef)
		} else {

			beta2=solve(A2)%*%(b-A1.star%*%relax.fit$coef)
		}
		fits.7.relaxed[[i]]$constrained.fit$coefs[,k][beta1.index]=relax.fit$coef
	}

	fits.7.relaxed[[i]]$constrained.fit$coefs[,k][beta2.index]=beta2

}

relax.index=NULL

for (j in 1:p) {

   if (lars.coefs.7[[i]][k,][j]!=0) {

	relax.index=c(relax.index,j)

   }

}

if(!is.null(relax.index)) {

x1.lars=fits.7[[i]]$data$x[,relax.index]
relax.fit=lm(y.data~x1.lars-1)
lars.coefs.7.relaxed[[i]][k,][relax.index]=relax.fit$coef

}

error[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],fits.7[[i]]$constrained.fit$coefs[,k])
error.lars[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],lars.coefs.7[[i]][k,])
error.relaxed[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],fits.7.relaxed[[i]]$constrained.fit$coefs[,k])
error.lars.relaxed[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],lars.coefs.7.relaxed[[i]][k,])

}

valid.place[i]=which.min(error)
valid.lam[i]=fits.7[[i]]$constrained.fit$lambda[valid.place[i]]
valid.place.lars[i]=which.min(error.lars)
valid.lam.lars[i]=fits.7[[i]]$constrained.fit$lambda[valid.place.lars[i]]
valid.place.relaxed[i]=which.min(error.relaxed)
valid.lam.relaxed[i]=fits.7[[i]]$constrained.fit$lambda[valid.place.relaxed[i]]
valid.place.lars.relaxed[i]=which.min(error.lars.relaxed)
valid.lam.lars.relaxed[i]=fits.7[[i]]$constrained.fit$lambda[valid.place.lars.relaxed[i]]


error.constrained.valid.7[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],fits.7[[i]]$constrained.fit$coefs[,valid.place[i]])/10000)
error.lars.valid.7[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],lars.coefs.7[[i]][valid.place.lars[i],])/10000)
error.constrained.valid.7.relaxed[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],fits.7.relaxed[[i]]$constrained.fit$coefs[,valid.place.relaxed[i]])/10000)
error.lars.valid.7.relaxed[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],lars.coefs.7.relaxed[[i]][valid.place.lars.relaxed[i],])/10000)

}

###################################################################################
###################    Calculating Error Rates     ################################
###################           (Table 2)            ################################
###################################################################################

## Table 2 Results

# Avg RMSE and SE, Setting 1, No Correlation
mean(error.lars.valid.4)
sd(error.lars.valid.4)/sqrt(iter)
mean(error.constrained.valid.4)
sd(error.constrained.valid.4)/sqrt(iter)
mean(error.lars.valid.4.relaxed)
sd(error.lars.valid.4.relaxed)/sqrt(iter)
mean(error.constrained.valid.4.relaxed)
sd(error.constrained.valid.4.relaxed)/sqrt(iter)

# Avg RMSE and SE, Setting 1, Correlation
mean(error.lars.valid.4.corr.dying)
sd(error.lars.valid.4.corr.dying)/sqrt(iter)
mean(error.constrained.valid.4.corr.dying)
sd(error.constrained.valid.4.corr.dying)/sqrt(iter)
mean(error.lars.valid.4.corr.dying.relaxed)
sd(error.lars.valid.4.corr.dying.relaxed)/sqrt(iter)
mean(error.constrained.valid.4.corr.dying.relaxed)
sd(error.constrained.valid.4.corr.dying.relaxed)/sqrt(iter)

# Avg RMSE and SE, Setting 2, No Correlation
mean(error.lars.valid.5)
sd(error.lars.valid.5)/sqrt(iter)
mean(error.constrained.valid.5)
sd(error.constrained.valid.5)/sqrt(iter)
mean(error.lars.valid.5.relaxed)
sd(error.lars.valid.5.relaxed)/sqrt(iter)
mean(error.constrained.valid.5.relaxed)
sd(error.constrained.valid.5.relaxed)/sqrt(iter)

# Avg RMSE and SE, Setting 2, Correlation
mean(error.lars.valid.5.corr.dying)
sd(error.lars.valid.5.corr.dying)/sqrt(iter)
mean(error.constrained.valid.5.corr.dying)
sd(error.constrained.valid.5.corr.dying)/sqrt(iter)
mean(error.lars.valid.5.corr.dying.relaxed)
sd(error.lars.valid.5.corr.dying.relaxed)/sqrt(iter)
mean(error.constrained.valid.5.corr.dying.relaxed)
sd(error.constrained.valid.5.corr.dying.relaxed)/sqrt(iter)

# Avg RMSE and SE, Setting 3, No Correlation
mean(error.lars.valid.7)
sd(error.lars.valid.7)/sqrt(iter)
mean(error.constrained.valid.7)
sd(error.constrained.valid.7)/sqrt(iter)
mean(error.lars.valid.7.relaxed)
sd(error.lars.valid.7.relaxed)/sqrt(iter)
mean(error.constrained.valid.7.relaxed)
sd(error.constrained.valid.7.relaxed)/sqrt(iter)

# Avg RMSE and SE, Setting 3, Correlation
mean(error.lars.valid.7.corr.dying)
sd(error.lars.valid.7.corr.dying)/sqrt(iter)
mean(error.constrained.valid.7.corr.dying)
sd(error.constrained.valid.7.corr.dying)/sqrt(iter)
mean(error.lars.valid.7.corr.dying.relaxed)
sd(error.lars.valid.7.corr.dying.relaxed)/sqrt(iter)
mean(error.constrained.valid.7.corr.dying.relaxed)
sd(error.constrained.valid.7.corr.dying.relaxed)/sqrt(iter)



###################################################################################
###################    Calculating Error Rates     ################################
###################       and Relaxed Fits         ################################
###################           (Table 3)            ################################
###################################################################################


### Section 5.2: Violations of Constraints
### Evaluating Results: PAC vs. Lasso (both relaxed and unrelaxed)

## Note this code follows from above (just repeated with error included)

error.constrained.valid.test1.25=rep(0,iter)
error.lars.valid.test1.25=rep(0,iter)
error.constrained.valid.test1.25.relaxed=rep(0,iter)
error.lars.valid.test1.25.relaxed=rep(0,iter)
error.constrained.valid.test1.5=rep(0,iter)
error.lars.valid.test1.5=rep(0,iter)
error.constrained.valid.test1.5.relaxed=rep(0,iter)
error.lars.valid.test1.5.relaxed=rep(0,iter)
error.constrained.valid.test1.75=rep(0,iter)
error.lars.valid.test1.75=rep(0,iter)
error.constrained.valid.test1.75.relaxed=rep(0,iter)
error.lars.valid.test1.75.relaxed=rep(0,iter)
error.constrained.valid.test11=rep(0,iter)
error.lars.valid.test11=rep(0,iter)
error.constrained.valid.test11.relaxed=rep(0,iter)
error.lars.valid.test11.relaxed=rep(0,iter)

error.constrained.valid.test1.25=rep(0,iter)
error.lars.valid.test1.25=rep(0,iter)
error.constrained.valid.test1.25.relaxed=rep(0,iter)
error.lars.valid.test1.25.relaxed=rep(0,iter)
error.constrained.valid.test1.5=rep(0,iter)
error.lars.valid.test1.5=rep(0,iter)
error.constrained.valid.test1.5.relaxed=rep(0,iter)
error.lars.valid.test1.5.relaxed=rep(0,iter)
error.constrained.valid.test1.75=rep(0,iter)
error.lars.valid.test1.75=rep(0,iter)
error.constrained.valid.test1.75.relaxed=rep(0,iter)
error.lars.valid.test1.75.relaxed=rep(0,iter)
error.constrained.valid.test11=rep(0,iter)
error.lars.valid.test11=rep(0,iter)
error.constrained.valid.test11.relaxed=rep(0,iter)
error.lars.valid.test11.relaxed=rep(0,iter)

error.constrained.valid.test2.25=rep(0,iter)
error.lars.valid.test2.25=rep(0,iter)
error.constrained.valid.test2.25.relaxed=rep(0,iter)
error.lars.valid.test2.25.relaxed=rep(0,iter)
error.constrained.valid.test2.5=rep(0,iter)
error.lars.valid.test2.5=rep(0,iter)
error.constrained.valid.test2.5.relaxed=rep(0,iter)
error.lars.valid.test2.5.relaxed=rep(0,iter)
error.constrained.valid.test2.75=rep(0,iter)
error.lars.valid.test2.75=rep(0,iter)
error.constrained.valid.test2.75.relaxed=rep(0,iter)
error.lars.valid.test2.75.relaxed=rep(0,iter)
error.constrained.valid.test21=rep(0,iter)
error.lars.valid.test21=rep(0,iter)
error.constrained.valid.test21.relaxed=rep(0,iter)
error.lars.valid.test21.relaxed=rep(0,iter)

error.constrained.valid.test3.25=rep(0,iter)
error.lars.valid.test3.25=rep(0,iter)
error.constrained.valid.test3.25.relaxed=rep(0,iter)
error.lars.valid.test3.25.relaxed=rep(0,iter)
error.constrained.valid.test3.5=rep(0,iter)
error.lars.valid.test3.5=rep(0,iter)
error.constrained.valid.test3.5.relaxed=rep(0,iter)
error.lars.valid.test3.5.relaxed=rep(0,iter)
error.constrained.valid.test3.75=rep(0,iter)
error.lars.valid.test3.75=rep(0,iter)
error.constrained.valid.test3.75.relaxed=rep(0,iter)
error.lars.valid.test3.75.relaxed=rep(0,iter)
error.constrained.valid.test31=rep(0,iter)
error.lars.valid.test31=rep(0,iter)
error.constrained.valid.test31.relaxed=rep(0,iter)
error.lars.valid.test31.relaxed=rep(0,iter)

lars.coefs.test1.25=vector("list",iter)
lars.coefs.test1.5=vector("list",iter)
lars.coefs.test1.75=vector("list",iter)
lars.coefs.test11=vector("list",iter)
lars.coefs.test2.25=vector("list",iter)
lars.coefs.test2.5=vector("list",iter)
lars.coefs.test2.75=vector("list",iter)
lars.coefs.test21=vector("list",iter)
lars.coefs.test3.25=vector("list",iter)
lars.coefs.test3.5=vector("list",iter)
lars.coefs.test3.75=vector("list",iter)
lars.coefs.test31=vector("list",iter)



## Setting 1

n=100
p=50
m=5

## a=0.25

#LASSO
for (i in 1:iter) {


lars.trial=lars(fits.test1.25[[i]]$data$x,fits.test1.25[[i]]$data$y,normalize=T,intercept=F)
lars.coefs.test1.25[[i]]=predict(lars.trial,s=fits.test1.25[[i]]$constrained.fit$lambda,type="coefficients",mode="lambda")$coef

}

lars.coefs.test1.25.relaxed=lars.coefs.test1.25
fits.test1.25.relaxed=fits.test1.25

valid.set.x=vector("list",iter)
valid.set.y=matrix(0,n,iter)
test.set.x=vector("list",iter)
test.set.y=matrix(0,n.test,iter)

for (i in 1:iter) {

 set.seed(seed[i]+1)
 valid.set.x[[i]]=matrix(rnorm(n*p),n,p)
 valid.set.y[,i]=as.vector(valid.set.x[[i]]%*%fits.test1.25[[i]]$data$beta+sigma*rnorm(n))
 test.set.x[[i]]=matrix(rnorm(n.test*p),n.test,p)
 test.set.y[,i]=as.vector(test.set.x[[i]]%*%fits.test1.25[[i]]$data$beta)
 
}

for (i in 1:iter) {

b=fits.test1.25[[i]]$data$b
y.data=fits.test1.25[[i]]$data$y


n.lam=length(fits.test1.25[[i]]$constrained.fit$lambda)
error=rep(0,n.lam)
error.lars=rep(0,n.lam)
error.relaxed=rep(0,n.lam)
error.lars.relaxed=rep(0,n.lam)


for (k in 1:n.lam) {


beta2.index=fits.test1.25[[i]]$constrained.fit$b2index[,k]
x2.data=fits.test1.25[[i]]$data$x[,beta2.index]
A2=fits.test1.25[[i]]$data$C.full[,beta2.index]
Y.til=y.data-x2.data%*%solve(A2)%*%b

total=c(1:p)

beta1.index.adj=NULL

for (j in 1:p) {
   if (fits.test1.25[[i]]$constrained.fit$coefs[,k][j]==0) {
   beta1.index.adj=c(beta1.index.adj,j)
   }
  }

beta1.index=total[-c(beta2.index,beta1.index.adj)]

if(length(beta1.index)==0) {

	fits.test1.25.relaxed[[i]]$constrained.fit$coefs[,k]==fits.test1.25[[i]]$constrained.fit$coefs[,k]

} else {

	x1.star=fits.test1.25[[i]]$data$x[,beta1.index]
	A1.star=fits.test1.25[[i]]$data$C.full[,beta1.index]
	X.til=x1.star-x2.data%*%solve(A2)%*%A1.star

	if (n<length(beta1.index)) {

   		lars.trial=lars(X.til,Y.til,normalize=T,intercept=F)
   		relax.fit=predict(lars.trial,s=fits.test1.25[[i]]$constrained.fit$lambda,type="coefficients",mode="lambda")$coef
   		fits.test1.25.relaxed[[i]]$constrained.fit$coefs[,k][beta1.index]=relax.fit[k,]
   		beta2=solve(A2)%*%(b-A1.star%*%relax.fit[k,])

	} else {
		relax.fit=lm(Y.til~X.til-1)
		if(length(beta1.index)==1) {
			beta2=solve(A2)%*%(b-A1.star*relax.fit$coef)
		} else {

			beta2=solve(A2)%*%(b-A1.star%*%relax.fit$coef)
		}
		fits.test1.25.relaxed[[i]]$constrained.fit$coefs[,k][beta1.index]=relax.fit$coef
	}

	fits.test1.25.relaxed[[i]]$constrained.fit$coefs[,k][beta2.index]=beta2

}

relax.index=NULL

for (j in 1:p) {

   if (lars.coefs.test1.25[[i]][k,][j]!=0) {

	relax.index=c(relax.index,j)

   }

}

if(!is.null(relax.index)) {

x1.lars=fits.test1.25[[i]]$data$x[,relax.index]
relax.fit=lm(y.data~x1.lars-1)
lars.coefs.test1.25.relaxed[[i]][k,][relax.index]=relax.fit$coef

}

error[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],fits.test1.25[[i]]$constrained.fit$coefs[,k])
error.lars[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],lars.coefs.test1.25[[i]][k,])
error.relaxed[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],fits.test1.25.relaxed[[i]]$constrained.fit$coefs[,k])
error.lars.relaxed[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],lars.coefs.test1.25.relaxed[[i]][k,])

}

valid.place[i]=which.min(error)
valid.lam[i]=fits.test1.25[[i]]$constrained.fit$lambda[valid.place[i]]
valid.place.lars[i]=which.min(error.lars)
valid.lam.lars[i]=fits.test1.25[[i]]$constrained.fit$lambda[valid.place.lars[i]]
valid.place.relaxed[i]=which.min(error.relaxed)
valid.lam.relaxed[i]=fits.test1.25[[i]]$constrained.fit$lambda[valid.place.relaxed[i]]
valid.place.lars.relaxed[i]=which.min(error.lars.relaxed)
valid.lam.lars.relaxed[i]=fits.test1.25[[i]]$constrained.fit$lambda[valid.place.lars.relaxed[i]]


error.constrained.valid.test1.25[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],fits.test1.25[[i]]$constrained.fit$coefs[,valid.place[i]])/10000)
error.lars.valid.test1.25[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],lars.coefs.test1.25[[i]][valid.place.lars[i],])/10000)
error.constrained.valid.test1.25.relaxed[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],fits.test1.25.relaxed[[i]]$constrained.fit$coefs[,valid.place.relaxed[i]])/10000)
error.lars.valid.test1.25.relaxed[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],lars.coefs.test1.25.relaxed[[i]][valid.place.lars.relaxed[i],])/10000)

}

## a=0.5


#LASSO
for (i in 1:iter) {


lars.trial=lars(fits.test1.5[[i]]$data$x,fits.test1.5[[i]]$data$y,normalize=T,intercept=F)
lars.coefs.test1.5[[i]]=predict(lars.trial,s=fits.test1.5[[i]]$constrained.fit$lambda,type="coefficients",mode="lambda")$coef

}

lars.coefs.test1.5.relaxed=lars.coefs.test1.5
fits.test1.5.relaxed=fits.test1.5

valid.set.x=vector("list",iter)
valid.set.y=matrix(0,n,iter)
test.set.x=vector("list",iter)
test.set.y=matrix(0,n.test,iter)

for (i in 1:iter) {

 set.seed(seed[i]+1)
 valid.set.x[[i]]=matrix(rnorm(n*p),n,p)
 valid.set.y[,i]=as.vector(valid.set.x[[i]]%*%fits.test1.5[[i]]$data$beta+sigma*rnorm(n))
 test.set.x[[i]]=matrix(rnorm(n.test*p),n.test,p)
 test.set.y[,i]=as.vector(test.set.x[[i]]%*%fits.test1.5[[i]]$data$beta)
 
}

for (i in 1:iter) {

b=fits.test1.5[[i]]$data$b
y.data=fits.test1.5[[i]]$data$y


n.lam=length(fits.test1.5[[i]]$constrained.fit$lambda)
error=rep(0,n.lam)
error.lars=rep(0,n.lam)
error.relaxed=rep(0,n.lam)
error.lars.relaxed=rep(0,n.lam)


for (k in 1:n.lam) {


beta2.index=fits.test1.5[[i]]$constrained.fit$b2index[,k]
x2.data=fits.test1.5[[i]]$data$x[,beta2.index]
A2=fits.test1.5[[i]]$data$C.full[,beta2.index]
Y.til=y.data-x2.data%*%solve(A2)%*%b

total=c(1:p)

beta1.index.adj=NULL

for (j in 1:p) {
   if (fits.test1.5[[i]]$constrained.fit$coefs[,k][j]==0) {
   beta1.index.adj=c(beta1.index.adj,j)
   }
  }

beta1.index=total[-c(beta2.index,beta1.index.adj)]

if(length(beta1.index)==0) {

	fits.test1.5.relaxed[[i]]$constrained.fit$coefs[,k]==fits.test1.5[[i]]$constrained.fit$coefs[,k]

} else {

	x1.star=fits.test1.5[[i]]$data$x[,beta1.index]
	A1.star=fits.test1.5[[i]]$data$C.full[,beta1.index]
	X.til=x1.star-x2.data%*%solve(A2)%*%A1.star

	if (n<length(beta1.index)) {

   		lars.trial=lars(X.til,Y.til,normalize=T,intercept=F)
   		relax.fit=predict(lars.trial,s=fits.test1.5[[i]]$constrained.fit$lambda,type="coefficients",mode="lambda")$coef
   		fits.test1.5.relaxed[[i]]$constrained.fit$coefs[,k][beta1.index]=relax.fit[k,]
   		beta2=solve(A2)%*%(b-A1.star%*%relax.fit[k,])

	} else {
		relax.fit=lm(Y.til~X.til-1)
		if(length(beta1.index)==1) {
			beta2=solve(A2)%*%(b-A1.star*relax.fit$coef)
		} else {

			beta2=solve(A2)%*%(b-A1.star%*%relax.fit$coef)
		}
		fits.test1.5.relaxed[[i]]$constrained.fit$coefs[,k][beta1.index]=relax.fit$coef
	}

	fits.test1.5.relaxed[[i]]$constrained.fit$coefs[,k][beta2.index]=beta2

}

relax.index=NULL

for (j in 1:p) {

   if (lars.coefs.test1.5[[i]][k,][j]!=0) {

	relax.index=c(relax.index,j)

   }

}

if(!is.null(relax.index)) {

x1.lars=fits.test1.5[[i]]$data$x[,relax.index]
relax.fit=lm(y.data~x1.lars-1)
lars.coefs.test1.5.relaxed[[i]][k,][relax.index]=relax.fit$coef

}

error[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],fits.test1.5[[i]]$constrained.fit$coefs[,k])
error.lars[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],lars.coefs.test1.5[[i]][k,])
error.relaxed[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],fits.test1.5.relaxed[[i]]$constrained.fit$coefs[,k])
error.lars.relaxed[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],lars.coefs.test1.5.relaxed[[i]][k,])

}

valid.place[i]=which.min(error)
valid.lam[i]=fits.test1.5[[i]]$constrained.fit$lambda[valid.place[i]]
valid.place.lars[i]=which.min(error.lars)
valid.lam.lars[i]=fits.test1.5[[i]]$constrained.fit$lambda[valid.place.lars[i]]
valid.place.relaxed[i]=which.min(error.relaxed)
valid.lam.relaxed[i]=fits.test1.5[[i]]$constrained.fit$lambda[valid.place.relaxed[i]]
valid.place.lars.relaxed[i]=which.min(error.lars.relaxed)
valid.lam.lars.relaxed[i]=fits.test1.5[[i]]$constrained.fit$lambda[valid.place.lars.relaxed[i]]


error.constrained.valid.test1.5[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],fits.test1.5[[i]]$constrained.fit$coefs[,valid.place[i]])/10000)
error.lars.valid.test1.5[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],lars.coefs.test1.5[[i]][valid.place.lars[i],])/10000)
error.constrained.valid.test1.5.relaxed[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],fits.test1.5.relaxed[[i]]$constrained.fit$coefs[,valid.place.relaxed[i]])/10000)
error.lars.valid.test1.5.relaxed[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],lars.coefs.test1.5.relaxed[[i]][valid.place.lars.relaxed[i],])/10000)

}

#a=0.75


#LASSO
for (i in 1:iter) {


lars.trial=lars(fits.test1.75[[i]]$data$x,fits.test1.75[[i]]$data$y,normalize=T,intercept=F)
lars.coefs.test1.75[[i]]=predict(lars.trial,s=fits.test1.75[[i]]$constrained.fit$lambda,type="coefficients",mode="lambda")$coef

}

lars.coefs.test1.75.relaxed=lars.coefs.test1.75
fits.test1.75.relaxed=fits.test1.75

valid.set.x=vector("list",iter)
valid.set.y=matrix(0,n,iter)
test.set.x=vector("list",iter)
test.set.y=matrix(0,n.test,iter)

for (i in 1:iter) {

 set.seed(seed[i]+1)
 valid.set.x[[i]]=matrix(rnorm(n*p),n,p)
 valid.set.y[,i]=as.vector(valid.set.x[[i]]%*%fits.test1.75[[i]]$data$beta+sigma*rnorm(n))
 test.set.x[[i]]=matrix(rnorm(n.test*p),n.test,p)
 test.set.y[,i]=as.vector(test.set.x[[i]]%*%fits.test1.75[[i]]$data$beta)
 
}

for (i in 1:iter) {

b=fits.test1.75[[i]]$data$b
y.data=fits.test1.75[[i]]$data$y


n.lam=length(fits.test1.75[[i]]$constrained.fit$lambda)
error=rep(0,n.lam)
error.lars=rep(0,n.lam)
error.relaxed=rep(0,n.lam)
error.lars.relaxed=rep(0,n.lam)


for (k in 1:n.lam) {


beta2.index=fits.test1.75[[i]]$constrained.fit$b2index[,k]
x2.data=fits.test1.75[[i]]$data$x[,beta2.index]
A2=fits.test1.75[[i]]$data$C.full[,beta2.index]
Y.til=y.data-x2.data%*%solve(A2)%*%b

total=c(1:p)

beta1.index.adj=NULL

for (j in 1:p) {
   if (fits.test1.75[[i]]$constrained.fit$coefs[,k][j]==0) {
   beta1.index.adj=c(beta1.index.adj,j)
   }
  }

beta1.index=total[-c(beta2.index,beta1.index.adj)]

if(length(beta1.index)==0) {

	fits.test1.75.relaxed[[i]]$constrained.fit$coefs[,k]==fits.test1.75[[i]]$constrained.fit$coefs[,k]

} else {

	x1.star=fits.test1.75[[i]]$data$x[,beta1.index]
	A1.star=fits.test1.75[[i]]$data$C.full[,beta1.index]
	X.til=x1.star-x2.data%*%solve(A2)%*%A1.star

	if (n<length(beta1.index)) {

   		lars.trial=lars(X.til,Y.til,normalize=T,intercept=F)
   		relax.fit=predict(lars.trial,s=fits.test1.75[[i]]$constrained.fit$lambda,type="coefficients",mode="lambda")$coef
   		fits.test1.75.relaxed[[i]]$constrained.fit$coefs[,k][beta1.index]=relax.fit[k,]
   		beta2=solve(A2)%*%(b-A1.star%*%relax.fit[k,])

	} else {
		relax.fit=lm(Y.til~X.til-1)
		if(length(beta1.index)==1) {
			beta2=solve(A2)%*%(b-A1.star*relax.fit$coef)
		} else {

			beta2=solve(A2)%*%(b-A1.star%*%relax.fit$coef)
		}
		fits.test1.75.relaxed[[i]]$constrained.fit$coefs[,k][beta1.index]=relax.fit$coef
	}

	fits.test1.75.relaxed[[i]]$constrained.fit$coefs[,k][beta2.index]=beta2

}

relax.index=NULL

for (j in 1:p) {

   if (lars.coefs.test1.75[[i]][k,][j]!=0) {

	relax.index=c(relax.index,j)

   }

}

if(!is.null(relax.index)) {

x1.lars=fits.test1.75[[i]]$data$x[,relax.index]
relax.fit=lm(y.data~x1.lars-1)
lars.coefs.test1.75.relaxed[[i]][k,][relax.index]=relax.fit$coef

}

error[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],fits.test1.75[[i]]$constrained.fit$coefs[,k])
error.lars[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],lars.coefs.test1.75[[i]][k,])
error.relaxed[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],fits.test1.75.relaxed[[i]]$constrained.fit$coefs[,k])
error.lars.relaxed[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],lars.coefs.test1.75.relaxed[[i]][k,])

}

valid.place[i]=which.min(error)
valid.lam[i]=fits.test1.75[[i]]$constrained.fit$lambda[valid.place[i]]
valid.place.lars[i]=which.min(error.lars)
valid.lam.lars[i]=fits.test1.75[[i]]$constrained.fit$lambda[valid.place.lars[i]]
valid.place.relaxed[i]=which.min(error.relaxed)
valid.lam.relaxed[i]=fits.test1.75[[i]]$constrained.fit$lambda[valid.place.relaxed[i]]
valid.place.lars.relaxed[i]=which.min(error.lars.relaxed)
valid.lam.lars.relaxed[i]=fits.test1.75[[i]]$constrained.fit$lambda[valid.place.lars.relaxed[i]]


error.constrained.valid.test1.75[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],fits.test1.75[[i]]$constrained.fit$coefs[,valid.place[i]])/10000)
error.lars.valid.test1.75[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],lars.coefs.test1.75[[i]][valid.place.lars[i],])/10000)
error.constrained.valid.test1.75.relaxed[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],fits.test1.75.relaxed[[i]]$constrained.fit$coefs[,valid.place.relaxed[i]])/10000)
error.lars.valid.test1.75.relaxed[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],lars.coefs.test1.75.relaxed[[i]][valid.place.lars.relaxed[i],])/10000)

}

#a=1


#LASSO
for (i in 1:iter) {


lars.trial=lars(fits.test11[[i]]$data$x,fits.test11[[i]]$data$y,normalize=T,intercept=F)
lars.coefs.test11[[i]]=predict(lars.trial,s=fits.test11[[i]]$constrained.fit$lambda,type="coefficients",mode="lambda")$coef

}

lars.coefs.test11.relaxed=lars.coefs.test11
fits.test11.relaxed=fits.test11

valid.set.x=vector("list",iter)
valid.set.y=matrix(0,n,iter)
test.set.x=vector("list",iter)
test.set.y=matrix(0,n.test,iter)

for (i in 1:iter) {

 set.seed(seed[i]+1)
 valid.set.x[[i]]=matrix(rnorm(n*p),n,p)
 valid.set.y[,i]=as.vector(valid.set.x[[i]]%*%fits.test11[[i]]$data$beta+sigma*rnorm(n))
 test.set.x[[i]]=matrix(rnorm(n.test*p),n.test,p)
 test.set.y[,i]=as.vector(test.set.x[[i]]%*%fits.test11[[i]]$data$beta)
 
}

for (i in 1:iter) {

b=fits.test11[[i]]$data$b
y.data=fits.test11[[i]]$data$y


n.lam=length(fits.test11[[i]]$constrained.fit$lambda)
error=rep(0,n.lam)
error.lars=rep(0,n.lam)
error.relaxed=rep(0,n.lam)
error.lars.relaxed=rep(0,n.lam)


for (k in 1:n.lam) {


beta2.index=fits.test11[[i]]$constrained.fit$b2index[,k]
x2.data=fits.test11[[i]]$data$x[,beta2.index]
A2=fits.test11[[i]]$data$C.full[,beta2.index]
Y.til=y.data-x2.data%*%solve(A2)%*%b

total=c(1:p)

beta1.index.adj=NULL

for (j in 1:p) {
   if (fits.test11[[i]]$constrained.fit$coefs[,k][j]==0) {
   beta1.index.adj=c(beta1.index.adj,j)
   }
  }

beta1.index=total[-c(beta2.index,beta1.index.adj)]

if(length(beta1.index)==0) {

	fits.test11.relaxed[[i]]$constrained.fit$coefs[,k]==fits.test11[[i]]$constrained.fit$coefs[,k]

} else {

	x1.star=fits.test11[[i]]$data$x[,beta1.index]
	A1.star=fits.test11[[i]]$data$C.full[,beta1.index]
	X.til=x1.star-x2.data%*%solve(A2)%*%A1.star

	if (n<length(beta1.index)) {

   		lars.trial=lars(X.til,Y.til,normalize=T,intercept=F)
   		relax.fit=predict(lars.trial,s=fits.test11[[i]]$constrained.fit$lambda,type="coefficients",mode="lambda")$coef
   		fits.test11.relaxed[[i]]$constrained.fit$coefs[,k][beta1.index]=relax.fit[k,]
   		beta2=solve(A2)%*%(b-A1.star%*%relax.fit[k,])

	} else {
		relax.fit=lm(Y.til~X.til-1)
		if(length(beta1.index)==1) {
			beta2=solve(A2)%*%(b-A1.star*relax.fit$coef)
		} else {

			beta2=solve(A2)%*%(b-A1.star%*%relax.fit$coef)
		}
		fits.test11.relaxed[[i]]$constrained.fit$coefs[,k][beta1.index]=relax.fit$coef
	}

	fits.test11.relaxed[[i]]$constrained.fit$coefs[,k][beta2.index]=beta2

}

relax.index=NULL

for (j in 1:p) {

   if (lars.coefs.test11[[i]][k,][j]!=0) {

	relax.index=c(relax.index,j)

   }

}

if(!is.null(relax.index)) {

x1.lars=fits.test11[[i]]$data$x[,relax.index]
relax.fit=lm(y.data~x1.lars-1)
lars.coefs.test11.relaxed[[i]][k,][relax.index]=relax.fit$coef

}

error[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],fits.test11[[i]]$constrained.fit$coefs[,k])
error.lars[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],lars.coefs.test11[[i]][k,])
error.relaxed[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],fits.test11.relaxed[[i]]$constrained.fit$coefs[,k])
error.lars.relaxed[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],lars.coefs.test11.relaxed[[i]][k,])

}

valid.place[i]=which.min(error)
valid.lam[i]=fits.test11[[i]]$constrained.fit$lambda[valid.place[i]]
valid.place.lars[i]=which.min(error.lars)
valid.lam.lars[i]=fits.test11[[i]]$constrained.fit$lambda[valid.place.lars[i]]
valid.place.relaxed[i]=which.min(error.relaxed)
valid.lam.relaxed[i]=fits.test11[[i]]$constrained.fit$lambda[valid.place.relaxed[i]]
valid.place.lars.relaxed[i]=which.min(error.lars.relaxed)
valid.lam.lars.relaxed[i]=fits.test11[[i]]$constrained.fit$lambda[valid.place.lars.relaxed[i]]


error.constrained.valid.test11[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],fits.test11[[i]]$constrained.fit$coefs[,valid.place[i]])/10000)
error.lars.valid.test11[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],lars.coefs.test11[[i]][valid.place.lars[i],])/10000)
error.constrained.valid.test11.relaxed[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],fits.test11.relaxed[[i]]$constrained.fit$coefs[,valid.place.relaxed[i]])/10000)
error.lars.valid.test11.relaxed[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],lars.coefs.test11.relaxed[[i]][valid.place.lars.relaxed[i],])/10000)

}



### Setting 2

n=50
p=500
m=10

## a=0.25


#LASSO
for (i in 1:iter) {


lars.trial=lars(fits.test2.25[[i]]$data$x,fits.test2.25[[i]]$data$y,normalize=T,intercept=F)
lars.coefs.test2.25[[i]]=predict(lars.trial,s=fits.test2.25[[i]]$constrained.fit$lambda,type="coefficients",mode="lambda")$coef

}

lars.coefs.test2.25.relaxed=lars.coefs.test2.25
fits.test2.25.relaxed=fits.test2.25

valid.set.x=vector("list",iter)
valid.set.y=matrix(0,n,iter)
test.set.x=vector("list",iter)
test.set.y=matrix(0,n.test,iter)

for (i in 1:iter) {

 set.seed(seed[i]+1)
 valid.set.x[[i]]=matrix(rnorm(n*p),n,p)
 valid.set.y[,i]=as.vector(valid.set.x[[i]]%*%fits.test2.25[[i]]$data$beta+sigma*rnorm(n))
 test.set.x[[i]]=matrix(rnorm(n.test*p),n.test,p)
 test.set.y[,i]=as.vector(test.set.x[[i]]%*%fits.test2.25[[i]]$data$beta)
 
}

for (i in 1:iter) {

b=fits.test2.25[[i]]$data$b
y.data=fits.test2.25[[i]]$data$y


n.lam=length(fits.test2.25[[i]]$constrained.fit$lambda)
error=rep(0,n.lam)
error.lars=rep(0,n.lam)
error.relaxed=rep(0,n.lam)
error.lars.relaxed=rep(0,n.lam)


for (k in 1:n.lam) {


beta2.index=fits.test2.25[[i]]$constrained.fit$b2index[,k]
x2.data=fits.test2.25[[i]]$data$x[,beta2.index]
A2=fits.test2.25[[i]]$data$C.full[,beta2.index]
Y.til=y.data-x2.data%*%solve(A2)%*%b

total=c(1:p)

beta1.index.adj=NULL

for (j in 1:p) {
   if (fits.test2.25[[i]]$constrained.fit$coefs[,k][j]==0) {
   beta1.index.adj=c(beta1.index.adj,j)
   }
  }

beta1.index=total[-c(beta2.index,beta1.index.adj)]

if(length(beta1.index)==0) {

	fits.test2.25.relaxed[[i]]$constrained.fit$coefs[,k]==fits.test2.25[[i]]$constrained.fit$coefs[,k]

} else {

	x1.star=fits.test2.25[[i]]$data$x[,beta1.index]
	A1.star=fits.test2.25[[i]]$data$C.full[,beta1.index]
	X.til=x1.star-x2.data%*%solve(A2)%*%A1.star

	if (n<length(beta1.index)) {

   		lars.trial=lars(X.til,Y.til,normalize=T,intercept=F)
   		relax.fit=predict(lars.trial,s=fits.test2.25[[i]]$constrained.fit$lambda,type="coefficients",mode="lambda")$coef
   		fits.test2.25.relaxed[[i]]$constrained.fit$coefs[,k][beta1.index]=relax.fit[k,]
   		beta2=solve(A2)%*%(b-A1.star%*%relax.fit[k,])

	} else {
		relax.fit=lm(Y.til~X.til-1)
		if(length(beta1.index)==1) {
			beta2=solve(A2)%*%(b-A1.star*relax.fit$coef)
		} else {

			beta2=solve(A2)%*%(b-A1.star%*%relax.fit$coef)
		}
		fits.test2.25.relaxed[[i]]$constrained.fit$coefs[,k][beta1.index]=relax.fit$coef
	}

	fits.test2.25.relaxed[[i]]$constrained.fit$coefs[,k][beta2.index]=beta2

}

relax.index=NULL

for (j in 1:p) {

   if (lars.coefs.test2.25[[i]][k,][j]!=0) {

	relax.index=c(relax.index,j)

   }

}

if(!is.null(relax.index)) {

x1.lars=fits.test2.25[[i]]$data$x[,relax.index]
relax.fit=lm(y.data~x1.lars-1)
lars.coefs.test2.25.relaxed[[i]][k,][relax.index]=relax.fit$coef

}

error[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],fits.test2.25[[i]]$constrained.fit$coefs[,k])
error.lars[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],lars.coefs.test2.25[[i]][k,])
error.relaxed[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],fits.test2.25.relaxed[[i]]$constrained.fit$coefs[,k])
error.lars.relaxed[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],lars.coefs.test2.25.relaxed[[i]][k,])

}

valid.place[i]=which.min(error)
valid.lam[i]=fits.test2.25[[i]]$constrained.fit$lambda[valid.place[i]]
valid.place.lars[i]=which.min(error.lars)
valid.lam.lars[i]=fits.test2.25[[i]]$constrained.fit$lambda[valid.place.lars[i]]
valid.place.relaxed[i]=which.min(error.relaxed)
valid.lam.relaxed[i]=fits.test2.25[[i]]$constrained.fit$lambda[valid.place.relaxed[i]]
valid.place.lars.relaxed[i]=which.min(error.lars.relaxed)
valid.lam.lars.relaxed[i]=fits.test2.25[[i]]$constrained.fit$lambda[valid.place.lars.relaxed[i]]


error.constrained.valid.test2.25[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],fits.test2.25[[i]]$constrained.fit$coefs[,valid.place[i]])/10000)
error.lars.valid.test2.25[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],lars.coefs.test2.25[[i]][valid.place.lars[i],])/10000)
error.constrained.valid.test2.25.relaxed[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],fits.test2.25.relaxed[[i]]$constrained.fit$coefs[,valid.place.relaxed[i]])/10000)
error.lars.valid.test2.25.relaxed[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],lars.coefs.test2.25.relaxed[[i]][valid.place.lars.relaxed[i],])/10000)

}


#a=0.5


#LASSO
for (i in 1:iter) {


lars.trial=lars(fits.test2.5[[i]]$data$x,fits.test2.5[[i]]$data$y,normalize=T,intercept=F)
lars.coefs.test2.5[[i]]=predict(lars.trial,s=fits.test2.5[[i]]$constrained.fit$lambda,type="coefficients",mode="lambda")$coef

}

lars.coefs.test2.5.relaxed=lars.coefs.test2.5
fits.test2.5.relaxed=fits.test2.5

valid.set.x=vector("list",iter)
valid.set.y=matrix(0,n,iter)
test.set.x=vector("list",iter)
test.set.y=matrix(0,n.test,iter)

for (i in 1:iter) {

 set.seed(seed[i]+1)
 valid.set.x[[i]]=matrix(rnorm(n*p),n,p)
 valid.set.y[,i]=as.vector(valid.set.x[[i]]%*%fits.test2.5[[i]]$data$beta+sigma*rnorm(n))
 test.set.x[[i]]=matrix(rnorm(n.test*p),n.test,p)
 test.set.y[,i]=as.vector(test.set.x[[i]]%*%fits.test2.5[[i]]$data$beta)
 
}

for (i in 1:iter) {

b=fits.test2.5[[i]]$data$b
y.data=fits.test2.5[[i]]$data$y


n.lam=length(fits.test2.5[[i]]$constrained.fit$lambda)
error=rep(0,n.lam)
error.lars=rep(0,n.lam)
error.relaxed=rep(0,n.lam)
error.lars.relaxed=rep(0,n.lam)


for (k in 1:n.lam) {


beta2.index=fits.test2.5[[i]]$constrained.fit$b2index[,k]
x2.data=fits.test2.5[[i]]$data$x[,beta2.index]
A2=fits.test2.5[[i]]$data$C.full[,beta2.index]
Y.til=y.data-x2.data%*%solve(A2)%*%b

total=c(1:p)

beta1.index.adj=NULL

for (j in 1:p) {
   if (fits.test2.5[[i]]$constrained.fit$coefs[,k][j]==0) {
   beta1.index.adj=c(beta1.index.adj,j)
   }
  }

beta1.index=total[-c(beta2.index,beta1.index.adj)]

if(length(beta1.index)==0) {

	fits.test2.5.relaxed[[i]]$constrained.fit$coefs[,k]==fits.test2.5[[i]]$constrained.fit$coefs[,k]

} else {

	x1.star=fits.test2.5[[i]]$data$x[,beta1.index]
	A1.star=fits.test2.5[[i]]$data$C.full[,beta1.index]
	X.til=x1.star-x2.data%*%solve(A2)%*%A1.star

	if (n<length(beta1.index)) {

   		lars.trial=lars(X.til,Y.til,normalize=T,intercept=F)
   		relax.fit=predict(lars.trial,s=fits.test2.5[[i]]$constrained.fit$lambda,type="coefficients",mode="lambda")$coef
   		fits.test2.5.relaxed[[i]]$constrained.fit$coefs[,k][beta1.index]=relax.fit[k,]
   		beta2=solve(A2)%*%(b-A1.star%*%relax.fit[k,])

	} else {
		relax.fit=lm(Y.til~X.til-1)
		if(length(beta1.index)==1) {
			beta2=solve(A2)%*%(b-A1.star*relax.fit$coef)
		} else {

			beta2=solve(A2)%*%(b-A1.star%*%relax.fit$coef)
		}
		fits.test2.5.relaxed[[i]]$constrained.fit$coefs[,k][beta1.index]=relax.fit$coef
	}

	fits.test2.5.relaxed[[i]]$constrained.fit$coefs[,k][beta2.index]=beta2

}

relax.index=NULL

for (j in 1:p) {

   if (lars.coefs.test2.5[[i]][k,][j]!=0) {

	relax.index=c(relax.index,j)

   }

}

if(!is.null(relax.index)) {

x1.lars=fits.test2.5[[i]]$data$x[,relax.index]
relax.fit=lm(y.data~x1.lars-1)
lars.coefs.test2.5.relaxed[[i]][k,][relax.index]=relax.fit$coef

}

error[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],fits.test2.5[[i]]$constrained.fit$coefs[,k])
error.lars[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],lars.coefs.test2.5[[i]][k,])
error.relaxed[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],fits.test2.5.relaxed[[i]]$constrained.fit$coefs[,k])
error.lars.relaxed[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],lars.coefs.test2.5.relaxed[[i]][k,])

}

valid.place[i]=which.min(error)
valid.lam[i]=fits.test2.5[[i]]$constrained.fit$lambda[valid.place[i]]
valid.place.lars[i]=which.min(error.lars)
valid.lam.lars[i]=fits.test2.5[[i]]$constrained.fit$lambda[valid.place.lars[i]]
valid.place.relaxed[i]=which.min(error.relaxed)
valid.lam.relaxed[i]=fits.test2.5[[i]]$constrained.fit$lambda[valid.place.relaxed[i]]
valid.place.lars.relaxed[i]=which.min(error.lars.relaxed)
valid.lam.lars.relaxed[i]=fits.test2.5[[i]]$constrained.fit$lambda[valid.place.lars.relaxed[i]]


error.constrained.valid.test2.5[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],fits.test2.5[[i]]$constrained.fit$coefs[,valid.place[i]])/10000)
error.lars.valid.test2.5[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],lars.coefs.test2.5[[i]][valid.place.lars[i],])/10000)
error.constrained.valid.test2.5.relaxed[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],fits.test2.5.relaxed[[i]]$constrained.fit$coefs[,valid.place.relaxed[i]])/10000)
error.lars.valid.test2.5.relaxed[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],lars.coefs.test2.5.relaxed[[i]][valid.place.lars.relaxed[i],])/10000)

}


#a=0.75


#LASSO
for (i in 1:iter) {


lars.trial=lars(fits.test2.75[[i]]$data$x,fits.test2.75[[i]]$data$y,normalize=T,intercept=F)
lars.coefs.test2.75[[i]]=predict(lars.trial,s=fits.test2.75[[i]]$constrained.fit$lambda,type="coefficients",mode="lambda")$coef

}

lars.coefs.test2.75.relaxed=lars.coefs.test2.75
fits.test2.75.relaxed=fits.test2.75

valid.set.x=vector("list",iter)
valid.set.y=matrix(0,n,iter)
test.set.x=vector("list",iter)
test.set.y=matrix(0,n.test,iter)

for (i in 1:iter) {

 set.seed(seed[i]+1)
 valid.set.x[[i]]=matrix(rnorm(n*p),n,p)
 valid.set.y[,i]=as.vector(valid.set.x[[i]]%*%fits.test2.75[[i]]$data$beta+sigma*rnorm(n))
 test.set.x[[i]]=matrix(rnorm(n.test*p),n.test,p)
 test.set.y[,i]=as.vector(test.set.x[[i]]%*%fits.test2.75[[i]]$data$beta)
 
}

for (i in 1:iter) {

b=fits.test2.75[[i]]$data$b
y.data=fits.test2.75[[i]]$data$y


n.lam=length(fits.test2.75[[i]]$constrained.fit$lambda)
error=rep(0,n.lam)
error.lars=rep(0,n.lam)
error.relaxed=rep(0,n.lam)
error.lars.relaxed=rep(0,n.lam)


for (k in 1:n.lam) {


beta2.index=fits.test2.75[[i]]$constrained.fit$b2index[,k]
x2.data=fits.test2.75[[i]]$data$x[,beta2.index]
A2=fits.test2.75[[i]]$data$C.full[,beta2.index]
Y.til=y.data-x2.data%*%solve(A2)%*%b

total=c(1:p)

beta1.index.adj=NULL

for (j in 1:p) {
   if (fits.test2.75[[i]]$constrained.fit$coefs[,k][j]==0) {
   beta1.index.adj=c(beta1.index.adj,j)
   }
  }

beta1.index=total[-c(beta2.index,beta1.index.adj)]

if(length(beta1.index)==0) {

	fits.test2.75.relaxed[[i]]$constrained.fit$coefs[,k]==fits.test2.75[[i]]$constrained.fit$coefs[,k]

} else {

	x1.star=fits.test2.75[[i]]$data$x[,beta1.index]
	A1.star=fits.test2.75[[i]]$data$C.full[,beta1.index]
	X.til=x1.star-x2.data%*%solve(A2)%*%A1.star

	if (n<length(beta1.index)) {

   		lars.trial=lars(X.til,Y.til,normalize=T,intercept=F)
   		relax.fit=predict(lars.trial,s=fits.test2.75[[i]]$constrained.fit$lambda,type="coefficients",mode="lambda")$coef
   		fits.test2.75.relaxed[[i]]$constrained.fit$coefs[,k][beta1.index]=relax.fit[k,]
   		beta2=solve(A2)%*%(b-A1.star%*%relax.fit[k,])

	} else {
		relax.fit=lm(Y.til~X.til-1)
		if(length(beta1.index)==1) {
			beta2=solve(A2)%*%(b-A1.star*relax.fit$coef)
		} else {

			beta2=solve(A2)%*%(b-A1.star%*%relax.fit$coef)
		}
		fits.test2.75.relaxed[[i]]$constrained.fit$coefs[,k][beta1.index]=relax.fit$coef
	}

	fits.test2.75.relaxed[[i]]$constrained.fit$coefs[,k][beta2.index]=beta2

}

relax.index=NULL

for (j in 1:p) {

   if (lars.coefs.test2.75[[i]][k,][j]!=0) {

	relax.index=c(relax.index,j)

   }

}

if(!is.null(relax.index)) {

x1.lars=fits.test2.75[[i]]$data$x[,relax.index]
relax.fit=lm(y.data~x1.lars-1)
lars.coefs.test2.75.relaxed[[i]][k,][relax.index]=relax.fit$coef

}

error[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],fits.test2.75[[i]]$constrained.fit$coefs[,k])
error.lars[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],lars.coefs.test2.75[[i]][k,])
error.relaxed[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],fits.test2.75.relaxed[[i]]$constrained.fit$coefs[,k])
error.lars.relaxed[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],lars.coefs.test2.75.relaxed[[i]][k,])

}

valid.place[i]=which.min(error)
valid.lam[i]=fits.test2.75[[i]]$constrained.fit$lambda[valid.place[i]]
valid.place.lars[i]=which.min(error.lars)
valid.lam.lars[i]=fits.test2.75[[i]]$constrained.fit$lambda[valid.place.lars[i]]
valid.place.relaxed[i]=which.min(error.relaxed)
valid.lam.relaxed[i]=fits.test2.75[[i]]$constrained.fit$lambda[valid.place.relaxed[i]]
valid.place.lars.relaxed[i]=which.min(error.lars.relaxed)
valid.lam.lars.relaxed[i]=fits.test2.75[[i]]$constrained.fit$lambda[valid.place.lars.relaxed[i]]


error.constrained.valid.test2.75[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],fits.test2.75[[i]]$constrained.fit$coefs[,valid.place[i]])/10000)
error.lars.valid.test2.75[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],lars.coefs.test2.75[[i]][valid.place.lars[i],])/10000)
error.constrained.valid.test2.75.relaxed[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],fits.test2.75.relaxed[[i]]$constrained.fit$coefs[,valid.place.relaxed[i]])/10000)
error.lars.valid.test2.75.relaxed[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],lars.coefs.test2.75.relaxed[[i]][valid.place.lars.relaxed[i],])/10000)

}


#a=1


#LASSO
for (i in 1:iter) {


lars.trial=lars(fits.test21[[i]]$data$x,fits.test21[[i]]$data$y,normalize=T,intercept=F)
lars.coefs.test21[[i]]=predict(lars.trial,s=fits.test21[[i]]$constrained.fit$lambda,type="coefficients",mode="lambda")$coef

}

lars.coefs.test21.relaxed=lars.coefs.test21
fits.test21.relaxed=fits.test21

valid.set.x=vector("list",iter)
valid.set.y=matrix(0,n,iter)
test.set.x=vector("list",iter)
test.set.y=matrix(0,n.test,iter)

for (i in 1:iter) {

 set.seed(seed[i]+1)
 valid.set.x[[i]]=matrix(rnorm(n*p),n,p)
 valid.set.y[,i]=as.vector(valid.set.x[[i]]%*%fits.test21[[i]]$data$beta+sigma*rnorm(n))
 test.set.x[[i]]=matrix(rnorm(n.test*p),n.test,p)
 test.set.y[,i]=as.vector(test.set.x[[i]]%*%fits.test21[[i]]$data$beta)
 
}

for (i in 1:iter) {

b=fits.test21[[i]]$data$b
y.data=fits.test21[[i]]$data$y


n.lam=length(fits.test21[[i]]$constrained.fit$lambda)
error=rep(0,n.lam)
error.lars=rep(0,n.lam)
error.relaxed=rep(0,n.lam)
error.lars.relaxed=rep(0,n.lam)


for (k in 1:n.lam) {


beta2.index=fits.test21[[i]]$constrained.fit$b2index[,k]
x2.data=fits.test21[[i]]$data$x[,beta2.index]
A2=fits.test21[[i]]$data$C.full[,beta2.index]
Y.til=y.data-x2.data%*%solve(A2)%*%b

total=c(1:p)

beta1.index.adj=NULL

for (j in 1:p) {
   if (fits.test21[[i]]$constrained.fit$coefs[,k][j]==0) {
   beta1.index.adj=c(beta1.index.adj,j)
   }
  }

beta1.index=total[-c(beta2.index,beta1.index.adj)]

if(length(beta1.index)==0) {

	fits.test21.relaxed[[i]]$constrained.fit$coefs[,k]==fits.test21[[i]]$constrained.fit$coefs[,k]

} else {

	x1.star=fits.test21[[i]]$data$x[,beta1.index]
	A1.star=fits.test21[[i]]$data$C.full[,beta1.index]
	X.til=x1.star-x2.data%*%solve(A2)%*%A1.star

	if (n<length(beta1.index)) {

   		lars.trial=lars(X.til,Y.til,normalize=T,intercept=F)
   		relax.fit=predict(lars.trial,s=fits.test21[[i]]$constrained.fit$lambda,type="coefficients",mode="lambda")$coef
   		fits.test21.relaxed[[i]]$constrained.fit$coefs[,k][beta1.index]=relax.fit[k,]
   		beta2=solve(A2)%*%(b-A1.star%*%relax.fit[k,])

	} else {
		relax.fit=lm(Y.til~X.til-1)
		if(length(beta1.index)==1) {
			beta2=solve(A2)%*%(b-A1.star*relax.fit$coef)
		} else {

			beta2=solve(A2)%*%(b-A1.star%*%relax.fit$coef)
		}
		fits.test21.relaxed[[i]]$constrained.fit$coefs[,k][beta1.index]=relax.fit$coef
	}

	fits.test21.relaxed[[i]]$constrained.fit$coefs[,k][beta2.index]=beta2

}

relax.index=NULL

for (j in 1:p) {

   if (lars.coefs.test21[[i]][k,][j]!=0) {

	relax.index=c(relax.index,j)

   }

}

if(!is.null(relax.index)) {

x1.lars=fits.test21[[i]]$data$x[,relax.index]
relax.fit=lm(y.data~x1.lars-1)
lars.coefs.test21.relaxed[[i]][k,][relax.index]=relax.fit$coef

}

error[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],fits.test21[[i]]$constrained.fit$coefs[,k])
error.lars[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],lars.coefs.test21[[i]][k,])
error.relaxed[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],fits.test21.relaxed[[i]]$constrained.fit$coefs[,k])
error.lars.relaxed[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],lars.coefs.test21.relaxed[[i]][k,])

}

valid.place[i]=which.min(error)
valid.lam[i]=fits.test21[[i]]$constrained.fit$lambda[valid.place[i]]
valid.place.lars[i]=which.min(error.lars)
valid.lam.lars[i]=fits.test21[[i]]$constrained.fit$lambda[valid.place.lars[i]]
valid.place.relaxed[i]=which.min(error.relaxed)
valid.lam.relaxed[i]=fits.test21[[i]]$constrained.fit$lambda[valid.place.relaxed[i]]
valid.place.lars.relaxed[i]=which.min(error.lars.relaxed)
valid.lam.lars.relaxed[i]=fits.test21[[i]]$constrained.fit$lambda[valid.place.lars.relaxed[i]]


error.constrained.valid.test21[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],fits.test21[[i]]$constrained.fit$coefs[,valid.place[i]])/10000)
error.lars.valid.test21[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],lars.coefs.test21[[i]][valid.place.lars[i],])/10000)
error.constrained.valid.test21.relaxed[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],fits.test21.relaxed[[i]]$constrained.fit$coefs[,valid.place.relaxed[i]])/10000)
error.lars.valid.test21.relaxed[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],lars.coefs.test21.relaxed[[i]][valid.place.lars.relaxed[i],])/10000)

}


### Setting 3

n=50
p=100
m=60

#a=0.25


#LASSO
for (i in 1:iter) {


lars.trial=lars(fits.test3.25[[i]]$data$x,fits.test3.25[[i]]$data$y,normalize=T,intercept=F)
lars.coefs.test3.25[[i]]=predict(lars.trial,s=fits.test3.25[[i]]$constrained.fit$lambda,type="coefficients",mode="lambda")$coef

}

lars.coefs.test3.25.relaxed=lars.coefs.test3.25
fits.test3.25.relaxed=fits.test3.25

valid.set.x=vector("list",iter)
valid.set.y=matrix(0,n,iter)
test.set.x=vector("list",iter)
test.set.y=matrix(0,n.test,iter)

for (i in 1:iter) {

 set.seed(seed[i]+1)
 valid.set.x[[i]]=matrix(rnorm(n*p),n,p)
 valid.set.y[,i]=as.vector(valid.set.x[[i]]%*%fits.test3.25[[i]]$data$beta+sigma*rnorm(n))
 test.set.x[[i]]=matrix(rnorm(n.test*p),n.test,p)
 test.set.y[,i]=as.vector(test.set.x[[i]]%*%fits.test3.25[[i]]$data$beta)
 
}

for (i in 1:iter) {

b=fits.test3.25[[i]]$data$b
y.data=fits.test3.25[[i]]$data$y


n.lam=length(fits.test3.25[[i]]$constrained.fit$lambda)
error=rep(0,n.lam)
error.lars=rep(0,n.lam)
error.relaxed=rep(0,n.lam)
error.lars.relaxed=rep(0,n.lam)


for (k in 1:n.lam) {


beta2.index=fits.test3.25[[i]]$constrained.fit$b2index[,k]
x2.data=fits.test3.25[[i]]$data$x[,beta2.index]
A2=fits.test3.25[[i]]$data$C.full[,beta2.index]
Y.til=y.data-x2.data%*%solve(A2)%*%b

total=c(1:p)

beta1.index.adj=NULL

for (j in 1:p) {
   if (fits.test3.25[[i]]$constrained.fit$coefs[,k][j]==0) {
   beta1.index.adj=c(beta1.index.adj,j)
   }
  }

beta1.index=total[-c(beta2.index,beta1.index.adj)]

if(length(beta1.index)==0) {

	fits.test3.25.relaxed[[i]]$constrained.fit$coefs[,k]==fits.test3.25[[i]]$constrained.fit$coefs[,k]

} else {

	x1.star=fits.test3.25[[i]]$data$x[,beta1.index]
	A1.star=fits.test3.25[[i]]$data$C.full[,beta1.index]
	X.til=x1.star-x2.data%*%solve(A2)%*%A1.star

	if (n<length(beta1.index)) {

   		lars.trial=lars(X.til,Y.til,normalize=T,intercept=F)
   		relax.fit=predict(lars.trial,s=fits.test3.25[[i]]$constrained.fit$lambda,type="coefficients",mode="lambda")$coef
   		fits.test3.25.relaxed[[i]]$constrained.fit$coefs[,k][beta1.index]=relax.fit[k,]
   		beta2=solve(A2)%*%(b-A1.star%*%relax.fit[k,])

	} else {
		relax.fit=lm(Y.til~X.til-1)
		if(length(beta1.index)==1) {
			beta2=solve(A2)%*%(b-A1.star*relax.fit$coef)
		} else {

			beta2=solve(A2)%*%(b-A1.star%*%relax.fit$coef)
		}
		fits.test3.25.relaxed[[i]]$constrained.fit$coefs[,k][beta1.index]=relax.fit$coef
	}

	fits.test3.25.relaxed[[i]]$constrained.fit$coefs[,k][beta2.index]=beta2

}

relax.index=NULL

for (j in 1:p) {

   if (lars.coefs.test3.25[[i]][k,][j]!=0) {

	relax.index=c(relax.index,j)

   }

}

if(!is.null(relax.index)) {

x1.lars=fits.test3.25[[i]]$data$x[,relax.index]
relax.fit=lm(y.data~x1.lars-1)
lars.coefs.test3.25.relaxed[[i]][k,][relax.index]=relax.fit$coef

}

error[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],fits.test3.25[[i]]$constrained.fit$coefs[,k])
error.lars[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],lars.coefs.test3.25[[i]][k,])
error.relaxed[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],fits.test3.25.relaxed[[i]]$constrained.fit$coefs[,k])
error.lars.relaxed[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],lars.coefs.test3.25.relaxed[[i]][k,])

}

valid.place[i]=which.min(error)
valid.lam[i]=fits.test3.25[[i]]$constrained.fit$lambda[valid.place[i]]
valid.place.lars[i]=which.min(error.lars)
valid.lam.lars[i]=fits.test3.25[[i]]$constrained.fit$lambda[valid.place.lars[i]]
valid.place.relaxed[i]=which.min(error.relaxed)
valid.lam.relaxed[i]=fits.test3.25[[i]]$constrained.fit$lambda[valid.place.relaxed[i]]
valid.place.lars.relaxed[i]=which.min(error.lars.relaxed)
valid.lam.lars.relaxed[i]=fits.test3.25[[i]]$constrained.fit$lambda[valid.place.lars.relaxed[i]]


error.constrained.valid.test3.25[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],fits.test3.25[[i]]$constrained.fit$coefs[,valid.place[i]])/10000)
error.lars.valid.test3.25[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],lars.coefs.test3.25[[i]][valid.place.lars[i],])/10000)
error.constrained.valid.test3.25.relaxed[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],fits.test3.25.relaxed[[i]]$constrained.fit$coefs[,valid.place.relaxed[i]])/10000)
error.lars.valid.test3.25.relaxed[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],lars.coefs.test3.25.relaxed[[i]][valid.place.lars.relaxed[i],])/10000)

}


#a=0.5


#LASSO
for (i in 1:iter) {


lars.trial=lars(fits.test3.5[[i]]$data$x,fits.test3.5[[i]]$data$y,normalize=T,intercept=F)
lars.coefs.test3.5[[i]]=predict(lars.trial,s=fits.test3.5[[i]]$constrained.fit$lambda,type="coefficients",mode="lambda")$coef

}

lars.coefs.test3.5.relaxed=lars.coefs.test3.5
fits.test3.5.relaxed=fits.test3.5

valid.set.x=vector("list",iter)
valid.set.y=matrix(0,n,iter)
test.set.x=vector("list",iter)
test.set.y=matrix(0,n.test,iter)

for (i in 1:iter) {

 set.seed(seed[i]+1)
 valid.set.x[[i]]=matrix(rnorm(n*p),n,p)
 valid.set.y[,i]=as.vector(valid.set.x[[i]]%*%fits.test3.5[[i]]$data$beta+sigma*rnorm(n))
 test.set.x[[i]]=matrix(rnorm(n.test*p),n.test,p)
 test.set.y[,i]=as.vector(test.set.x[[i]]%*%fits.test3.5[[i]]$data$beta)
 
}

for (i in 1:iter) {

b=fits.test3.5[[i]]$data$b
y.data=fits.test3.5[[i]]$data$y


n.lam=length(fits.test3.5[[i]]$constrained.fit$lambda)
error=rep(0,n.lam)
error.lars=rep(0,n.lam)
error.relaxed=rep(0,n.lam)
error.lars.relaxed=rep(0,n.lam)


for (k in 1:n.lam) {


beta2.index=fits.test3.5[[i]]$constrained.fit$b2index[,k]
x2.data=fits.test3.5[[i]]$data$x[,beta2.index]
A2=fits.test3.5[[i]]$data$C.full[,beta2.index]
Y.til=y.data-x2.data%*%solve(A2)%*%b

total=c(1:p)

beta1.index.adj=NULL

for (j in 1:p) {
   if (fits.test3.5[[i]]$constrained.fit$coefs[,k][j]==0) {
   beta1.index.adj=c(beta1.index.adj,j)
   }
  }

beta1.index=total[-c(beta2.index,beta1.index.adj)]

if(length(beta1.index)==0) {

	fits.test3.5.relaxed[[i]]$constrained.fit$coefs[,k]==fits.test3.5[[i]]$constrained.fit$coefs[,k]

} else {

	x1.star=fits.test3.5[[i]]$data$x[,beta1.index]
	A1.star=fits.test3.5[[i]]$data$C.full[,beta1.index]
	X.til=x1.star-x2.data%*%solve(A2)%*%A1.star

	if (n<length(beta1.index)) {

   		lars.trial=lars(X.til,Y.til,normalize=T,intercept=F)
   		relax.fit=predict(lars.trial,s=fits.test3.5[[i]]$constrained.fit$lambda,type="coefficients",mode="lambda")$coef
   		fits.test3.5.relaxed[[i]]$constrained.fit$coefs[,k][beta1.index]=relax.fit[k,]
   		beta2=solve(A2)%*%(b-A1.star%*%relax.fit[k,])

	} else {
		relax.fit=lm(Y.til~X.til-1)
		if(length(beta1.index)==1) {
			beta2=solve(A2)%*%(b-A1.star*relax.fit$coef)
		} else {

			beta2=solve(A2)%*%(b-A1.star%*%relax.fit$coef)
		}
		fits.test3.5.relaxed[[i]]$constrained.fit$coefs[,k][beta1.index]=relax.fit$coef
	}

	fits.test3.5.relaxed[[i]]$constrained.fit$coefs[,k][beta2.index]=beta2

}

relax.index=NULL

for (j in 1:p) {

   if (lars.coefs.test3.5[[i]][k,][j]!=0) {

	relax.index=c(relax.index,j)

   }

}

if(!is.null(relax.index)) {

x1.lars=fits.test3.5[[i]]$data$x[,relax.index]
relax.fit=lm(y.data~x1.lars-1)
lars.coefs.test3.5.relaxed[[i]][k,][relax.index]=relax.fit$coef

}

error[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],fits.test3.5[[i]]$constrained.fit$coefs[,k])
error.lars[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],lars.coefs.test3.5[[i]][k,])
error.relaxed[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],fits.test3.5.relaxed[[i]]$constrained.fit$coefs[,k])
error.lars.relaxed[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],lars.coefs.test3.5.relaxed[[i]][k,])

}

valid.place[i]=which.min(error)
valid.lam[i]=fits.test3.5[[i]]$constrained.fit$lambda[valid.place[i]]
valid.place.lars[i]=which.min(error.lars)
valid.lam.lars[i]=fits.test3.5[[i]]$constrained.fit$lambda[valid.place.lars[i]]
valid.place.relaxed[i]=which.min(error.relaxed)
valid.lam.relaxed[i]=fits.test3.5[[i]]$constrained.fit$lambda[valid.place.relaxed[i]]
valid.place.lars.relaxed[i]=which.min(error.lars.relaxed)
valid.lam.lars.relaxed[i]=fits.test3.5[[i]]$constrained.fit$lambda[valid.place.lars.relaxed[i]]


error.constrained.valid.test3.5[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],fits.test3.5[[i]]$constrained.fit$coefs[,valid.place[i]])/10000)
error.lars.valid.test3.5[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],lars.coefs.test3.5[[i]][valid.place.lars[i],])/10000)
error.constrained.valid.test3.5.relaxed[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],fits.test3.5.relaxed[[i]]$constrained.fit$coefs[,valid.place.relaxed[i]])/10000)
error.lars.valid.test3.5.relaxed[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],lars.coefs.test3.5.relaxed[[i]][valid.place.lars.relaxed[i],])/10000)

}


#a=0.75


#LASSO
for (i in 1:iter) {


lars.trial=lars(fits.test3.75[[i]]$data$x,fits.test3.75[[i]]$data$y,normalize=T,intercept=F)
lars.coefs.test3.75[[i]]=predict(lars.trial,s=fits.test3.75[[i]]$constrained.fit$lambda,type="coefficients",mode="lambda")$coef

}

lars.coefs.test3.75.relaxed=lars.coefs.test3.75
fits.test3.75.relaxed=fits.test3.75

valid.set.x=vector("list",iter)
valid.set.y=matrix(0,n,iter)
test.set.x=vector("list",iter)
test.set.y=matrix(0,n.test,iter)

for (i in 1:iter) {

 set.seed(seed[i]+1)
 valid.set.x[[i]]=matrix(rnorm(n*p),n,p)
 valid.set.y[,i]=as.vector(valid.set.x[[i]]%*%fits.test3.75[[i]]$data$beta+sigma*rnorm(n))
 test.set.x[[i]]=matrix(rnorm(n.test*p),n.test,p)
 test.set.y[,i]=as.vector(test.set.x[[i]]%*%fits.test3.75[[i]]$data$beta)
 
}

for (i in 1:iter) {

b=fits.test3.75[[i]]$data$b
y.data=fits.test3.75[[i]]$data$y


n.lam=length(fits.test3.75[[i]]$constrained.fit$lambda)
error=rep(0,n.lam)
error.lars=rep(0,n.lam)
error.relaxed=rep(0,n.lam)
error.lars.relaxed=rep(0,n.lam)


for (k in 1:n.lam) {


beta2.index=fits.test3.75[[i]]$constrained.fit$b2index[,k]
x2.data=fits.test3.75[[i]]$data$x[,beta2.index]
A2=fits.test3.75[[i]]$data$C.full[,beta2.index]
Y.til=y.data-x2.data%*%solve(A2)%*%b

total=c(1:p)

beta1.index.adj=NULL

for (j in 1:p) {
   if (fits.test3.75[[i]]$constrained.fit$coefs[,k][j]==0) {
   beta1.index.adj=c(beta1.index.adj,j)
   }
  }

beta1.index=total[-c(beta2.index,beta1.index.adj)]

if(length(beta1.index)==0) {

	fits.test3.75.relaxed[[i]]$constrained.fit$coefs[,k]==fits.test3.75[[i]]$constrained.fit$coefs[,k]

} else {

	x1.star=fits.test3.75[[i]]$data$x[,beta1.index]
	A1.star=fits.test3.75[[i]]$data$C.full[,beta1.index]
	X.til=x1.star-x2.data%*%solve(A2)%*%A1.star

	if (n<length(beta1.index)) {

   		lars.trial=lars(X.til,Y.til,normalize=T,intercept=F)
   		relax.fit=predict(lars.trial,s=fits.test3.75[[i]]$constrained.fit$lambda,type="coefficients",mode="lambda")$coef
   		fits.test3.75.relaxed[[i]]$constrained.fit$coefs[,k][beta1.index]=relax.fit[k,]
   		beta2=solve(A2)%*%(b-A1.star%*%relax.fit[k,])

	} else {
		relax.fit=lm(Y.til~X.til-1)
		if(length(beta1.index)==1) {
			beta2=solve(A2)%*%(b-A1.star*relax.fit$coef)
		} else {

			beta2=solve(A2)%*%(b-A1.star%*%relax.fit$coef)
		}
		fits.test3.75.relaxed[[i]]$constrained.fit$coefs[,k][beta1.index]=relax.fit$coef
	}

	fits.test3.75.relaxed[[i]]$constrained.fit$coefs[,k][beta2.index]=beta2

}

relax.index=NULL

for (j in 1:p) {

   if (lars.coefs.test3.75[[i]][k,][j]!=0) {

	relax.index=c(relax.index,j)

   }

}

if(!is.null(relax.index)) {

x1.lars=fits.test3.75[[i]]$data$x[,relax.index]
relax.fit=lm(y.data~x1.lars-1)
lars.coefs.test3.75.relaxed[[i]][k,][relax.index]=relax.fit$coef

}

error[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],fits.test3.75[[i]]$constrained.fit$coefs[,k])
error.lars[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],lars.coefs.test3.75[[i]][k,])
error.relaxed[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],fits.test3.75.relaxed[[i]]$constrained.fit$coefs[,k])
error.lars.relaxed[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],lars.coefs.test3.75.relaxed[[i]][k,])

}

valid.place[i]=which.min(error)
valid.lam[i]=fits.test3.75[[i]]$constrained.fit$lambda[valid.place[i]]
valid.place.lars[i]=which.min(error.lars)
valid.lam.lars[i]=fits.test3.75[[i]]$constrained.fit$lambda[valid.place.lars[i]]
valid.place.relaxed[i]=which.min(error.relaxed)
valid.lam.relaxed[i]=fits.test3.75[[i]]$constrained.fit$lambda[valid.place.relaxed[i]]
valid.place.lars.relaxed[i]=which.min(error.lars.relaxed)
valid.lam.lars.relaxed[i]=fits.test3.75[[i]]$constrained.fit$lambda[valid.place.lars.relaxed[i]]


error.constrained.valid.test3.75[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],fits.test3.75[[i]]$constrained.fit$coefs[,valid.place[i]])/10000)
error.lars.valid.test3.75[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],lars.coefs.test3.75[[i]][valid.place.lars[i],])/10000)
error.constrained.valid.test3.75.relaxed[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],fits.test3.75.relaxed[[i]]$constrained.fit$coefs[,valid.place.relaxed[i]])/10000)
error.lars.valid.test3.75.relaxed[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],lars.coefs.test3.75.relaxed[[i]][valid.place.lars.relaxed[i],])/10000)

}


#a=1


#LASSO
for (i in 1:iter) {


lars.trial=lars(fits.test31[[i]]$data$x,fits.test31[[i]]$data$y,normalize=T,intercept=F)
lars.coefs.test31[[i]]=predict(lars.trial,s=fits.test31[[i]]$constrained.fit$lambda,type="coefficients",mode="lambda")$coef

}

lars.coefs.test31.relaxed=lars.coefs.test31
fits.test31.relaxed=fits.test31

valid.set.x=vector("list",iter)
valid.set.y=matrix(0,n,iter)
test.set.x=vector("list",iter)
test.set.y=matrix(0,n.test,iter)

for (i in 1:iter) {

 set.seed(seed[i]+1)
 valid.set.x[[i]]=matrix(rnorm(n*p),n,p)
 valid.set.y[,i]=as.vector(valid.set.x[[i]]%*%fits.test31[[i]]$data$beta+sigma*rnorm(n))
 test.set.x[[i]]=matrix(rnorm(n.test*p),n.test,p)
 test.set.y[,i]=as.vector(test.set.x[[i]]%*%fits.test31[[i]]$data$beta)
 
}

for (i in 1:iter) {

b=fits.test31[[i]]$data$b
y.data=fits.test31[[i]]$data$y


n.lam=length(fits.test31[[i]]$constrained.fit$lambda)
error=rep(0,n.lam)
error.lars=rep(0,n.lam)
error.relaxed=rep(0,n.lam)
error.lars.relaxed=rep(0,n.lam)


for (k in 1:n.lam) {


beta2.index=fits.test31[[i]]$constrained.fit$b2index[,k]
x2.data=fits.test31[[i]]$data$x[,beta2.index]
A2=fits.test31[[i]]$data$C.full[,beta2.index]
Y.til=y.data-x2.data%*%solve(A2)%*%b

total=c(1:p)

beta1.index.adj=NULL

for (j in 1:p) {
   if (fits.test31[[i]]$constrained.fit$coefs[,k][j]==0) {
   beta1.index.adj=c(beta1.index.adj,j)
   }
  }

beta1.index=total[-c(beta2.index,beta1.index.adj)]

if(length(beta1.index)==0) {

	fits.test31.relaxed[[i]]$constrained.fit$coefs[,k]==fits.test31[[i]]$constrained.fit$coefs[,k]

} else {

	x1.star=fits.test31[[i]]$data$x[,beta1.index]
	A1.star=fits.test31[[i]]$data$C.full[,beta1.index]
	X.til=x1.star-x2.data%*%solve(A2)%*%A1.star

	if (n<length(beta1.index)) {

   		lars.trial=lars(X.til,Y.til,normalize=T,intercept=F)
   		relax.fit=predict(lars.trial,s=fits.test31[[i]]$constrained.fit$lambda,type="coefficients",mode="lambda")$coef
   		fits.test31.relaxed[[i]]$constrained.fit$coefs[,k][beta1.index]=relax.fit[k,]
   		beta2=solve(A2)%*%(b-A1.star%*%relax.fit[k,])

	} else {
		relax.fit=lm(Y.til~X.til-1)
		if(length(beta1.index)==1) {
			beta2=solve(A2)%*%(b-A1.star*relax.fit$coef)
		} else {

			beta2=solve(A2)%*%(b-A1.star%*%relax.fit$coef)
		}
		fits.test31.relaxed[[i]]$constrained.fit$coefs[,k][beta1.index]=relax.fit$coef
	}

	fits.test31.relaxed[[i]]$constrained.fit$coefs[,k][beta2.index]=beta2

}

relax.index=NULL

for (j in 1:p) {

   if (lars.coefs.test31[[i]][k,][j]!=0) {

	relax.index=c(relax.index,j)

   }

}

if(!is.null(relax.index)) {

x1.lars=fits.test31[[i]]$data$x[,relax.index]
relax.fit=lm(y.data~x1.lars-1)
lars.coefs.test31.relaxed[[i]][k,][relax.index]=relax.fit$coef

}

error[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],fits.test31[[i]]$constrained.fit$coefs[,k])
error.lars[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],lars.coefs.test31[[i]][k,])
error.relaxed[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],fits.test31.relaxed[[i]]$constrained.fit$coefs[,k])
error.lars.relaxed[k]=pred.error(valid.set.x[[i]],valid.set.y[,i],lars.coefs.test31.relaxed[[i]][k,])

}

valid.place[i]=which.min(error)
valid.lam[i]=fits.test31[[i]]$constrained.fit$lambda[valid.place[i]]
valid.place.lars[i]=which.min(error.lars)
valid.lam.lars[i]=fits.test31[[i]]$constrained.fit$lambda[valid.place.lars[i]]
valid.place.relaxed[i]=which.min(error.relaxed)
valid.lam.relaxed[i]=fits.test31[[i]]$constrained.fit$lambda[valid.place.relaxed[i]]
valid.place.lars.relaxed[i]=which.min(error.lars.relaxed)
valid.lam.lars.relaxed[i]=fits.test31[[i]]$constrained.fit$lambda[valid.place.lars.relaxed[i]]


error.constrained.valid.test31[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],fits.test31[[i]]$constrained.fit$coefs[,valid.place[i]])/10000)
error.lars.valid.test31[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],lars.coefs.test31[[i]][valid.place.lars[i],])/10000)
error.constrained.valid.test31.relaxed[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],fits.test31.relaxed[[i]]$constrained.fit$coefs[,valid.place.relaxed[i]])/10000)
error.lars.valid.test31.relaxed[[i]]=sqrt(pred.error(test.set.x[[i]],test.set.y[,i],lars.coefs.test31.relaxed[[i]][valid.place.lars.relaxed[i],])/10000)

}

###################################################################################
###################    Calculating Error Rates     ################################
###################           (Table 3)            ################################
###################################################################################


###### Table 3 Results

### Setting 1

## a=0.25

mean(error.lars.valid.test1.25)
sd(error.lars.valid.test1.25)/sqrt(iter)
mean(error.constrained.valid.test1.25)
sd(error.constrained.valid.test1.25)/sqrt(iter)
mean(error.lars.valid.test1.25.relaxed)
sd(error.lars.valid.test1.25.relaxed)/sqrt(iter)
mean(error.constrained.valid.test1.25.relaxed)
sd(error.constrained.valid.test1.25.relaxed)/sqrt(iter)

## a=0.5

mean(error.lars.valid.test1.5)
sd(error.lars.valid.test1.5)/sqrt(iter)
mean(error.constrained.valid.test1.5)
sd(error.constrained.valid.test1.5)/sqrt(iter)
mean(error.lars.valid.test1.5.relaxed)
sd(error.lars.valid.test1.5.relaxed)/sqrt(iter)
mean(error.constrained.valid.test1.5.relaxed)
sd(error.constrained.valid.test1.5.relaxed)/sqrt(iter)

## a=0.75

mean(error.lars.valid.test1.75)
sd(error.lars.valid.test1.75)/sqrt(iter)
mean(error.constrained.valid.test1.75)
sd(error.constrained.valid.test1.75)/sqrt(iter)
mean(error.lars.valid.test1.75.relaxed)
sd(error.lars.valid.test1.75.relaxed)/sqrt(iter)
mean(error.constrained.valid.test1.75.relaxed)
sd(error.constrained.valid.test1.75.relaxed)/sqrt(iter)

## a=1

mean(error.lars.valid.test11)
sd(error.lars.valid.test11)/sqrt(iter)
mean(error.constrained.valid.test11)
sd(error.constrained.valid.test11)/sqrt(iter)
mean(error.lars.valid.test11.relaxed)
sd(error.lars.valid.test11.relaxed)/sqrt(iter)
mean(error.constrained.valid.test11.relaxed)
sd(error.constrained.valid.test11.relaxed)/sqrt(iter)


### Setting 2

## a=0.25

mean(error.lars.valid.test2.25)
sd(error.lars.valid.test2.25)/sqrt(iter)
mean(error.constrained.valid.test2.25)
sd(error.constrained.valid.test2.25)/sqrt(iter)
mean(error.lars.valid.test2.25.relaxed)
sd(error.lars.valid.test2.25.relaxed)/sqrt(iter)
mean(error.constrained.valid.test2.25.relaxed)
sd(error.constrained.valid.test2.25.relaxed)/sqrt(iter)

## a=0.5

mean(error.lars.valid.test2.5)
sd(error.lars.valid.test2.5)/sqrt(iter)
mean(error.constrained.valid.test2.5)
sd(error.constrained.valid.test2.5)/sqrt(iter)
mean(error.lars.valid.test2.5.relaxed)
sd(error.lars.valid.test2.5.relaxed)/sqrt(iter)
mean(error.constrained.valid.test2.5.relaxed)
sd(error.constrained.valid.test2.5.relaxed)/sqrt(iter)

## a=0.75

mean(error.lars.valid.test2.75)
sd(error.lars.valid.test2.75)/sqrt(iter)
mean(error.constrained.valid.test2.75)
sd(error.constrained.valid.test2.75)/sqrt(iter)
mean(error.lars.valid.test2.75.relaxed)
sd(error.lars.valid.test2.75.relaxed)/sqrt(iter)
mean(error.constrained.valid.test2.75.relaxed)
sd(error.constrained.valid.test2.75.relaxed)/sqrt(iter)

## a=1

mean(error.lars.valid.test21)
sd(error.lars.valid.test21)/sqrt(iter)
mean(error.constrained.valid.test21)
sd(error.constrained.valid.test21)/sqrt(iter)
mean(error.lars.valid.test21.relaxed)
sd(error.lars.valid.test21.relaxed)/sqrt(iter)
mean(error.constrained.valid.test21.relaxed)
sd(error.constrained.valid.test21.relaxed)/sqrt(iter)


### Setting 3

## a=0.25

mean(error.lars.valid.test3.25)
sd(error.lars.valid.test3.25)/sqrt(iter)
mean(error.constrained.valid.test3.25)
sd(error.constrained.valid.test3.25)/sqrt(iter)
mean(error.lars.valid.test3.25.relaxed)
sd(error.lars.valid.test3.25.relaxed)/sqrt(iter)
mean(error.constrained.valid.test3.25.relaxed)
sd(error.constrained.valid.test3.25.relaxed)/sqrt(iter)

## a=0.5

mean(error.lars.valid.test3.5)
sd(error.lars.valid.test3.5)/sqrt(iter)
mean(error.constrained.valid.test3.5)
sd(error.constrained.valid.test3.5)/sqrt(iter)
mean(error.lars.valid.test3.5.relaxed)
sd(error.lars.valid.test3.5.relaxed)/sqrt(iter)
mean(error.constrained.valid.test3.5.relaxed)
sd(error.constrained.valid.test3.5.relaxed)/sqrt(iter)

## a=0.75

mean(error.lars.valid.test3.75)
sd(error.lars.valid.test3.75)/sqrt(iter)
mean(error.constrained.valid.test3.75)
sd(error.constrained.valid.test3.75)/sqrt(iter)
mean(error.lars.valid.test3.75.relaxed)
sd(error.lars.valid.test3.75.relaxed)/sqrt(iter)
mean(error.constrained.valid.test3.75.relaxed)
sd(error.constrained.valid.test3.75.relaxed)/sqrt(iter)

## a=1

mean(error.lars.valid.test31)
sd(error.lars.valid.test31)/sqrt(iter)
mean(error.constrained.valid.test31)
sd(error.constrained.valid.test31)/sqrt(iter)
mean(error.lars.valid.test31.relaxed)
sd(error.lars.valid.test31.relaxed)/sqrt(iter)
mean(error.constrained.valid.test31.relaxed)
sd(error.constrained.valid.test31.relaxed)/sqrt(iter)







