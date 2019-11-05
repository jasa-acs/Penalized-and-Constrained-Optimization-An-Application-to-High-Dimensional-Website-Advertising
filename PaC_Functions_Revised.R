#######################################

### compare.to.lasso() function ###

### Function overview: This function is only used for reproducibility.
###   This is the main function to call in order to compare the PaC and
###   the unconstrained LASSO for simulated data. To do so, it calls 
###   generate.data to create the simulated data set according to given
###   values of observations (n), variables (p), and constraints (m).
###   It then passes this data to the lasso.c function to calculate the 
###   PaC coefficient paths and the lars function (part of the lars 
###   package) to calculate corresponding LASSO coefficient paths.

#######################################
### Function Arguments 
###
### n: number of rows in randomly-generated data set (default is 1000)
### p: number of variables in randomly-generated data set (default is 10)
### m: number of constraints in randomly-generated constraint matrix (default is 5)
### s: number of true non-zero elements in coefficient vector beta1 (default is 5)
### sigma: standard deviation of noise in response (default is 1, indicating standard normal)
### cov.mat a covariance matrix applied in the generation of data to impose a correlation structure. Default is NULL (no correlation) 
### err: error to be introduced in random generation of coefficient values. Default is no error (err = 0)
### trace: should function print output as algorithm runs. Default is FALSE
### seed: option to choose a number to reproduce results given by the function. Default is NULL
### l.min: lowest value of lambda to consider (used as 10^l.min). Default is -2
### l.max: largest value of lambda to consider (used as 10^l.max). Default is 6
### backwards: which direction should algorithm go, backwards from lambda = 10^l.max (TRUE)
###             or forwards from 10^l.max and then backwards if algorithm gets stuck (FALSE).
###             Default is FALSE.
### glasso: should the generalized Lasso be used (TRUE) or standard Lasso (FALSE). Default is FALSE
### plots: should plots of the coefficent paths and other method results be produced (TRUE) or not (FALSE). Default is FALSE
#######################################


#######################################
### Function Output
###
### data: randomly-generated data
### constrained.fit: fit of model (output from lasso.c)
### lars.error: sum of squared errors for standard Lasso at each lambda
### constrained.error: sum of squared errors for constrained Lasso at each lambda
#######################################

compare.to.lasso <-
  function(n = 1000, p = 10, m = 5, s = 5, sigma = 1, cov.mat = NULL, 
           err = 0, trace = F, seed = NULL, l.min = -2, l.max = 6,
           backwards = F, glasso = F, plots = F){
    
    if (!is.null(seed))
      set.seed(seed)
    temp.data = generate.data(n = n, p = p, m = m, cov.mat = cov.mat,
                              s = s, sigma = sigma, glasso, err = err)
    constrained.fit = lasso.c(temp.data$x, temp.data$y, temp.data$C.full,
                              temp.data$b, l.min = l.min, l.max = l.max,
                              verbose = trace,intercept = F, normalize = T,
                              backwards = backwards)
    tmp.lars = lars(temp.data$x, temp.data$y, normalize = T, intercept = F)
    tmp.lars.coefs = predict(tmp.lars, s = constrained.fit$lambda,
                             type = "coefficients", mode = "lambda")$coef
    constrained.error = apply((constrained.fit$coefs-temp.data$beta)^2,2,sum)
    lars.error = apply((t(tmp.lars.coefs)-temp.data$beta)^2,2,sum)
    
    if(plots==T) {
      par(mfrow = c(1,3))
      plot(log10(constrained.fit$lambda), lars.error, type = 'b',
           ylim = range(c(lars.error,constrained.error)))
      lines(log10(constrained.fit$lambda), constrained.error, col = 2,type = 'b')
      c.L1 = apply(abs(constrained.fit$coefs),2,sum)
      matplot(c.L1,t(constrained.fit$coefs), type = 'l')
      lars.L1 = apply(abs(tmp.lars.coefs), 1, sum)
      matplot(lars.L1, tmp.lars.coefs, type = 'l')
    }
    list(data = temp.data, constrained.fit = constrained.fit,
         lars.error = lars.error, constrained.error = constrained.error)
  }

#######################################

### generate.data() function ###

### Function overview: This function is only used for reproducibility.
###   This function is called by compare.to.lasso to generate random data.

#######################################
### Function Arguments
###
### n: number of rows in randomly-generated data set (default is 1000)
### p: number of variables in randomly-generated data set (default is 10)
### m: number of constraints in randomly-generated constraint matrix (default is 5)
### s: number of true non-zero elements in coefficient vector beta1 (default is 5)
### sigma: standard deviation of noise in response (default is 1, indicating standard normal)
### cov.mat a covariance matrix applied in the generation of data to impose a correlation structure. Default is NULL (no correlation) 
### err: error to be introduced in random generation of coefficient values. Default is no error (err = 0)
### glasso: should the generalized Lasso be used (TRUE) or standard Lasso (FALSE). Default is FALSE
#######################################


#######################################
### Function Output
###
### x: generated x data
### y: generated response y vector
### C.full: generated full constraint matrix
### b: generated constraint vector b
### b.run: if error was included, the error-adjusted value of b
### beta: the complete beta vector, including generated beta1 and beta2
#######################################

generate.data <-
  function(n = 1000, p = 10, m = 5, cov.mat = cov.mat, s = 5, sigma = 1,
           glasso = F, err = 0){
    
    if(!is.null(cov.mat)){
      x<-mvrnorm(n = n, rep(0, p), cov.mat)
    } 
    else {x <- matrix(rnorm(n*p), n, p)}
    C.full = matrix(rnorm(m*p), m, p)
    b = rnorm(m)
    b.run = b*(1+err)
    if (glasso)
      b = rep(0,M)
    beta1 = c(runif(s)+1,rep(0,p-m-s))
    index = (1:(p-m))
    C1 = C.full[,index]
    C2 = C.full[,-index]
    beta2 = as.vector(solve(C2)%*%(b.run-C1%*%beta1))
    s.beta = sd(beta2)
    beta2 = beta2/s.beta
    C.full[,-index] = C.full[,-index]*s.beta
    beta = c(beta1,beta2)
    y = as.vector(x%*%beta+sigma*rnorm(n))
    list(x = x, y = y, C.full = C.full, b = b, b.run = b.run, beta = beta)}

#######################################

### lars.c() function ###

### Function overview: This function computes the PaC constrained LASSO
###   coefficient paths following the methodology laid out in the PaC 
###   paper. This function could be called directly as a standalone 
###   function, but the authors recommend using lasso.c for any 
###   implementation. This is because lasso.c has additional checks for
###   errors across the coefficient paths and allows for users to go 
###   forwards and backwards through the paths if the paths are unable
###   to compute in a particular direction for a particular run.

#######################################
### Function Arguments
###
### x: independent variable matrix of data to be used in calculating PaC coefficient paths
### y: response vector of data to be used in calculating PaC coefficient paths
### C.full: complete constraint matrix C
### b: constraint vector b
### l.min: lowest value of lambda to consider (used as 10^l.min). Default is -2
### l.max: largest value of lambda to consider (used as 10^l.max). Default is 6
### step: step size increase in lambda attempted at each iteration (by a factor of 10^step). Default is 0.2
### beta0: initial guess for beta coefficient vector. Default is NULL (indicating
###         initial vector should be calculated by algorithm)
### verbose: should function print output at each iteration (TRUE) or not (FALSE). Default is FALSE
### max.it: maximum number of times step size is halved before the algorithm terminates and gives a warning. Default is 12
### intercept: should intercept be included in modeling (TRUE) or not (FALSE). Default is TRUE.
### normalize: should X data be normalized. Default is TRUE
### forwards: if forwards = F, then the algorithm starts at 10^l.max and
###             moves backwards (without the forward step). If forwards = T,
###             algorithm starts at 10^l.min and works forward. Default is FALSE
#######################################


#######################################
### Function Output
###
### coefs: A p by length(lambda) matrix with each column corresponding to the beta estimate for that lambda
### lambda: the grid of lambdas used to calculate the coefficients on the coefficient path
### intercept: the intercept value
### error: did the algorithm terminate due to too many iterations (TRUE or FALSE)
### b2index: the index of the beta2 values identified by the algorithm at each lambda
#######################################

lars.c <-
  function(x, y, C.full, b, l.min = -2, l.max = 6, step = 0.2,
           beta0 = NULL, verbose = F, max.it = 12, intercept = T,
           normalize = T, forwards = T){
    p = ncol(x)
    n = nrow(x)
    m = nrow(C.full)
    beta.new = rep(0,p)
    one <- rep(1, n)
    if (intercept) {
      meanx <- drop(one %*% x)/n
      x <- scale(x, meanx, FALSE)
      mu <- mean(y)
      y <- drop(y - mu)
    }
    else {
      meanx <- rep(0, p)
      mu <- 0
      y <- drop(y)
    }
    normx <- rep(1, p)
    if (normalize) {
      normx <- sqrt(n)*apply(x,2,sd,na.rm=T)
      x <- scale(x, FALSE, normx)}
    C.full = t(t(C.full)/normx)
    if (!forwards){
      lambda = lambda.old = 10^l.max
      step = -step
      if (is.null(beta0))
        beta0 = lin.int(C.full,b)
    }
    else{
      lambda = lambda.old = 10^l.min
      if (is.null(beta0))
        beta0 = quad.int(x,y,C.full,b,lambda)
    }
    step.orig = step
    coefs = grid = b2index = NULL
    t.data=transform.data(x, y, C.full, b, lambda, beta0)
    beta1.old = t.data$beta1
    beta2.old = t.data$beta2
    iterations = 1
    end.loop = F
    while (!end.loop & (iterations <= max.it)){
      iterations = 1
      loop = T
      
      while (loop & (iterations <= max.it)){
        
        t.data$y = t.data$y+(lambda-lambda.old)*t.data$C
        beta1.new = rep(0,length(beta1.old))
        fit = lars(t.data$x[,t.data$active], t.data$y, normalize = F,
                   intercept = F)
        beta1.new[t.data$active] = predict(fit, s = lambda,
                                           type = "coefficients",
                                           mode = "lambda")$coef
        beta2.new = beta2.old + t.data$a2 %*% (beta1.old - beta1.new)
        bad.beta2 = (sum(abs(sign(t.data$beta2)-sign(beta2.new))) != 0)
        X_star = t.data$x
        derivs = abs(as.vector(t(X_star)%*%(X_star%*%beta1.new))-t(X_star)%*%t.data$Y_star-lambda*t.data$C2)
        bad.active = F
        if (n<(p-m))
          bad.active = (max(derivs[-t.data$active])>lambda)
        if (bad.beta2 | bad.active){
          t.data$y = t.data$y-(lambda-lambda.old)*t.data$C
          step=step/2
          lambda = lambda.old*10^step
          iterations = iterations+1}
        else
          loop = F
      }
      if (iterations <= max.it){
        if (verbose==T){
          print(paste("Lambda =",round(lambda,3)))
          if (abs(step) < abs(step.orig))
            print(paste("Step size reduced to ",step))
        }
        step = step.orig
        beta.new[t.data$beta2.index] = beta2.new
        beta.new[-t.data$beta2.index] = beta1.new
        coefs = cbind(coefs,beta.new)
        b2index = cbind(b2index,t.data$beta2.index)
        change.beta = (min(abs(beta2.new))<max(abs(beta1.new)))
        change.active = F
        if (n<(p-m))
          change.active = (min(derivs[t.data$active]) < max(derivs[-t.data$active]))
        if (change.beta | change.active){
          t.data = transform.data(x, y, C.full, b, lambda, beta.new)
          beta1.new = t.data$beta1
          beta2.new = t.data$beta2}
        beta1.old = beta1.new
        beta2.old = beta2.new
        grid = c(grid,lambda)
        lambda.old = lambda
        lambda = lambda*10^step}
      if ((forwards & (lambda>10^l.max)) | (!forwards & (lambda<10^l.min)))
        end.loop = T
    }
    if (iterations > max.it)
      print(paste("Warning: Algorithm terminated at lambda =",round(lambda.old,1),": Maximum iterations exceeded."))
    colnames(coefs) = intercept = NULL
    if (!is.null(grid)){
      coefs = coefs/normx
      coefs = coefs[,order(grid)]
      b2index = b2index[,order(grid)]
      grid = sort(grid)
      intercept = mu-drop(t(coefs)%*%meanx)}
    list(coefs = coefs, lambda = grid, intercept = intercept,
         error = (iterations>max.it), b2index = b2index)
  }

#######################################

### lasso.c() function ###

### Function overview: This is a wrapper function for the lars.c PaC 
###   constrained Lasso function. lasso.c controls the overall path,
###   providing checks for the path and allowing the user to control
###   how the path is computed (and what to do in the case of a stopped path).

#######################################
### Function Arguments
###
### x: independent variable matrix of data to be used in calculating PaC coefficient paths
### y: response vector of data to be used in calculating PaC coefficient paths
### C.full: complete constraint matrix C
### b: constraint vector b
### l.min: lowest value of lambda to consider (used as 10^l.min). Default is -2
### l.max: largest value of lambda to consider (used as 10^l.max). Default is 6
### step: step size increase in lambda attempted at each iteration (by a factor of 10^step). Default is 0.2
### beta0: initial guess for beta coefficient vector. Default is NULL (indicating
###         initial vector should be calculated by algorithm)
### verbose: should function print output at each iteration (TRUE) or not (FALSE). Default is FALSE
### max.it: maximum number of times step size is halved before the algorithm terminates and gives a warning. Default is 12
### intercept: should intercept be included in modeling (TRUE) or not (FALSE). Default is TRUE.
### normalize: should X data be normalized. Default is TRUE
### backwards: which direction should algorithm go, backwards from lambda = 10^l.max (TRUE)
###             or forwards from 10^l.max and then backwards if algorithm gets stuck (FALSE).
###             Default is FALSE.
#######################################


#######################################
### Function Output
###
### coefs: A p by length(lambda) matrix with each column corresponding to the beta estimate for that lambda
### lambda: vector of values of lambda that were fit
### intercept: vector with each element corresponding to beta 0 for corresponding lambda
### error: Indicator of whether the algorithm terminated early because max.it was reached
#######################################

lasso.c <-
  function(x, y, C.full, b, l.min = -2, l.max = 6, step = 0.2,
           beta0 = NULL, verbose = F, max.it = 12, intercept = T,
           normalize = T, backwards = F){
    
    if (!backwards){
      fit = lars.c(x, y, C.full, b, l.min = l.min, l.max = l.max,
                   step = step, beta0 = beta0, verbose = verbose,
                   max.it = max.it, intercept = intercept,
                   normalize = normalize, forwards = T)
      if (fit$error | is.null(fit$coefs)){
        if (is.null(fit$lambda))
          fit$lambda = 10^l.min
        fit2 = lars.c(x, y, C.full, b, l.min = log10(max(fit$lambda)),
                      l.max = l.max, step = step, beta0 = beta0,
                      verbose = verbose, max.it = max.it,
                      intercept = intercept, normalize = normalize,
                      forwards = F)
        if (is.null(fit$coefs))
          fit$lambda = NULL
        fit$coefs = cbind(fit$coefs,fit2$coefs)
        fit$lambda = c(fit$lambda,fit2$lambda)
        fit$intercept = c(fit$intercept,fit2$intercept)
        fit$b2index = cbind(fit$b2index,fit2$b2index)
      }}
    else
      fit = lars.c(x, y, C.full, b, l.min = l.min, l.max = l.max,
                   step = step, beta0 = beta0, verbose = verbose,
                   max.it = max.it, intercept = intercept,
                   normalize = normalize, forwards = F)
    fit}

#######################################

### lin.int() function ###

### Function overview: This function is called internally by lars.c
###   to get the linear programming initial fit if the user requests
###   implementation of the algorithm starting at the largest lambda
###   value and proceeding backwards.

#######################################
### Function Arguments 
### C.full: complete constraint matrix C
### b: constraint vector b
#######################################


#######################################
### Function Output
###
### beta: the initial beta vector of coefficients to use for the PaC algorithm
#######################################

lin.int <-
  function(C.full, b){
    
    p = ncol(C.full)
    Cmat = cbind(C.full,-C.full)
    cvec = rep(1,2*p)
    temp = linp(E = Cmat,F = b,Cost = cvec)
    beta = temp$X[1:p]-temp$X[(p+1):(2*p)]
    beta
    
  }


#######################################

### quad.int() function ###

### Function overview: This function is called internally by lars.c
###   to get the quadratic programming fit if the user requests
###   implementation of the algorithm starting at the smallest lambda
###   value and proceeding forwards.

#######################################
### Function Arguments 
###
### X: independent variable matrix of data to be used in calculating PaC coefficient paths
### Y: response vector of data to be used in calculating PaC coefficient paths
### C.full: complete constraint matrix C
### b: constraint vector b
### lambda: value of lambda
### d: value close to zero used for SVD matrix decomposition. Default is 10^-7
#######################################


#######################################
### Function Output
###
### beta: the initial beta vector of coefficients to use for the PaC algorithm
#######################################

quad.int <-
  function(X, Y, C.full, b, lambda, d=10^-7){
    
    Xstar = cbind(X,-X)
    p = ncol(X)
    Dmat = t(Xstar)%*%Xstar+d*diag(2*p)
    dvec = as.vector(t(Xstar)%*%Y-lambda*rep(1,2*p))
    Cmat = cbind(t(cbind(C.full,-C.full)),diag(2*p))
    bvec = c(b,rep(0,2*p))
    temp = solve.QP(Dmat, dvec, Cmat, bvec, meq = length(b))
    beta = temp$sol[1:p]-temp$sol[(p+1):(2*p)]
    beta
    
  }

#######################################

### transform.data() function ###

### Function overview: This function is called internally by lars.c
###   to compute the transformed versions of the X, Y, and constraint
###   matrix data, as shown in the PaC paper.

#######################################
### Function Arguments 
###
### X: independent variable matrix of data to be used in calculating PaC coefficient paths
### Y: response vector of data to be used in calculating PaC coefficient paths
### C.full: complete constraint matrix C
### b: constraint vector b
### lambda: value of lambda
### beta0: initial guess for beta coefficient vector
### eps: value close to zero used to verify SVD decomposition. Default is 10^-8
#######################################


#######################################
### Function Output
###
### x: transformed x data to be used in the PaC algorithm
### y: transformed y data to be used in the PaC algorithm
### Y_star: transformed Y* value to be used in the PaC algorithm
### a2: index of A used in the calculation of beta2 (the non-zero coefficients)
### beta1: beta1 values
### beta2: beta2 values
### C: constraint matrix
### C2: subset of constraint matrix corresponding to non-zero coefficients
### active.beta: index of non-zero coefficient values
### beta2.index: index of non-zero coefficient values
#######################################

transform.data <-
  function(x, y, C.full, b, lambda, beta0, eps = 10^-8){
    
    p = ncol(x)
    m = nrow(C.full)
    n = nrow(x)
    beta.order = order(-abs(beta0))
    s = svd(C.full[,beta.order[1:m]])
    if (min(s$d)<eps){
      i = 1
      while (i <= m){
        s = svd(C.full[,beta.order[1:i]])
        if (min(s$d)<eps)
          beta.order = beta.order[-i]
        else
          i = i+1
      }
    }
    beta2.index = sort(beta.order[1:m])
    beta2 = beta0[beta2.index]
    beta1 = beta0[-beta2.index]
    s2 = sign(beta2)
    # Find C (split by A) and X split matrices
    
    A2 = C.full[,beta2.index]
    A1 = C.full[,-beta2.index]
    if(length(beta2.index)==1){A1 = t(A1)}
    A2.inv = solve(A2)
    a2 = A2.inv %*% A1
    X2 = x[,beta2.index]
    X1 = x[,-beta2.index]
    
    # Compute X_star, X_dash, Y_star, and Y_tilde as referenced in writeup
    
    X_star = X1 - (X2 %*% A2.inv %*% A1)
    Y_star = as.vector(y - (X2 %*% (A2.inv %*% b)))
    C2 = (t(A1)%*%(t(A2.inv)%*%s2))
    active.beta = sort(order(-abs(as.vector(t(X_star)%*%(X_star%*%beta1))-t(X_star)%*%Y_star-lambda*C2))[1:min(c(p-m,n))])
    s = svd(X_star[,active.beta])
    X_dash = s$u %*% diag(1/s$d) %*% t(s$v)
    C = X_dash %*% C2[active.beta]
    Y_til = Y_star + lambda*C
    
    list(x = X_star, y = Y_til, Y_star = Y_star, a2 = a2, beta1 = beta1,
         beta2 = beta2, C = C, C2 = C2, active.beta = active.beta,
         beta2.index = beta2.index)
  }

#######################################

### pred.error() function ###

### Function overview: This function is only used for reproducibility.
###   Its purpose is to compute the prediction SSE for a given
###   combination of data and estimated coefficients.

#######################################
### Function Arguments 
### x: matrix of independent x data
### y: true response vector y
### betas: vector of coefficient values
#######################################


#######################################
### Function Output
###
### error: calculated SSE
#######################################

pred.error<-
  function(x, y, betas){
    
    y.est = x %*% betas
    error = sum((y-y.est)^2)
    error
    
  }

#######################################

### ELMSO() function ###

### Function overview: This function is used in the calculation of the case
###   study data of Section 6. This function implements the calculation
###   of the ELMSO method given in Paulson et al. (2018) for calculating
###   reach/click rate on page view data for Internet website advertising problems.

#######################################
### Function Arguments 
###
### z: matrix of views to each advertising channel (e.g. page view matrix to Internet websites)
### beta: parameter in ELMSO only to control the shape of budget allocation curve.
###           Default value is 1, indicating each allocation channel has a linear slope for budget
###           allocation (each additional dollar of revenue has a constant buying effect).
### gamma: vector of gamma values (the fraction of all ads purchased at website j for every dollar spent by a campaign)
### step: step size in lambda for calculation of coefficient path. Default is 0.05
### size: total number of budget points that coefficients should be calculated for. Default is 100
### q: vector of clickthrough rates. Default is NULL, indicating reach calculation.
#######################################


#######################################
### Function Output
###
### w: A p by length(lambda) matrix with each column corresponding to the coefficient estimate for that lambda
### w.sum: a total sum of all coefficients at each lambda (used as budget calculation for advertising example)
### lam: grid of lambda values used in the estimation
### click.est: vector of click rates if q is not NULL
#######################################

ELMSO <-
  function(z, beta = 1, gamma, step = 0.05, size = 100, q = NULL){
    
    if(is.null(q)) {q = rep(1,dim(z)[1])}
    
    gamma = gamma*q
    
    B = 0
    n = nrow(z)
    p = ncol(z)
    
    az = z
    azw = z
    
    w.til = rep(0,p)
    w = rep(0,p)
    
    w.mat = w.sum = lam.grid = click.est = NULL
    Q = rep(1,n)
    
    for(j in 1:p){
      az[,j] = gamma[j]*z[,j]
      azw[,j] = gamma[j]*z[,j]*((w.til[j]+1)^beta-1)
    }
 
    
    Q = exp(-1*as.matrix(azw)%*%rep(1,p))
    
    s = 0
    
    for(j in 1:p){
      
      L1.temp = -az[,j]*beta*(w.til[j]+1)^(beta-1)*Q
      L2.temp = -L1.temp/((w.til[j]+1))*(1-beta+beta*az[,j]*(w.til[j]+1)^beta)
      
      S.temp = sum(w.til[j]*L2.temp-L1.temp)
      
      if(s < abs(S.temp)){s = abs(S.temp)}
      
      
    }
    
    lam = (2*s)^0.1
    lam.min = ((2*s)^0.1-step*size)
    
    
    while(lam > lam.min){
      converge = 1
      iter = 0
      while(converge>10^-3 & iter<200){
        
        iter = iter+1
        
        for(j in 1:p){
          
          L1 = -az[,j]*beta*(w.til[j]+1)^(beta-1)*Q
          L2 = -L1/((w.til[j]+1))*(1-beta+beta*az[,j]*(w.til[j]+1)^beta)
          
          S = sum(w.til[j]*L2-L1)
          sum.L = sum(L2)
          
          if(abs(S)>10^lam){
            w[j] = w.til[j]-(sum(L1)+10^lam*sign(S))/sum.L
          }
          
          else {
            w[j] = 0
          }
          
          Q = Q*exp(az[,j]*((w.til[j]+1)^beta-(w[j]+1)^beta))
          
        }
        
        if(sum(w)==0){
          converge = 0
        }
        
        else{
          converge = sum(abs(w.til-w))/sum(w)
        }
        
        w.til = w
        B = sum(w.til)
        
      }
      
      if(iter>200){
        
        print("Max iterations reached. Converge value: ",converge)
        
      }
      
      w.mat = cbind(w.mat,w)
      w.sum = cbind(w.sum,sum(w))
      lam.grid = c(lam.grid,10^lam)
      lam = lam-step
      click.est = c(click.est,1/n*sum(Q))
      
    }
    
    list(w = w.mat, w.sum = w.sum, lam = lam.grid, click.est = click.est)
    
  }


#######################################

### reach.calc() function ###

### Function overview: This function calculates the reach of an
###     advertising model given a set of allocation coefficients and
###     a matrix of views across a set of unique advertising
###     channels/opportunities. In the case study presented in the
###     PaC paper, these are Internet websites and views by users to
###     those sites.

#######################################
### Function Arguments 
###
### z: matrix of views to each advertising channel (e.g. page view matrix to Internet websites)
### gamma: vector of gamma values (the fraction of all ads purchased at website j for every dollar spent by a campaign)
### step: step size in lambda for calculation of coefficient path. Default is 0.05
### size: total number of budget points that coefficients should be calculated for. Default is 100
### q: vector of clickthrough rates. Default is NULL, indicating reach calculation.
#######################################


#######################################
### Function Output
###
### w: A p by length(lambda) matrix with each column corresponding to the coefficient estimate for that lambda
### w.sum: a total sum of all coefficients at each lambda (used as budget calculation for advertising example)
### lam: grid of lambda values used in the estimation
### click.est: vector of click rates if q is not NULL
#######################################

reach.calc<-function(w, z, gamma, q = NULL) {
  
  if(is.null(q)) {q = rep(1,dim(z)[1])}
  
  gamma = gamma*q
  
  p = ncol(z)
  azw = z
  
  for(j in 1:p){
    azw[,j] = gamma[j]*z[,j]*w[j]
  }
  
  sum1 = exp(-1*as.matrix(azw)%*%rep(1,p))
  
  1-1/nrow(z)*sum(sum1)
  
}

#######################################

### compare.to.ELMSO() function ###

### Function overview: This function is used in the calculation of the
###     case study data of Section 6. This function implements the PaC
###     calculation in comparison to the ELMSO method given in Paulson
###     et al. (2018) for calculating reach on page view data for
###     Internet website advertising problems.

#######################################
### Function Arguments 
###
### z: matrix of views to each advertising channel (e.g. page view matrix to Internet websites)
### gamma: vector of gamma values (the fraction of all ads purchased at website j for every dollar spent by a campaign)
### step: step size in lambda for calculation of coefficient path. Default is 0.05
### size: total number of budget points that coefficients should be calculated for. Default is 100
### q: vector of clickthrough rates. Default is NULL, indicating reach calculation.
#######################################


#######################################
### Function Output
###
### w: A p by length(lambda) matrix with each column corresponding to the coefficient estimate for that lambda
### w.sum: a total sum of all coefficients at each lambda (used as budget calculation for advertising example)
### lam: grid of lambda values used in the estimation
### click.est: vector of click rates if q is not NULL
#######################################

 compare.to.ELMSO<-
  function(z, gamma, step = 0.05, size = 100, q = NULL){
    
    if(is.null(q)) {q = rep(1,dim(z)[1])}
    
    gamma = gamma*q
    
    empty = 0
    
    n = nrow(z)
    p = ncol(z)
    
    az = z
    a2z = z
    azw = z
    
    w.til = rep(0,p)
    w = rep(0,p)
    
    w.mat = w.sum = lam.grid = click.est = NULL
    Q = rep(1,n)
    
    for(j in 1:p){
      az[,j] = gamma[j]*z[,j]
      azw[,j] = gamma[j]*z[,j]*w.til[j]
    }

    
    s = 0
    
    for(j in 1:p){
      
      L1.temp = -az[,j]*(1-gamma[j]*w.til[j])^(z[,j]-1)
      L2.temp = az[,j]^2*(1-gamma[j]*w.til[j])^(z[,j]-2)
      
      S.temp = sum(w.til[j]*L2.temp-L1.temp)
      
      if(s < abs(S.temp)){s = abs(S.temp)}
      
      
    }
    
    lam = (2*s)^0.1
    lam.min = ((2*s)^0.1-step*size)
    
    
    while(lam > lam.min){
      converge = 1
      iter = 0
      while(converge>10^-3 & iter<200){
        
        iter = iter+1
        
        for(j in 1:p){
          
          L1 = -az[,j]*(1-gamma[j]*w.til[j])^(z[,j]-1)
          L2 = az[,j]^2*(1-gamma[j]*w.til[j])^(z[,j]-2)
          
          S = sum(w.til[j]*L2-L1)
          sum.L = sum(L2)
          
          
          if(S=="NaN"){
            
            w.til[j] = w.til[j]+1
            L1 = -az[,j]*(1-gamma[j]*w.til[j])^(z[,j]-1)
            L2= az[,j]^2*(1-gamma[j]*w.til[j])^(z[,j]-2)
            
            S = sum(w.til[j]*L2-L1)
            sum.L = sum(L2)
            
          }
          
          
          if(abs(S) > 10^lam){
            w[j] = w.til[j]-(sum(L1)+10^lam*sign(S))/sum.L
            empty = empty+1

          }
          
          else {
            w[j] = 0
          }
          
          
        }
        
        if(sum(w)==0){
          converge = 0
        }
        
        else{
          converge = sum(abs(w.til-w))/sum(w)
        }
        
        w.til = w
        
      }
      
      if(iter>200){
        
        print("Max iterations reached. Converge value: ",converge)
        
      }
      
      w.mat = cbind(w.mat,w)
      w.sum = cbind(w.sum,sum(w))
      lam.grid = c(lam.grid,10^lam)
      lam = lam-step
      click.est = c(click.est,1/n*sum(Q))
      
    }
    
    list(w = w.mat, w.sum = w.sum, lam = lam.grid, click.est = click.est)
    
  }
 