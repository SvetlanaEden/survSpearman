############################################### estimating equations for  
############################################### parametric regressions
mestimation.survregPSRs = function(object, inverseA = TRUE, beta = NULL, ...){
  ################### for objectects of type survreg.scores
  ################### returns the M-estimation stuff
  ### object - is the fitted survival model
    
  ### similar to mean(.) function from sandwich library
  dl.dtheta = t(sandwich::estfun(object))
  
  ######################################## d2l.dtheta.dtheta (matrix A)
  if(inverseA){
    d2l.dtheta.dtheta = object$var
  }else{
    d2l.dtheta.dtheta = solve(object$var)
  }
  
  ######################################## prob. scale residuals (PSR)
  ### When presidTMP.survreg() is merged with presid.survreg() in PResiduals
  ### Replace presid = presidTMP.survreg(object) with presid = presid(object)
  #presid = presid(object)
  presid = presidTMP.survreg(object)
  
  ######################################## prepare stuff for
  ######################################## PSR derivatives by theta (dpresid.dtheta)
  X1 = model.matrix(object)
  if(is.null(beta)){
    Xbeta =  object$linear.predictors
  }else{
    Xbeta =  X1 %*% matrix(beta, ncol=1)
  }
  nSubj = length(Xbeta)
  nParam = dim(object$var)[1]
  delta = object$y[,2]
  dpresid.dtheta = matrix(NA, nrow = nParam, ncol = nSubj)

  switch(object$dist,
    ######################################## exponential
    exponential = {
      Y = exp(object$y[,1])
      dpresid.dtheta[1,] = -(1+delta)*Y*exp(-exp(-Xbeta)*Y - Xbeta)
      if(nParam>=2){
        for (i in 2:nParam){
          dpresid.dtheta[i, ] = dpresid.dtheta[1,] * X1[, i]
        }
      }
    },
 
    ######################################## weibull
    weibull = {
      scale = object$scale
      Y = exp(object$y[,1])
      dpresid.dtheta[1,] = -(1+delta)*(Y^(1/scale))*exp(-(Y*exp(-Xbeta))^(1/scale) - Xbeta/scale)/scale
      if(nParam-1 >= 2){
        for (i in 2:(nParam-1)){
          dpresid.dtheta[i, ] = dpresid.dtheta[1,] * X1[, i]
        }
      }
      dpresid.dtheta[nParam, ] = -(1+delta)*(log(Y) - Xbeta)*(Y^(1/scale))*exp(-(Y*exp(-Xbeta))^(1/scale) - Xbeta/scale - log(scale))
    },
 
    ######################################## log-logistic
    loglogistic = {
      Y = exp(object$y[,1])
      gamma = (1/object$scale)
      monsterTerm = (Y^gamma)*exp(-Xbeta*gamma)
      dpresid.dtheta[1,] = -gamma*(1+delta)*monsterTerm/((1+monsterTerm)^2)
      if(nParam-1 >= 2){
        for (i in 2:(nParam-1)){
          dpresid.dtheta[i, ] = dpresid.dtheta[1,] * X1[, i]
        }
      }
      dpresid.dtheta[nParam, ] = -dpresid.dtheta[1,] * (log(Y) - Xbeta)
    },

    ######################################## log-normal
    lognormal = {
      logY = object$y[,1]
      scale = summary(object)$scale
      dpresid.dtheta[1,] = -(1+delta)*dnorm(logY, mean = Xbeta, sd = scale)
      if(nParam-1 >= 2){
        for (i in 2:(nParam-1)){
          dpresid.dtheta[i, ] = dpresid.dtheta[1,] * X1[, i]
        }
      }
      dpresid.dtheta[nParam, ] = (1+delta)* (-(logY - Xbeta)*dnorm(logY, mean = Xbeta, sd = scale)  +  pnorm(logY, mean = Xbeta, sd = scale)  -  plnorm(exp(logY), meanlog=Xbeta, sdlog = scale, log.p=FALSE))
      },
    ### Stop when unknown distribution:
    stop("Unhandled dist", object$dist)
  ) ### end of switch
         
  res = list(dl.dtheta = dl.dtheta,
       d2l.dtheta.dtheta = d2l.dtheta.dtheta,
       resid = NULL,
       dresid.dtheta = NULL,
       presid = presid,
       presid.k= NULL,
       dpresid.dtheta = dpresid.dtheta,
       dpresid.dtheta.k = NULL)
}



