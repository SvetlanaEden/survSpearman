###############################################
############################################### M-estimation for Cox, full likelihood
###############################################
mestimation.coxph.full.likelihood = function(object, time, continuous = TRUE, ...){
  ### in general, argument time can be extracted from object: time = object$y[, 1]
  ### but sometimes, numeric issues get in the way
  ################################# M-estimation for Cox  
  Y = object$y[,1]
  nSubj = length(Y)
  delta = object$y[,2]
  resid <- residuals(object, type="martingale")
  orderedY = Y[order(Y)]
  orderedDelta = delta[order(Y)]
  orderedUniqueFailureTimes = unique(orderedY[orderedDelta == 1])
  orderedFailureTimesWithZero = c(0, orderedUniqueFailureTimes)
  nHazardParam = length(orderedUniqueFailureTimes)
  orderedMartRes = resid[order(Y)]
  
  orderMatrix = matrix(c(1:nSubj, (1:nSubj)[order(Y)]), ncol=2)
  colnames(orderMatrix) = c("new", "old")
  originalOrder = orderMatrix[order(orderMatrix[,2]), 1]

  X1 = model.matrix(object)
  if(dim(X1)[2] == 1){
    orderedX1 = as.matrix(X1[order(Y),], ncol = 1)
  }else{
    orderedX1 = X1[order(Y),]
  }
  Xbeta =  X1 %*% matrix(object$coefficients, nrow = length(object$coefficients))
  expXbeta = exp(Xbeta)
  orderedExpXbeta = expXbeta[order(Y)]

  nParam = dim(object$var)[1]
  if(is.null(nParam)){
    nParam = 1
  }
  
  H0 = survival::basehaz(object, centered=FALSE)
  H0_minus = hazardStretcher_take_two(H0, Y)
  H0_minus$delta = orderedDelta
  uniqueBlHaz = H0_minus[H0_minus$delta == 1,]
  uniqueBlHaz = unique(uniqueBlHaz[, c("hazard", "time")])
  timeDiff = diff(c(0, uniqueBlHaz[, "time"]))
  lambdas = diff(c(0, uniqueBlHaz[, "hazard"]))/timeDiff
  
  ######################################## dl.dtheta
  ######################################## dl.dtheta
  ######################################## dl.dtheta
  dl.dtheta = matrix(NA, nrow = nParam + nHazardParam, ncol = nSubj)
  ######################################## dl.dtheta   for beta:
  for (j in 1:nParam){
    dl.dtheta[j, ] = (orderedDelta - H0_minus[,"hazard"]*orderedExpXbeta)*orderedX1[,j]
  }  
  ######################################## dl.dtheta   for lambda:
  dl.dlambda = matrix(0, nrow=nHazardParam, ncol = nSubj)
  for (k in 1:nHazardParam){
    hadEventAtT = (orderedDelta == 1 & orderedY == orderedFailureTimesWithZero[k+1])
    hadEventAfterT = (orderedY > orderedFailureTimesWithZero[k+1] | (orderedDelta == 0 & orderedY == orderedFailureTimesWithZero[k+1]))
    dl.dlambda[k, hadEventAtT] = orderedDelta[hadEventAtT]/lambdas[k] - timeDiff[k] * orderedExpXbeta[hadEventAtT]
    dl.dlambda[k, hadEventAfterT] = -timeDiff[k] * orderedExpXbeta[hadEventAfterT]
  }
  ######################################## merge dl.dtheta for betand and for lambda:
  dl.dtheta[(nParam+1):(nParam+nHazardParam), ] = dl.dlambda
  
  ######################################## d2l.dtheta.dtheta (matrix A)
  ######################################## d2l.dtheta.dtheta (matrix A)
  ######################################## d2l.dtheta.dtheta (matrix A)
  d2l.dtheta.dtheta = diag(nParam + nHazardParam)
  ######################################## d2l.dtheta.dtheta  for beta/beta:
  for (i in 1:nParam){
    for (j in 1:i){
      element = H0_minus[,"hazard"] * orderedExpXbeta * orderedX1[,i] * orderedX1[,j]
      d2l.dtheta.dtheta[i, j] = d2l.dtheta.dtheta[j, i] = -sum(element)
    }
  }  
  ######################################## d2l.dtheta.dtheta  for lambda/lambda:
  for (k in (1:nHazardParam)){
    hadEventAtT = (orderedY == orderedFailureTimesWithZero[k+1])
    ############################################
    element = orderedDelta/(lambdas[k]^2)
    d2l.dtheta.dtheta[nParam+k, nParam+k] = -sum(element[hadEventAtT])
  }

  ######################################## d2l.dtheta.dtheta  for beta/lambda:
  for (j in 1:nParam){
    for (k in 1:nHazardParam){
      greaterOrEqual = (orderedY >= orderedFailureTimesWithZero[k+1])
      element = timeDiff[k] * orderedExpXbeta * orderedX1[, j]
      d2l.dtheta.dtheta[j, nParam + k] = d2l.dtheta.dtheta[nParam + k, j] = -sum(element[greaterOrEqual])
    }
  }  
  
  ######################################## prob. scale residuals (PSR)
  presid = PResiduals::presid(object, continuous = continuous)
  orderedPresid = presid[order(Y)]         #### ORDERING BY Y
  
  ######################################## PSR derivatives by beta and lambda_k
  ######################################## PSR derivatives by beta and lambda_k
  dpresid.dtheta = matrix(0, nrow = nParam + nHazardParam, ncol = nSubj)
  if(continuous){
    ######################################## PSR derivatives by beta  
    expression1 = (1 + orderedDelta)*exp(orderedMartRes - orderedDelta)*(orderedDelta-orderedMartRes)
    dpresid.dtheta[1:nParam, ] = matrix(expression1, nrow = nParam, ncol = nSubj, byrow = TRUE) * t(orderedX1)

    ######################################## PSR derivatives by lambda_k
    expression2 = (1 + orderedDelta)*exp(orderedMartRes - orderedDelta) * orderedExpXbeta
    for(k in 1:nHazardParam){
      condGreaterOrEqual = (orderedY >= orderedFailureTimesWithZero[k+1])
      dpresid.dtheta[nParam + k, condGreaterOrEqual]  = timeDiff[k] * expression2[condGreaterOrEqual]
    }
  }else{
    element1 = exp(-H0_minus[, "hazard"] * orderedExpXbeta)
    element2 = exp(-H0_minus[, "Hminus"] * orderedExpXbeta)
    element3 = element1 * orderedExpXbeta
    element4 = element2 * orderedExpXbeta

    ######################################## PSR derivatives by beta  
    if(nParam>=1){
      for (j in 1:nParam){
        dpresid.dtheta[j, ] = (element1*H0_minus[, "hazard"] + orderedDelta * element2 * H0_minus[, "Hminus"]) * orderedExpXbeta * orderedX1[, j]
      }
    }

    ######################################## PSR derivatives by lambda_k
    for(k in 1:nHazardParam){
      condGreaterOrEqual = (orderedY >= orderedFailureTimesWithZero[k+1])
      dpresid.dtheta[nParam + k, condGreaterOrEqual]  = (timeDiff[k] * (element3 + orderedDelta * element4) )[condGreaterOrEqual]
    }
  }
  
  res = list(dl.dtheta = dl.dtheta[, originalOrder],
       d2l.dtheta.dtheta = d2l.dtheta.dtheta,
       resid = NULL,
       dresid.dtheta = NULL,
       presid = orderedPresid[originalOrder],
       presid.k= NULL,
       dpresid.dtheta = dpresid.dtheta[, originalOrder],
       dpresid.dtheta.k = NULL)

  res

}

