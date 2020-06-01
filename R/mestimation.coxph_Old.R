#' @importFrom stats model.matrix residuals
###############################################
############################################### OLD M-estimation for Cox regression
###############################################
mestimation.coxph_Old = function(object, ...){
  ### for non-time dependent covariates
  
  #stopifnot(require("sandwich"))
    
if (requireNamespace("sandwich", quietly = TRUE)) {
  dl.dtheta = t(sandwich::estfun(object))   ### similar to mean(.) function from sandwich library
  
  ######################################## d2l.dtheta.dtheta (matrix A)
  d2l.dtheta.dtheta = solve(object$var)
  
  ######################################## prob. scale residuals (PSR)
#  presid = presid.coxph(object)
  presid = PResiduals::presid(object)
  
  ######################################## PSR derivatives by theta (dpresid.dtheta)
  Y = object$y[,1]
  nSubj = length(Y)
  delta = object$y[,2]
  orderedY = Y[order(Y)]
  orderedDelta = delta[order(Y)]
  
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
  dpresid.dtheta = matrix(NA, nrow = nParam, ncol = nSubj)
  
  ######################################## PSR derivatives by beta  
  H0 = survival::basehaz(object, centered=FALSE)
  H0_minus = hazardStretcher_take_two(H0, orderedY)
  H0_minus$delta = orderedDelta
  uniqueBlHaz = H0_minus[H0_minus$delta == 1,]
  uniqueBlHaz = unique(uniqueBlHaz[, c("hazard", "time")])
  timeDiff = diff(c(0, uniqueBlHaz[, "time"]))
  lambdas = diff(c(0, uniqueBlHaz[, "hazard"]))/timeDiff
  element1 = exp(-H0_minus[, "hazard"] * orderedExpXbeta)
  element2 = exp(-H0_minus[, "Hminus"] * orderedExpXbeta)
  if(nParam>=1){
    for (j in 1:nParam){
      dpresid.dtheta[j, ] = (element1*H0_minus[, "hazard"] + orderedDelta * element2 * H0_minus[, "Hminus"]) * orderedExpXbeta * orderedX1[, j]
    }
  }   
  orderMatrix = matrix(c(1:nSubj, (1:nSubj)[order(Y)]), ncol=2)
  colnames(orderMatrix) = c("new", "old")
  originalOrder = orderMatrix[order(orderMatrix[,2]), 1]
  dpresid.dtheta = dpresid.dtheta[, originalOrder]
  
  res = list(dl.dtheta = dl.dtheta,
       d2l.dtheta.dtheta = d2l.dtheta.dtheta,
       resid = NULL,
       dresid.dtheta = NULL,
       presid = presid,
       presid.k= NULL,
       dpresid.dtheta = dpresid.dtheta,
       dpresid.dtheta.k = NULL)
  res
  }else{
    cat("Package 'sandwich' is  required.\n")
    return(NULL)
  }
}

