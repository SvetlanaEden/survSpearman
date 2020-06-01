#' @importFrom stats quantile density lm
###################################################################################
##### function for obtaining the standard error of conditional estimate
#####    using large sample approximation, author: Qi Liu and Svetlana Eden
#####    Qi originally wrote this function. Svetlana modified it 
#####    to be able to use splines with z and also
#####    to get predicted values and CIs right away
###################################################################################
prodTS.new.SKE = function(x.resid, y.resid, z, numKnots = 3, newZ = z,
                       xres2.method=c("emp", "constant", "model"),
                       yres2.method=c("emp", "constant", "model"),
                       xz.dl.dtheta, yz.dl.dtheta,
                       xz.d2l.dtheta.dtheta, yz.d2l.dtheta.dtheta,
                       dxresid.dtheta, dyresid.dtheta, debugMode = FALSE, ci = 0.95){

  if (!any(numKnots == c(-1, 0, 3, 4, 5, 6, 7))) stop("The number of knots should be either 0, 3, 4, or 5")
  if(numKnots == -1){
    designMatrOfZ = matrixForZ = matrix(1, ncol = 1, nrow = length(z))
    newDesignMatrOfZ = matrixForNewZ = matrix(1, ncol = 1, nrow = length(newZ))
  }else{
    if(numKnots == 0){
      matrixForZ = matrix(z, ncol = 1)
      matrixForNewZ = matrix(newZ, ncol = 1)
    }else{
      splinePercentileList = list("3" = c(.1, .5, .9), "4" = c(.05, .35, .65, .95), "5" = c(.05, .275, .5, .725, .95), "6" = c(.05, .23, .41, .59, .77, .95), "7" = c(.05, .1833, .3417, .5, .6583, .8167, .975)) ### see Regression Modeling Strategies by F. Harrell.
      knots = quantile(z, splinePercentileList[[as.character(numKnots)]], na.rm = TRUE)
      matrixForZ = componentsOfSplines(z, knots)
      matrixForNewZ = componentsOfSplines(newZ, knots)
    }
    designMatrOfZ <- cbind(1, matrixForZ)
    newDesignMatrOfZ <- cbind(1, matrixForNewZ)
  }  
  npar.xz <- dim(xz.dl.dtheta)[2]
  npar.yz <- dim(yz.dl.dtheta)[2]

  score.prod <- lm.scoresForBivarSurvival(x.resid*y.resid, matrixForZ)

  npar.prod <- dim(score.prod$dl.dtheta)[2]
  prod.dl.dtheta <- score.prod$dl.dtheta
  prod.d2l.dtheta.dtheta <- score.prod$d2l.dtheta.dtheta
  N <- dim(score.prod$dl.dtheta)[1]
  
  if(xres2.method[1]=="model"){
    score.xres <- lm.scoresForBivarSurvival(x.resid, matrixForZ)
    score.xres2 <- lm.scoresForBivarSurvival(x.resid^2, matrixForZ)
    npar.xres2 <- dim(score.xres2$dl.dtheta)[2]
    xres2.dl.dtheta <- score.xres2$dl.dtheta
    xres2.d2l.dtheta.dtheta <- score.xres2$d2l.dtheta.dtheta
    
  } else if (xres2.method[1]=="emp"){
    xres2.dl.dtheta <- as.matrix(x.resid^2-mean(x.resid^2) )
    npar.xres2 <- 1
    xres2.d2l.dtheta.dtheta <- -N
    
  } else if (xres2.method[1]=="constant"){
    npar.xres2 <- 0
    xres2.dl.dtheta <- NULL
    xres2.d2l.dtheta.dtheta <- NULL
  }
  
  if(yres2.method[1]=="model"){
    score.yres <- lm.scoresForBivarSurvival(y.resid, matrixForZ)
    score.yres2 <- lm.scoresForBivarSurvival(y.resid^2, matrixForZ)
    npar.yres2 <- dim(score.yres2$dl.dtheta)[2]
    yres2.dl.dtheta <- score.yres2$dl.dtheta
    yres2.d2l.dtheta.dtheta <- score.yres2$d2l.dtheta.dtheta
    
  } else if (yres2.method[1]=="emp"){
    yres2.dl.dtheta <- as.matrix( y.resid^2-mean(y.resid^2)  )
    npar.yres2 <- 1
    yres2.d2l.dtheta.dtheta <- -N
    
  } else if (yres2.method[1]=="constant"){
    npar.yres2 <- 0
    yres2.dl.dtheta <- NULL
    yres2.d2l.dtheta.dtheta <- NULL
  }
  
  Ntheta <- npar.xz + npar.yz + npar.prod + npar.xres2 + npar.yres2
  bigphi <- cbind(xz.dl.dtheta,
                  yz.dl.dtheta,
                  prod.dl.dtheta,
                  xres2.dl.dtheta,
                  yres2.dl.dtheta)
  
  A <- matrix(0, Ntheta, Ntheta)
  ### initiate the diagonal of A 
  A[(1:npar.xz), (1:npar.xz)] <- as.matrix(xz.d2l.dtheta.dtheta)
  A[npar.xz + (1:npar.yz), npar.xz + (1:npar.yz)] <- as.matrix(yz.d2l.dtheta.dtheta)
  A[npar.xz + npar.yz + (1: npar.prod), npar.xz + npar.yz + (1: npar.prod)] <- prod.d2l.dtheta.dtheta
  if (npar.xres2>0){
    A[npar.xz + npar.yz + npar.prod+ (1: npar.xres2), npar.xz + npar.yz + npar.prod+ (1: npar.xres2)] <- xres2.d2l.dtheta.dtheta
  }
  if(npar.yres2>0){
    A[npar.xz + npar.yz + npar.prod+ npar.xres2+(1: npar.yres2), npar.xz + npar.yz + npar.prod+ npar.xres2+(1: npar.yres2)] <- yres2.d2l.dtheta.dtheta
  }
  par.1 <- t(y.resid * designMatrOfZ) %*% t(dxresid.dtheta)
  par.2 <- t(x.resid * designMatrOfZ) %*% t(dyresid.dtheta)
  A[npar.xz + npar.yz + (1: npar.prod), 1:npar.xz] <- par.1
  A[npar.xz + npar.yz + (1: npar.prod), npar.xz + (1:npar.yz)] <- par.2
  
  if (xres2.method[1]=="model"){
    par.3 <- t(2*x.resid * designMatrOfZ) %*% t(dxresid.dtheta)
  } else if (xres2.method[1]=="emp"){
    par.3 <- t(2*x.resid ) %*% t(dxresid.dtheta)
  } else if(xres2.method[1]=="constant"){
    par.3 <- NULL
  }
  
  if(npar.xres2>0){
    A[npar.xz + npar.yz + npar.prod+ (1: npar.xres2), 1:npar.xz] <- par.3
  }
  
  if (yres2.method[1]=="model"){
    par.4 <- t(2*y.resid * designMatrOfZ) %*% t(dyresid.dtheta)
  } else if(yres2.method[1]=="emp"){
    par.4 <- t(2*y.resid ) %*% t(dyresid.dtheta)
  } else if(yres2.method[1]=="constant"){
    par.4 <- NULL
  }
  
  if(npar.yres2>0){
    A[npar.xz + npar.yz + npar.prod+ npar.xres2+(1:npar.yres2), npar.xz + (1:npar.yz)] <- par.4
  }
  
  if(debugMode){
    # A = diag(1, nrow(A))
    #require(matlib)
    matrB = t(bigphi) %*% bigphi
    SS = MASS::ginv(A) %*% matrB %*% t(MASS::ginv(A))
    warning("Used generalized inverse")
  }else{
    SS <- solve(A, t(bigphi))
  }
  var.theta <- tcrossprod(SS, SS)

  prod.coef <- score.prod$mod$coef
  names(prod.coef) <- paste("prod:", names(prod.coef))
  if(xres2.method[1]=="model"){
    xres2.coef <- score.xres2$mod$coef
  } else if (xres2.method[1]=="emp"){
    xres2.coef <- mean(x.resid^2)
  } else if (xres2.method[1]=="constant"){
    xres2.coef <- 1/3
  }
  names(xres2.coef) <- paste("xres2:", names(xres2.coef))
  
  if(yres2.method[1]=="model"){
    yres2.coef <- score.yres2$mod$coef
  } else if (yres2.method[1]=="emp"){
    yres2.coef <- mean(y.resid^2)
  } else if (yres2.method[1]=="constant"){
    yres2.coef <- 1/3
  }
  names(yres2.coef) <- paste("yres2:", names(yres2.coef))
  coef <- rbind(prod.coef, xres2.coef, yres2.coef)
  
  var.coef <- var.theta[(npar.xz+npar.yz+1): Ntheta , (npar.xz+npar.yz+1): Ntheta]
  
  ####### compute predicted values
  predYX = matrix(prod.coef, nrow = 1) %*% t(newDesignMatrOfZ)
  predX2 = matrix(xres2.coef, nrow = 1) %*% t(newDesignMatrOfZ)
  predY2 = matrix(yres2.coef, nrow = 1) %*% t(newDesignMatrOfZ)
  predRho = predYX/sqrt(predX2*predY2)
  
  ####### compute derivatives for delta method:
  ### npar.prod + npar.xres2 + npar.yres2
  derivForDeltaMethod = cbind(newDesignMatrOfZ, newDesignMatrOfZ, newDesignMatrOfZ)
  ### divide by the denominator
  derivForDeltaMethod = derivForDeltaMethod/matrix(rep(sqrt(predX2 * predY2), npar.prod + npar.xres2 + npar.yres2), nrow = nrow(derivForDeltaMethod), ncol = npar.prod + npar.xres2 + npar.yres2)

  ### multiply by -(1/2) and by predYX
  derivForDeltaMethod[, npar.prod + (1 : (npar.xres2 + npar.yres2))] = -(1/2)*derivForDeltaMethod[, npar.prod + (1 : (npar.xres2 + npar.yres2))] * matrix(rep(predYX, npar.xres2 + npar.yres2), nrow = nrow(derivForDeltaMethod), ncol = npar.xres2 + npar.yres2)
  
  ### divide by predX2
  derivForDeltaMethod[, npar.prod + (1 : npar.xres2)] = derivForDeltaMethod[, npar.prod + (1 : npar.xres2)] / matrix(rep(predX2, npar.xres2), nrow = nrow(derivForDeltaMethod), ncol = npar.xres2) 
  
  ### divide by predY2
  derivForDeltaMethod[, npar.prod + npar.xres2 + (1 : npar.yres2)] = derivForDeltaMethod[, npar.prod + npar.xres2 + (1 : npar.yres2)] / matrix(rep(predY2, npar.yres2), nrow = nrow(derivForDeltaMethod), ncol = npar.yres2) 
  
  resultingVar = diag(derivForDeltaMethod %*% var.coef %*% t(derivForDeltaMethod))
  
  
  pointEstAndCIs = data.frame(pointEst = as.vector(predRho), lower = as.vector(predRho - abs(qnorm(0.5*(1-ci)))*sqrt(resultingVar)), upper = as.vector(predRho + abs(qnorm(0.5*(1-ci)))*sqrt(resultingVar)))
    
  result <- list(prod=prod.coef,
                 xres2=xres2.coef,
                 yres2=yres2.coef,
                 var.coef=var.coef,
                 pointEstAndCIs = pointEstAndCIs)
                 
  return(result)
}

