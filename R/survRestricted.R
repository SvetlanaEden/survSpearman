#' @export
#################################################################
######### Conditional survival surface (final version)
#################################################################
survRestricted = function(bivarSurf, tauX = Inf, tauY = Inf){
  Xs = as.numeric(rownames(bivarSurf))
  Ys = as.numeric(colnames(bivarSurf))
  if (tauX < min(Xs)) stop("tauY is out of bounds.")
  if (tauY < min(Ys)) stop("tauY is out of bounds.")

  deltaDiff = min(abs(c(diff(Xs), diff(Ys))))/2
  if(is.null(tauX)){tauX = max(Xs) + deltaDiff}
  if(is.null(tauY)){tauY = max(Ys) + deltaDiff}

  dimX = nrow(bivarSurf)
  dimY = ncol(bivarSurf)
  
  Fxy = bivarSurf*0
  Fxy[1,1] = 0
  Fxy[1, 2:dimY] = 1 - bivarSurf[1, 2:dimY]
  Fxy[2:dimX, 1] = 1 - bivarSurf[2:dimX, 1]
  Fxy[2:dimX, 2:dimY] = 1 - matrix(bivarSurf[1, 2:dimY], nrow = dimX-1, ncol = dimY-1, byrow = TRUE) - matrix(bivarSurf[2:dimX, 1], nrow = dimX-1, ncol = dimY-1) + bivarSurf[2:dimX, 2:dimY]
  
  transition = Fxy[as.numeric(rownames(Fxy)) < tauX, ]
  transition = transition[, as.numeric(colnames(transition)) < tauY]
  if(is.null(dim(transition))){
    stop("Not enough data points or the restricted region is too small.")
  }
  if(all(abs(transition[2:nrow(transition), 2:ncol(transition)]) < 10^(-13))){
    stop("Not enough data points or the restricted region is too small.")
  }  
  FxyCond = transition

  FxyCond[2:nrow(FxyCond), 2:ncol(FxyCond)] = FxyCond[2:nrow(FxyCond), 2:ncol(FxyCond)]/transition[nrow(transition), ncol(transition)]
  
  if(transition[nrow(transition), ncol(transition)] <= 0){
    cat("The probability of being in the restricted region is less or equal to zero, so the conditional survival probability cannot be computed.")
    res = list(Sxy = NA, SxMyM = NA,  Sx = NA, Sy = NA, Sdx = NA, Sdy = NA, SxM = NA, SyM = NA, SxM_y = NA, Sx_yM = NA, Sdx_y = NA, Sx_dy = NA, Sdx_yM = NA, SxM_dy = NA, Sdxdy = NA)
    return(res)
  }
  
  ### Computing the joint probability mass:
  #################################################################
  FxyCond[, 1] = 0
  FxyCond[1, ] = 0
  massMatr = FxyCond
  for(i in 2:nrow(FxyCond)){
    for(j in 2:ncol(FxyCond)){
      massMatr[i, j] = FxyCond[i, j] - FxyCond[i-1, j] - FxyCond[i, j-1] + FxyCond[i-1, j-1]
    }
  }
  
  ### Computing the marginal probability mass:
  #################################################################
  reducedRows = massMatr[2:nrow(massMatr), ]
  if(is.null(dim(reducedRows))){reducedRows = matrix(reducedRows, nrow = nrow(massMatr)-1, ncol = ncol(massMatr))}
  massMatr[1, ] = apply(reducedRows, 2, sum)
  reducedCols = massMatr[, 2:ncol(massMatr)]
  if(is.null(dim(reducedCols))){reducedCols = matrix(reducedCols, nrow = nrow(massMatr), ncol = ncol(massMatr) - 1)}
  massMatr[, 1] = apply(reducedCols, 1, sum)
  massMatr[1, 1] = 0
  
  ### Finishing conditional CDF:
  #################################################################
  FxyCond[1, ] = cumsum(massMatr[1, ])
  FxyCond[, 1] = cumsum(massMatr[, 1])
  
  ### Computing conditional survival surface:
  #################################################################
  SurvCond = FxyCond*0
  SurvCond[1,1] = 1
  SurvCond[2:nrow(SurvCond), 1] = 1 - FxyCond[2:nrow(FxyCond), 1]
  SurvCond[1, 2:ncol(SurvCond)] = 1 - FxyCond[1, 2:ncol(FxyCond)]
  SurvCond[2:nrow(SurvCond), 2:ncol(SurvCond)] = 1 - matrix(FxyCond[1, 2:ncol(SurvCond)], nrow = nrow(SurvCond)-1, ncol = ncol(SurvCond)-1, byrow = TRUE)- matrix(FxyCond[2:nrow(SurvCond), 1], nrow = nrow(SurvCond)-1, ncol = ncol(SurvCond)-1) + FxyCond[2:nrow(SurvCond), 2:ncol(SurvCond)]

  ### Prepare the rest of the surfaces
  #################################################################
  newDimX = nrow(SurvCond); newDimY = ncol(SurvCond)
  baseMatr = matrix(0, nrow = newDimX - 1, ncol = newDimY - 1)
  rownames(baseMatr) = rownames(SurvCond)[2:newDimX]
  colnames(baseMatr) = colnames(SurvCond)[2:newDimY]
  unitVecX = matrix(1, ncol = newDimX - 1)
  unitVecY = matrix(1, ncol = newDimY - 1)
  Sxy = Sx = Sy = Sdx = Sdy = SxM = SyM = SxM_y = Sx_yM = Sdx_y = Sx_dy = Sdx_yM = SxM_dy = SxMyM = Fxy = baseMatr
  
  Sxy = SurvCond[2:nrow(SurvCond), 2:ncol(SurvCond)]
  Sx = SurvCond[2:nrow(SurvCond), 1] %*% unitVecY
  Sy = t(SurvCond[1, 2:ncol(SurvCond)] %*% unitVecX)

  SxMyM = SurvCond[1:(newDimX-1), 1:(newDimY-1)]

  Sdx = -diff(c(1, Sx[, 1])) %*% unitVecY
  Sdy = t(-diff(c(1, Sy[1, ])) %*% unitVecX)
  SxM = Sdx + Sx; SyM = Sdy + Sy
  
  SxM_y = rbind(Sy[1, ], Sxy[1:(newDimX-2), ])
  Sx_yM = cbind(Sx[, 1], Sxy[, 1:(newDimY-2)])
  Sdx_y = SxM_y - Sxy; Sx_dy = Sx_yM - Sxy 
  Sdx_yM = SxMyM - Sx_yM; SxM_dy = SxMyM - SxM_y
  
  Sdxdy = Sdx_yM - Sdx_y
  
  res = list(Sxy = SurvCond, SxMyM = SxMyM,  Sx = Sx, Sy = Sy, Sdx = Sdx, Sdy = Sdy, SxM = SxM, SyM = SyM, SxM_y = SxM_y, Sx_yM = Sx_yM, Sdx_y = Sdx_y, Sx_dy = Sx_dy, Sdx_yM = Sdx_yM, SxM_dy = SxM_dy, Sdxdy = Sdxdy)
  
  for(n_i in setdiff(names(res), "Sxy")){
    rownames(res[[n_i]]) = NULL
    colnames(res[[n_i]]) = NULL
  }
  
  res
}

