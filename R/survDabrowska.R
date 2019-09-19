#' @export
###################################################################
######### Dabrowska estimator  ####################################
###################################################################
survDabrowska = function(X, Y, deltaX, deltaY){
  data = data.frame(X = X, Y = Y, deltaX = deltaX, deltaY = deltaY)

  if(any(!c("X", "Y", "deltaX", "deltaY") %in% names(data))) stop("data.frame data has to contain the following variables: X, Y, deltaX, deltaY")
  if(nrow(data)<2) stop("Not enough data")
  X = data[,"X"]
  Y = data[,"Y"]
  deltaX = data[,"deltaX"]
  deltaY = data[,"deltaY"]
  
  ### Compute marginal survival curves:
  ###########################################################
  Sx = myOwnKM(time = X, delta = data[,"deltaX"])
  Sy = myOwnKM(time = Y, delta = data[,"deltaY"])
  SxUnique = unique(Sx[order(Sx$time), c("time", "KM")])
  SyUnique = unique(Sy[order(Sy$time), c("time", "KM")])
  
  ### Compute bivariate and univariate hazards:
  ###########################################################
  nX = nY = length(X)
  uniqueTimeX = unique(unlist(X))
  uniqueTimeY = unique(unlist(Y))
  numXTimes = length(uniqueTimeX)
  numYTimes = length(uniqueTimeY)
  Lam = matrix(NA, nrow=numXTimes, ncol=numYTimes)
  LamX = matrix(NA, nrow=numXTimes, ncol=numYTimes)
  LamY = matrix(NA, nrow=numXTimes, ncol=numYTimes)
  XEventMatr = matrix(NA, nrow=numXTimes, ncol=nX)
  YEventMatr = matrix(NA, nrow=numYTimes, ncol=nY)
  XAtRiskMatr =  matrix(NA, nrow=numXTimes, ncol=nX)
  YAtRiskMatr =  matrix(NA, nrow=numYTimes, ncol=nY)
  for(i in 1:nX){
    XAtRiskMatr[,i] = (X[i] >= uniqueTimeX)
    YAtRiskMatr[,i] = (Y[i] >= uniqueTimeY)
    XEventMatr[,i] = (X[i] == uniqueTimeX)*deltaX[i]
    YEventMatr[,i] = (Y[i] == uniqueTimeY)*deltaY[i]
  }
  atRiskMatr = XAtRiskMatr %*% t(YAtRiskMatr)
  M11 = XEventMatr %*% t(YEventMatr)
  M10 = XEventMatr %*% t(YAtRiskMatr)
  M01 = XAtRiskMatr %*% t(YEventMatr)
  
  ######## double hazard
  ###########################################################
  Lam = M11/atRiskMatr
  
  ######## single hazards
  ###########################################################
  LamX = M10 /atRiskMatr
  LamY = M01 /atRiskMatr

  ######## Dabrowska estimator:
  ###########################################################
  bigL = (LamX*LamY - Lam)/((1 - LamX)*(1 - LamY))
  indicesX = data.frame(rowNames = uniqueTimeX, rowIndex = 1:length(uniqueTimeX))
  indicesY = data.frame(colNames = uniqueTimeY, colIndex = 1:length(uniqueTimeY))
  indicesX = indicesX[order(uniqueTimeX), ]
  indicesY = indicesY[order(uniqueTimeY), ]
  DabrowskaEstNew = matrix(NA, nrow=numXTimes+1, ncol=numYTimes+1)
  DabrowskaEstNew[1, 1] = 1
  DabrowskaEstNew[2:nrow(DabrowskaEstNew), 1] = SxUnique["KM"][,1]
  DabrowskaEstNew[1, 2:ncol(DabrowskaEstNew)] = SyUnique["KM"][,1]
  colnames(DabrowskaEstNew) = c("0", indicesY$colNames)
  rownames(DabrowskaEstNew) = c("0", indicesX$rowNames)
  DabrowskaCDF = 1 - DabrowskaEstNew
  subBigLNonMiss = 1 - bigL
  subBigLNonMiss[is.na(subBigLNonMiss)] = 1
  rowInd = indicesX$rowIndex
  colInd = indicesY$colIndex
  for (i in 1:numXTimes){
    for (j in 1:numYTimes){
      DabrowskaEstNew[i+1, j+1] = (DabrowskaEstNew[i, j+1] * DabrowskaEstNew[i+1, j] / DabrowskaEstNew[i, j]) * subBigLNonMiss[rowInd[i], colInd[j]]
    }
  }
  ### Some values will be NaN because of division by zero (because we DabrowskaEstNew[i, j]
  ### = 0 at some point), so we will assign NA's to 0 
  ###########################################################
  DabrowskaEstNew[is.na(DabrowskaEstNew)] = 0
  
  ### Compute CDF
  ###########################################################
  DabrowskaCDF = 1 + DabrowskaEstNew - matrix(DabrowskaEstNew[1, ], nrow = nrow(DabrowskaEstNew), ncol = ncol(DabrowskaEstNew), byrow = TRUE) - 
  matrix(DabrowskaEstNew[, 1], nrow = nrow(DabrowskaEstNew), ncol = ncol(DabrowskaEstNew))
  list(DabrowskaEst = DabrowskaEstNew, DabrowskaCDF = DabrowskaCDF)
}
