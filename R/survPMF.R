#' @export
###################################################################
######### Get probability mass function from survival surface
###################################################################
survPMF = function(bivarSurf){  
  if(bivarSurf[1,1] != 1)
    {stop("Element bivarSurf[1,1] should be equal to 1")}
  numericX = as.numeric(rownames(bivarSurf))
  if(any(is.na(numericX)) | any(numericX != numericX[order(numericX)]) | numericX[1] != 0)
    {stop("Please check the row names of 'bivarSurf': they should be ORDERED NUMERIC values with FIRST element equal to '0'")}
  numericY = as.numeric(colnames(bivarSurf))
  if(any(is.na(numericY)) | any(numericY != numericY[order(numericY)]) | numericY[1] != 0)
    {stop("Please check the column names of 'bivarSurf': they should be ORDERED NUMERIC values with FIRST element equal to '0'")}
  
  nX = nrow(bivarSurf) - 1
  nY = ncol(bivarSurf) - 1
  unitVecX = matrix(1, ncol = nX)
  unitVecY = matrix(1, ncol = nY)

  Sx = matrix(as.numeric(bivarSurf[(1:nX)+1,1]), nrow = nX) %*% unitVecY
  Sy = t(matrix(as.numeric(bivarSurf[1,(1:nY)+1]), nrow = nY) %*% unitVecX)
  SxM = matrix(as.numeric(bivarSurf[1:nX,1]), nrow = nX) %*% unitVecY
  SyM = t(matrix(as.numeric(bivarSurf[1,1:nY]), nrow = nY) %*% unitVecX)
  Sdx = SxM - Sx
  Sdy = SyM - Sy

  SxMyM = rbind(rep(NA, nY + 1), cbind(rep(NA, nX), bivarSurf[1:nX, 1:nY]))
  SxM_y = rbind(rep(NA, nY+1), bivarSurf[1:nX, ])
  Sx_yM = cbind(rep(NA, nX+1), bivarSurf[, 1:nY])
  Sdx_y = SxM_y - bivarSurf
  Sx_dy = Sx_yM - bivarSurf
  Sdx_yM = SxMyM - Sx_yM
  SxM_dy = SxMyM - SxM_y
  Sdxdy = SxM_dy - Sx_dy
  rangeX = 1+(1:nX)
  rangeY = 1+(1:nY)
  
  pmfWithMarginalPMFs = Sdxdy
  rownames(pmfWithMarginalPMFs) = rownames(bivarSurf)
  colnames(pmfWithMarginalPMFs) = colnames(bivarSurf)
  pmfWithMarginalPMFs[,1] = c(0, Sdx[,1])
  pmfWithMarginalPMFs[1,] = c(0, Sdy[1,])
  
  returnRes = list(Sdxdy = pmfWithMarginalPMFs, Sxy = bivarSurf[rangeX, rangeY], SxMyM = SxMyM[rangeX, rangeY],  Sx = Sx, Sy = Sy, Sdx = Sdx, Sdy = Sdy, SxM = SxM, SyM = SyM, SxM_y = SxM_y[rangeX, rangeY], Sx_yM = Sx_yM[rangeX, rangeY], Sdx_y = Sdx_y[rangeX, rangeY], Sx_dy = Sx_dy[rangeX, rangeY], Sdx_yM = Sdx_yM[rangeX, rangeY], SxM_dy = SxM_dy[rangeX, rangeY])

  returnRes
}
