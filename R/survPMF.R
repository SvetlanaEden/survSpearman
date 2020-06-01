#' @name survPMF
#' @aliases survPMF
#' @title Computes marginal and joint probability mass function from marginal and joint survival probabilities.
#' @description The function computes marginal and joint probability mass functions from marginal and joint survival probabilities.
#' 
#' @usage survPMF(bivarSurf)
#' 
#' @param bivarSurf A matrix containing the marginal and joint survival probabilities. The first column is the marginal survival probability corresponding to variable \code{X}. The first row is the marginal survival probability corresponding to variable \code{Y}. The rest of the matrix contains the joint survival probabilities. The row names of \code{bivarSurf} are ordered \code{X} values. The column names of \code{bivarSurf} are ordered \code{Y} values. Element \code{bivarSurf[1,1]} equals 1. Its row and column name is \code{'0'} (see the documentation for the return value \code{DabrowskaEst} in function \code{survDabrowska}).
#' @return The following list of survival surfaces and their differentials is returned. \code{Sdxdy} is the marginal and joint probability mass functions in the same format as argument \code{bivarSurf}; \code{Sxy} is the joint survival probability; \code{SxMyM} is \code{Sxy} at point \code{(x-, y-)}, where \code{x-} is the left limit of \code{x}; \code{Sx} is the marginal survival probability function for variable X; \code{Sy} is the marginal survival probability function for variable Y; \code{Sdx} is the marginal probability mass function for variable X; \code{Sdy} is the marginal probability mass function for variable Y; \code{SxM} is the marginal survival probability function for X at point \code{x-}; \code{SyM} is the marginal survival probability function for Y at point \code{y-}; \code{SxM_y} is the joint survival probability function at point \code{(x-, y)}; \code{Sx_yM} is the joint survival probability function at point \code{(x, y-)}; \code{Sdx_y} is \code{SxM_y - Sxy}; \code{Sx_dy} is \code{Sx_yM - Sxy}; \code{Sdx_yM} is \code{SxMyM - Sx_yM}; \code{SxM_dy} is \code{SxMyM - SxM_y}.

#' @details The function returns a list of survival surfaces and their differentials. Element \code{Sdxdy} of this list is the marginal and joint probability mass function in the same format as argument \code{bivarSurf}. The rest of the returned list elements are matrices in the same format as \code{bivarSurf} except that they do not contain marginal values and row/column names.
#' 
#' @examples
#' X = c(0.5, 0.6, 0.8)
#' Y = c(0.44, 0.77, 0.99)
#' deltaX = c(1, 0, 1)
#' deltaY = c(1, 1, 1)
#' 
#' bivarSurf = survDabrowska(X, Y, deltaX, deltaY)$DabrowskaEst
#' bivarSurf
#' 
#' bivarPMF = survPMF(bivarSurf)$Sdxdy
#' bivarPMF
#' 
#' @keywords bivariate survival PMF
#' @references Eden, S., Li, C., Shepherd B. (2020). Non-parametric Estimation of Spearman's Rank Correlation with Bivariate Survival Data.
#' @author Svetlana K Eden, \email{svetlanaeden@gmail.com}
#' 
#' @export
###################################################################
######### Get probability mass function from survival surface
###################################################################
survPMF = function(bivarSurf){  
  if(bivarSurf[1,1] != 1)
    {stop("Element bivarSurf[1,1] should be equal to 1")}
  numericX = as.numeric(rownames(bivarSurf))
  if(any(is.na(numericX)) | any(numericX != numericX[order(numericX)]) | numericX[1] != 0)
    {stop("Please check the row names of 'bivarSurf': they should be ORDERED NUMERIC POSITIVE values with FIRST element equal to '0'")}
  numericY = as.numeric(colnames(bivarSurf))
  if(any(is.na(numericY)) | any(numericY != numericY[order(numericY)]) | numericY[1] != 0)
    {stop("Please check the column names of 'bivarSurf': they should be ORDERED NUMERIC POSITIVE values with FIRST element equal to '0'")}
  
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
