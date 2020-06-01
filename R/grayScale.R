###################################################################
######### Create gray scale with lighter tone being smaller positive values
#########      and darker tone being higher positive values
######### It is used by function survPMFPlot()
###################################################################
grayScale <- function(x, X0, X1, Y0, Y1, method, lowerThreshold = -Inf, YValueForLowerThreshold = NA){
  ### there are points (X0,Y0) and (X1, Y1) that we aim for
  if(X0 == X1){
    stop("X0 should not be equal to X1.")
  }else{    
    if(!(method %in% c(1, 2, 3))){
      stop("Argument 'method' can equal only to 1, 2, or 3.")
    }
    if(method == 1){
      a = (Y0 - Y1)/(X0 - X1)
      b = Y1 - a*X1
      res = a*x  + b
    }
    if(method == 2){
      a = (Y1 - Y0)/((X1 - X0)^2)
      b = -2*a*X0
      c = Y0 + a*X0^2
      res = a*x^2  + b*x + c
    }
    if(method == 3){
      a = -(Y1 - Y0)/((X1 - X0)^2)
      b = -2*a*X1
      c = Y0 - a*X0^2 + 2*a*X1*X0
      res = a*x^2  + b*x + c
    }
  }
  res[x < lowerThreshold] = YValueForLowerThreshold
  res
}
