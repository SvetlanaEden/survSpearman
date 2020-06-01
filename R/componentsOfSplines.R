###########################################
########################################### creates splines compondents according to Dr. Frank Harrell
###########################################
componentsOfSplines = function(x, knots){
### gives sline components
  positive = function(x){
  ### used by spline functions
    x[x <= 0]=0
    x
  }
  k = length(knots)
  res = matrix(NA, nrow=length(x), ncol=k-1)
  res[, 1] = x
  kd = (knots[k] - knots[1])^(2/3)
  for (j in 1:(k-2)){
    res[, j+1] = positive(((x-knots[j])/kd)^3) - positive(((x-knots[k-1])/kd)^3) * (knots[k]-knots[j])/(knots[k]-knots[k-1]) + positive(((x-knots[k])/kd)^3) * (knots[k-1]-knots[j])/(knots[k]-knots[k-1])
  }
  res
}
