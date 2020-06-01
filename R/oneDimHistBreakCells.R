###################################################################
######### Prepares data for function plot1DHistWithRotation()
###################################################################
oneDimHistBreakCells <- function(Sdx, divPointsX){
  ### Takes a one dimensional function (as a vector)
  ### this vector has to have names as points
  ### Takes division points: divPointsX
  ### and sums up the pieces of the vector based on these division points.
  if(FALSE){
    Sdx = c(1, 2, 3, 1, 2, 5, 9); names(Sdx) = as.character(1:length(Sdx))
    divPointsX = c(0, 3.5, 5, 5.1, 5.2, 8, 9)
    oneDimHistBreakCells(Sdx, divPointsX)
  }
  # browser()
  arrangeIndicesForHist = function(X, divPointsX){
    indicesX = matrix(NA, ncol = 2, nrow = length(divPointsX) - 1)
    # indicesX[1, 1] = 1
    i = iterX = 1
    while(iterX <= nrow(indicesX) & i <= length(X)){
      anythingInside = divPointsX[iterX] < X[i] & X[i] <= divPointsX[iterX+1]
      if(anythingInside){
        indicesX[iterX, 1] = i
        while( divPointsX[iterX] < X[i] & X[i] <= divPointsX[iterX+1] & i <= length(X)){
          i = i + 1
        }
        indicesX[iterX, 2] = i - 1
      }
      iterX = iterX + 1
    }
    divPX = cbind(divPointsX[1:(length(divPointsX) - 1)], divPointsX[2:(length(divPointsX))])
    list(divisP = divPX, indices = indicesX)
  }
  
  resX = arrangeIndicesForHist(as.numeric(names(Sdx)), divPointsX)
  indicesX = resX$indices
  
  newVec = rep(NA, nrow(indicesX))
  for(i in 1:nrow(indicesX)){
    if(all(is.na(indicesX[i, ]))){
      newVec[i] = 0
    }else{
      newVec[i] = sum(Sdx[indicesX[i, 1]:indicesX[i, 2]])
    }
  }
  list(newVec = newVec, divPX = resX$divisP)
}
