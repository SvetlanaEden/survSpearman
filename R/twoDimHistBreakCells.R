###################################################################
######### Prepares data for function survPMFPlot()
###################################################################
twoDimHistBreakCells <- function(Sdxdy, divPointsX, divPointsY){
  ### Takes a two dimensional function (as a matrix)
  ### this matrix has to have col and row names as points
  ### rows are Y's and cols are X's
  ### Takes division points:
  ### divPointsX are for colnames
  ### divPointsY are for rownames
  ### and sums up the pieces of the matrix based on these division points.
  
  arrangeIndicesForHist = function(X, divPointsX){
    indicesX = matrix(NA, ncol = 2, nrow = length(divPointsX) - 1)
    indicesX[1, 1] = 1
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
  
  resX = arrangeIndicesForHist(as.numeric(rownames(Sdxdy)), divPointsX)
  resY = arrangeIndicesForHist(as.numeric(colnames(Sdxdy)), divPointsY)
  indicesX = resX$indices
  indicesY = resY$indices
  
  newMatr = matrix(NA, nrow = nrow(indicesX), ncol = nrow(indicesY))
  for(i in 1:nrow(indicesX)){
    for(j in 1:nrow(indicesY)){
      if(any(is.na(indicesX[i, ])) | any(is.na(indicesY[j, ]))){
        newMatr[i, j] = 0
      }else{
        newMatr[i, j] = sum(Sdxdy[indicesX[i, 1]:indicesX[i, 2], indicesY[j, 1]:indicesY[j, 2]])
      }
    }
  }
  list(newMatr = newMatr, divPX = resX$divisP, divPY = resY$divisP)
}

