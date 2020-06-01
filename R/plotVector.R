#' @name plotVector
#' @aliases plotVector
#' @title Plots a set of rectangles that represent a numeric vector.
#' @description Plots a vector of numeric values as a set of rectangles with user specified height, width, color, and border color. The bases of the rectangles are aligned like in a histogram.
#' 
#' @usage plotVector(vectValues, width, coordX, coordY, rotationRadians, vectColor, bordColor)
#' 
#' @param vectValues The heights of the rectangles. If \code{vectValues[i]}>=0 then the base of rectangle \code{i} is plotted at \code{coordY} and the top is plotted at \code{coordY} + (\code{abs(vectValues[i])}).  If the \code{vectValues[i]}<0 then the base of rectangle \code{i} is plotted at \code{coordY} and the bottom is plotted at \code{coordY} - \code{abs(vectValues[i])}.
#' @param width The width (positive number) of the rectangles. It can be either one number or a vector of numbers with the same length as \code{vectValues}.
#' @param coordX X coordinate of the left base corner of the first rectangle, \code{vectValues[1]}.
#' @param coordY Y coordinate of the left base corner of the first rectangle, \code{vectValues[1]}.
#' @param rotationRadians By how much (in radians) to rotate the aligned rectangles around \code{(coordX, coordY)}. The rotation is performed counter clock-wise.
#' @param vectColor The colors of the rectangles in \code{R} color format. It can be either one value or a vector of values with the same length as \code{vectValues}.
#' @param bordColor The color of the border of the rectangles. It can be either one value or a vector of values with the same length as \code{vectValues}.

#' @return None
#' @details \code{R}-function \code{plot} has to be called first. The rectangles are aligned in the following manner. If \code{rotationRadians} is zero, then the base of all rectangles is at \code{coordY}. If \code{vectValues[i]>=0} then the top of the rectangle is plotted at \code{coordY} + \code{abs(vectValues[i])}. If the \code{vectValues[i]}<0 then the bottom of the rectangle is plotted at \code{coordY} - \code{abs(vectValues[i])}. If \code{rotationRadians} is not zero, then the aligned rectangles are rotated counter clock-wise around point (\code{coordX},\code{coordY}).
#' 
#' @examples
#' ### Histogram-like plot
#' plot(c(-1, 2), c(-1, 1), type = "n")
#' vectValues = c(0.057, 0.336, 0.260, 0.362, 0.209)
#' plotVector(vectValues, width = 0.2, coordX = 0, coordY = 0, rotationRadians = 0,
#'    vectColor = "gray", bordColor = "white")
#'    
#' ### Rotated and flipped sequence of rectangles.
#' width = c(0.10, 0.20, 0.10, 0.80, 0.12)
#' vectColor = c("orange", "green", "orchid", "blue", "goldenrod1")
#' plot(c(-1, 2), c(-1, 2), type = "n")
#' plotVector(vectValues = vectValues, width, coordX = 0, coordY = 0,
#'    rotationRadians = pi/2, vectColor, bordColor = "white")
#' plotVector(vectValues = -vectValues, width, coordX = 0, coordY = 0,
#'    rotationRadians = 0, vectColor, bordColor = "white")
#'    
#' ### Histogram-like plot with positive and negative values.
#' vectValues = c(0.057, -0.336, 0.260, -0.222, 0.209)
#' plot(c(-1, 1), c(-1, 1), type = "n")
#' vectColor = rep("goldenrod1", length(vectValues))
#' vectColor[vectValues<0] = "royalblue1"
#' bordColor = rep("red", length(vectValues))
#' bordColor[vectValues<0] = "darkblue"
#' plotVector(vectValues, width = 0.4, coordX = -1, coordY = 0, rotationRadians = 0,
#'     vectColor, bordColor = bordColor)
#' 
#' @keywords color rotated histogram
#' @author Svetlana K Eden, \email{svetlanaeden@gmail.com}
#' 
#' @importFrom grDevices dev.size gray
#' @importFrom graphics mtext par plot polygon text
#' @export
###################################################################
######### Plots a histogram starting at a given point and with rotation option.
###################################################################
plotVector <- function(vectValues, width, coordX, coordY, rotationRadians = 0, vectColor = NULL, bordColor = NULL){
  if (length(width) == 1){
    width = rep(width, length(vectValues))
  }else{
    if (length(width) < length(vectValues) | length(width) <= 0) stop("The length of 'width' should be either 1 or the same as the length of 'vectValues'")
  }
  
  if(is.null(vectColor)){
    vectColor = gray(rep(0.2, length(vectValues)))
  }
  if (length(vectColor) == 1){
    vectColor = rep(vectColor, length(vectValues))
  }else{
    if (length(vectColor) < length(vectValues) | length(vectColor) <= 0) stop("The length of 'vectColor' should be either 1 or the same as the length of 'vectValues'")
  }
  
  if(is.null(bordColor)){
    bordColor = vectColor
  }
  if (length(bordColor) == 1){
    bordColor = rep(bordColor, length(vectValues))
  }else{
    if (length(bordColor) < length(vectValues) | length(bordColor) <= 0) stop("The length of 'bordColor' should be either 1 or the same as the length of 'vectValues'")
  }
  
  divPX = cbind(c(0, cumsum(width)[1:(length(width) - 1)]), cumsum(width))
  rotaionMatrix = matrix(c(cos(rotationRadians), - sin(rotationRadians), sin(rotationRadians), cos(rotationRadians)), nrow = 2, ncol = 2, , byrow = TRUE)
  
  ### checking the range of the plotted vector:
  xyRange = rotaionMatrix %*% rbind(c(coordX, coordX + sum(width)), c(coordY, coordY + max(vectValues)))
  usr <- par('usr')
  xr <- (usr[2] - usr[1]) / 27 # 27 = (100 + 2*4) / 4
  yr <- (usr[4] - usr[3]) / 27
  xlim <- c(usr[1] + xr, usr[2] - xr)
  ylim <- c(usr[3] + yr, usr[4] - yr)
  
  if((all(xyRange[,1] <= xlim[1])  |  all(xyRange[,1] >= xlim[2]))  &  (all(xyRange[,2] <= ylim[1])  |  all(xyRange[,2] >= ylim[2]))){warning("The plotted area in function 'plotVector':  [", paste(xyRange[,1], collapse = ", "), "]x[",  paste(xyRange[,2], collapse = ", "), "] is out of the following plotted range: [", paste(xlim, collapse = ", "), "]x[",  paste(ylim, collapse = ", "), "], and might not be displayed in the plot.")}
  
  #cat("xlim", xlim, "ylim", ylim, "\n")
  
  for(i in 1:length(vectValues)){
    x = c(divPX[i, 1], divPX[i, 1], divPX[i, 2], divPX[i, 2])
    y = c(0, vectValues[i], vectValues[i], 0)
    newCoord = rotaionMatrix %*% rbind(x, y)
    polygon(newCoord[1,] + coordX, newCoord[2,] + coordY, border = bordColor[i], col = vectColor[i])
  }
}
