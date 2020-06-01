#' @name plotMatrix
#' @aliases plotMatrix
#' @title Plots a color matrix.
#' @description Plots a matrix of colors at a given point. Each element of the color matrix is plotted as a rectangle with user specified side lengths and border color.
#' @usage plotMatrix(colorMatr, borderCol, coordX, coordY, widthX, widthY)
#' 
#' @param colorMatr A matrix of colors in R format.
#' @param borderCol The border color of plotted rectangles.
#' @param coordX X coordinate of the lower left corner of \code{colorMatr}.
#' @param coordY Y coordinate of the lower left corner of \code{colorMatr}.
#' @param widthX The vector of lengths of the rectangle sides parallel to the X-axis. The length of \code{widthX} can be either 1 (all rectangles have the same side length for X-axis) or the number of columns of \code{colorMatr}
#' @param widthY The vector of lengths of the rectangle sides parallel to the Y-axis. The length of \code{widthY} can be either 1 (all rectangles have the same side length for Y-axis) or the number of rows of \code{colorMatr}
#' @return None
#' @details \code{R}-function \code{plot()} has to be called first. The function plots a color matrix at a given coordinate. Each element of the color matrix is plotted as a rectangle with user specified border color and side lengths on X- and Y-axes. Element \code{colorMatr[nrow(colMatr), 1]} is displayed as the bottom left rectangle. Element \code{colorMatr[1, ncol(colMatr)]} is displayed as the top right rectangle.
#' 
#' @examples
#' ### Uneven rectangles
#' colorMatr = matrix(c("goldenrod1", "mediumpurple3", "palegreen3",
#'    "royalblue1", "orchid", "firebrick1"), nrow = 2, byrow = TRUE)
#' plot(c(1, 4), c(2, 4), type = "n", xlab = "", ylab = "")
#' plotMatrix(colorMatr, borderCol = "white", coordX = 1, coordY = 2,
#'    widthX = c(1/4, 1, 1/6), widthY = c(1, 1/2))
#'    
#' ### Plotting the legend:
#' reshapeColMatr = matrix(t(colorMatr), ncol = 1, nrow = nrow(colorMatr)*ncol(colorMatr),
#'    byrow = TRUE)
#' plotMatrix(reshapeColMatr, borderCol = "white", coordX = 2.8, coordY = 2,
#'    widthX = c(1/4), widthY = c(1/4))
#' text(x = rep(3.03, nrow(reshapeColMatr)),
#'    y = 2.12+c(0,cumsum(rep(1/4, nrow(reshapeColMatr)-1))),
#'    labels = reshapeColMatr[nrow(reshapeColMatr):1], pos = 4, cex = 0.7)
#'    
#' ### Same length and width rectangles on the black background
#' plot(c(1, 4), c(2, 4), type = "n", xlab = "", ylab = "")
#' polygon(x = c(0, 5, 5, 0, 0), y = c(0, 0, 5, 5, 0), col = "black")
#' plotMatrix(colorMatr, borderCol = "black", coordX = 1, coordY = 2, widthX = 1, widthY = 1)
#' 
#' @keywords color matrix
#' @author Svetlana K Eden, \email{svetlanaeden@gmail.com}
#' 
#' @importFrom grDevices dev.size gray
#' @importFrom graphics mtext par plot polygon text
#' @export
###################################################################
######### Plots matrix of colored rectangles.
###################################################################
plotMatrix <- function(colorMatr, borderCol = "white", coordX, coordY, widthX, widthY){
  ####### Plots a matrix as a set of rectangles with colors defined by the colorMatr
  ####### The columns ???
  ### colorMatr - matrix with color
  ### coordX, coordY - where to place the bottom left corder of the matrix
  ### widthX - the width of the plotted rectangles on X axis
  ### widthY - the width of the plotted rectangles on Y axis
  ### borderCol - the color of the rectanges border
  if(FALSE){
  if(class(colorMatr) != "matrix") stop("'colorMatr' is a matrix of character strings representing colors")
  if(lwX != nrow(colorMatr)) {stop("The number of rows in 'colorMatr' should be exactly the length of 'widthX'")}
  if(lwY != ncol(colorMatr)) {stop("The number of columns in 'colorMatr' should be exactly the length of 'widthY'")}
  
  xValues = coordX + c(0, cumsum(widthX))
  yValues = coordY + c(0, cumsum(widthY))
  for(i in 1:nrow(colorMatr)){
    for(j in 1:ncol(colorMatr)){
      polygon(x = c(xValues[i], xValues[i], xValues[i+1], xValues[i+1]), y = c(yValues[c(j, j+1)], yValues[c(j+1, j)]), border = borderCol, col = colorMatr[i,j])
    }
  }
  }
  
  if(class(colorMatr) != "matrix") stop("'colorMatr' is a matrix of character strings representing colors")
  lwX = length(widthX)
  lwY = length(widthY)
  if(lwX == 1){
    widthX = rep(widthX, ncol(colorMatr))
    lwX = length(widthX)
  }
  if(lwY == 1){
    widthY = rep(widthY, nrow(colorMatr))
    lwY = length(widthY)
  }
  if(lwX != ncol(colorMatr)) {stop("The length of 'widthX' should be 1 or the number of columns in 'colorMatr'")}
  if(lwY != nrow(colorMatr)) {stop("The length of 'widthY' should be 1 or the number of rows in 'colorMatr'")}

  xValues = coordX + c(0, cumsum(widthX))
  yValues = coordY + c(0, cumsum(widthY))
  yValues = yValues[length(yValues):1]
  for(i in 1:ncol(colorMatr)){
    for(j in 1:nrow(colorMatr)){
      polygon(x = c(xValues[i], xValues[i], xValues[i+1], xValues[i+1]), y = c(yValues[c(j, j+1)], yValues[c(j+1, j)]), border = borderCol, col = colorMatr[j,i])
    }
  }
  
}
