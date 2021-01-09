#' @name visualBivarTimeToEvent
#' @aliases visualBivarTimeToEvent
#' @title Plots bivariate right-censored data.
#' @description Plots bivariate right-censored data distinguishing between uncensored, singly censored, and doubly censored observations.
#' 
#' @usage visualBivarTimeToEvent(X, Y, deltaX, deltaY, labelX, labelY,
#'    xlim = NULL, ylim = NULL, dotSize = 0.7,
#'    segLength=abs(mean(diff(c(X, Y)))), scaleLegendGap = 1,
#'    legendCex = 1, labCex = 1, axisCex = 1)
#' 
#' @param X Time to event or censoring for variable \code{X}. It indicates time to event if argument \code{deltaX}=1 and time to censoring if argument \code{deltaX}=0.
#' @param Y Time to event or censoring for variable \code{Y}. It indicates time to event if argument \code{deltaY}=1 and time to censoring if argument \code{deltaY}=0.
#' @param deltaX Event indicator for variable \code{X}. \code{deltaX=1} if the event is observed and \code{0} if it is censored.
#' @param deltaY Event indicator for variable \code{Y}. \code{deltaY=1} if the event is observed and \code{0} if it is censored.
#' @param labelX Label for the X event.
#' @param labelY Label for the Y event.
#' @param xlim The range for the X-axis (the same as parameter \code{xlim} in function \code{plot()}).
#' @param ylim The range on the X-axis (the same as parameter \code{ylim} in function \code{plot()}).
#' @param dotSize The size of the points (the same as parameter \code{cex} in function \code{points()}).
#' @param segLength The length of the segment representing censored observations.
#' @param scaleLegendGap Increases (if > 1) or decreases (if < 1) the distance between the labels in the legend; 
#' @param legendCex The size of the legend font (the same as parameter \code{cex} in function \code{text()}).
#' @param labCex The size of \code{xlab} and \code{ylab} (the same as parameter \code{cex.lab} in function \code{plot()}).
#' @param axisCex The size of axis labels (the same as parameter \code{cex.axis} in function \code{plot()}).
#' @return None
#' @details Plots bivariate right-censored data distinguishing between uncensored, singly censored,
#' and doubly censored observations. The singly and doubly censored observations are plotted
#'  as diamonds and squares, respectively, with a short segment on the right,
#'  which length is the same for X and Y.
#'  The legend is always plotted in the left upper corner.
#' 
#' @examples
#' X = c(0.5, 0.61, 0.6, 0.8, 0.78, 0.7, 0.9)
#' Y = c(0.44, 0.15, 0.77, 0.88, 0.22, 0.99, .33)
#' deltaX = c(1, 0, 1, 1, 0, 1, 0)
#' deltaY = c(1, 0, 1, 0, 1, 1, 1)
#' 
#' visualBivarTimeToEvent(X, Y, deltaX, deltaY, xlim = c(0, 1), ylim = c(0, 1),
#'                        labelX = "X", labelY = "Y", segLength = 0.05)
#' 
#' @keywords bivariate data plot
#' @importFrom graphics plot segments points text
#' @author Svetlana K Eden, \email{svetlanaeden@gmail.com}
#' @export
###################################################################
### Visualize bivariate right-cenosored data
###################################################################
visualBivarTimeToEvent = function(X, Y, deltaX, deltaY, labelX, labelY, xlim = NULL, ylim = NULL, dotSize = 0.7, segLength=abs(mean(diff(c(X, Y)))), scaleLegendGap = 1, legendCex = 1, labCex = 1, axisCex = 1){
  cex = dotSize
  censoredX = deltaX==0 & deltaY==1 
  censoredY = deltaX==1 & deltaY==0
  bothCensored = deltaX==0 & deltaY==0
  notCensored = deltaX==1 & deltaY==1
  dbCenCol = gray(0.8)
  siCenCol = gray(0.5)
  noCenCol = "black"
  
  oldMarPar = par("mar")
  par(mar=c(4, 5, 2, 1))
  on.exit(par(mar = oldMarPar))
  
  ##################################  Not censored
  if(is.null(xlim) & is.null(ylim)){
    xlim = ylim = range(c(X, Y))
  }
  if(is.null(xlim)){
    xlim = range(X)
  }
  if(is.null(ylim)){
    ylim = range(Y)
  }
  
  plot(range(X), range(Y), type="n", xlab=labelX, ylab=labelY, xlim = xlim, ylim = ylim,
       cex.lab = labCex, cex.axis = axisCex)
  ##################################  X censored
  segments(X[censoredX], Y[censoredX], X[censoredX]+segLength, Y[censoredX], col=siCenCol)
  points(X[censoredX], Y[censoredX], pch=23, col = siCenCol, bg = siCenCol, cex = cex)
  ##################################  Y censored
  segments(X[censoredY], Y[censoredY], X[censoredY], Y[censoredY]+segLength, col=siCenCol)
  points(X[censoredY], Y[censoredY], pch=23, col = siCenCol, bg = siCenCol, cex = cex)
  set.seed(1)
  newX = X[bothCensored]
  newY = Y[bothCensored]
  ##################################  Both censored
  points(newX, newY, pch=22, col=dbCenCol, bg = dbCenCol, cex = cex)
  segments(newX, newY, newX+segLength/sqrt(2), newY+segLength/sqrt(2), col=dbCenCol)
  ##################################  Not censored
  points((X[notCensored]), Y[notCensored], pch=21, col=noCenCol, bg = noCenCol, cex = cex)

  ### legend
  legendY = seq(ylim[2], by =-(diff(ylim)/20)*scaleLegendGap, length.out = 4)
  legendX = rep(xlim[1], 4)
  points(legendX, legendY, pch=c(21, 23, 23, 22), col=c(noCenCol, siCenCol, siCenCol, dbCenCol), bg = c("black", siCenCol, siCenCol, dbCenCol), cex = rep(cex, 4))
  text(legendX + segLength * 1.2, legendY, c("Uncensored", paste0(labelX, " is censored"), paste0(labelY, " is censored"), "Doubly Censored"), cex = legendCex, pos = 4)
  segments(legendX, legendY, legendX + c(0, segLength, 0, segLength/sqrt(2)), legendY + c(0, 0, segLength, segLength/sqrt(2)), lty = 1, col = c(noCenCol, siCenCol, siCenCol, dbCenCol))  
}
