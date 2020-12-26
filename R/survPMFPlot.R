#' @name survPMFPlot
#' @aliases survPMFPlot
#' @title Plots marginal and joint probability mass functions.
#' @description Plots marginal and joint probability mass functions for bivariate survival data.
#' 
#' @usage survPMFPlot(bivarSurf, gridWidthX, gridWidthY,
#' scaleGapX, scaleGapY, scaleHistX, scaleHistY,
#' XAxisLabel, YAxisLabel, lineXAxisLabel, lineYAxisLabel,
#' timeLabelScale, axisLabelScale,
#' labelSkipX, labelSkipY, roundX, roundY,
#' plotLabel, linePlotLabel)
#' 
#' @param bivarSurf A matrix containing marginal and joint survival probabilities.  The first column is the marginal survival probability corresponding to variable \code{X}. The first row is the marginal survival probability corresponding to variable \code{Y}. The rest of the matrix contains the joint survival probabilities. The row names of \code{bivarSurf} are ordered \code{X} values. The column names of \code{bivarSurf} are ordered \code{Y} values. Element \code{bivarSurf[1,1]} equals 1. Its row and column name is \code{'0'} (see the documentation for the return value \code{DabrowskaEst} in function \code{survDabrowska})
#' @param gridWidthX Grid size on \code{X}-axis.
#' @param gridWidthY Grid size on \code{Y}-axis.
#' @param scaleGapX The proportion by which the gap between the marginal and joint histograms along the X-axis is increased (if >1) or decreased (if >0 and <1).
#' @param scaleGapY The proportion by which the gap between the marginal and joint histograms along the Y-axis is increased (if >1) or decreased (if >0 and <1).
#' @param scaleHistX The proportion by which to increase (if >1) or decrease (if >0 and <1) the height of the marginal histogram along the X-axis.
#' @param scaleHistY The proportion by which to increase (if >1) or decrease (if >0 and <1) the height of the marginal histogram along the Y-axis.
#' @param XAxisLabel \code{X}-axis label.
#' @param YAxisLabel \code{Y}-axis label.
#' @param lineXAxisLabel Line where to place the \code{X}-axis label (used by function \code{mtext}).
#' @param lineYAxisLabel Line where to place the \code{Y}-axis label (used by function \code{mtext}).
#' @param timeLabelScale The proportion by which the axis labels is increased (if >1) or decreased (if >0 and <1).
#' @param axisLabelScale The proportion by which the time grid labels are increased (if >1) or decreased (if >0 and <1).
#' @param labelSkipX If the time grid looks too busy, not all \code{X}-label time grid labels have to be printed. For example, if \code{labelSkipX = 1} then every other label is printed starting from the first (time 0). If \code{labelSkipX = 2} then the 1st, 4rd, 7th, 10th, ... and so on, labels are printed. The default is \code{0}, which means that all time grid labels are printed.
#' @param labelSkipY If the time grid looks too busy, not all \code{Y}-label time grid labels have to be printed. For example, if \code{labelSkipY = 1} then every other label is printed starting from the first (time 0). If \code{labelSkipX = 2} then the 1st, 4rd, 7th, 10th, ... and so on, labels are printed. The default is \code{0}, which means that all time grid labels are printed.
#' @param roundX Number of decimal places in time grid labels for \code{X}-label.
#' @param roundY Number of decimal places in time grid labels for \code{Y}-label.
#' @param plotLabel Plot label.
#' @param linePlotLabel Line where to place the plot label (used by function \code{mtext}).
#' @return None
#' @details Plots marginal and joint probability mass functions (PMFs) from marginal and joint survival probabilities. The probability mass gets aggregated into cells according to the user-specified arguments \code{gridWidthX} and \code{gridWidthY}. After this aggregation, the negative values (if any) are set to zero. Marginal probability mass functions are displayed as histograms. Joint probability mass function is displayed as a matrix with darker cells indicating larger probability mass aggregated in this cell. Zero mass is denoted with a very faint gray shade. Because the shading is relative, with greater range of probability mass more gray shades are observed. A future version of the function will allow users to choose their own shading/color function.
#' 
#' @examples
#' X = c(0.5, 0.61, 0.6, 0.8, 0.78, 0.7, 0.9)
#' Y = c(0.44, 0.15, 0.77, 0.88, 0.22, 0.99, .33)
#' deltaX = c(1, 0, 1, 1, 0, 1, 0)
#' deltaY = c(1, 0, 1, 0, 1, 1, 1)
#' 
#' dabrSurf = survDabrowska(X, Y, deltaX, deltaY)$DabrowskaEst
#' 
#' grid = 0.1
#' survPMFPlot(bivarSurf = dabrSurf, gridWidthX = grid, gridWidthY = grid,
#'             scaleGapX = 1, scaleGapY = 1,
#'             XAxisLabel = "X", YAxisLabel = "Y", timeLabelScale = 1,
#'             axisLabelScale = 1,
#'             labelSkipX = 0, labelSkipY = 0, roundX = 2, roundY = 2,
#'             plotLabel = "Bivariate PMF")
#' 
#' @keywords bivariate PMF plot
#' @references Eden, S., Li, C., Shepherd B. (2021). Non-parametric Estimation of Spearman's Rank Correlation with Bivariate Survival Data. Biometrics (under revision).
#' @author Svetlana K Eden, \email{svetlanaeden@gmail.com}
#' 
#' @importFrom grDevices dev.size gray
#' @importFrom graphics mtext par plot polygon text
#' @export
###################################################################
######### Plots bivariate probability mass function and historgrams of each variable
###################################################################
survPMFPlot <- function(bivarSurf, gridWidthX, gridWidthY, scaleGapX = 1, scaleGapY = 1, scaleHistX = 1, scaleHistY = 1, XAxisLabel = "X", YAxisLabel = "Y", lineXAxisLabel = -0.5, lineYAxisLabel = -0.5, timeLabelScale = 1, axisLabelScale = 1, labelSkipX = 0, labelSkipY = 0, roundX = 2, roundY = 2, plotLabel = "Bivariate PMF", linePlotLabel = 0.5){
  
  whitestCol = 0.99; grayMaxCol = 0.00
  #maxX = max(sdata$X)
  #maxY = max(sdata$Y)
  maxX = as.numeric(rownames(bivarSurf))
  maxY = as.numeric(colnames(bivarSurf))
  totalMax = max(maxX, maxY)
  
  KMX = unlist(bivarSurf[2:nrow(bivarSurf), 1])
  KMY = unlist(bivarSurf[1, 2:ncol(bivarSurf)])
  tmpSurf = survPMF(bivarSurf)$Sdxdy
  Sdxdy = tmpSurf[2:nrow(tmpSurf), 2:ncol(tmpSurf)]
  Sdx = tmpSurf[2:nrow(tmpSurf), 1]
  Sdy = tmpSurf[1, 2:ncol(tmpSurf)]
  Sdx[Sdx < 0] = 0 
  Sdy[Sdy < 0] = 0 
  #################################
  ### get rid of negative mass:
  
  divPointsX = seq(0, max(gridWidthX*ceiling(maxX/gridWidthX), maxX), gridWidthX)
  divPointsX = divPointsX[order(divPointsX)]
  divPointsY = seq(0, max(gridWidthY*ceiling(maxY/gridWidthY), maxY), gridWidthY)
  divPointsY = divPointsY[order(divPointsY)]
  timeRangeX = divPointsX[seq(1, length(divPointsX), (labelSkipX+1))]
  timeRangeY = divPointsY[seq(1, length(divPointsY), (labelSkipY+1))]
  timeRangeX = timeRangeX[2:length(timeRangeX)]
  timeRangeY = timeRangeY[2:length(timeRangeY)]
  timeRangeLabelX = fillUpStr(fillUpDecimals(round(timeRangeX, roundX)))
  timeRangeLabelY = fillUpStr(fillUpDecimals(round(timeRangeY, roundY)))
  
  oneDimKMX = oneDimHistBreakCells(Sdx = Sdx, divPointsX)
  oneDimKMY = oneDimHistBreakCells(Sdx = Sdy, divPointsY)
  oneDimKMY$newVec = oneDimKMY$newVec/sum(oneDimKMY$newVec)
  oneDimKMX$newVec = oneDimKMX$newVec/sum(oneDimKMX$newVec)
  twoDStuff = twoDimHistBreakCells(Sdxdy, divPointsX, divPointsY)
  newMatr = twoDStuff$newMatr; divPX = twoDStuff$divPX; divPY = twoDStuff$divPY
  newMatr[newMatr < 0] = 0
  newMatr = newMatr/sum(newMatr)
  colMatr = grayScale(newMatr, X0 = min(newMatr), X1 = max(newMatr), Y0 = 0.9, Y1 = -20, method = 3,
                      lowerThreshold = 10^(-13), YValueForLowerThreshold = whitestCol)
  colMatr[colMatr < 0] = 0
  
  totalHistMax = max(oneDimKMX$newVec, oneDimKMY$newVec)
  ### note that the scales and gaps are switched because ....
  # gapX = 3 * gridWidthX * scaleGapY
  # gapY = 3 * gridWidthY * scaleGapX
  gapX = 3*min(gridWidthX, gridWidthY) * scaleGapY
  gapY = 3*min(gridWidthX, gridWidthY) * scaleGapX
  startY = 0 - totalHistMax - gapX
  startX = 0 - totalHistMax - gapY
  # startY = 0 - maxY - gapX
  # startX = 0 - maxX - gapY
  
  par(mar = rep(4, 4))
  xlim = c(startX, max(divPointsX))
  ylim = c(startY, max(divPointsY))
  plot(0, 0, xlim = xlim, ylim = ylim, type = "n", main = "", xlab = "", ylab = "", axes = FALSE)
  # plot(0, 0, xlim = , ylim = c(min(startX, startY), totalMax), type = "n", main = "", xlab = "", ylab = "", axes = FALSE)
  widowSize = dev.size("in")
  cexNumbers = (0.7/7)*min(widowSize)*timeLabelScale
  axisLabelCex = (0.9/7)*min(widowSize)*axisLabelScale
  
  # plot1DHistWithRotation(startX = -gapX, startY = 0, rotationRadians = pi/2, newVec = oneDimKMY$newVec, divPX = oneDimKMY$divPX, whitestCol = whitestCol, grayMaxCol = grayMaxCol)
  # plot1DHistWithRotation (startX = 0, startY = -gapY, rotationRadians = 0, oneDimKMX$newVec * (-1), oneDimKMX$divPX, xlab = "X", whitestCol = whitestCol, grayMaxCol = grayMaxCol)
  ### plot(c(-1, 2), c(-1, 2), type = "n")

  # plotVector(vectValues = oneDimKMY$newVec, width = oneDimKMY$divPX[,2] - oneDimKMY$divPX[,1],
  #                       coordX = -gapX, coordY = 0, rotationRadians = pi/2, bordColor = "white")
  # plotVector(vectValues = oneDimKMX$newVec * (-1), width = oneDimKMX$divPX[,2] - oneDimKMX$divPX[,1],
  #                       coordX = 0, coordY = -gapY, rotationRadians = 0, bordColor = "white")
  
  plotVector(vectValues = scaleHistX*oneDimKMY$newVec, width = oneDimKMY$divPX[,2] - oneDimKMY$divPX[,1],
             coordX = -gapX, coordY = 0, rotationRadians = pi/2, bordColor = "white")
  plotVector(vectValues = scaleHistY*oneDimKMX$newVec * (-1), width = oneDimKMX$divPX[,2] - oneDimKMX$divPX[,1],
             coordX = 0, coordY = -gapY, rotationRadians = 0, bordColor = "white")

  # plot2DHist(startX = 0, startY = 0, newMatr, divPX, divPY, whitestCol = whitestCol, grayMaxCol = grayMaxCol, colMatr = colMatr)
  # plotMatrix(colorMatr = colorMatr, borderCol = "white", coordX = 0, coordY = 0, widthX = divPX[,2] - divPX[,1], widthY = divPY[,2] - divPY[,1])

  colorMatr = matrix(gray(colMatr), nrow = nrow(divPX), ncol = nrow(divPY))
  colorMatr = t(colorMatr)
  colorMatr = colorMatr[nrow(colorMatr):1, ]
  plotMatrix(colorMatr = colorMatr, borderCol = "white", coordX = 0, coordY = 0, widthY = divPY[,2] - divPY[,1], widthX = divPX[,2] - divPX[,1])
  
  # text(x = rep(0, length(timeRangeX)) - gapX/2, y = timeRangeX, label = timeRangeLabelX, srt = 0, cex = cexNumbers)
  # text(y = rep(0, length(timeRangeY)) - gapY/2, x = timeRangeY, label = timeRangeLabelY, srt = 90, cex = cexNumbers)
  text(x = timeRangeX, y = rep(0, length(timeRangeX)) - gapY/2, label = timeRangeLabelX, srt = 90, cex = cexNumbers)
  text(x = rep(0, length(timeRangeY)) - gapX/2, y = timeRangeY, label = timeRangeLabelY, srt = 0, cex = cexNumbers)
  text(x = - gapX/2, y = - gapY/2, label = 0, srt = 135, cex = cexNumbers, font = 1)
  
  mtext(text = XAxisLabel, side = 1, at = mean(timeRangeX[c(1,length(timeRangeX))]), line = lineXAxisLabel, cex = axisLabelCex) ### -4
  mtext(text = YAxisLabel, side = 2, at = mean(timeRangeY[c(1,length(timeRangeY))]), line = lineYAxisLabel, cex = axisLabelCex)
  mtext(text = plotLabel, side = 3,  at = mean(timeRangeX[c(1,length(timeRangeX))]), line = linePlotLabel, cex = axisLabelCex, font = 2)
}
