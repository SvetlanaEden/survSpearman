#' @name unadjusted.CorPSRs
#' @aliases unadjusted.CorPSRs
#' @title Unadjusted correlation of probability scale residuals (PSRs)
#' @description Unadjusted correlation of probability scale residuals (PSRs).
#' 
#' @usage unadjusted.CorPSRs(X, Y, deltaX, deltaY)
#' 
#' @param X Time to event or censoring for variable \code{X}. It indicates time to event if argument \code{deltaX}=1 and time to censoring if argument \code{deltaX}=0.
#' @param Y Time to event or censoring for variable \code{Y}. It indicates time to event if argument \code{deltaY}=1 and time to censoring if argument \code{deltaY}=0.
#' @param deltaX Event indicator for variable \code{X}. \code{deltaX=1} if the event is observed and \code{0} if it is censored.
#' @param deltaY Event indicator for variable \code{Y}. \code{deltaY=1} if the event is observed and \code{0} if it is censored.
#'
#' @return The point estimate, standard error, P-value, and confidence interval based on the \code{0.95\%} significance level.
#' @details Computes the unadjusted correlation of probability scale residuals. The variance is computed using the Kaplan-Meier estimating equations proposed by Stute [1995].
#' 
#' @examples
#' X = c(0.5, 0.61, 0.6, 0.8, 0.78, 0.7, 0.9)
#' Y = c(0.44, 0.15, 0.77, 0.88, 0.22, 0.99, .33)
#' deltaX = c(1, 0, 1, 1, 0, 1, 0)
#' deltaY = c(1, 0, 1, 0, 1, 1, 1)
#' unadjusted.CorPSRs(X, Y, deltaX, deltaY)
#'
#' @references Winfried Stute. The central limit theorem under random censorship. The Annals of Statistics, pages 422–439, 1995.
#' @references Svetlana K Eden, Chun Li, and Bryan E Shepherd. Spearman’s rank correlation adjusting for covariates in bivariate survival data. In preparation, 2020.
#'
#' @author Svetlana K Eden, \email{svetlanaeden@gmail.com}
#'
#' @keywords probability scale residuals
#' @importFrom stats cor pnorm qnorm
#' @export
##############################################################################################
############################################### compute unadjusted PSRs correlation
############################################### using Stute's estim. equations [1995]
##############################################################################################
unadjusted.CorPSRs = function(X, Y, deltaX, deltaY){
  l1 = length(X)
  l2 = length(Y)
  l3 = length(deltaX)
  l4 = length(deltaY)
  if (l1 != l2 | l1 != l3 | l1 != l4) 
      stop("Arguments 'X', 'Y', 'deltaX', and 'deltaY' are numeric vectors of equal length.\n'")

  resX = computePSRs(time = X, delta = deltaX)
  resY = computePSRs(time = Y, delta = deltaY)

  xz.d2l.dtheta.dtheta = diag(length(X), length(unique(X[deltaX == 1])))
  yz.d2l.dtheta.dtheta = diag(length(Y), length(unique(Y[deltaY == 1])))

  bundle1 = make_estimating_equations_stute(myKaplanMeier = resX$KM)
  bundle2 = make_estimating_equations_stute(myKaplanMeier = resY$KM)

  psrWithStute = corTS.survreg.modified(xresid = resX$PSRs, yresid = resY$PSRs,
    xz.dl.dtheta = t(bundle1$dl.dtheta), yz.dl.dtheta = t(bundle2$dl.dtheta),
    xz.d2l.dtheta.dtheta = xz.d2l.dtheta.dtheta, yz.d2l.dtheta.dtheta = yz.d2l.dtheta.dtheta,
    dxresid.dthetax = bundle1$dpresid.dtheta, dyresid.dthetay = bundle2$dpresid.dtheta,
    fisher=FALSE, inverseA = FALSE)
  
  res = c(est = psrWithStute$TS, strerr = sqrt(psrWithStute$varTS/length(X)), P = psrWithStute$pvalTS, lower.CI = psrWithStute$CI_TS[1], upper.CI = psrWithStute$CI_TS[2])
  res
}

