#' @name survSpearman
#' @aliases survSpearman
#' @title Computes Spearman's Correlation for Bivariate Survival Data.
#' @description Computes non-parametric estimates of Spearman's rank correlation for bivariate survival data. Two correlations are returned: a highest rank correlation that can be interpreted as Spearman's correlation after assigning a highest rank to observations beyond a specified region, and a restricted correlation that estimates Spearman's correlation within the specified region.
#' 
#' @usage survSpearman(X = NULL, Y = NULL, deltaX = NULL, deltaY = NULL, data = NULL,
#' tauX = Inf, tauY = Inf, bivarSurf = NULL)
#' 
#' @param X Time to event or censoring for variable \code{X}. It indicates time to event if argument \code{deltaX}=1 and time to censoring if argument \code{deltaX}=0.
#' @param Y Time to event or censoring for variable \code{Y}. It indicates time to event if argument \code{deltaY}=1 and time to censoring if argument \code{deltaY}=0.
#' @param deltaX Event indicator for variable \code{X}. \code{deltaX=1} if the event is observed and \code{0} if it is censored.
#' @param deltaY Event indicator for variable \code{Y}. \code{deltaY=1} if the event is observed and \code{0} if it is censored.
#' @param data Data frame containing variables (arguments) \code{X}, \code{Y}, \code{deltaX}, and \code{deltaY}.
#' @param tauX The \code{X} value that defines the restricted region for \code{X} variable.
#' @param tauY The \code{Y} value that defines the restricted region for \code{Y} variable.
#' @param bivarSurf A matrix containing the marginal and joint survival probabilities. The first column is the marginal survival probability corresponding to variable \code{X}. The first row is the marginal survival probability corresponding to variable \code{Y}. The rest of the matrix contains the joint survival probabilities. The row names of \code{bivarSurf} are ordered \code{X}-values. The column names of \code{bivarSurf} are ordered \code{Y}-values. Element \code{bivarSurf[1,1]} equals 1. Its row and column name is \code{'0'} (see the documentation for the return value \code{DabrowskaEst} in function \code{survDabrowska}).
#' 
#' @return The function returns the following list of values. \code{'Restricted region set by user'} is a character vector of two user-specified restricted region values, \code{tauX} and \code{tauY}. \code{'Effective restricted region'} is character vector of two values that define the effective restricted region, the values that are just above the latest observed event times within the user-specified restricted region. \code{'Correlation'} is a numeric vector of two correlation values: the highest rank Spearman's correlation (\code{'HighestRank'}) and the restricted region Spearman's correlation (\code{'Restricted'}), where the restricted region is defined by the values in \code{'Effective restricted region'}.
#' @details The function computes the highest rank and restricted Spearman's correlations with bivariate survival data. The data can be supplied in three ways: 1) as vectors \code{X}, \code{Y}, \code{deltaX}, and \code{deltaY}; 2) as data frame \code{data} that contains the variables mentioned in 1); and 3) as matrix \code{bivarSurf} containing marginal and joint survival probabilities. If \code{bivarSurf} is not \code{NULL} then 1) and 2) are ignored. If \code{bivarSurf} is \code{NULL} and \code{data} is not then 2) is used. If \code{bivarSurf} and \code{data} are \code{NULL} then 1) is used. The highest rank correlation is the Spearman's correlation that can be interpreted as Spearman's rank correlation computed after assigning the highest rank to the events outside of \code{tauX} and \code{tauY}. The restricted Spearman's correlation is Spearman's correlation computed within the restricted region defined by \code{tauX} and \code{tauY}. Note that given \code{tauX} and \code{tauY} the survival probability is estimated using the values that are just above the latest observed event times within that region, what we call an effective restricted region. This means that if, for example, \code{tauX} is greater than the latest observed event time for \code{X} variable and \code{tauY} is greater than the latest observed event time for \code{Y} variable, then \code{tauX} and \code{tauY} do not affect the correlation values since the effective restricted region remains the same (as defined by the maximum observed event or censoring event times). The method of Dabrowska can result in negative probability mass for some points. This may result in zero or negative probability of failure in the restricted region, in which case the restricted Spearman's correlation cannot be computed and NA value is returned. This only happens when the sample size is small and censoring is heavy.
#' 
#' @examples
#' ### Compute correlation from data
#' X <- c(0.5, 0.6, 0.7, 0.8)
#' Y <- c(0.44, 0.77, 0.88, 0.99)
#' deltaX <- c(1, 0, 1, 1)
#' deltaY <- c(1, 1, 1, 1)
#' survSpearman(X, Y, deltaX, deltaY)
#' survSpearman(X, Y, deltaX, deltaY, tauX = 100, tauY = 100)
#' survSpearman(X, Y, deltaX, deltaY, tauX = 100, tauY = 0.99)
#' survSpearman(X, Y, deltaX, deltaY, tauX = 0.8, tauY = 0.99)
#' 
#' ### Compute correlation from survival surface
#' someSurf <- survDabrowska(X, Y, deltaX, deltaY)$DabrowskaEst
#' survSpearman(tauY = 0.9, bivarSurf = someSurf)
#' 
#' @keywords bivariate survival Spearman correlation
#' @references Dabrowska, D. M. (1988) Kaplan–Meier estimate on the plane. The Annals of Statistics 16, 1475–1489.
#' @references Eden, S., Li, C., Shepherd B. (2019). Non-parametric Estimation of Spearman's Rank Correlation with Bivariate Survival Data.
#' @author Svetlana K Eden, \email{svetlanaeden@gmail.com}
#' 
#' @export
################################################################### Highest Rank and Restrcited Spearman's Rho in the Restricted
################################################################### Region
survSpearman = function(X = NULL, Y = NULL, deltaX = NULL, deltaY = NULL, data = NULL, tauX = Inf, tauY = Inf, bivarSurf = NULL) {
    ### Arguments' check
    # if(method == "Dabrowska"){
    if(is.null(bivarSurf)){
      if(is.null(data)){
        l1 = length(X); l2 = length(Y); l3 = length(deltaX); l4 = length(deltaY)
        if(l1 != l2 | l1 != l3 | l1 != l4) stop("Arguments 'X', 'Y', 'deltaX', and 'deltaY' are numeric vectors of equal length when 'data' and 'bivarSurf' are NULL.\n'")
        data = data.frame(X = X, Y = Y, deltaX = deltaX, deltaY = deltaY)
      }
      if(nrow(data) == 0) stop("When bivarSurf' is not supplied either 'data' or 'X', 'Y', 'deltaX', and 'deltaY' should be provided with non-zero number of records")
      if(any(is.na(data))) stop("Missing values are not allowed in 'X', 'Y', 'deltaX', and 'deltaY'.\n")
      uniqueDeltas = unique(unlist(data[, c("deltaX", "deltaY")]))
      if(!(setequal(uniqueDeltas, c(0, 1)) | all(uniqueDeltas == 1) | all(uniqueDeltas == 0))) stop("Arguments 'deltaX', and 'deltaY' are binary.\n")
      if(all(data$deltaX == 0))stop("There are no events for X, all values of 'deltaX' = 0.\n")
      if(all(data$deltaY == 0))stop("There are no events for Y, all values of 'deltaY' = 0.\n")
      bivarSurf = survDabrowska(data$X, data$Y, data$deltaX, data$deltaY)$DabrowskaEst
    }else{
      #if(method != "Custom"){stop("Argument 'method' can be either 'Dabrowska' or 'Custom'.\n")}
      #if(is.null(bivarSurf)){stop("When 'method' is 'Custom', the user should supply a bivariate survival surface, argument 'bivarSurf', in the format described in the help file.\n")}
      if(any(is.na(bivarSurf))) stop("NA values are not allowed in 'bivarSurf'.\n")
      colN = as.numeric(colnames(bivarSurf))
      rowN = as.numeric(rownames(bivarSurf))
      if(length(colN)==0 | any(colN[order(colN)] != colN)) stop("Column names of 'bivarSurf' should be ordered X values.")
      if(length(rowN)==0 | any(rowN[order(rowN)] != rowN)) stop("Row names of 'bivarSurf' should be ordered Y values.")
      if(colN[1]!=0) stop("The first element of column names of 'bivarSurf' should be equal to 0.")
      if(rowN[1]!=0) stop("The first element of row names of 'bivarSurf' should be equal to 0.")
      if(bivarSurf[1,1] != 1) stop("Element 'bivarSurf[1,1]' should be equal to 1.")
    }
  
    ### Find the last observed events and censoring events
    lastX = rownames(bivarSurf)[nrow(bivarSurf)]
    lastY = colnames(bivarSurf)[ncol(bivarSurf)]
    lastEventOnX = names(bivarSurf[, 1] == min(bivarSurf[, 1]))[bivarSurf[, 
        1] == min(bivarSurf[, 1])][1]
    lastEventOnY = names(bivarSurf[1, ] == min(bivarSurf[1, ]))[bivarSurf[1, 
        ] == min(bivarSurf[1, ])][1]
    
    ### Find the last observed events and censoring events in the
    ### restricted region
    if(sum(as.numeric(rownames(bivarSurf)) < tauX) < 2){stop("'tauX' is too small: it excludes all points X.")}
    if(sum(as.numeric(colnames(bivarSurf)) < tauY) < 2){stop("'tauY' is too small: it excludes all points Y.")}
    restrBivarSurf = bivarSurf
    restrBivarSurf = restrBivarSurf[as.numeric(rownames(restrBivarSurf)) < 
        tauX, ]
    restrBivarSurf = restrBivarSurf[, as.numeric(colnames(restrBivarSurf)) < 
        tauY]
    restrLastX = rownames(restrBivarSurf)[nrow(restrBivarSurf)]
    restrLastY = colnames(restrBivarSurf)[ncol(restrBivarSurf)]
    restrLastEventOnX = names(restrBivarSurf[, 1] == min(restrBivarSurf[, 
        1]))[restrBivarSurf[, 1] == min(restrBivarSurf[, 1])][1]
    restrLastEventOnY = names(restrBivarSurf[1, ] == min(restrBivarSurf[1, 
        ]))[restrBivarSurf[1, ] == min(restrBivarSurf[1, ])][1]
    
    actualTauX = paste(restrLastEventOnX, "+")
    actualTauY = paste(restrLastEventOnY, "+")
    
    ### compute restricted correlation
    res1 = HighestRankAndRestrictedSpearman(bivarSurf, tauX = tauX, 
        tauY = tauY)
    bivarSurfForHR = bivarSurf
    bivarSurfForHR[as.numeric(rownames(bivarSurfForHR)) >= tauX, ] = 0
    bivarSurfForHR[, as.numeric(colnames(bivarSurfForHR)) >= tauY] = 0
    
    ### compute highest rank correlation in the restricted region
    res2 = HighestRankAndRestrictedSpearman(bivarSurfForHR, tauX = Inf, 
        tauY = Inf)
    
    resCor = c(HighestRank = res2["HighestRank", ], Restricted = res1["Restricted", 
        ])
    resRegUser = c(tauX = as.character(tauX), tauY = as.character(tauY))
    resReg = c(tauX = actualTauX, tauY = actualTauY)
    #list(RestrictedRegionSetByUser = resRegUser, RestrictedRegionEffective = resReg, 
    #    Correlation = resCor)
    list("Restricted region set by user" = resRegUser, "Effective restricted region" = resReg, 
         Correlation = resCor)
}
