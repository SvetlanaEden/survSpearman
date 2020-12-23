#' CCASAnetData
#' 
#' Bivariate survival data simulated from Frank's copula family with parameter 4.426265, which is equivalent to 0.6 Spearman's correlation.
#' In the first 100 rows, the follow-up time for all observations is restricted to 2, and some observations may be censored before the follow-up time.
#' In rows 101:200, the follow-up time is not restricted, but thirty percent of observations are independently censored.
#' In rows 201:300, there is no censoring.
#' 
#' @format A data frame with 6691 rows and 5 variables:
#' \describe{
#'   \item{X}{Time to event X.}
#'   \item{deltaX}{Event indicator for event X (1 - event; 0 - censoring event).}
#'   \item{Y}{Time to event Y}
#'   \item{deltaY}{Event indicator for event Y (1 - event; 0 - censoring event).}
#' }
#' @source Simulated data.