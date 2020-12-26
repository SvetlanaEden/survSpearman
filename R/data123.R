#' data123
#' 
#' Bivariate survival data simulated from Frank's copula family with parameter 4.426265, which is equivalent to Spearman's correlation of 0.6.
#' The marginal distributions of time to event and time to censoring are exponential with mean 1.
#' In the first 100 rows, the follow-up time is restricted to time 2, and the observations can also be independently censored before time 2.
#' In rows 101:200, the follow-up time is not restricted, but some observations are independently censored.
#' In rows 201:300, there is no censoring.
#' 
#' @format A data frame with 300 rows and 4 variables:
#' \describe{
#'   \item{X}{Time to event X.}
#'   \item{deltaX}{Event indicator for event X (1 - event; 0 - censoring event).}
#'   \item{Y}{Time to event Y}
#'   \item{deltaY}{Event indicator for event Y (1 - event; 0 - censoring event).}
#' }
#' @source Simulated data.
"data123"