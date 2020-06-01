#' @importFrom stats pweibull pexp qnorm plogis plnorm
##########################################################
### Eventually this function has to be merged with an existing function in PResiduals
###   which is called presid.survreg()
### The current version of presidTMP.survreg() in Presiduals does not support log-logistic
###   and log-normal distributions
##########################################################
presidTMP.survreg <- function(object, ...){  
    time <- object$y[,1]
    delta <- object$y[,2]
    
    switch(object$dist,
           weibull = {
               prob <- pweibull(exp(time), shape=1/summary(object)$scale,
                                scale=exp(object$linear.predictors),
                                lower.tail=TRUE, log.p=FALSE)
               prob + delta*(prob - 1)
           },
           
           exponential = {
               ### should time be exp(time)?  I am pretty sure about this.
               # prob <- pexp(exp(time), rate=1/exp(object$linear.predictors),
               prob <- pexp(exp(time), rate=1/exp(object$linear.predictors),
                            lower.tail=TRUE, log.p=FALSE)
               prob + delta*(prob - 1)
           },
           
           gaussian = {
               prob <- pnorm(time, mean=object$linear.predictors,
                             sd=summary(object)$scale, lower.tail=TRUE,
                             log.p=FALSE)
               prob + delta*(prob - 1)
           },
           
           logistic = {
               prob <- plogis(time, location=object$linear.predictors,
                              scale=summary(object)$scale, lower.tail=TRUE,
                              log.p=FALSE)
               prob + delta*(prob - 1)
           },
         
           ######### Svetlana's update
           loglogistic = {
                 gamma = (1/object$scale)
                 monsterTerm = ((exp(time))^gamma)*exp(-object$linear.predictors*gamma)
                 ### presid = 1 - (1+delta)*1/(1+monsterTerm)   ### PSR using survival probability
                 prob = 1 - 1/(1+monsterTerm)
                 prob + delta*(prob - 1)
           },
                    
           lognormal = {
               ### should time be exp(time)?
               # prob <- plnorm(exp(time), meanlog=object$linear.predictors,
               prob <- plnorm(exp(time), meanlog=object$linear.predictors,
                              sdlog=summary(object)$scale, lower.tail=TRUE,
                              log.p=FALSE)
               prob + delta*(prob - 1)
           },
           stop("Unhandled dist", object$dist))
}
