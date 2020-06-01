###################################################################
### Qi's code that is necessary for conditinal correlation
### It is part of the package PResiduals, but it is "invisible" to users
### So I duplicated it, but changed the name, of course
### It has to replace with lm.scores() eventually
### I could not do it because I do not maintain PResiduals
###################################################################
lm.scoresForBivarSurvival = function(y, X){
  N = length(y)  
  mod = lm(y~X)
  smod = summary(mod)
  resid = smod$residuals

  d2l.dtheta.dtheta = -crossprod(cbind(1, X))

  dl.dtheta <- resid*cbind(1, X)
  presid = 2*pnorm((y - mod$fitted.values)/smod$sigma) -1
  dresid.dtheta = t(cbind(-1, -X))
  dpresid.dtheta = t(cbind(-2*dnorm((y - mod$fitted.values)/smod$sigma)/smod$sigma,
                           -2*dnorm((y - mod$fitted.values)/smod$sigma)/smod$sigma *
                             X))
  
  f.y<-density(resid)
  fy.ry <- NULL
  presid.k <- NULL
  for (i in 1:length(resid)){
    fy.ry[i] <- f.y$y[which(abs(f.y$x-resid[i])==min(abs(f.y$x-resid[i])))]
    presid.k[i] <- sum(resid<resid[i])/length(resid) - sum(resid>resid[i])/length(resid)
  }
  dpresid.dtheta.k <- t(cbind(-2*fy.ry,
                              -2*fy.ry*X))
  list(mod = mod, 
       dl.dtheta = dl.dtheta,
       d2l.dtheta.dtheta = d2l.dtheta.dtheta,
       resid = resid,
       dresid.dtheta = dresid.dtheta,
       presid = presid,
       presid.k= presid.k,
       dpresid.dtheta = dpresid.dtheta,
       dpresid.dtheta.k = dpresid.dtheta.k)
}

