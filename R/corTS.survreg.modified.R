###################################################################
######### Partial correlation usingM-estimation the modification of 
######### Qi's code that have correct sign for matrix A and that
######### that is built in a way that one does not have to invert
######### the part of matrix A that is obtained from the models
######### because it may sometime not invert and the models 
######### provide its inverse anyway ...
###################################################################
corTS.survreg.modified = function(xresid, yresid,
                 xz.dl.dtheta, yz.dl.dtheta,
                 xz.d2l.dtheta.dtheta, yz.d2l.dtheta.dtheta,
                 dxresid.dthetax, dyresid.dthetay, inverseA = TRUE,
                 fisher=FALSE, confLevel = 0.95){
  
  TS = cor(xresid, yresid)
  
  xresid2 = xresid^2
  yresid2 = yresid^2
  xbyyresid = xresid * yresid
  mean.xresid = mean(xresid)
  mean.yresid = mean(yresid)
  mean.xbyyresid = mean(xbyyresid)
  
  bigphi = cbind(xz.dl.dtheta,
                 yz.dl.dtheta,
                 mean.xresid - xresid,
                 mean.yresid - yresid,
                 mean.xbyyresid - xbyyresid,
                 mean(xresid2)-xresid2,
                 mean(yresid2)-yresid2,
                 0)
  
  npar.xz = dim(xz.dl.dtheta)[2]
  npar.yz = dim(yz.dl.dtheta)[2]
  Ntheta = npar.xz + npar.yz + 6
  N = dim(xz.dl.dtheta)[1]
  
  A = matrix(0,Ntheta,Ntheta)
  A[1:npar.xz, 1:npar.xz] = xz.d2l.dtheta.dtheta
  A[npar.xz+(1:npar.yz), npar.xz+(1:npar.yz)] = yz.d2l.dtheta.dtheta

  A[Ntheta-6+(1:6), Ntheta-6+(1:6)] = diag(N, 6)
  
  bigpartial = rbind(c(dxresid.dthetax %*% rep(1, N), rep(0, npar.yz)),
                     c(rep(0, npar.xz), dyresid.dthetay %*% rep(1, N)),
                     c(dxresid.dthetax %*% yresid, dyresid.dthetay %*% xresid),
                     c(dxresid.dthetax %*% (2*xresid), rep(0, npar.yz)),
                     c(rep(0, npar.xz), dyresid.dthetay %*% (2*yresid)))
  
  A[Ntheta-6+(1:5), 1:(npar.xz+npar.yz)] = bigpartial
  
  ## TS also equals numTS / sqrt(varprod) = numTS * revsvp
  numTS = mean.xbyyresid - mean.xresid * mean.yresid
  var.xresid = mean(xresid2) - mean.xresid^2
  var.yresid = mean(yresid2) - mean.yresid^2
  varprod = var.xresid * var.yresid
  revsvp = 1/sqrt(varprod)
  dTS.dvarprod = numTS * (-0.5) * revsvp^3
  
  smallpartial = N *
    c(-mean.yresid * revsvp + dTS.dvarprod * (-2*mean.xresid*var.yresid),
      -mean.xresid * revsvp + dTS.dvarprod * (-2*mean.yresid*var.xresid),
      revsvp,
      dTS.dvarprod * var.yresid,
      dTS.dvarprod * var.xresid)

  A[Ntheta, Ntheta-6+(1:5)] = -smallpartial
  
  if (inverseA){ ### this option means that the part of matrix A (A11) is supplied in inverse form
    A11_inv = A[1:(npar.xz+npar.yz), 1:(npar.xz+npar.yz)]
    A21 = A[Ntheta-6+(1:6), 1:(npar.xz+npar.yz)]
    B_inv = solve(A[Ntheta-6+(1:6), Ntheta-6+(1:6)])
    A_inv = A
    A_inv[Ntheta-6+(1:6), 1:(npar.xz+npar.yz)] = - B_inv %*% A21 %*% A11_inv
    A_inv[Ntheta-6+(1:6), Ntheta-6+(1:6)] = B_inv
    var.theta = A_inv %*% t(bigphi) %*% bigphi %*% t(A_inv)
  }else{
    SS = solve(A, t(bigphi))
    var.theta = tcrossprod (SS, SS)
  }
  varTS = var.theta[Ntheta, Ntheta]  ### this is actually a squared standard error
  pvalTS = 2 * pnorm( -abs(TS)/sqrt(varTS))
  addToValue = abs(qnorm((1-confLevel)/2))*sqrt(varTS)
  CI_TS = TS + c(-1, 1)*addToValue
  

  ####Fisher's transformation
  TS_f <- log( (1+TS)/(1-TS) )
  varTS_f <- varTS*(2/(1-TS^2))^2
  pvalTS_f <- 2 * pnorm( -abs(TS_f)/sqrt(varTS_f))
  addToValue_f = abs(qnorm((1-confLevel)/2))*sqrt(varTS_f)
  CI_TS_f = TS_f + c(-1, 1)*addToValue_f
    
  list(TS=TS, varTS=varTS, pvalTS=pvalTS, CI_TS = CI_TS, pvalTS=pvalTS_f, CI_TS_f = CI_TS_f, var.thetaPSR = diag(var.theta)[dim(var.theta)[1] + c(-5 : 0)])
}


