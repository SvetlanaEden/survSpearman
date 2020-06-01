###################################################################
######### compute PSRs for censored data
###################################################################
computePSRs = function(time, delta){
  dataWithKM = myOwnKM(time = time, delta = delta)
  PSRs = dataWithKM$CDF - dataWithKM$delta*(1 - dataWithKM$CDF_M)
  list(PSRs = PSRs, KM = dataWithKM)
}

