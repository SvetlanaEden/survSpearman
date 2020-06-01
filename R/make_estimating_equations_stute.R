############################################### estimating equations proposed 
############################################### by Stute's [1995]
############################################### and explained in Shepherd's [2007]
make_estimating_equations_stute = function(myKaplanMeier){
  orderedTime = myKaplanMeier$time[order(myKaplanMeier$time)]
  orderedDelta =  myKaplanMeier$delta[order(myKaplanMeier$time)]
  orderedUniqueKM = unique(myKaplanMeier$KM[order(myKaplanMeier$time)][orderedDelta == 1])
  
  uniqueTime = unique(myKaplanMeier$time)
  orderedUniqueTime = uniqueTime[order(uniqueTime)]
  orderedUniqueFailureTime = unique(orderedTime[orderedDelta == 1])

  N = length(orderedTime)
  uniqueN = length(orderedUniqueTime)
  uniqueFailureN = length(orderedUniqueFailureTime)
  
  orderMatrix = matrix(c(1:N, (1:N)[order(myKaplanMeier$time)]), ncol=2)
  colnames(orderMatrix) = c("new", "old")
  originalOrder = orderMatrix[order(orderMatrix[,2]), 1]
  
  dpresid.dtheta = matrix(0, nrow = uniqueFailureN, ncol = N)
  colnames(dpresid.dtheta) = orderedTime
  rownames(dpresid.dtheta) = orderedUniqueFailureTime
  condition1 = (orderedTime == orderedUniqueFailureTime[1]) | ((orderedTime <= orderedUniqueFailureTime[1]) & (orderedDelta == 0))
  dpresid.dtheta[1, (1:N)[condition1]] = 1

  phi_ji = matrix(0, nrow = uniqueFailureN, ncol = N)
  colnames(phi_ji) = orderedTime
  rownames(phi_ji) = orderedUniqueFailureTime
  indices = (1:N)[orderedTime <= orderedUniqueFailureTime[1]]
  phi_ji[1, indices] = 1

  for (j in 2:uniqueFailureN){
    condition2 = (orderedTime == orderedUniqueFailureTime[j] & (orderedDelta == 1))
    condition3 = ((orderedTime > orderedUniqueFailureTime[j-1]) & (orderedTime <= orderedUniqueFailureTime[j]) & (orderedDelta == 0))
    dpresid.dtheta[j, (1:N)[condition2 | condition3]] = 1
    dpresid.dtheta[j-1, (1:N)[condition2]] = 1
    indices = (1:N)[orderedTime <= orderedUniqueFailureTime[j]]
    phi_ji[j, indices] = 1
  }
  condition4 = (orderedTime > orderedUniqueFailureTime[uniqueFailureN] & (orderedDelta == 0))
  dpresid.dtheta[uniqueFailureN, (1:N)[condition4]] = 1
  orderedDeltaMatr = matrix(orderedDelta, nrow = uniqueFailureN, ncol = N, byrow=TRUE)
  
  H = H0 = H1 = rep(NA, N)
  for (i in 1:N){
    H[i] = mean(orderedTime <= orderedTime[i])
    H0[i] = mean((orderedTime <= orderedTime[i])*(1 - orderedDelta))
    H1[i] = mean((orderedTime <= orderedTime[i])*orderedDelta)
  }
  
  H0dv = diff(c(0, H0))
  H1dw = diff(c(0, H1))
  H1dwPerY = H0dvPerY = rep(0, N)
  H0dvPerY[orderedTime * (1 - orderedDelta) == orderedTime] = H0dv[orderedTime * (1 - orderedDelta) == orderedTime]
  H1dwPerY[orderedTime * orderedDelta == orderedTime] = H1dw[orderedTime * orderedDelta == orderedTime]
  
  ################ making sure that 1-H is never zero: Bryan's suggestion
  Hadj = H
  Hadj[Hadj == 1] = .999999
  
  multiplier = 1/(1 - Hadj) 
  multiplier[H==1] = 0   ### this porbably fixes the the problem of division by 0 anyway
  
  gamma0 = exp(cumsum(c(0, (H0dvPerY * multiplier)  )))[1:N]
  Vji = gamma_j2 = gamma_j1 = phi_ji*0
  vValue = wValue  = orderedTime
  
  for(j in 1:nrow(gamma_j1)){    
    for(i in 1:ncol(gamma_j1)){  
      indicatorForGamma1 = as.numeric((orderedTime[i] < wValue) & phi_ji[j, ]) 
      gamma_j1[j, i] =  multiplier[i]  *  sum(indicatorForGamma1  *  gamma0  *  H1dwPerY)
    }
  } 
  
  for(j in 1:nrow(gamma_j1)){    
    for(i in 1:ncol(gamma_j1)){  
      indicatorForGamma2 = as.numeric((vValue < orderedTime[i]))
      gamma_j2[j, i] =   sum( multiplier * indicatorForGamma2 * gamma_j1[j, ] * H0dvPerY)
      Vji[j, i] = phi_ji[j, i] * gamma0[i] * orderedDelta[i] + gamma_j1[j, i] * (1-orderedDelta[i]) - gamma_j2[j, i]
    }
  }   
  res = Vji[, originalOrder] - (1 - matrix(orderedUniqueKM, nrow = uniqueFailureN, ncol = N))
  
  list(dl.dtheta = res, dpresid.dtheta = dpresid.dtheta[, originalOrder])
}


