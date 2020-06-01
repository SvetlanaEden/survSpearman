############################################### Hazard Stretcher:
### it turns out that coxph does not stretch hazard for multiple events per time point
###############################################
hazardStretcher_take_two = function(data, time){
  ### data - data with hazard, time  what is supplied by cox model
  ### time - is the original time variable (with duplicates and everything...)
  
  ### data = data.frame(hazard = c(0, 0, .11, .11, .11, .22, .33, .33), time = c(1.1, 1.2, 2.2, 2.22, 2.3, 3, 4, 5.5)); time = c(1.1, 1.2, 2.2, 2.22, 2.22, 2.3, 3, 4, 4, 4, 5.5)
  ### data = data.frame(hazard = c(0, 0, .11, .11, .11, .22, .33, .44), delta = c(0, 0, 1, 0, 0, 1, 1, 1), time = c(1.1, 1.2, 2, 2.2, 2.3, 3, 4, 5.5))
  ### data = data.frame(hazard = c(.11, .11, .22, .33, .44), delta = c(1, 0, 1, 1, 1), time = c(2, 2.2, 3, 4, 5.5))
  ### data = data.frame(hazard = c(.11, .11, .22, .33, .44), delta = c(1, 0, 1, 1, 1), time = c(2, 2.2, 3, 4, 5.5))
  
  ######################### important note:
  ### when this funciton is applied to coxph baseline hazard (model1), 
  ### it is important that time comes from object$y[,1] because of issue with
  ### number representation (see argument timefix in coxph.control and Terri Therneau's email)

  if(!all(c("hazard", "time") %in% names(data))){
    stop("The data has to have the following names: hazard", "time")
  }
  
  orderedTime = time[order(time)]
  data = data[order(data$time),]
  uniqueHazard = unique(data[, c("hazard", "time")])
  #uniqueHazard = rbind(c(0, 1), uniqueHazard)
  uniqueHazard = rbind(c(0, 0), uniqueHazard)
  data$Hminus = NA
  
  ##################### stretch the data:
  ##################### define hazard for every point of the original time:
  orig_i = 1
  stretchedData = data.frame(time = orderedTime, hazard = NA)

  for(i in 1:nrow(uniqueHazard)){
    while((uniqueHazard$time[i] == orderedTime[orig_i]) & (orig_i <= length(orderedTime))){
      stretchedData$hazard[orig_i] = uniqueHazard$hazard[i]
      orig_i = orig_i + 1
    }
  }
  
  i_minus = 1
  i_unique = 2
  while(stretchedData$time[i_minus] <= uniqueHazard$time[nrow(uniqueHazard)] & i_minus <= nrow(stretchedData) & i_unique < nrow(uniqueHazard) ){
    ### fill in censored observations 
    if(stretchedData$time[i_minus] < uniqueHazard$time[i_unique+1]){
      stretchedData$Hminus[i_minus] = uniqueHazard$hazard[i_unique-1]
      i_minus = i_minus + 1
    }else{
      stretchedData$Hminus[i_minus] = uniqueHazard$hazard[i_unique]
      i_minus = i_minus + 1
      i_unique = i_unique + 1
    }
  }
  
  ### if there are any observations that are censored after the last failure point.
  while(i_minus <= nrow(stretchedData)){
    stretchedData$Hminus[i_minus] = uniqueHazard$hazard[nrow(uniqueHazard) - 1]
    i_minus = i_minus + 1
  }
  
  stretchedData
}

