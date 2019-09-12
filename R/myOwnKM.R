################################################################### computes KM
myOwnKM = function(time, delta, returnToOriginalOrder = TRUE) {
    ###-----------------------------------------------
    ### time - time to event delta - event indicator
    ### returnToOriginalOrder = TRUE - the order of time values is the
    ### same as in original data it also returns the delta (event
    ### indicator) returnToOriginalOrder = FALSE - the time is unique
    ### and ordered and the resulting data frame does not contain delta
    ### (event indicator).
    uniqueAndOrderedTime = unique(time)[order(unique(time))]
    nEvents = tapply(delta, time, sum)
    nDrop = tapply(delta, time, length)
    atRisk = c(length(time), length(time) - cumsum(nDrop))[1:length(nDrop)]
    probForEachTime = (1 - nEvents/atRisk)
    dataKM = data.frame(time = uniqueAndOrderedTime, nEvents = nEvents, 
        atRisk = atRisk, KM = cumprod(probForEachTime))
    dataKM$CDF = 1 - dataKM$KM
    dataKM$CDF_M = c(0, dataKM$CDF[1:(nrow(dataKM) - 1)])
    if (returnToOriginalOrder) {
        rownames(dataKM) = uniqueAndOrderedTime
        ### let's order the output according to the original order of time
        dataKM = dataKM[as.character(time), ]
        dataKM$delta = delta
    }
    rownames(dataKM) = 1:nrow(dataKM)
    dataKM
}

