
################################################## Random and type-I censoring simulation
simulateRealisticData = function(nSubj = 10, thetaPar, family, restrictedTimeX = Inf, 
    restrictedTimeY = Inf, censoringProb1 = 0, censoringProb2 = 0, 
    dependCensProb = 0, independentCensoring = TRUE) {
    firstStep = dataSim(nSim = nSubj, thetaPar = thetaPar, lambda1 = 1, 
        lambda2 = 1, censoringProb1 = censoringProb1, censoringProb2 = censoringProb2, 
        independentCensoring = independentCensoring, censProbForDependCens = dependCensProb, 
        family = family)
    data = data.frame(X = firstStep[, "t1"], deltaX = firstStep[, 
        "delta1"], Y = firstStep[, "t2"], deltaY = firstStep[, "delta2"])
    moreXCens = data$X > restrictedTimeX
    moreYCens = data$Y > restrictedTimeY
    data$X[moreXCens] = restrictedTimeX
    data$deltaX[moreXCens] = 0
    data$Y[moreYCens] = restrictedTimeY
    data$deltaY[moreYCens] = 0
    data
}

