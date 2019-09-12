
################################################### generates clayton's or frank's copula data
dataSim = function(nSim, thetaPar = 0.1432, lambda1 = 1, lambda2 = 1, 
    censoringProb1 = 0, censoringProb2 = 0, independentCensoring = TRUE, 
    censProbForDependCens = 0, family = "clayton") {
    ### Allows two types of censoring: independent (b/w C1 and C2) and
    ### dependent (C1 = C2)
    if (family == "clayton") {
        if (thetaPar == 0) {
            ###### trying to decide which one is better
            t1t2 = rarchi(n = nSim, "clayton", d = 2, theta = 0)
        } else {
            cop12 = archmCopula(family = "clayton", param = thetaPar)
            t1t2 <- rCopula(nSim, cop12)
        }
    } else {
        if (family != "frank") 
            stop("Only Clayton's and Frank's families were implemented\n")
        if (thetaPar == 0) {
            t1t2 = rarchi(n = nSim, "clayton", d = 2, theta = 0)
        } else {
            cop12 = archmCopula(family = "frank", param = thetaPar)
            t1t2 <- rCopula(nSim, cop12)
        }
    }
    t1 = -log(1 - t1t2[, 1])/lambda1
    t2 = -log(1 - t1t2[, 2])/lambda2
    delta1 = delta2 = rep(1, nSim)
    if (independentCensoring) {
        if (censoringProb1 > 0 & censoringProb1 < 1) {
            beta1 = lambda1/censoringProb1 - lambda1
            censoringTime1 = rexp(n = nSim, rate = 1/beta1)
            delta1 = as.numeric(t1 <= censoringTime1)
            t1[delta1 == 0] = censoringTime1[delta1 == 0]
        }
        if (censoringProb2 > 0 & censoringProb2 < 1) {
            beta2 = lambda2/censoringProb2 - lambda2
            censoringTime2 = rexp(n = nSim, rate = 1/beta2)
            delta2 = as.numeric(t2 <= censoringTime2)
            t2[delta2 == 0] = censoringTime2[delta2 == 0]
        }
    } else {
        if (censProbForDependCens > 0 & censProbForDependCens < 1) {
            beta3 = 1/censProbForDependCens - 1
            censoringTime3 = rexp(n = nSim, rate = 1/beta3)
            delta1 = as.numeric(t1 <= censoringTime3)
            delta2 = as.numeric(t2 <= censoringTime3)
            t1[delta1 == 0] = censoringTime3[delta1 == 0]
            t2[delta2 == 0] = censoringTime3[delta2 == 0]
        }
    }
    res = cbind(t1, delta1, t2, delta2)
    colnames(res) = c("t1", "delta1", "t2", "delta2")
    res
}

