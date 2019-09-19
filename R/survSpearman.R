#' @export
################################################################### Highest Rank and Restrcited Spearman's Rho in the Restricted
################################################################### Region
survSpearman = function(bivarSurf, tauX = Inf, tauY = Inf) {
    ### The function returns two correlation values: 1)
    ### 'HighestRankInRestrictedRegion' is computed using Highest Rank
    ### approach in the restricted region 2) 'RestrictedRegion' is
    ### computed in the restricted region
    
    ### Find the last observed events and censoring events
    lastX = rownames(bivarSurf)[nrow(bivarSurf)]
    lastY = colnames(bivarSurf)[ncol(bivarSurf)]
    lastEventOnX = names(bivarSurf[, 1] == min(bivarSurf[, 1]))[bivarSurf[, 
        1] == min(bivarSurf[, 1])][1]
    lastEventOnY = names(bivarSurf[1, ] == min(bivarSurf[1, ]))[bivarSurf[1, 
        ] == min(bivarSurf[1, ])][1]
    
    ### Find the last observed events and censoring events in the
    ### restricted region
    restrBivarSurf = bivarSurf
    restrBivarSurf = restrBivarSurf[as.numeric(rownames(restrBivarSurf)) < 
        tauX, ]
    restrBivarSurf = restrBivarSurf[, as.numeric(colnames(restrBivarSurf)) < 
        tauY]
    restrLastX = rownames(restrBivarSurf)[nrow(restrBivarSurf)]
    restrLastY = colnames(restrBivarSurf)[ncol(restrBivarSurf)]
    restrLastEventOnX = names(restrBivarSurf[, 1] == min(restrBivarSurf[, 
        1]))[restrBivarSurf[, 1] == min(restrBivarSurf[, 1])][1]
    restrLastEventOnY = names(restrBivarSurf[1, ] == min(restrBivarSurf[1, 
        ]))[restrBivarSurf[1, ] == min(restrBivarSurf[1, ])][1]
    
    actualTauX = paste(restrLastEventOnX, "+")
    actualTauY = paste(restrLastEventOnY, "+")
    
    ### compute restricted correlation
    res1 = HighestRankAndRestrictedSpearman(bivarSurf, tauX = tauX, 
        tauY = tauY)
    bivarSurfForHR = bivarSurf
    bivarSurfForHR[as.numeric(rownames(bivarSurfForHR)) >= tauX, 
        ] = 0
    bivarSurfForHR[, as.numeric(colnames(bivarSurfForHR)) >= tauY] = 0
    
    ### compute highest rank correlation in the restricted region
    res2 = HighestRankAndRestrictedSpearman(bivarSurfForHR, tauX = Inf, 
        tauY = Inf)
    
    resCor = c(HighestRank = res2["HighestRank", ], Restricted = res1["Restricted", 
        ])
    resRegUser = c(tauX = as.character(tauX), tauY = as.character(tauY))
    resReg = c(tauX = actualTauX, tauY = actualTauY)
    list(`Restricted region set by user` = resRegUser, `Effective restricted region` = resReg, 
        Correlation = resCor)
}
