################################################################### Overall and Restrcited Spearman's Rho
HighestRankAndRestrictedSpearman = function(bivarSurf, tauX = Inf, 
    tauY = Inf) {
    ######## computes 1) Overall Spearmans correlation using highest rank
    ######## approach 2) Restricted Spearmans correlation 'bivarSurf' is a
    ######## matrix representing a bivariate surface fit with the following
    ######## format: 1) the 1st row is the marginal survival function for Y
    ######## 2) the 1nd column is the marginal survival function for X 3)
    ######## bivarSurf[1,1] equals to 1 4) row names of 'bivarSurf' are
    ######## ordered X values including 0 5) column names of 'bivarSurf' are
    ######## ordered Y values including 0 Note: the covariance (the
    ######## numerator of the correlation) of the highest rank estimator can
    ######## be also computed using formula: Numerator = \sum_i \sum_j
    ######## S(x_i,y_i)S_X(dx_i)S_Y(dy_i) - S(x_i,dy_i)S_X(dx_i)S_Y(y_i) -
    ######## S(dx_i,y_i)S_X(x_i)S_Y(dy_i) + S(dx_i,dy_i)S_X(x_i)S_Y(y_i) The
    ######## above is equivalent to the highest rank estimator Numerator =
    ######## \sum _{i^*}\sum _{j^*} \left\{1 - \widehat{S}^H_X(x_{i^*})
    ######## - \widehat{S}^H_X(x_{i^*}^-) \right\} \left\{1 -
    ######## \widehat{S}^H_Y(y_{j^*}) - \widehat{S}^H_Y(y_{j^*}^-)
    ######## \right\} \widehat{S}^H(dx_{i^*}, dy_{j^*})
    
    ######## By setting the survival proability to zero outside of the data
    ######## support we effectively assume a left-over probability mass
    ######## which is what highest rank estimator is doing This ensures that
    ######## the above two estimators are the same
    oldrownames = rownames(bivarSurf)
    oldcolnames = colnames(bivarSurf)
    bivarSurfnew = bivarSurf
    bivarSurfnew = cbind(bivarSurfnew, rep(0, nrow(bivarSurfnew)))
    bivarSurfnew = rbind(bivarSurfnew, rep(0, ncol(bivarSurfnew)))
    rownames(bivarSurfnew) = c(oldrownames, "Inf")
    colnames(bivarSurfnew) = c(oldcolnames, "Inf")
    bivarSurf = bivarSurfnew
    ######## End of the part that ensures that the above two estimators are
    ######## the same
    
    Xs = as.numeric(rownames(bivarSurf))
    Ys = as.numeric(colnames(bivarSurf))
    rP1 = rP2 = c(NA, NA)
    
    ######## Identify the resticted region for X:
    maxX = max(Xs[Xs < tauX])
    maxY = max(Ys[Ys < tauY])
    if ((tauX <= Inf) & (sum(Xs < tauX) <= nrow(bivarSurf))) {
        rowsInRegion = Xs <= maxX
        rX1 = 1:(sum(rowsInRegion) - 1)
    } else {
        rX1 = 1:(nrow(bivarSurf) - 2)
    }
    rP1[1] = rX1[length(rX1)] + 1
    rP2[1] = rP1[1] + 1
    rX2 = rX1 + 1
    
    ######## Identify the resticted region for Y:
    if ((tauY <= Inf) & (sum(Ys < tauY) <= ncol(bivarSurf))) {
        colsInRegion = Ys <= maxY
        rY1 = 1:(sum(colsInRegion) - 1)
    } else {
        rY1 = 1:(ncol(bivarSurf) - 2)
    }
    rP1[2] = rY1[length(rY1)] + 1
    rP2[2] = rP1[2] + 1
    rY2 = rY1 + 1
    nX = nrow(bivarSurf) - 1
    nY = ncol(bivarSurf) - 1
    rangeX = 1 + (1:nX)
    rangeY = 1 + (1:nY)
    
    unitVecX = matrix(1, ncol = nX)
    unitVecY = matrix(1, ncol = nY)
    Sx = matrix(as.numeric(bivarSurf[(1:nX) + 1, 1]), nrow = nX) %*% 
        unitVecY
    Sy = t(matrix(as.numeric(bivarSurf[1, (1:nY) + 1]), nrow = nY) %*% 
        unitVecX)
    SxM = matrix(as.numeric(bivarSurf[1:nX, 1]), nrow = nX) %*% unitVecY
    SyM = t(matrix(as.numeric(bivarSurf[1, 1:nY]), nrow = nY) %*% 
        unitVecX)
    Sdx = SxM - Sx
    Sdy = SyM - Sy
    SxMyM = rbind(rep(NA, nY + 1), cbind(rep(NA, nX), bivarSurf[1:nX, 
        1:nY]))
    SxM_y = rbind(rep(NA, nY + 1), bivarSurf[1:nX, ])
    Sx_yM = cbind(rep(NA, nX + 1), bivarSurf[, 1:nY])
    Sdx_y = SxM_y - bivarSurf
    Sx_dy = Sx_yM - bivarSurf
    Sdx_yM = SxMyM - Sx_yM
    SxM_dy = SxMyM - SxM_y
    Sdxdy = SxM_dy - Sx_dy
    
    ### my favorite estimator
    numerator = sum(SxMyM[rangeX, rangeY] * Sdx * Sdy - SxM_dy[rangeX, 
        rangeY] * Sdx * SyM - Sdx_yM[rangeX, rangeY] * SxM * Sdy + 
        Sdxdy[rangeX, rangeY] * SxM * SyM)
    
    ### Highest rank estimator, which is equivalent to my favorite
    ### estimator when compu
    numerator0 = sum((2 * SxM - Sdx - 1) * (2 * SyM - Sdy - 1) * 
        Sdxdy[rangeX, rangeY], na.rm = TRUE)
    
    Sx_Stuff_Tail = c(Sx[, 1], 0)
    SxM_Stuff_Tail = c(1, Sx[, 1])
    Sdx_Stuff_Tail = SxM_Stuff_Tail - Sx_Stuff_Tail
    Sy_Stuff_Tail = c(Sy[1, ], 0)
    SyM_Stuff_Tail = c(1, Sy[1, ])
    Sdy_Stuff_Tail = SyM_Stuff_Tail - Sy_Stuff_Tail
    PSR_X_like = (1 - Sx_Stuff_Tail - SxM_Stuff_Tail)
    PSR_Y_like = (1 - Sy_Stuff_Tail - SyM_Stuff_Tail)
    PSR_X_mean = sum(PSR_X_like * Sdx_Stuff_Tail)
    PSR_Y_mean = sum(PSR_Y_like * Sdy_Stuff_Tail)
    PSR_X_var = sum(Sdx_Stuff_Tail * (PSR_X_like - PSR_X_mean)^2)
    PSR_Y_var = sum(Sdy_Stuff_Tail * (PSR_Y_like - PSR_Y_mean)^2)
    rhoConst = 1/sqrt(PSR_X_var * PSR_Y_var)
    
    ####### Compute overall rho
    localRhoRhoConst = rhoConst * numerator
    prsLikeNoTails = rhoConst * numerator0
    
    ####### Restricted Estimator
    ProbOmegaRMatr = 1 - Sx[rX1, rY1] - Sy[rX1, rY1] + bivarSurf[rX2, 
        rY2]
    if (length(rX1) == 1 | length(rY1) == 1) {
        ProbOmegaRMatr = matrix(ProbOmegaRMatr, nrow = length(rX1), 
            ncol = length(rY1))
        ProbOmegaR = ProbOmegaRMatr[nrow(ProbOmegaRMatr), ncol(ProbOmegaRMatr)]
    } else {
        ProbOmegaR = ProbOmegaRMatr[nrow(ProbOmegaRMatr), ncol(ProbOmegaRMatr)]
    }
    SX_tauXM = Sx[rX1[length(rX1)], 1]
    SY_tauYM = Sy[1, rY1[length(rY1)]]
    
    ####### Compute S(x, tauY-), S(x-, tauY-), S(tauX-, y), S(tauX-, y)
    S_x_TauYM = matrix(bivarSurf[rX2, rY2[length(rY2)]], nrow = length(rX2), 
        ncol = length(rY2))
    S_xM_TauYM = matrix(SxM_y[rX2, rY2[length(rY2)]], nrow = length(rX2), 
        ncol = length(rY2))
    S_TauXM_y = matrix(bivarSurf[rX2[length(rX2)], rY2], nrow = length(rX2), 
        ncol = length(rY2), byrow = TRUE)
    S_TauXM_yM = matrix(Sx_yM[rX2[length(rX2)], rY2], nrow = length(rX2), 
        ncol = length(rY2), byrow = TRUE)
    g_X_No_OmegaR = (2 - Sx[rX1, rY1] - SxM[rX1, rY1] - 2 * SY_tauYM + 
        S_x_TauYM + S_xM_TauYM - ProbOmegaR)
    g_Y_No_OmegaR = (2 - Sy[rX1, rY1] - SyM[rX1, rY1] - 2 * SX_tauXM + 
        S_TauXM_y + S_TauXM_yM - ProbOmegaR)
    
    ####### Compute SX(dx|OmegaR), SY(dy|OmegaR) for covariance According
    ####### to our paper, it is: SX(dx) - S(dx, tauYM)
    SX_dx_No_OmR = (Sdx[rX1, 1] - (S_xM_TauYM[, 1] - S_x_TauYM[, 
        1]))
    SY_dy_No_OmR = (Sdy[1, rY1] - (S_TauXM_yM[1, ] - S_TauXM_y[1, 
        ]))
    SX_dx_OmR = SX_dx_No_OmR/ProbOmegaR
    SY_dy_OmR = SY_dy_No_OmR/ProbOmegaR
    
    ####### Compute SX(dx|OmegaR), SY(dy|OmegaR) for rho (the normalizing
    ####### const) take care of negative variance
    SX_dx_No_OmR_const = pmax(SX_dx_No_OmR, 0)
    SY_dy_No_OmR_const = pmax(SY_dy_No_OmR, 0)
    SX_dx_No_OmR_const = ProbOmegaR * SX_dx_No_OmR_const/sum(SX_dx_No_OmR_const)
    SY_dy_No_OmR_const = ProbOmegaR * SY_dy_No_OmR_const/sum(SY_dy_No_OmR_const)
    
    ####### Compute terms related to variance
    var_g_X_No_OmR = sum(((g_X_No_OmegaR[, 1] - sum(g_X_No_OmegaR[, 
        1] * SX_dx_No_OmR)/ProbOmegaR)^2) * SX_dx_No_OmR)
    var_g_Y_No_OmR = sum(((g_Y_No_OmegaR[1, ] - sum(g_Y_No_OmegaR[1, 
        ] * SY_dy_No_OmR)/ProbOmegaR)^2) * SY_dy_No_OmR)
    
    ####### If the probability of being in the region is negative -> the
    ####### restricted correlation is not defined
    if (ProbOmegaR <= 0) {
        warning("Restricted Spearman's correlation is not defined because the probability of being in the restricted region is less or equal to zero.")
        rhoX0Y0_restr_No_OmR = NA
    } else {
        if (var_g_X_No_OmR < 0 | var_g_Y_No_OmR < 0) {
            cat("Restricted Spearman's correlation was corrected.\n")
            warning("Restricted Spearman's correlation was corrected.")
        }
        if (var_g_X_No_OmR < 0) {
            var_g_X_No_OmR = sum(((g_X_No_OmegaR[, 1] - sum(g_X_No_OmegaR[, 
                1] * SX_dx_No_OmR_const)/ProbOmegaR)^2) * SX_dx_No_OmR_const)
        }
        if (var_g_Y_No_OmR < 0) {
            var_g_Y_No_OmR = sum(((g_Y_No_OmegaR[1, ] - sum(g_Y_No_OmegaR[1, 
                ] * SY_dy_No_OmR_const)/ProbOmegaR)^2) * SY_dy_No_OmR_const)
        }
        
        ####### Compute the constant for restricted correlation:
        inv_rhoCons_cond_form_No_OmR = sqrt(var_g_X_No_OmR * var_g_Y_No_OmR)
        ####### Compute restricted correlation:
        rhoX0Y0_restr_No_OmR = sum(g_X_No_OmegaR * g_Y_No_OmegaR * 
            Sdxdy[rX2, rY2], na.rm = TRUE)/inv_rhoCons_cond_form_No_OmR
    }
    # res = matrix(c(localRhoRhoConst, prsLikeNoTails,
    # rhoX0Y0_restr_No_OmR), ncol = 1) rownames(res) =
    # c('HighestRank', 'HighestRankNoTails', 'Restricted')
    res = matrix(c(prsLikeNoTails, rhoX0Y0_restr_No_OmR), ncol = 1)
    rownames(res) = c("HighestRank", "Restricted")
    colnames(res) = "Spearman's Correlation"
    
    ####### Make sure that all values are in the proper range
    for (n_i in rownames(res)) {
        res[n_i, ] = min(abs(res[n_i, ]), 1) * sign(res[n_i, ])
    }
    res
}
