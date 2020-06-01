###############################################
############################################### M-estimation for survival regressoin
###############################################
mestimation.coxph = function(object1, object2, method = "full"){
  if(method == "old"){
    scoreRes1 = mestimation.coxph_Old(object = object1)
    scoreRes2 = mestimation.coxph_Old(object = object2)
  }
  if(method == "full"){
    scoreRes1 = mestimation.coxph.full.likelihood(object = object1, continuos = FALSE)
    scoreRes2 = mestimation.coxph.full.likelihood(object = object2, continuos = FALSE)
    if(det(scoreRes1$d2l.dtheta.dtheta) == 0){
      scoreRes1 = mestimation.coxph_Old(object = object1)
      cat("FUNCTION mestimation.coxph: Reverted to M-estimation assuming known hazard for object1\n")
      warning("FUNCTION mestimation.coxph: Reverted to M-estimation assuming known hazard for object1")
    }
    if(det(scoreRes2$d2l.dtheta.dtheta) == 0){
      scoreRes2 = mestimation.coxph_Old(object = object2)
      cat("FUNCTION mestimation.coxph: Reverted to M-estimation assuming known hazard for object2\n")
      warning("FUNCTION mestimation.coxph: Reverted to M-estimation assuming known hazard for object2")
    }
  }
  cor1 = corTS.survreg.modified(xresid = scoreRes1$presid, yresid = scoreRes2$presid,
               xz.dl.dtheta = t(scoreRes1$dl.dtheta),
               yz.dl.dtheta = t(scoreRes2$dl.dtheta),
               xz.d2l.dtheta.dtheta = scoreRes1$d2l.dtheta.dtheta,
               yz.d2l.dtheta.dtheta = scoreRes2$d2l.dtheta.dtheta,
               dxresid.dthetax = scoreRes1$dpresid.dtheta,
               dyresid.dthetay = scoreRes2$dpresid.dtheta, inverseA = FALSE)
  cor1
}

