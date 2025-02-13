# This function will calculate the matrix GammaEjperp, P, Q, R
# detailed math can be found in manuscript Page 5

PQR = function(joint_distcal, select_E) {
  # Input:
  # joint_distcal: the result of function joint_dist
  # select_E: the result of variable_selection (function variable_selection_PY)
  
  betaEjperp = joint_distcal[["betaEjperp_cal"]][["betaEjperp"]]
  betaEperp = joint_distcal[["betaEperp_cal"]][["betaEperp"]]
  HEE = joint_distcal[["H"]][["HEE"]]
  HNEE = joint_distcal[["H"]][["HNEE"]]
  sigmaTT = joint_distcal[["Sigma"]][["sigmaTT"]]
  sigmaST = joint_distcal[["Sigma"]][["sigmaST"]]
  E = joint_distcal[["E"]]
  NE = joint_distcal[["NE"]]
  n = joint_distcal[["n"]]
  ej = joint_distcal[["ej"]]
  
  Enum = length(E)+1 # add one because E doesn't include intercept
  NEnum = length(NE)
  pnum = Enum + NEnum
  
  ###### Below code need to edit after you build randomized lasso ######
  lam = select_E[["lam"]]
  if(length(lam) > 1){lam = lam(which(lam != 0)[1])} # need to edit it if do group lasso
  se = select_E[["sign_soln"]]
  if(length(se) < Enum) {se = c(sign(joint_distcal[["betaEM"]][["coefficients"]][["(Intercept)"]]),se)} 
  
  # Gamma Ej perp
  GammaEjPerp = rbind(betaEjperp, betaEperp)
  
  # matrix P
  p1 = HEE %*% sigmaTT %*% ej/c(t(ej) %*% sigmaTT %*% ej)
  p2 = HEE
  p3 = matrix(0, nrow = Enum, ncol = NEnum)
  p4 = (sigmaST %*% ej + HNEE %*% sigmaTT %*% ej)/c(t(ej) %*% sigmaTT %*% ej)
  p5 = sigmaST %*% solve(sigmaTT) + HNEE
  p6 = diag(rep(1,NEnum))
  
  Pup = cbind(p1, p2, p3)
  Pdown = cbind(p4, p5, p6)
  
  P = rbind(Pup, Pdown) * sqrt(n)
  
  # matrix Q
  Q1 = HEE * c(-sqrt(n))
  Q2 = matrix(0, nrow = Enum, ncol = NEnum)
  Q3 = HNEE * c(-sqrt(n))
  Q4 = diag(rep(lam, NEnum)) 
  
  Qup = cbind(Q1,Q2)
  Qdown = cbind(Q3,Q4)
  
  Q = rbind(Qup, Qdown)
  
  # matrix R
  Rup = matrix(se*c(lam), ncol = 1)
  R = rbind(Rup, matrix(rep(0, NEnum), ncol = 1))
  
  return(list(GammaEjPerp = GammaEjPerp,
              P = P,
              Q = Q,
              R = R,
              Enum = Enum,
              NEnum = NEnum, 
              pnum = pnum, 
              p1 = p1* sqrt(n),
              p4 = p4* sqrt(n),
              se = se))
  
  # Output:
  # Enum: number of selected parameters
  # NEnum: number of unselected parameters
  # pnum: number of all potential parameters
  # p1: HEE %*% sigmaTT %*% ej/c(t(ej) %*% sigmaTT %*% ej)
  # p4: (sigmaST %*% ej + HNEE %*% sigmaTT %*% ej)/c(t(ej) %*% sigmaTT %*% ej)
  # se: Signs of the estimated coefficients for selected variables.
}

############################################################
# PQR allow intercept being penalized
PQR_Pint = function(joint_distcal, select_E) {
  # Input:
  # joint_distcal: the result of function joint_dist
  # select_E: the result of variable_selection (function variable_selection_PY)
  
  betaEjperp = joint_distcal[["betaEjperp_cal"]][["betaEjperp"]]
  betaEperp = joint_distcal[["betaEperp_cal"]][["betaEperp"]]
  HEE = joint_distcal[["H"]][["HEE"]]
  HNEE = joint_distcal[["H"]][["HNEE"]]
  sigmaTT = joint_distcal[["Sigma"]][["sigmaTT"]]
  sigmaST = joint_distcal[["Sigma"]][["sigmaST"]]
  E = joint_distcal[["E"]]
  NE = joint_distcal[["NE"]]
  n = joint_distcal[["n"]]
  ej = joint_distcal[["ej"]]
  
  Enum = length(E)
  NEnum = length(NE)
  pnum = Enum + NEnum
  
  ###### Below code need to edit after you build randomized lasso ######
  lam = select_E[["lam"]]
  if(length(lam) > 1){lam = lam(which(lam != 0)[1])} # need to edit it if do group lasso
  se = select_E[["sign_soln"]]
  
  # Gamma Ej perp
  GammaEjPerp = rbind(betaEjperp, betaEperp)
  
  # matrix P
  p1 = HEE %*% sigmaTT %*% ej/c(t(ej) %*% sigmaTT %*% ej)
  p2 = HEE
  p3 = matrix(0, nrow = Enum, ncol = NEnum)
  p4 = (sigmaST %*% ej + HNEE %*% sigmaTT %*% ej)/c(t(ej) %*% sigmaTT %*% ej)
  p5 = sigmaST %*% solve(sigmaTT) + HNEE
  p6 = diag(rep(1,NEnum))
  
  Pup = cbind(p1, p2, p3)
  Pdown = cbind(p4, p5, p6)
  
  P = rbind(Pup, Pdown) * sqrt(n)
  
  # matrix Q
  Q1 = HEE * c(-sqrt(n))
  Q2 = matrix(0, nrow = Enum, ncol = NEnum)
  Q3 = HNEE * c(-sqrt(n))
  Q4 = diag(rep(lam, NEnum)) 
  
  Qup = cbind(Q1,Q2)
  Qdown = cbind(Q3,Q4)
  
  Q = rbind(Qup, Qdown)
  
  # matrix R
  Rup = matrix(se*c(lam), ncol = 1)
  R = rbind(Rup, matrix(rep(0, NEnum), ncol = 1))
  
  return(list(GammaEjPerp = GammaEjPerp,
              P = P,
              Q = Q,
              R = R,
              Enum = Enum,
              NEnum = NEnum, 
              pnum = pnum, 
              p1 = p1* sqrt(n),
              p4 = p4* sqrt(n),
              se = se))
  
  # Output:
  # Enum: number of selected parameters
  # NEnum: number of unselected parameters
  # pnum: number of all potential parameters
  # p1: HEE %*% sigmaTT %*% ej/c(t(ej) %*% sigmaTT %*% ej)
  # p4: (sigmaST %*% ej + HNEE %*% sigmaTT %*% ej)/c(t(ej) %*% sigmaTT %*% ej)
  # se: Signs of the estimated coefficients for selected variables.

}






