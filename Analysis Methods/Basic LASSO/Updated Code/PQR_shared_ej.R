PQR_Pint_shared = function(joint_distcal_shared, select_E) {
  # joint_distcal_shared: the result of function joint_dist
  # select_E: the result of function variable_selection_PY_penal_int
  
  HEE = joint_distcal_shared[["H"]][["HEE"]]
  HNEE = joint_distcal_shared[["H"]][["HNEE"]]
  sigmaTT = joint_distcal_shared[["Sigma"]][["sigmaTT"]]
  sigmaST = joint_distcal_shared[["Sigma"]][["sigmaST"]]
  E = joint_distcal_shared[["E"]]
  NE = joint_distcal_shared[["NE"]]
  n = joint_distcal_shared[["n"]]
  
  Enum = length(E)
  NEnum = length(NE)
  pnum = Enum + NEnum
  
  ###### Below code need to edit after you build randomized lasso ######
  lam = select_E[["lam"]]
  if(length(lam) > 1){lam = lam(which(lam != 0)[1])} # need to edit it if do group lasso
  se = select_E[["sign_soln"]]
  
  # matrix P shared part
  p2 = HEE
  p3 = matrix(0, nrow = Enum, ncol = NEnum)
  p5 = sigmaST %*% solve(sigmaTT) + HNEE
  p6 = diag(rep(1,NEnum))
  
  Pshareup = cbind(p2, p3)
  Psharedown = cbind(p5, p6)
  
  P_share = rbind(Pshareup, Psharedown) * sqrt(n)
  
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
  
  return(list(P_share = P_share,
              Q = Q,
              R = R,
              Enum = Enum,# number of selected parameters
              NEnum = NEnum, # number of unselected parameters
              pnum = pnum, # number of all potential parameters
              se = se))
  
  # Output:
  # P_share: the part of P matrix that will not change with ej value
  # Q: Q matrix
  # R: R matrix
  # Enum: number of selected parameters
  # NEnum: number of unselected parameters
  # pnum: number of all potential parameters
  # se: Signs of the estimated coefficients for selected variables.

}

PQR_Pint_ej = function(PQR_Pint_shared, joint_dist_Penal_Int_ej, joint_distcal_shared) {
  
  betaEjperp = joint_dist_Penal_Int_ej[["betaEjperp_cal"]][["betaEjperp"]]
  ej = joint_dist_Penal_Int_ej[["ej"]]
  betaEperp = joint_distcal_shared[["betaEperp_cal"]][["betaEperp"]]
  sigmaTT = joint_distcal_shared[["Sigma"]][["sigmaTT"]]
  P_shared = PQR_Pint_shared[["P_share"]]
  HEE = joint_distcal_shared[["H"]][["HEE"]]
  sigmaST = joint_distcal_shared[["Sigma"]][["sigmaST"]]
  HNEE = joint_distcal_shared[["H"]][["HNEE"]]
  n = joint_distcal_shared[["n"]]
  
  # Gamma Ej perp
  GammaEjPerp = rbind(betaEjperp, betaEperp)
  
  # matrix P depends on ej part
  p1 = HEE %*% sigmaTT %*% ej/c(t(ej) %*% sigmaTT %*% ej)
  p4 = (sigmaST %*% ej + HNEE %*% sigmaTT %*% ej)/c(t(ej) %*% sigmaTT %*% ej)
  P = cbind(rbind(p1 * sqrt(n), p4 * sqrt(n)), P_shared)
  
  
  return(list(GammaEjPerp = GammaEjPerp,
              p1 = p1* sqrt(n),
              p4 = p4* sqrt(n),
              P = P))
  # Output:
  # GammaEjPerp: the GammaEjPerp matrix
  # p1: the 1st block of P matrix. This value will be used to calculate eta
  # p4: the 4st block of P matrix. This value will be used to calculate eta
  # P: the P matrix
}

