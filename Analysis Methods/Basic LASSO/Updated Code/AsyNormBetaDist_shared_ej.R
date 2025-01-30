joint_dist_Penal_Int_shared = function(E, NE, pes_outcome, data, id, time) {
  # E: vector of selected predictors (don't include intercept)
  # NE: vector of unselected predictors
  # pes_outcome: column name for pesudo-outcome
  # data: the output of pesudo_outcomecal function
  # id: column name, a vector which identifies individuals
  # time: column name, a vector that records the decision points for each individual
  
  require(dplyr)
  if("(Intercept)" %in% E) {ftStE = t(cbind(1,as.matrix(data[,E[E!="(Intercept)"]])))} else {ftStE = t(as.matrix(data[,E]))}
  if("(Intercept)" %in% NE) {ftStNE = t(cbind(1,as.matrix(data[,NE[NE!="(Intercept)"]])))} else {ftStNE = t(as.matrix(data[,NE]))}
  
  id = data[,id]
  n = n_distinct(id)
  
  # get betaE point estiamtes
  wt = data$ptSt * (1-data$ptSt)
  idf = as.factor(id)
  time = data[,time]
  
  if("(Intercept)" %in% E) {
    formula = as.formula(paste(pes_outcome, "~", paste(E[which(E != "(Intercept)")], collapse = "+")))
  } else {
    formula = as.formula(paste(pes_outcome, "~-1+", paste(E, collapse = "+")))  
  }
  
  betaEM = geepack::geeglm(formula, data = data, weights = wt/sqrt(n), corstr = "independence", id = idf,
                           waves = time)
  
  # S
  res = betaEM[["residuals"]] # this is unweighted residuals
  S = ftStNE %*% matrix(res*wt, ncol = 1)/(c(sqrt(n)))
  
  #sigmaTT
  HEE = H(ftStE, ftStNE, wt, "EE", n)
  KEE = K(ftStE, ftStNE, wt, "EE", betaEM, n)
  sigmaTT = solve(HEE) %*% KEE %*% solve(HEE)
  
  #sigmaST
  KNEE = K(ftStE, ftStNE, wt, "-EE", betaEM, n)
  HNEE = H(ftStE, ftStNE, wt, "-EE", n)
  sigmaST = KNEE %*% solve(HEE) - HNEE %*% sigmaTT
  
  #sigmaSS
  KNENE = K(ftStE, ftStNE, wt, "-E-E", betaEM, n)
  HENE = H(ftStE, ftStNE, wt, "E-E", n)
  sigmaSS = KNENE + HNEE %*% sigmaTT %*% HENE - KNEE %*% solve(HEE) %*% HENE - HNEE %*% solve(HEE) %*% t(KNEE)

  
  # betaEperp
  betaEperp = S/(c(sqrt(n))) - sigmaST %*% solve(sigmaTT) %*% betaE
  sigma2 = sigmaSS - sigmaST %*% solve(sigmaTT) %*% t(sigmaST)
  

  return(list(wt = wt,
              betaEM = betaEM,
              S = S,
              H = list(HEE = HEE, HNEE = HNEE, HENE = HENE),
              K = list(KEE = KEE, KNEE = KNEE, KNENE = KNENE),
              Sigma = list(sigmaTT = sigmaTT, sigmaST = sigmaST, sigmaSS = sigmaSS),
              betaEperp_cal = list(betaEperp = betaEperp, sigma2 = sigma2),
              n = n,
              E = E,
              NE = NE
  ))
}

joint_dist_Penal_Int_ej = function(joint_dist_joint, ej) {
  # joint_dist_joint: the results of function joint_dist_Penal_Int_shared
  # ej: jth standard basis vector. used to get betaE,j, don't need in matrix form
  
  betaEM = joint_dist_joint[["betaEM"]]
  sigmaTT = joint_dist_joint[["Sigma"]][["sigmaTT"]]
  
  # point estimate betaE
  betaE = matrix(betaEM[["coefficients"]], ncol = 1)
  
  # betaEj
  ej = matrix(ej, ncol = 1)
  betaEj = t(ej) %*% betaE
  sigmasq1 = t(ej) %*% sigmaTT %*% ej
  
  # betaEjperp
  betaEjperp = betaE - sigmaTT %*% ej %*% solve(t(ej) %*% sigmaTT %*% ej) * c(betaEj)
  sigma3 = sigmaTT - sigmaTT %*% ej %*% t(ej) %*% sigmaTT * c(solve(t(ej) %*% sigmaTT %*% ej))
  
  return(list(ej = ej,
              betaEj_cal = list(betaEj = betaEj, sigmasq1 = sigmasq1),
              betaEjperp_cal = list(betaEjperp = betaEjperp, sigma3 = sigma3)))

}


