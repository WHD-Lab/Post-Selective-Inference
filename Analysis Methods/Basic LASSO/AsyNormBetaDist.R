# This code calculate the joint multivariate normal distribution for betaEj, betaEperp, betaEjPerp
# for notation, refer the math manual Page 1 to 4

joint_dist = function(E, NE, pes_outcome, data, ej, id, time) {
  # E: vector of selected predictors (don't include intercept)
  # NE: vector of unselected predictors
  # pes_outcome: column name for pesudo-outcome
  # data: the output of pesudo_outcomecal function
  # ej: jth standard basis vector. used to get betaE,j, don't need in matrix form
  # id: a vector which identifies individuals
  # time: a vector that records the decision points for each individual
  
  require(dplyr)
  ftStE = t(cbind(1,as.matrix(data[,E]))) # include intercept as default
  ftStNE = t(as.matrix(data[,NE]))
  id = data[,id]
  n = n_distinct(id)
  
  # get betaE point estiamtes
  wt = data$ptSt * (1-data$ptSt)
  idf = as.factor(id)
  time = data[,time]
  formula = as.formula(paste(pes_outcome, "~", paste(E, collapse = "+")))
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
  
  # point estimate betaE
  betaE = matrix(betaEM[["coefficients"]], ncol = 1)
  
  # betaEj
  ej = matrix(ej, ncol = 1)
  betaEj = t(ej) %*% betaE
  sigmasq1 = t(ej) %*% sigmaTT %*% ej
  
  # betaEperp
  betaEperp = S/(c(sqrt(n))) - sigmaST %*% solve(sigmaTT) %*% betaE
  sigma2 = sigmaSS - sigmaST %*% solve(sigmaTT) %*% t(sigmaST)
  
  # betaEjperp
  betaEjperp = betaE - sigmaTT %*% ej %*% solve(t(ej) %*% sigmaTT %*% ej) * c(betaEj)
  sigma3 = sigmaTT - sigmaTT %*% ej %*% t(ej) %*% sigmaTT * c(solve(t(ej) %*% sigmaTT %*% ej))
  
  return(list(wt = wt,
              betaEM = betaEM,
              S = S,
              ej = ej,
              H = list(HEE = HEE, HNEE = HNEE, HENE = HENE),
              K = list(KEE = KEE, KNEE = KNEE, KNENE = KNENE),
              Sigma = list(sigmaTT = sigmaTT, sigmaST = sigmaST, sigmaSS = sigmaSS),
              betaEj_cal = list(betaEj = betaEj, sigmasq1 = sigmasq1),
              betaEperp_cal = list(betaEperp = betaEperp, sigma2 = sigma2),
              betaEjperp_cal = list(betaEjperp = betaEjperp, sigma3 = sigma3),
              n = n,
              E = E,
              NE = NE
  ))
}

#############################################################################
# Now intercept is allowed to be penalized
joint_dist_Penal_Int = function(E, NE, pes_outcome, data, ej, id, time) {
  # E: vector of selected predictors (don't include intercept)
  # NE: vector of unselected predictors
  # pes_outcome: column name for pesudo-outcome
  # data: the output of pesudo_outcomecal function
  # ej: jth standard basis vector. used to get betaE,j, don't need in matrix form
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
  
  # point estimate betaE
  betaE = matrix(betaEM[["coefficients"]], ncol = 1)
  
  # betaEj
  ej = matrix(ej, ncol = 1)
  betaEj = t(ej) %*% betaE
  sigmasq1 = t(ej) %*% sigmaTT %*% ej
  
  # betaEperp
  betaEperp = S/(c(sqrt(n))) - sigmaST %*% solve(sigmaTT) %*% betaE
  sigma2 = sigmaSS - sigmaST %*% solve(sigmaTT) %*% t(sigmaST)
  
  # betaEjperp
  betaEjperp = betaE - sigmaTT %*% ej %*% solve(t(ej) %*% sigmaTT %*% ej) * c(betaEj)
  sigma3 = sigmaTT - sigmaTT %*% ej %*% t(ej) %*% sigmaTT * c(solve(t(ej) %*% sigmaTT %*% ej))
  
  return(list(wt = wt,
              betaEM = betaEM,
              S = S,
              ej = ej,
              H = list(HEE = HEE, HNEE = HNEE, HENE = HENE),
              K = list(KEE = KEE, KNEE = KNEE, KNENE = KNENE),
              Sigma = list(sigmaTT = sigmaTT, sigmaST = sigmaST, sigmaSS = sigmaSS),
              betaEj_cal = list(betaEj = betaEj, sigmasq1 = sigmasq1),
              betaEperp_cal = list(betaEperp = betaEperp, sigma2 = sigma2),
              betaEjperp_cal = list(betaEjperp = betaEjperp, sigma3 = sigma3),
              n = n,
              E = E,
              NE = NE
  ))
}


#############################################################################
H = function(ftStE, ftStNE, wt, type, n) {
  # type": EE, -EE, E-E
  
  if(type == "EE") {return(t(t(ftStE) * c(wt)) %*% t(ftStE)/c(n))}
  if(type == "-EE") {return(t(t(ftStNE) * c(wt)) %*% t(ftStE)/c(n))}
  if(type == "E-E") {return(t(t(ftStE) * c(wt)) %*% t(ftStNE)/c(n))}
}

K = function(ftStE, ftStNE, wt, type, betaEM, n) {
  # type": EE, -EE, -E-E
  # betaEM: fitted model for selected predictors
  
  require(dplyr)
  
  res = betaEM[["residuals"]]
  wtres_prod = wt * res
  id = betaEM[["id"]]
  uniqueID = unique(id)
  
  
  if (type == "EE") {
    # initialize K
    K = matrix(0, nrow = dim(ftStE)[1], ncol = dim(ftStE)[1])
    for(i in uniqueID) {
      range = which(i == id)
      select_wtres_prod = wtres_prod[range]
      ftStE_select = ftStE[,range]
      innerprod = ftStE_select %*% matrix(select_wtres_prod, ncol = 1) 
      K = K + innerprod %*% t(innerprod)
    }
  }
  
  if (type == "-EE") {
    # initialize K
    K = matrix(0, nrow = dim(ftStNE)[1], ncol = dim(ftStE)[1])
    for(i in uniqueID) {
      range = which(i == id)
      select_wtres_prod = wtres_prod[range]
      ftStE_select = ftStE[,range]
      ftStNE_select = ftStNE[,range]
      inner1 = ftStNE_select %*% matrix(select_wtres_prod, ncol = 1)
      inner2 = ftStE_select %*% matrix(select_wtres_prod, ncol = 1)
      K = K + (inner1 %*% t(inner2))
    }
  }
  
  if (type == "-E-E") {
    # initialize K
    K = matrix(0, nrow = dim(ftStNE)[1], ncol = dim(ftStNE)[1])
    for(i in uniqueID) {
      range = which(i == id)
      select_wtres_prod = wtres_prod[range]
      ftStNE_select = ftStNE[,range]
      inner = ftStNE_select %*% matrix(select_wtres_prod, ncol = 1)
      K = K + inner %*% t(inner)
    }
  }
  
  return(K/n)
  
}


