# This calculate the I-n, I+n. Details can be found at the beginning of page 9
support_update = function(se, eta, LAMBDA, Aeta, Qn) {
  
  diagSe = diag(se)
  
  lowerbond = -Inf
  upperbond = Inf
  
  for(i in 1:length(se)) {
    # warning: I manually set i starts from 2 to skip intercept
    # Because we lack se and point estimate for intercept from randomized lasso
    # to avoid manipulated value impact later result I skip it
    
    frac = (t(diagSe[,i]) %*% Aeta)/(-t(diagSe[,i]) %*% Qn)
    
    if(t(diagSe[,i]) %*% LAMBDA %*% eta > 0) {
      lowerbond = max(frac, lowerbond)
    }
    
    if(t(diagSe[,i]) %*% LAMBDA %*% eta < 0) {
      upperbond = min(frac, upperbond)
    }
  }
  # what if exactly equal to 0
  
  return(list(upperbond = upperbond,
              lowerbond = lowerbond))
}

#####################################################################################
# when intercept is penalized
pivot_split_update = function(PQR, cond_dist, joint_distcal, select_E, level = 0.9,pes_outcome, data,
                      id, time) {
  # temp_support: result of function support
  # PQR: result of function PQR
  # cond_dist: result of function conditional_dist
  # joint_distcal: result of function joint_dist
  # level: the significant level
  
  E = select_E[["E"]]
  NE = select_E[["NE"]]
  ej = joint_distcal[["ej"]]
  HEE = joint_distcal[["H"]][["HEE"]]
  HNEE = joint_distcal[["H"]][["HNEE"]]
  HENE = joint_distcal[["H"]][["HENE"]]
  S = joint_distcal[["S"]]
  theta22 = cond_dist[["theta22"]]
  omega = cond_dist[["omega"]]
  ZNE = matrix(cond_dist[["subgradient_unselected"]], ncol = 1)
  n = cond_dist[["n"]]
  GammaEjPerp = PQR[["GammaEjPerp"]]
  lam = select_E[["lam"]]
  R = PQR[["R"]]
  LAMBDA = cond_dist[["LAMBDA"]]
  pnum = PQR[["pnum"]]
  Enum = PQR[["Enum"]]
  se = cond_dist[["se"]]
  hat_betaE_lambda = select_E[["soln"]]
  betaEj = joint_distcal[["betaEj_cal"]][["betaEj"]]
  sigmasq1 = joint_distcal[["betaEj_cal"]][["sigmasq1"]]

  
  lowp = (1-level)/2
  upp = 1 - (1-level)/2
 
  
  zeroEPNE = matrix(0, Enum, pnum-Enum)
  lamPNE = diag(rep(lam, pnum - Enum))
  
  
  # recreate self-define function to construct pivot
  my_function = function(b, mu_final,  sigmasq_final, eta,p1,p4,Pleft,Var_etabeta,HE,
                         lower, upper) {
    
    mu_etabeta = sapply(b, function(b) {t(eta) %*% solve(t(HE) %*% solve(omega) %*% HE) %*% (t(HE)%*% solve(omega)) %*% 
        (b*rbind(p1,p4) + Pleft %*% GammaEjPerp + rbind(zeroEPNE, lamPNE) %*% ZNE + R)/c(sqrt(n))})
    
    out = dnorm(b, mean = mu_final, sqrt(sigmasq_final)) * 
      (pnorm(upper,mu_etabeta, sqrt(Var_etabeta)) - pnorm(lower,mu_etabeta, sqrt(Var_etabeta))) 
    
    return(out)
  }
  
  # function that search betaE,j value using half split
  betaEj_search = function(initial, range, max_tol, max_iterate) {
    # initial: it is the GEE point estimate
    # range: from how far away we start the first search
    # max_tol: the maximum acceptable difference between target values and prop from pivot
    # max_iterate: the maximum iteration number
    
    ##########################################################
    # find the lower bound that p value will greater than upp#
    ##########################################################

    # update joint beta distribution
    update_beta_dist = Asynorm_beta_under_hy(E, NE, pes_outcome, data, ej, id, time, null_value = initial,
                                             HEE, HNEE, HENE, S)
    # update P matrix
    update_P = P_modify(update_beta_dist)
    p1 = update_P[["p1"]]
    p4 = update_P[["p4"]]
    P = update_P[["P"]]
    Pleft = P[,-1]
    
    # update mu_final and sigmasq_final
    update_mu_sigmasq_final = mu_sigmasq_final_modify(null_value = initial, omega, P_modify = update_P,
                                                      PQR_org = PQR, theta22,
                                                      Asynorm_beta_under_hy = update_beta_dist, ZNE)
    mu_final = update_mu_sigmasq_final[["mu_final"]]
    sigmasq_final = update_mu_sigmasq_final[["sigmasq_final"]]
    
    HE = rbind(HEE, HNEE)
    eta = c(sqrt(n)) * t(HE) %*% solve(omega) %*% rbind(p1, p4)
    Var_etabeta  = t(eta) %*% LAMBDA %*% eta
    Qn = LAMBDA %*% eta/ c(t(eta) %*% LAMBDA %*% eta)
    hat_betaE_lambda = hat_betaE_lambda[which(hat_betaE_lambda != 0)]
    Aeta = hat_betaE_lambda - LAMBDA %*% eta %*% t(eta) %*% hat_betaE_lambda/c(t(eta) %*% LAMBDA %*% eta)
    
    # calculate upper and lower bound
    temp_support = support_update(se = se, eta = eta, LAMBDA = LAMBDA, Aeta = Aeta, Qn = Qn)
    lower = temp_support[["lowerbond"]]
    upper = temp_support[["upperbond"]]
    
    denominator = integrate(my_function,
                            lower = min(mu_final, initial) - 15*sqrt(sigmasq_final),
                            upper = max(mu_final, initial) + 15*sqrt(sigmasq_final),
                            mu_final = mu_final, sigmasq_final, eta,p1,p4,Pleft,Var_etabeta,HE,
                            lower, upper)
    numerator = integrate(my_function,
                          lower = min(mu_final, initial) - 15*sqrt(sigmasq_final),
                          upper = initial,
                          mu_final = mu_final, sigmasq_final, eta,p1,p4,Pleft,Var_etabeta,HE,
                          lower, upper)
    prop_initial = numerator$value/denominator$value
    
    if(prop_initial > upp & prop_initial < 1) {
      bottom = initial
      top = initial + range
    } else {
      bottom = initial - range
      top = initial
    }
    
    i = 0
    prop = 0
    while((abs(upp - prop) > max_tol) & (i < max_iterate)) {
        if(i == 0) {
          cal = (top + bottom)/2
        } else {
          if(prop > upp) {
            bottom = cal
            cal = (bottom + top)/2
          } else {
            top = cal
            cal = (bottom + top)/2
          }
        }
          
      # update joint beta distribution
      update_beta_dist = Asynorm_beta_under_hy(E, NE, pes_outcome, data, ej, id, time, null_value = cal,
                                               HEE, HNEE, HENE, S)
      # update P matrix
      update_P = P_modify(update_beta_dist)
      p1 = update_P[["p1"]]
      p4 = update_P[["p4"]]
      P = update_P[["P"]]
      Pleft = P[,-1]
      
      # update mu_final and sigmasq_final
      update_mu_sigmasq_final = mu_sigmasq_final_modify(null_value = cal, omega, P_modify = update_P,
                                                        PQR_org = PQR, theta22,
                                                        Asynorm_beta_under_hy = update_beta_dist, ZNE)
      mu_final = update_mu_sigmasq_final[["mu_final"]]
      sigmasq_final = update_mu_sigmasq_final[["sigmasq_final"]]
      
      HE = rbind(HEE, HNEE)
      eta = c(sqrt(n)) * t(HE) %*% solve(omega) %*% rbind(p1, p4)
      Var_etabeta  = t(eta) %*% LAMBDA %*% eta
      Qn = LAMBDA %*% eta/ c(t(eta) %*% LAMBDA %*% eta)
      hat_betaE_lambda = hat_betaE_lambda[which(hat_betaE_lambda != 0)]
      Aeta = hat_betaE_lambda - LAMBDA %*% eta %*% t(eta) %*% hat_betaE_lambda/c(t(eta) %*% LAMBDA %*% eta)
      
      # calculate upper and lower bound
      temp_support = support_update(se = se, eta = eta, LAMBDA = LAMBDA, Aeta = Aeta, Qn = Qn)
      lower = temp_support[["lowerbond"]]
      upper = temp_support[["upperbond"]]
      
      denominator = integrate(my_function,
                              lower = min(mu_final, initial) - 15*sqrt(sigmasq_final),
                              upper = max(mu_final, initial) + 15*sqrt(sigmasq_final),
                              mu_final = mu_final, sigmasq_final, eta,p1,p4,Pleft,Var_etabeta,HE,
                              lower, upper)
      numerator = integrate(my_function,
                            lower = min(mu_final, initial) - 15*sqrt(sigmasq_final),
                            upper = initial,
                            mu_final = mu_final, sigmasq_final, eta,p1,p4,Pleft,Var_etabeta,HE,
                            lower, upper)
      prop = numerator$value/denominator$value
      
      i = i+ 1
    }
    # the lower bound value
    low_bound = cal
    prop_up = prop
    
    ########################################################
    # find the upper bound that p value will less than lowp#
    ########################################################
    # now find the upper bound for CI
    ## reset top and bottom
    prop = 0
    if(prop_initial < lowp) {
      bottom = initial - range
      top = initial
    } else {
      bottom = initial
      top = initial + range
    }
    
    i = 0
    while((abs(lowp - prop) > max_tol) & (i < max_iterate)) {
      if(i == 0) {
        cal = (top + bottom)/2
      } else {
        if(prop > lowp) {
          bottom = cal
          cal = (bottom + top)/2
        } else {
          top = cal
          cal = (bottom + top)/2
        }
      }
      
      # update joint beta distribution
      update_beta_dist = Asynorm_beta_under_hy(E, NE, pes_outcome, data, ej, id, time, null_value = cal,
                                               HEE, HNEE, HENE, S)
      # update P matrix
      update_P = P_modify(update_beta_dist)
      p1 = update_P[["p1"]]
      p4 = update_P[["p4"]]
      P = update_P[["P"]]
      Pleft = P[,-1]
      
      # update mu_final and sigmasq_final
      update_mu_sigmasq_final = mu_sigmasq_final_modify(null_value = cal, omega, P_modify = update_P,
                                                        PQR_org = PQR, theta22,
                                                        Asynorm_beta_under_hy = update_beta_dist, ZNE)
      mu_final = update_mu_sigmasq_final[["mu_final"]]
      sigmasq_final = update_mu_sigmasq_final[["sigmasq_final"]]
      
      HE = rbind(HEE, HNEE)
      eta = c(sqrt(n)) * t(HE) %*% solve(omega) %*% rbind(p1, p4)
      Var_etabeta  = t(eta) %*% LAMBDA %*% eta
      Qn = LAMBDA %*% eta/ c(t(eta) %*% LAMBDA %*% eta)
      hat_betaE_lambda = hat_betaE_lambda[which(hat_betaE_lambda != 0)]
      Aeta = hat_betaE_lambda - LAMBDA %*% eta %*% t(eta) %*% hat_betaE_lambda/c(t(eta) %*% LAMBDA %*% eta)
      
      # calculate upper and lower bound
      temp_support = support_update(se = se, eta = eta, LAMBDA = LAMBDA, Aeta = Aeta, Qn = Qn)
      lower = temp_support[["lowerbond"]]
      upper = temp_support[["upperbond"]]
      
      denominator = integrate(my_function,
                              lower = min(mu_final, initial) - 15*sqrt(sigmasq_final),
                              upper = max(mu_final, initial) + 15*sqrt(sigmasq_final),
                              mu_final = mu_final, sigmasq_final, eta,p1,p4,Pleft,Var_etabeta,HE,
                              lower, upper)
      numerator = integrate(my_function,
                            lower = min(mu_final, initial) - 15*sqrt(sigmasq_final),
                            upper = initial,
                            mu_final = mu_final, sigmasq_final, eta,p1,p4,Pleft,Var_etabeta,HE,
                            lower, upper)
      prop = numerator$value/denominator$value
      
      i = i+ 1
    }
    
    # the upper bound value
    up_bound = cal
    prop_low = prop
    
    return(list(low_bound = low_bound,
                up_bound = up_bound,
                prop_up = prop_up,
                prop_low = prop_low))
  }
  
  CI = betaEj_search(betaEj, sqrt(sigmasq1/n)*5, 10^{-5}, 10^{6})
  low = CI$low_bound
  up = CI$up_bound
  
  ############################
  # now calculate the p value#
  ############################
  
  # under the null hypothesis betaEj = 0
  
  # update joint beta distribution
  update_beta_dist = Asynorm_beta_under_hy(E, NE, pes_outcome, data, ej, id, time, null_value = 0,
                                           HEE, HNEE, HENE, S)
  # update P matrix
  update_P = P_modify(update_beta_dist)
  p1 = update_P[["p1"]]
  p4 = update_P[["p4"]]
  P = update_P[["P"]]
  Pleft = P[,-1]
  
  # update mu_final and sigmasq_final
  update_mu_sigmasq_final = mu_sigmasq_final_modify(null_value = 0, omega, P_modify = update_P,
                                                    PQR_org = PQR, theta22,
                                                    Asynorm_beta_under_hy = update_beta_dist, ZNE)
  mu_final = update_mu_sigmasq_final[["mu_final"]]
  sigmasq_final = update_mu_sigmasq_final[["sigmasq_final"]]
  
  HE = rbind(HEE, HNEE)
  eta = c(sqrt(n)) * t(HE) %*% solve(omega) %*% rbind(p1, p4)
  Var_etabeta  = t(eta) %*% LAMBDA %*% eta
  Qn = LAMBDA %*% eta/ c(t(eta) %*% LAMBDA %*% eta)
  hat_betaE_lambda = hat_betaE_lambda[which(hat_betaE_lambda != 0)]
  Aeta = hat_betaE_lambda - LAMBDA %*% eta %*% t(eta) %*% hat_betaE_lambda/c(t(eta) %*% LAMBDA %*% eta)
  
  # calculate upper and lower bound
  temp_support = support_update(se = se, eta = eta, LAMBDA = LAMBDA, Aeta = Aeta, Qn = Qn)
  lower = temp_support[["lowerbond"]]
  upper = temp_support[["upperbond"]]
  
  denominator = integrate(my_function,
                          lower = min(mu_final, betaEj) - 15*sqrt(sigmasq_final),
                          upper = max(mu_final, betaEj) + 15*sqrt(sigmasq_final),
                          mu_final = mu_final, sigmasq_final, eta,p1,p4,Pleft,Var_etabeta,HE,
                          lower, upper)
  numerator = integrate(my_function,
                        lower = min(mu_final, betaEj) - 15*sqrt(sigmasq_final),
                        upper = betaEj,
                        mu_final = mu_final, sigmasq_final, eta,p1,p4,Pleft,Var_etabeta,HE,
                        lower, upper)
  if(numerator$value > denominator$value) {
    p_value = NA
  } else {
    p_value = numerator$value/denominator$value
    p_value = min(1 - p_value, p_value)*2
  }
  
  ############################################################################
  # Under the true value for mu_final, changes hat_beta_E,j get pivot values #
  ############################################################################
  
  # update joint beta distribution
  #update_beta_dist = Asynorm_beta_under_hy(E, NE, pes_outcome, data, ej, id, time, null_value = select_E[["postbeta"]][which(ej == 1)],
  #                                         HEE, HNEE, HENE, S)
  # update P matrix
  #update_P = P_modify(update_beta_dist)
  #p1 = update_P[["p1"]]
  #p4 = update_P[["p4"]]
  #P = update_P[["P"]]
  #Pleft = P[,-1]
  
  # update mu_final and sigmasq_final
  #update_mu_sigmasq_final = mu_sigmasq_final_modify(null_value = select_E[["postbeta"]][which(ej == 1)], omega, P_modify = update_P,
  #                                                  PQR_org = PQR, theta22,
  #                                                  Asynorm_beta_under_hy = update_beta_dist, ZNE)
  #mu_final = update_mu_sigmasq_final[["mu_final"]]
  #sigmasq_final = update_mu_sigmasq_final[["sigmasq_final"]]
  
  #HE = rbind(HEE, HNEE)
  #eta = c(sqrt(n)) * t(HE) %*% solve(omega) %*% rbind(p1, p4)
  #Var_etabeta  = t(eta) %*% LAMBDA %*% eta
  #Qn = LAMBDA %*% eta/ c(t(eta) %*% LAMBDA %*% eta)
  #hat_betaE_lambda = hat_betaE_lambda[which(hat_betaE_lambda != 0)]
  #Aeta = hat_betaE_lambda - LAMBDA %*% eta %*% t(eta) %*% hat_betaE_lambda/c(t(eta) %*% LAMBDA %*% eta)
  
  # calculate upper and lower bound
  #temp_support = support_update(se = se, eta = eta, LAMBDA = LAMBDA, Aeta = Aeta, Qn = Qn)
  #lower = temp_support[["lowerbond"]]
  #upper = temp_support[["upperbond"]]
  
  #denominator = integrate(my_function,
  #                        lower = min(mu_final, betaEj) - 100*sqrt(sigmasq_final),
  #                        upper = max(mu_final, betaEj) + 100*sqrt(sigmasq_final),
  #                        mu_final = mu_final, sigmasq_final, eta,p1,p4,Pleft,Var_etabeta,HE,
  #                        lower, upper)
  
  # a sequence of made-up hat_beta_E,j values
  #fake_betaEj = seq(from= mu_final-100*sigmasq_final, to = mu_final+100*sigmasq_final, length.out = 200)
  #pivot_values = lapply(fake_betaEj, FUN = function(x) {
  #  numerator = integrate(my_function,
  #                        lower = min(mu_final, betaEj) - 100*sqrt(sigmasq_final),
  #                        upper = x,
  #                        mu_final = mu_final, sigmasq_final, eta,p1,p4,Pleft,Var_etabeta,HE,
  #                        lower, upper)
  #  return(numerator$value/denominator$value)
  #})
  
  
  return(list(E = E[which(ej == 1)],
              GEE_est = betaEj,
              post_beta = select_E[["postbeta"]][which(ej == 1)],
              pvalue = p_value,
              lowCI = low,
              upperCI = up,
              prop_low = CI$prop_low,
              prop_up = CI$prop_up
              #pivot_values = unlist(pivot_values)
              ))
  
}

#######################################################
# If I don't simply consider under hypothesis beta_E,j = b, will
# only change the way I calculate mu_final. Instead I recalculate
# all previous values under this hypothesis value.

# step 1. Get an estimate for the theoretical beta value under
# the hypothesis

Asynorm_beta_under_hy = function(E, NE, pes_outcome, data, ej, id, time, null_value,
                             HEE, HNEE, HENE, S) {
  # E: vector of selected predictors (don't include intercept)
  # NE: vector of unselected predictors
  # pes_outcome: column name for pesudo-outcome
  # data: the output of pesudo_outcomecal function
  # ej: jth standard basis vector. used to get betaE,j, don't need in matrix form
  # id: a vector which identifies individuals
  # time: a vector that records the decision points for each individual
  # null_value: the value of betaE,j under the null hypothesis
  # HEE: can recycle from previous calculation
  # HNEE: can recycle from previous calculation
  # HENE: can recycle from previous calculation
  # S: can recycle from previous calculation
  
  require(dplyr)
  
  id = data[,id]
  n = n_distinct(id)
  
  if("(Intercept)" %in% E) {ftStE = t(cbind(1,as.matrix(data[,E[E!="(Intercept)"]])))} else {ftStE = t(as.matrix(data[,E]))}
  if("(Intercept)" %in% NE) {ftStNE = t(cbind(1,as.matrix(data[,NE[NE!="(Intercept)"]])))} else {ftStNE = t(as.matrix(data[,NE]))}
  
  # modify the pesudo-outcome to remove the impact of betaEj
  data$yDR_modify = data[,pes_outcome] - t(ftStE) %*% matrix(ej, ncol = 1) * c(null_value)
  # modify E to remove betaE,j from selection
  E_modify = E[which(ej != 1)]
  if("(Intercept)" %in% E_modify) {
    formula = as.formula(paste("yDR_modify ~", paste(E_modify[which(E_modify != "(Intercept)")], collapse = "+")))
  } else {
    formula = as.formula(paste("yDR_modify~-1+", paste(E_modify, collapse = "+")))  
  }
  
  # get modified betaE point estiamtes
  wt = data$ptSt * (1-data$ptSt)
  idf = as.factor(id)
  time = data[,time]
  betaEM_modify = geepack::geeglm(formula, data = data, weights = wt/sqrt(n), corstr = "independence", id = idf,
                           waves = time)
  
  #sigmaTT
  KEE = K(ftStE, ftStNE, wt, "EE", betaEM_modify, n)
  sigmaTT = solve(HEE) %*% KEE %*% solve(HEE)
  
  #sigmaST
  KNEE = K(ftStE, ftStNE, wt, "-EE", betaEM_modify, n)
  sigmaST = KNEE %*% solve(HEE) - HNEE %*% sigmaTT
  
  #sigmaSS
  KNENE = K(ftStE, ftStNE, wt, "-E-E", betaEM_modify, n)
  sigmaSS = KNENE + HNEE %*% sigmaTT %*% HENE - KNEE %*% solve(HEE) %*% HENE - HNEE %*% solve(HEE) %*% t(KNEE)
  
  # point estimate betaE
  matrix_pre = diag(c(1-ej))
  matrix_pre = matrix_pre[,-which(ej == 1)]
  betaE_rest = matrix(betaEM_modify[["coefficients"]], ncol = 1)
  ej = matrix(ej, ncol = 1)
  betaE = ej*c(null_value) + matrix_pre %*% betaE_rest
  
  # betaEj
  betaEj = null_value
  sigmasq1 = t(ej) %*% sigmaTT %*% ej
  
  # betaEperp
  betaEperp = S/(c(sqrt(n))) - sigmaST %*% solve(sigmaTT) %*% betaE
  sigma2 = sigmaSS - sigmaST %*% solve(sigmaTT) %*% t(sigmaST)
  
  # betaEjperp
  betaEjperp = betaE - sigmaTT %*% ej %*% solve(t(ej) %*% sigmaTT %*% ej) * c(betaEj)
  sigma3 = sigmaTT - sigmaTT %*% ej %*% t(ej) %*% sigmaTT * c(solve(t(ej) %*% sigmaTT %*% ej))
  
  # return those updated value
  return(list(betaEM_modify = betaEM_modify,
              H = list(HEE = HEE, HNEE = HNEE, HENE = HENE),
              K_modify = list(KEE = KEE, KNEE = KNEE, KNENE = KNENE),
              Sigma_modify = list(sigmaTT = sigmaTT, sigmaST = sigmaST, sigmaSS = sigmaSS),
              betaEj_cal_modify = list(betaEj = betaEj, sigmasq1 = sigmasq1),
              betaEperp_cal_modify = list(betaEperp = betaEperp, sigma2 = sigma2),
              betaEjperp_cal_modify = list(betaEjperp = betaEjperp, sigma3 = sigma3),
              n = n,
              ej = ej,
              Enum = length(E),
              NEnum = length(NE)
              ))
}

# now update PQR calculation
# only P matrix will change with the new betaE,j value
# due to SigmaTT and SigmaST is changed
P_modify = function(Asynorm_beta_under_hy){
  # Asynorm_beta_under_hy: the result from function Asynorm_beta_under_hy
  
  HEE = Asynorm_beta_under_hy[["H"]][["HEE"]]
  HNEE = Asynorm_beta_under_hy[["H"]][["HNEE"]]
  sigmaTT = Asynorm_beta_under_hy[["Sigma_modify"]][["sigmaTT"]]
  sigmaST = Asynorm_beta_under_hy[["Sigma_modify"]][["sigmaST"]]
  n = Asynorm_beta_under_hy[["n"]]
  ej = Asynorm_beta_under_hy[["ej"]]
  Enum = Asynorm_beta_under_hy[["Enum"]]
  NEnum = Asynorm_beta_under_hy[["NEnum"]]
  
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
  
  return(list(P = P,
         p1 = p1 * sqrt(n),
         p4 = p4*sqrt(n), n = n))
}

# Update mu_final and sigmasq_final
mu_sigmasq_final_modify = function(null_value, omega, P_modify,
                                   PQR_org, theta22,
                                   Asynorm_beta_under_hy,ZNE) {
  # null_value: the value of betaE,j under the null hypothesis
  # Omega: the var-covariance of random variable, omega
  # P_modify: the result from function P_modify
  # PQR_org: the original calculation for PQR matrix
  # theta22: recycle from the previous calculation
  # Asynorm_beta_under_hy: the results from function Asynorm_beta_under_hy
  # ZNE: recycle from the previous calculation
  
  Q = PQR_org[["Q"]]
  R = PQR_org[["R"]]
  P = P_modify[["P"]]
  GammaEjPerp = PQR_org[["GammaEjPerp"]]
  n = P_modify[["n"]]
  sigma1sq = Asynorm_beta_under_hy[["betaEj_cal_modify"]][["sigmasq1"]]
  Enum = Asynorm_beta_under_hy[["Enum"]]
  NEnum = Asynorm_beta_under_hy[["NEnum"]]
  pnum = Enum + NEnum
  
  # update matrix A and C
  zeros = matrix(0, nrow = NEnum, ncol = Enum)
  I = diag(rep(1, NEnum))
  shared = cbind(zeros, I)
  zero1p = matrix(0, nrow = 1, ncol = pnum)
  Ip = diag(rep(1, pnum))
  A = -shared %*% solve(t(Q) %*% solve(omega) %*% Q) %*% t(Q) %*% solve(omega) %*% P %*% matrix(c(1,rep(0,pnum)), ncol = 1)
  C = -shared %*% solve(t(Q) %*% solve(omega) %*% Q) %*% t(Q) %*% solve(omega) %*% (R + P %*% rbind(zero1p, Ip) %*% GammaEjPerp)
  
  mu_final = (n*null_value + c(sigma1sq) * (t(ZNE - C) %*% solve(theta22) %*% A))/(n + c(sigma1sq) * (t(A) %*% solve(theta22) %*% A))
  sigmasq_final = sigma1sq / (n + c(sigma1sq) * (t(A) %*% solve(theta22) %*% A))
  
  return(list(mu_final = mu_final,
              sigmasq_final = sigmasq_final))
}









