# this file is built upon file "pivot_half_split_updateValue.R"
# the only difference between this version and previous version
# is in this version upper and lower bound search process is replaced by a function
# so you will not see repeated code

# This calculate the I-n, I+n. Details can be found at the beginning of page 9
support_update = function(se, eta, LAMBDA, Aeta, Qn) {
  # Input:
  # se: the sign solution of selected variables
  # eta: the updated eta value
  # LAMBDA: the updated variance of conditional distribution hat{beta}^{lambda}_E | (hat{beta}_Ej, hat{Gamma}EjPerp, Z_{-E}).
  # Aeta: the updated value of Aeta
  # Qn: the updated value of Qn
  
  # Output:
  # the lower and upper bounds for inner integration
  
  diagSe = diag(se)
  
  lowerbond = -Inf
  upperbond = Inf
  
  for(i in 1:length(se)) {
    
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
pivot_split_update = function(PQR_shared, PQR_ej, cond_dist, joint_distcal_shared, joint_distcal_ej, select_E, level = 0.9, pes_outcome, data,
                              id, time) {
  # PQR_shared: result of function PQR_Pint_shared
  # PQR_ej: result of function PQR_Pint_ej
  # cond_dist: result of function conditional_dist
  # joint_distcal_shared: result of function joint_dist_Penal_Int_shared
  # joint_distcal_ej: result of function joint_dist_Penal_Int_ej
  # select_E: result of function variable_selection_PY_penal_int
  # level: the significant level of confidence interval
  # pes_outcome: the column names of pseudo outcome
  # data: the dataset with pseudo outcome
  # id: column names of participant id
  # time: column names of decidion time points
  
  
  E = select_E[["E"]]
  NE = select_E[["NE"]]
  ej = joint_distcal_ej[["ej"]]
  HEE = joint_distcal_shared[["H"]][["HEE"]]
  HNEE = joint_distcal_shared[["H"]][["HNEE"]]
  HENE = joint_distcal_shared[["H"]][["HENE"]]
  S = joint_distcal_shared[["S"]]
  theta22 = cond_dist[["theta22"]]
  omega = cond_dist[["omega"]]
  ZNE = matrix(cond_dist[["subgradient_unselected"]], ncol = 1)
  n = cond_dist[["n"]]
  GammaEjPerp = PQR_ej[["GammaEjPerp"]]
  lam = select_E[["lam"]]
  R = PQR_shared[["R"]]
  LAMBDA = cond_dist[["LAMBDA"]]
  pnum = PQR_shared[["pnum"]]
  Enum = PQR_shared[["Enum"]]
  se = cond_dist[["se"]]
  hat_betaE_lambda = select_E[["soln"]]
  betaEjhat = joint_distcal_ej[["betaEj_cal"]][["betaEj"]]
  sigmasq1 = joint_distcal_ej[["betaEj_cal"]][["sigmasq1"]]
  
  
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
  
  pivot_prop_value = function(null_value) {
    # update joint beta distribution
    update_beta_dist = Asynorm_beta_under_hy(E, NE, pes_outcome, data, ej, id, time, null_value = null_value,
                                             HEE, HNEE, HENE, S)
    # update P matrix
    update_P = P_modify(update_beta_dist)
    p1 = update_P[["p1"]]
    p4 = update_P[["p4"]]
    P = update_P[["P"]]
    Pleft = P[,-1]
    
    # update mu_final and sigmasq_final
    update_mu_sigmasq_final = mu_sigmasq_final_modify(null_value = null_value, omega, P_modify = update_P,
                                                      PQR_org_shared = PQR_shared, PQR_org_ej = PQR_ej, theta22,
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
                            lower = min(mu_final, betaEjhat) - 15*sqrt(sigmasq_final),
                            upper = max(mu_final, betaEjhat) + 15*sqrt(sigmasq_final),
                            mu_final = mu_final, sigmasq_final, eta,p1,p4,Pleft,Var_etabeta,HE,
                            lower, upper)
    numerator = integrate(my_function,
                          lower = min(mu_final, betaEjhat) - 15*sqrt(sigmasq_final),
                          upper = betaEjhat,
                          mu_final = mu_final, sigmasq_final, eta,p1,p4,Pleft,Var_etabeta,HE,
                          lower, upper)
    prop_null_value = numerator$value/denominator$value
    
    return(prop_null_value)
  }
  
  # function that search betaE,j value using half split
  betaEj_search = function(range, max_tol, max_iterate) {
    # check the pivot value under the GEE point estimate
    prop_initial = pivot_prop_value(betaEjhat)
    
    #######################
    # find the lower bound#
    #######################
    
    if(prop_initial > upp & prop_initial < 1) {
      bottom = betaEjhat
      top = betaEjhat + range
      
      # check the pivot value under top
      top_prop = pivot_prop_value(top)
      while(top_prop > upp & top_prop < 1){
        range = range*1.25
        top = betaEjhat + range
        top_prop = pivot_prop_value(top)
      }
    } else {
      bottom = betaEjhat - range
      top = betaEjhat
      
      # check the pivot value under bottom
      bottom_prop = pivot_prop_value(bottom)
      while(bottom_prop < upp & !is.na(bottom_prop)) {
        range = range*1.25
        bottom = betaEjhat - range
        bottom_prop = pivot_prop_value(bottom)
      }
      
      while(bottom_prop > upp + 0.05 | is.na(bottom_prop)) {
        range = range/1.15
        bottom = betaEjhat - range
        bottom_prop = pivot_prop_value(bottom)
      }
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
      
      prop = pivot_prop_value(cal)
      i = i + 1
    }
    
    # the lower bound value
    low_bound = cal
    prop_up = prop
    
    #######################
    # find the upper bound#
    #######################
    ## reset top and bottom
    prop = 0
    if(prop_initial < lowp) {
      bottom = betaEjhat - range
      top = betaEjhat
      
      # check the pivot value under bottom
      bottom_prop = pivot_prop_value(bottom)
      while(bottom_prop < lowp) {
        range = range*1.25
        bottom = betaEjhat - range
        bottom_prop = pivot_prop_value(bottom)
      }
      
      
    } else {
      bottom = betaEjhat
      top = betaEjhat + range
      
      # check the pivot value under top
      top_prop = pivot_prop_value(top)
      while(top_prop > lowp & top_prop < 1 & !is.na(top_prop)) {
        range = range * 1.25
        top = betaEjhat + range
        top_prop = pivot_prop_value(top)
      }
      while(is.na(top_prop) | top_prop < lowp - 0.05) {
        range = range/1.15
        top = betaEjhat + range
        top_prop = pivot_prop_value(top)
      }
      
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
      prop = pivot_prop_value(cal)
      i = i + 1
    }
    
    # the upper bound value
    up_bound = cal
    prop_low = prop
    
    
    return(list(low_bound = low_bound,
                up_bound = up_bound,
                prop_up = prop_up,
                prop_low = prop_low))
  }
  
  CI = betaEj_search(sqrt(sigmasq1/n)*5, 10^{-5}, 10^{6})
  low = CI$low_bound
  up = CI$up_bound
  
  ############################
  # now calculate the p value#
  ############################
  
  prop0 = pivot_prop_value(0)
  if(prop0 > 1|is.na(prop0)) {pvalue = NA} else {pvalue = min(prop0, 1 - prop0)*2}
  
  return(list(E = E[which(ej == 1)],
              GEE_est = betaEjhat,
              post_beta = select_E[["postbeta"]][which(ej == 1)],
              pvalue = pvalue,
              lowCI = low,
              upperCI = up,
              prop_low = CI$prop_low,
              prop_up = CI$prop_up
  ))
  
  # Output:
  # E: the selected variables for which CI is calculated
  # GEE_est: the GEE estimate for this predictor
  # post_beta: the post selection true value for this predictor if simulation is conducted
  # pvalue: the p value
  # lowCI: the lower bound of the confidence interval
  # upperCI: the upper bound of the confidence interval
  # prop_low: the true corresponding pivot value for the lower bound. For example, if it's 90% CI, this 
  #          value will close to 0.05
  # prop_up: the true corresponding pivot value for the upper bound. For example, if it's 90% CI, this 
  #          value will close to 0.95
}

#######################################################
# I recalculate all previous values under this hypothesis value beta_E,j = b.

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
                                   PQR_org_shared, PQR_org_ej, theta22,
                                   Asynorm_beta_under_hy,ZNE) {
  # null_value: the value of betaE,j under the null hypothesis
  # Omega: the var-covariance of random variable, omega
  # P_modify: the result from function P_modify
  # PQR_org_shared: the original calculation for PQR matrix shared part
  # PQR_org_ej: the original calculation for PQR matrix will change with ej part
  # theta22: recycle from the previous calculation
  # Asynorm_beta_under_hy: the results from function Asynorm_beta_under_hy
  # ZNE: recycle from the previous calculation
  
  Q = PQR_org_shared[["Q"]]
  R = PQR_org_shared[["R"]]
  P = P_modify[["P"]]
  GammaEjPerp = PQR_org_ej[["GammaEjPerp"]]
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






















