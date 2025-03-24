# This function calculate all conditional distribution for betaEj, GammaPerp which will be used 
# to built later pivot

conditional_dist = function(PQR_shared, PQR_ej, joint_distcal_shared, joint_distcal_ej, select_E) {
  # PQR_shared: output of function PQR_Pint_shared
  # PQR_ej: output of function PQR_Pint_ej
  # joint_distcal_shared: output of function joint_dist_Penal_Int_shared
  # joint_distcal_ej: output of function joint_dist_Penal_Int_ej
  # select_E: output of function variable_selection_PY_penal_int
  
  OMEGA = select_E[["OMEGA"]] # need to change it when get randomized lasso down
  pnum = PQR_shared[["pnum"]]
  Enum = PQR_shared[["Enum"]]
  P = PQR_ej[["P"]]
  Q = PQR_shared[["Q"]]
  R = PQR_shared[["R"]]
  GammaEjPerp = PQR_ej[["GammaEjPerp"]]
  betaEj = joint_distcal_ej[["betaEj_cal"]][["betaEj"]]
  HEE = joint_distcal_shared[["H"]][["HEE"]]
  HNEE = joint_distcal_shared[["H"]][["HNEE"]]
  lam = select_E[["lam"]] # need to change it when get randomized lasso down
  n = joint_distcal_shared[["n"]]
  p1 = PQR_ej[["p1"]]
  p4 = PQR_ej[["p4"]]
  ZNE = matrix(select_E[["Z"]], ncol = 1)
  hat_betaE_lambda = select_E[["soln"]]
  
  # create omega
  if(length(OMEGA) == 1) {
    omega = diag(rep(OMEGA, pnum)) 
  } else {
    omega = OMEGA
  }
  
  # conditional distribution of hat{beta}^{lambda}_E, Z_{-E}
  delta = -solve(t(Q) %*% solve(omega) %*% Q) %*% t(Q) %*% solve(omega) %*% (R + P %*% rbind(betaEj, GammaEjPerp))
  theta = solve(t(Q) %*% solve(omega) %*% Q)
  
  # conditional distribution of Z_{-E}
  delta2 = matrix(delta[(Enum+1):pnum,], ncol = 1)
  theta22 = theta[(Enum+1):pnum, (Enum+1):pnum]
  
  # conditional distribution of hat{beta}^{lambda}_E
  HE = rbind(HEE, HNEE)
  Zero = matrix(0,nrow = Enum, ncol = pnum - Enum)
  lammatrix = diag(rep(lam, pnum - Enum))
  mu = solve(t(HE) %*% solve(omega) %*% HE) %*% (t(HE) %*% solve(omega)) %*% 
    (P %*% rbind(betaEj, GammaEjPerp) + rbind(Zero, lammatrix) %*% ZNE + R) /c(sqrt(n))
  LAMBDA = solve(c(n)*t(HE) %*% solve(omega) %*% HE)
  
  # calculate values related to pivot
  eta = t(HE) %*% solve(omega) %*% rbind(p1, p4) * sqrt(n)
  Qn = LAMBDA %*% eta/ c(t(eta) %*% LAMBDA %*% eta)
  
  hat_betaE_lambda = hat_betaE_lambda[which(hat_betaE_lambda != 0)]
  # if(length(hat_betaE_lambda) < Enum) {hat_betaE_lambda = matrix(c(0.001,hat_betaE_lambda), ncol = 1)} 
  # WARNING: hat{beta}^lambda_E doesn't provide estimate for intercept. Here I use 0.001 for default
  # try to use fitted intercept from gee
  # if(length(hat_betaE_lambda) < Enum) {hat_betaE_lambda = matrix(c(joint_distcal[["betaEM"]][["coefficients"]][["(Intercept)"]]
  #                                                                 ,hat_betaE_lambda), ncol = 1)} 
  
  Aeta = hat_betaE_lambda - LAMBDA %*% eta %*% t(eta) %*% hat_betaE_lambda/c(t(eta) %*% LAMBDA %*% eta)
  
  return(list(omega = omega,
              delta = delta,
              theta = theta,
              delta2 = delta2,
              theta22 = theta22,
              HE = HE,
              subgradient_unselected = ZNE,
              mu = mu,
              LAMBDA = LAMBDA,
              eta = eta,
              Aeta = Aeta,
              Qn = Qn,
              se = PQR_shared[["se"]],
              n = n,hat_betaE_lambda = hat_betaE_lambda))
  # Output:
  # omega: Matrix version of OMEGA. It the variance of added random noised.
  # delta: the mean of conditional distribution (hat{beta}^{lambda}_E, Z_{-E}) | (hat{beta}_Ej, hat{Gamma}EjPerp).
  # theta: the variance of conditional distribution (hat{beta}^{lambda}_E, Z_{-E}) | (hat{beta}_Ej, hat{Gamma}EjPerp).
  # delta2: the mean of conditional distribution Z_{-E} | (hat{beta}_Ej, hat{Gamma}EjPerp).
  # theta22: the variance of conditional distribution Z_{-E} | (hat{beta}_Ej, hat{Gamma}EjPerp).
  # HE: the combined vector of HEE and HNEE.
  # subgradient_unselected: the subgradient vector of unselected variables 
  # mu: the mean of conditional distribution hat{beta}^{lambda}_E | (hat{beta}_Ej, hat{Gamma}EjPerp, Z_{-E}).
  # LAMBDA: the variance of conditional distribution hat{beta}^{lambda}_E | (hat{beta}_Ej, hat{Gamma}EjPerp, Z_{-E}).
  # eta: the specifically designed vector that will be used to reduce integration dimension when calculate pivot
  # se: Signs of the estimated coefficients for selected variables.
  # n: # of unique subjects in the dataset.
}