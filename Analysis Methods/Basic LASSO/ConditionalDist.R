# This function calculate all conditional distribution for betaEj, GammaPerp which will be used 
# to built later pivot

conditional_dist = function(PQR, joint_distcal, select_E) {
  # Input:
  # PQR: the result of function PQR
  # joint_distcal: the result of function joint_dist
  # select_E: the result of variable_selection (function variable_selection_PY)
  
  OMEGA = select_E[["OMEGA"]] # need to change it when get randomized lasso down
  pnum = PQR[["pnum"]]
  Enum = PQR[["Enum"]]
  P = PQR[["P"]]
  Q = PQR[["Q"]]
  R = PQR[["R"]]
  GammaEjPerp = PQR[["GammaEjPerp"]]
  betaEj = joint_distcal[["betaEj_cal"]][["betaEj"]]
  HEE = joint_distcal[["H"]][["HEE"]]
  HNEE = joint_distcal[["H"]][["HNEE"]]
  lam = select_E[["lam"]] # need to change it when get randomized lasso down
  n = joint_distcal[["n"]]
  p1 = PQR[["p1"]]
  p4 = PQR[["p4"]]
  ZNE = matrix(select_E[["Z"]], ncol = 1)
  hat_betaE_lambda = select_E[["soln"]]
  
  # create omega
  if(length(OMEGA) == 1) {
    omega = diag(rep(OMEGA, pnum)) 
  } else {
    omega = OMEGA
  }
  
  # conditional distribution of (hat{beta}^{lambda}_E, Z_{-E}) | (hat{beta}_Ej, hat{Gamma}EjPerp)
  delta = -solve(t(Q) %*% solve(omega) %*% Q) %*% t(Q) %*% solve(omega) %*% (R + P %*% rbind(betaEj, GammaEjPerp))
  theta = solve(t(Q) %*% solve(omega) %*% Q)
  
  # conditional distribution of Z_{-E} | (hat{beta}_Ej, hat{Gamma}EjPerp)
  delta2 = matrix(delta[(Enum+1):pnum,], ncol = 1)
  theta22 = theta[(Enum+1):pnum, (Enum+1):pnum]
  
  # conditional distribution of hat{beta}^{lambda}_E | (hat{beta}_Ej, hat{Gamma}EjPerp, Z_{-E})
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
              se = PQR[["se"]],
              n = n,hat_betaE_lambda = hat_betaE_lambda))
  # Output:
  # omega: Matrix version of OMEGA (output of function RandomLASSO).
  # delta: the mean of conditional distribution (hat{beta}^{lambda}_E, Z_{-E}) | (hat{beta}_Ej, hat{Gamma}EjPerp).
  # theta: the variance of conditional distribution (hat{beta}^{lambda}_E, Z_{-E}) | (hat{beta}_Ej, hat{Gamma}EjPerp).
  # delta2: the mean of conditional distribution Z_{-E} | (hat{beta}_Ej, hat{Gamma}EjPerp).
  # theta22: the variance of conditional distribution Z_{-E} | (hat{beta}_Ej, hat{Gamma}EjPerp).
  # HE: the combined vector of HEE and HNEE.
  # subgradient_unselected: the subgradient vector of unselected variables 
  # mu: the mean of conditional distribution hat{beta}^{lambda}_E | (hat{beta}_Ej, hat{Gamma}EjPerp, Z_{-E}).
  # LAMBDA: the variance of conditional distribution hat{beta}^{lambda}_E | (hat{beta}_Ej, hat{Gamma}EjPerp, Z_{-E}).
  # eta: t(HE) %*% solve(omega) %*% rbind(p1, p4) * sqrt(n) (?)
  # se: Signs of the estimated coefficients for selected variables.
  # n: # of unique subjects in the dataset.
  
}
