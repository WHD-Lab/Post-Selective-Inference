# Data simulation function
## The first function generate_individual simulates data for an individual participant
## generate_dataset: simulates data for all participants

generate_individual <- function(id=1, T, P, sigma_residual, sigma_randint, main_rand = 0.5, rho, beta_logit, model, beta, theta1) {
  ## Generate an individual with T time points
  ## Using parameter inputs above
  ## Order of generation: 
  ## (1) generate state S_t: VAR only depend on state
  ## (2) generate action A_t: depend on A_{t-1} and S_t, 
  ## (3) generate outcome Y_{t+1}
  all_states = matrix(ncol = P, nrow = T)
  for(p in 1:P) {
    all_states[,p] = arima.sim(list(ar = rho, ma = 1), n = T)
  }
  all_actions = rep(0,length = T)
  all_probabilities = vector(length = T)
  current_action = 0
  for(t in 1:T) {
    current_states = all_states[t,]
    if(t ==1) {current_action = 0} else{current_action = all_actions[t-1]}
    prob_action = 1/(1+exp(-c(current_action, current_states)%*%beta_logit)) ## Calculate probability of action at current time
    current_action = rbinom(n=1, size=1, prob = prob_action) ## Choose current action at time t
    all_probabilities[t] = prob_action
    all_actions[t] = current_action
  }
  
  # Set state1, state2, state3, and intercept are the true S_t
  # And assume an easy linear function to start
  temp = data.frame(all_states)
  colnames(temp) = paste("state",1:P, sep="")
  X = data.matrix(modelr::model_matrix(temp, model))
  treatment_effect = X%*% beta ## E[ Y_{t+1} | S_t ]
  main_effect = theta1 * rowSums(all_states) ## E[ Y_{t+1} | H_t ] 
  
  meanY = main_effect + (all_actions-all_probabilities)*treatment_effect ## Combines Main and Txt 
  #errorY = arima.sim(list(ar = rho), n = T, sd = main_rand)
  errorY = jmuOutlier::rlaplace(n=T, sd = main_rand)
  #txterrorY = arima.sim(list(ar = rho), n = T, sd = sigma_residual) + rnorm(n=1, mean =0, sd = sigma_randint)
  txterrorY = jmuOutlier::rlaplace(n=T, sd = sigma_residual) +  jmuOutlier::rlaplace(n=1, sd = sigma_randint) 
  obsY = meanY + errorY + txterrorY*(all_actions - all_probabilities) # Add error back
  
  noise_total = txterrorY*(all_actions - all_probabilities) + errorY
  
  df_individual = data.frame(cbind(id, 1:T, all_states, all_probabilities, all_actions, obsY, treatment_effect, noise_total))
  colnames(df_individual)[P+5+1 : 2] = c("signal", "noise")
  colnames(df_individual)[2+ P+1:3] = c("prob", "action", "outcome")
  colnames(df_individual)[2 + 1:P] =  paste("state",1:P, sep="")
  colnames(df_individual)[1:2] = c("id", "decision_point")
  
  return(df_individual)
}

generate_dataset <- function(N, T, P, sigma_residual, sigma_randint, main_rand = 0.5, rho, beta_logit, model, beta, theta1) {
  MRT_data = rep(0,0)
  for(n in 1:N) {
    fake_individual = generate_individual(n, T, P, sigma_residual, sigma_randint, main_rand, rho, beta_logit, model, beta, theta1)  
    MRT_data = rbind(MRT_data,fake_individual)
  }
  return(MRT_data)
}





