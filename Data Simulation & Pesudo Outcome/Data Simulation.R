# Data simulation function
## generate_individual: simulates data for an individual participant
## generate_dataset: simulates data for all participants

generate_individual <- function(id=1, T, P, sigma_residual, sigma_randint, main_rand = 0.5, rho, beta_logit, model, beta, theta1) {
  # id: the assigned id for this specific individual
  # T: the number of observations for this individual
  # P: the number of measurements or predictors
  # sigma_residual: the noise added to the moderator formula when generate pesudo outcome. This noise is 
  #                 generated using ARIMA
  # sigma_randint: the noise added to the moderator formula when generate pesudo outcome.
  # main_rand: the noise added to the control formula when generate pesudo outcome. This noise is generated using ARIMA
  # rho: this value is used to simulate state or predictors. We also use ARIMA when generating predictors, and
  #      rho defines how much of the past values influence the current value in an autoregressive model
  # beta_logit: the true value of beta when simulate the probability of assign treatment p(At = 1|Ht)
  # model: the true model for moderator when generate pesudo outcomes
  # beta: the true coefficients for moderator
  # theta1: the coefficients for the working model of control variables. For now, we just simply assume the working model
  #       contains the all predictors with coefficients equal to theta1
  
  ## Input: id, T, P (Dimension of state matrix), sigma_residual, sigma_randint, main_rand, rho (Correlation between time points), beta_logit, model, beta, theta1
  ## Function: Generates an individual with T time points, including state (S_t), action (A_t), and outcome (Y_{t+1}) using ARIMA for states and logistic regression for actions.
  ## Output: A data frame containing the generated states S_t, action A_t, outcomes Y_{t+1}, treatment effects, and noise components for a single participant.
  
  ## Order of generation: 
  ## (1) generate state S_t: VAR only depend on state
  ## (2) generate action A_t: depend on A_{t-1} and S_t
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
    if(t == 1) {current_action = 0} else{current_action = all_actions[t-1]}
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
  # N: the number of participants you hope to simulate
  # T: how many observations per participant
  # P: how many predictors to simulate
  # sigma_residual: the noise added to the moderator formula when generate pesudo outcome. This noise is 
  #                 generated using ARIMA
  # sigma_randint: the noise added to the moderator formula when generate pesudo outcome.
  # main_rand: the noise added to the control formula when generate pesudo outcome. This noise is generated using ARIMA
  # rho: this value is used to simulate state or predictors. We also use ARIMA when generating predictors, and
  #      rho defines how much of the past values influence the current value in an autoregressive model
  # beta_logit: the true value of beta when simulate the probability of assign treatment p(At = 1|Ht)
  # model: the true model for moderator when generate pesudo outcomes
  # beta: the true coefficients for moderator
  # theta1: the coefficients for the working model of control variables. For now, we just simply assume the working model
  #       contains the all predictors with coefficients equal to theta1
  
  
  MRT_data = rep(0,0)
  for(n in 1:N) {
    fake_individual = generate_individual(n, T, P, sigma_residual, sigma_randint, main_rand, rho, beta_logit, model, beta, theta1)  
    MRT_data = rbind(MRT_data,fake_individual)
  }
  return(MRT_data)
  # Output:
  # ID (1:N), Decision_point (1:T), state1 - stateP, prob, action, outcome (observed Y), signal (treatment_effect), noise
  # A large longitudinal dataset that contain records for N participants
}






