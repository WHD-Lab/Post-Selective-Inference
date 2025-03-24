# Run random forest to train working model
pesudo_outcome_generator_rf_v2 = function(fold, ID, data, Ht, St, At, outcome, core_num = NULL) {
  # fold: # of folds hope to split
  # ID: the name of column where participants' ID are stored
  # data: simulated dataset
  # Ht: a vector that contains column names of Ht variables
  # St: a vector that contains column names of St variables; St should be a subset of Ht
  # At: column names of treatment (At)
  # outcome: column names of outcome variable
  # core_num: number of cores will be used for calculation
  
  # Output
  # This function returns a dataset with pseudo outcome. It will call function split_data to generate folds.
  # Then call function ps_random_forest_v2 to train Random Forest on these folds to obtain estimates that will be used for 
  # pesudo outcome calculation. Then the last function will be called is pesudo_outcome_cal_rf_v2, and a 
  # column with name "yDR" will be generated.
  
  fold_ind = split_data(data[,ID], fold = fold)
  MRT_rf = ps_random_forest_v2(fold_indices = fold_ind, fold = fold, ID = ID, 
                               data = data, Ht = Ht, St = St, At = At, outcome = outcome, core_num)
  pesudo = pesudo_outcome_cal_rf_v2(MRT_rf)
  return(pesudo)
}

ps_random_forest_v2 = function(fold_indices, fold, ID, data, Ht, St, At, outcome, core_num = NULL) {
  # fold_indices: the result of function fold_indMRTSim
  # fold: # of folds hope to split
  # ID: the name of column where participants' ID are stored
  # data: simulated dataset
  # Ht: a vector that contains column names of Ht variables
  # St: a vector that contains column names of St variables; St should be a subset of Ht
  # At: column names of treatment (At)
  # outcome: column names of outcome variable
  # core_num: number of cores will be used for calculation
  
  # Output
  # this function use Random Forest to train model on folds and provide estimates that will be used for later 
  # pseudo outcome calculation
  
  
  data_withpred = data.frame()
  
  expectation_cal = function(i) {
    require(ranger)
    require(tidyverse)
    
    # print(i)
    reserve = data[(data[,ID] %in% fold_indices[[i]]), ]
    train_fold = data[!(data[,ID] %in% fold_indices[[i]]), ]
    
    ## gt(Ht,At)
    X_test = reserve[, c(Ht, At)]
    gt_rf = ranger(formula = as.formula(paste(outcome, "~", paste(c(Ht, At), collapse = " + "))), 
                   data = train_fold, 
                   num.trees = 500)
    reserve$gt_pred_rf = predict(gt_rf, data = reserve)$predictions
    
    ## gt(Ht, At = 1)
    X_testfixAt1 = X_test
    X_testfixAt1[,At] = 1
    reserve$gtAt1_pred_rf = predict(gt_rf, data = X_testfixAt1)$predictions
    
    ## gt(Ht, At = 0)
    X_testfixAt0 = X_test
    X_testfixAt0[,At] = 0
    reserve$gtAt0_pred_rf = predict(gt_rf, data = X_testfixAt0)$predictions
    
    ## E[At|Ht]
    ptHt_rf = ranger(formula = as.formula(paste(At, "~", paste(Ht, collapse = " + "))), 
                     data = train_fold, 
                     num.trees = 500, 
                     probability = TRUE)
    reserve$ptHt_pred_rf = predict(ptHt_rf, data = reserve)$predictions[, 2]  
    
    ## E[At|St]
    ptSt_rf = ranger(formula = as.formula(paste(At, "~", paste(St, collapse = " + "))), 
                     data = train_fold, 
                     num.trees = 500, 
                     probability = TRUE)
    
    reserve$ptSt_pred_rf = predict(ptSt_rf, data = reserve)$predictions[, 2]
    
    return(reserve)
  }
  
  ################## do it parallel
  require(parallel)
  
  folds_list = 1:fold
  
  if(is.null(core_num)) {cl = makeCluster(core_num)} else {cl = makeCluster(detectCores())}
  
  var_names = ls(envir= environment())
  
  clusterExport(cl, varlist = c(var_names, "expit"), envir = environment())
  
  results = parLapply(cl, folds_list, expectation_cal)
  
  stopCluster(cl)
  
  for(i in folds_list) {
    data_withpred = rbind(data_withpred, results[[i]])
  }
  ###########################
  
  colnames(data_withpred)[ncol(data_withpred)+1 - 5:1] = c("gt_pred_rf","gtAt1_pred_rf","gtAt0_pred_rf", "ptHt", "ptSt")
  data_withpred$ptHtobs = data_withpred[,At]*data_withpred$ptHt + (1-data_withpred[,At])*(1-data_withpred$ptHt)
  data_withpred$ptStobs = data_withpred[,At]*data_withpred$ptSt + (1-data_withpred[,At])*(1-data_withpred$ptSt)
  return(data_withpred)
}

pesudo_outcome_cal_rf_v2 = function(data_withpred) {
  # Input:
  # the outcome of function ps_random_forest_v2
  
  # Output:
  # produce dataset that contains a new column "yDR". This is the pseudo outcome we hope to get
  
  At = data_withpred$action
  ptSt = data_withpred$ptSt
  y = data_withpred$outcome
  gt = data_withpred$gt_pred_rf
  gtAt1 = data_withpred$gtAt1_pred_rf
  gtAt0 = data_withpred$gtAt0_pred_rf
  
  # wt
  # wt = data_withpred$ptStobs/data_withpred$ptHtobs
  wt = data_withpred$ptStobs/(data_withpred$prob*At + (1-data_withpred$prob)* (1-At))
  
  data_withpred$yDR = wt*(At - ptSt)*(y - gt)/(ptSt * (1-ptSt)) + (gtAt1 - gtAt0)
  # print(data_withpred$yDR)
  return(data_withpred)
}
