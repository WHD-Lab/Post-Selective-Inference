pseudo_outcome_generator_gbm = function(fold, ID, data, Ht, St, At, outcome, core_num = NULL) {
  fold_ind = split_data(data[[ID]], fold = fold)
  MRT_gbm = ps_gradient_boosting(fold_indices = fold_ind, fold = fold, ID = ID, 
                                 data = data, Ht = Ht, St = St, At = At, outcome = outcome, core_num)
  pseudo = pseudo_outcome_cal_gbm(MRT_gbm)
  return(pseudo)
}


ps_gradient_boosting = function(fold_indices, fold, ID, data, Ht, St, At, outcome, core_num = NULL) {   
  data_withpred = data.frame()    
  
  expectation_cal = function(i) {     
    require(xgboost)     
    require(parallel)     
    require(dplyr)          
    
    reserve = data[data[[ID]] %in% fold_indices[[i]], ]     
    train_fold = data[!(data[[ID]] %in% fold_indices[[i]]), ]          
    
    train_matrix = xgb.DMatrix(data = as.matrix(train_fold[, c(Ht, At)]), label = train_fold[[outcome]])     
    test_matrix = xgb.DMatrix(data = as.matrix(reserve[, c(Ht, At)]))          
    
    params_reg = list(objective = "reg:squarederror", booster = "gbtree", eta = 0.1, max_depth = 6, nrounds = 500)       
    params_class = list(objective = "binary:logistic", booster = "gbtree", eta = 0.1, max_depth = 6, nrounds = 500)  
    
    ## gt(Ht,At)     
    gt_gbm = xgb.train(params = params_reg, data = train_matrix, nrounds = params_reg$nrounds)     
    reserve$gt_pred_gbm = predict(gt_gbm, test_matrix)          
    
    ## gt(Ht, At = 1) and gt(Ht, At = 0)     
    X_testfixAt1 = reserve[, c(Ht, At)]; X_testfixAt1[[At]] = 1     
    X_testfixAt0 = reserve[, c(Ht, At)]; X_testfixAt0[[At]] = 0          
    
    reserve$gtAt1_pred_gbm = predict(gt_gbm, xgb.DMatrix(as.matrix(X_testfixAt1)))     
    reserve$gtAt0_pred_gbm = predict(gt_gbm, xgb.DMatrix(as.matrix(X_testfixAt0)))          
    
    ## E[At|Ht] and E[At|St]    
    ptHt_gbm = xgb.train(params = params_class, data = xgb.DMatrix(as.matrix(train_fold[, Ht]), label = train_fold[[At]]), nrounds = params_class$nrounds)     
    ptSt_gbm = xgb.train(params = params_class, data = xgb.DMatrix(as.matrix(train_fold[, St]), label = train_fold[[At]]), nrounds = params_class$nrounds)          
    
    reserve$ptHt_pred_gbm = predict(ptHt_gbm, xgb.DMatrix(as.matrix(reserve[, Ht])))     
    reserve$ptSt_pred_gbm = predict(ptSt_gbm, xgb.DMatrix(as.matrix(reserve[, St])))      
          
    
    return(reserve)   
  }    
  
  folds_list = 1:fold   
  if (is.null(core_num)) { cl = parallel::makeCluster(detectCores()) } else { cl = parallel::makeCluster(core_num) }   
  clusterExport(cl, varlist = ls(envir= environment()), envir = environment())   
  results = parLapply(cl, folds_list, expectation_cal)   
  stopCluster(cl)   
  data_withpred = dplyr::bind_rows(results)    
  
  colnames(data_withpred)[ncol(data_withpred)+1 - 5:1] = c("gt_pred_gbm","gtAt1_pred_gbm","gtAt0_pred_gbm", "ptHt", "ptSt")   
  data_withpred$ptHtobs = data_withpred[,At]*data_withpred$ptHt + (1-data_withpred[,At])*(1-data_withpred$ptHt)   
  data_withpred$ptStobs = data_withpred[,At]*data_withpred$ptSt + (1-data_withpred[,At])*(1-data_withpred$ptSt)    
  
  return(data_withpred) 
}

pseudo_outcome_cal_gbm = function(data_withpred) {
  At = data_withpred$action
  ptSt = data_withpred$ptSt
  y = data_withpred$outcome
  gt = data_withpred$gt_pred_gbm
  gtAt1 = data_withpred$gtAt1_pred_gbm
  gtAt0 = data_withpred$gtAt0_pred_gbm
  
  wt = data_withpred$ptStobs / (data_withpred$prob * At + (1 - data_withpred$prob) * (1 - At))
  data_withpred$yDR = wt * (At - ptSt) * (y - gt) / (ptSt * (1 - ptSt)) + (gtAt1 - gtAt0)
  
  return(data_withpred)
}
