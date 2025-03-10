# Run random forest to train working model
pesudo_outcome_generator_rf = function(fold, ID, data, Ht, St, At, outcome) {
  # fold: # of folds hope to split
  # ID: the name of column where participants' ID are stored
  # data: simulated dataset
  # Ht: a vector that contains column names of Ht variables
  # St: a vector that contains column names of St variables; St should be a subset of Ht
  # At: column names of treatment (At)
  # outcome: column names of outcome variable
  
  fold_ind = split_data(data[,ID], fold = fold)
  MRT_rf = ps_random_forest(fold_indices = fold_ind, fold = fold, ID = ID, 
                            data = data, Ht = Ht, St = St, At = At, outcome = outcome)
  pesudo = pesudo_outcome_cal_rf(MRT_rf)
  return(pesudo)
}

## Function split simulated data into K folds
split_data = function(id, fold) {
  # id: the id column in data set
  # fold: # of folds hope to split
  
  uniID <- unique(id)
  numPat <- length(uniID)
  fsize <- floor(numPat/fold)
  permute_ID <- sample(uniID, size = numPat, replace = F) 
  fsize_exact <- c(rep(fsize, fold - 1), numPat-(fold - 1)*fsize) 
  ### incorporating the reaminder if numPat/fold is not whole num
  fold_labels <- rep(1:fold, fsize_exact)
  ### split the shuffled ID accroding to the folds
  fold_indices <- split(permute_ID, fold_labels)
  return(fold_indices)
}

ps_random_forest = function(fold_indices, fold, ID, data, Ht, St, At, outcome) {
  # fold_indices: the result of function fold_indMRTSim
  # fold: # of folds hope to split
  # ID: the name of column where participants' ID are stored
  # data: simulated dataset
  # Ht: a vector that contains column names of Ht variables
  # St: a vector that contains column names of St variables; St should be a subset of Ht
  # At: column names of treatment (At)
  # outcome: column names of outcome variable
  
  require(randomForest)
  require(tidyverse)
  
  data_withpred = data.frame()
  
  for(i in 1:fold) {
    print(i)
    reserve = data[(data[,ID] %in% fold_indices[[i]]), ]
    train_fold = data[!(data[,ID] %in% fold_indices[[i]]), ]
    
    
    ## gt(Ht,At)
    X_test = reserve[, c(Ht, At)]
    X_train = train_fold[, c(Ht, At)]
    y_train = train_fold[, outcome]
    
    gt_rf = randomForest(x = X_train, y = y_train)
    reserve$gt_pred_rf = predict(gt_rf, newdata = X_test)
    
    ## gt(Ht, At = 1)
    X_testfixAt1 = X_test
    X_testfixAt1[,At] = 1
    reserve$gtAt1_pred_rf = predict(gt_rf, newdata = X_testfixAt1)
    
    ## gt(Ht, At = 0)
    X_testfixAt0 = X_test
    X_testfixAt0[,At] = 0
    reserve$gtAt0_pred_rf = predict(gt_rf, newdata = X_testfixAt0)
    
    ## E[At|Ht]
    X_testHt = X_test[,Ht] #remove true action
    X_trainHt = X_train[,Ht]
    y_trainAt = as.factor(train_fold[,At])
    ptHt_rf = randomForest(x = X_trainHt, y = y_trainAt)
    reserve$ptHt_pred_rf = predict(ptHt_rf, newdata = X_testHt, type = "prob")[, 2]
    
    ## E[At|St]
    X_testSt = X_testHt[,St]
    X_trainSt = X_trainHt[,St]
    
    ptSt_rf = randomForest(x = X_trainSt, y = y_trainAt)
    reserve$ptSt_pred_rf = predict(ptSt_rf, newdata = X_testSt, type = "prob")[, 2]
    
    data_withpred = rbind(data_withpred, reserve)
    print(data_withpred)
  }
  colnames(data_withpred)[ncol(data_withpred)+1 - 5:1] = c("gt_pred_rf","gtAt1_pred_rf","gtAt0_pred_rf", "ptHt", "ptSt")
  data_withpred$ptHtobs = data_withpred[,At]*data_withpred$ptHt + (1-data_withpred[,At])*(1-data_withpred$ptHt)
  data_withpred$ptStobs = data_withpred[,At]*data_withpred$ptSt + (1-data_withpred[,At])*(1-data_withpred$ptSt)
  return(data_withpred)
}

pesudo_outcome_cal_rf = function(data_withpred) {
  At = data_withpred$action
  ptSt = data_withpred$ptSt
  y = data_withpred$outcome
  gt = data_withpred$gt_pred_rf
  gtAt1 = data_withpred$gtAt1_pred_rf
  gtAt0 = data_withpred$gtAt0_pred_rf
  
  # wt
  # wt = data_withpred$ptStobs/data_withpred$ptHtobs
  wt = data_withpred$ptStobs/(data_withpred$prob*At + (1-data_withpred$prob)* (1-At))
  
  data_withpred$yDR_rf = wt*(At - ptSt)*(y - gt)/(ptSt * (1-ptSt)) + (gtAt1 - gtAt0)
  print(data_withpred$yDR_rf)
  return(data_withpred)
}
