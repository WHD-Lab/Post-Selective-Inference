# Generate Pesudo-outcome DR-WCLS version

################# Main Function #########################
pesudo_outcome_generator_CVlasso = function(fold, ID, data, Ht, St, At, outcome, core_num = NULL) {
  # fold: # of folds hope to split
  # ID: the name of column where participants' ID are stored
  # data: dataset name
  # Ht: a vector that contains column names of control variables
  # St: a vector that contains column names of moderator variables; St should be a subset of Ht
  # At: column names of treatment (At)
  # outcome: column names of outcome variable
  # core_num: number of cores will be used for calculation
  
  # Output
  # This function returns a dataset with pseudo outcome. It will call function split_data to generate folds.
  # Then call function simple_lasso to train CV LASSO on these folds to obtain estimates that will be used for 
  # pesudo outcome calculation. Then the last function will be called is pesudo_outcomecal, and a 
  # column with name "yDR" will be generated.
  
  fold_ind = split_data(data[,ID], fold = fold)
  MRT_sim_lasso = simple_lasso(fold_indices = fold_ind, fold = fold, ID = ID, data = data,
                               Ht = Ht, St = St, At = At, outcome = outcome, core_num)[[1]]
  pesudo = pesudo_outcomecal(MRT_sim_lasso)
  return(pesudo)
}


################# Supportive Function ###################

expit = function(x) {exp(x)/(1+exp(x))}

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

## test
# split_data(rep(1:30,40), fold = 8)
# split_data(rep(1:100, 100), fold = 10)

# Run simple LASSO to train working model
simple_lasso = function(fold_indices, fold, ID, data, Ht, St, At, outcome, core_num = NULL) {
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
  # this function use CV LASSO to train model on folds and provide estimates that will be used for later 
  # pseudo outcome calculation

  
  data_withpred <- data.frame()
  varselect_gt <- list()
  varselect_ptHt <- list()
  varselect_ptSt <- list()
  
  expectation_cal = function(i) {
    require(glmnet)
    require(tidyverse)
    
    reserve <- data[(data[,ID] %in% fold_indices[[i]]), ]
    train_fold <- data[!(data[,ID] %in% fold_indices[[i]]), ]  ## data used to train model
    
    ## gt(Ht,At)
    X_test <- reserve[, c(Ht, At)] # Add intercept
    X_train <- train_fold[, c(Ht, At)] # Add intercept
    y_train <- train_fold[, outcome]
    
    ### standardize data and store mean and var for y_train to convert prediction back to normal scale
    X_test = scale(as.matrix(X_test))
    X_train = scale(as.matrix(X_train))
    ymean = mean(y_train)
    ysd = sd(y_train)
    y_train = scale(y_train)
    
    # set.seed(1)
    cv.lasso.gt = cv.glmnet(as.matrix(X_train), y_train, alpha = 1)
    bestlambda.gt = cv.lasso.gt$lambda.min
    lasso.gt.mod = glmnet(as.matrix(X_train), y_train, alpha = 1, lambda = bestlambda.gt)
    varselect_gt[[i]] = lapply(coef(lasso.gt.mod), is.numeric) %>% unlist()
    reserve$lasso.pred.gt = predict(lasso.gt.mod, newx = as.matrix(X_test))*ysd + ymean
    # Add predicted gt(Ht, At = 1)
    X_testfixAt1 = X_test
    X_testfixAt1[,At] = 1
    reserve$lasso.pred_gtAt1 = predict(lasso.gt.mod, newx = as.matrix(X_testfixAt1))*ysd + ymean
    # Add predicted gt(Ht, At = 0)
    X_testfixAt0 = X_test
    X_testfixAt0[,At] = 0
    reserve$lasso.pred_gtAt0 = predict(lasso.gt.mod, newx = as.matrix(X_testfixAt0))*ysd + ymean
    
    
    ## E[At|Ht]
    X_testAt = X_test[,Ht] #remove true action
    X_trainAt = X_train[,Ht]
    y_testAt = reserve[,At]
    y_trainAt = train_fold[,At]
    
    # set.seed(1)
    cv.lasso.ptHt = cv.glmnet(as.matrix(X_trainAt), y_trainAt, alpha = 1, family = "binomial")
    bestlambda.ptHt = cv.lasso.ptHt$lambda.min
    lasso.ptHt.mod = glmnet(as.matrix(X_trainAt), y_trainAt, alpha = 1, family = "binomial", lambda = bestlambda.ptHt)
    varselect_ptHt[[i]] = lapply(coef(lasso.ptHt.mod), is.numeric) %>% unlist()
    reserve$lasso.pred.ptHt = expit(predict(lasso.ptHt.mod, newx = as.matrix(X_testAt)))
    
    ## E[At|St]
    X_testSt = X_testAt[,St]
    X_trainSt = X_trainAt[,St]
    
    # set.seed(1)
    cv.lasso.ptSt = cv.glmnet(as.matrix(X_trainSt), y_trainAt, alpha = 1, family = "binomial")
    bestlambda.ptSt = cv.lasso.ptSt$lambda.min
    lasso.ptSt.mod = glmnet(as.matrix(X_trainSt), y_trainAt, alpha = 1, family = "binomial", lambda = bestlambda.ptSt)
    varselect_ptSt[[i]] = lapply(coef(lasso.ptSt.mod), is.numeric) %>% unlist()
    reserve$lasso.pred.ptSt = expit(predict(lasso.ptSt.mod, newx = as.matrix(X_testSt)))
    
    return(reserve)
  }
  
  # do it parallel
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
  
  ##################################################
  
  colnames(data_withpred)[ncol(data_withpred)+1 - 5:1] = c("lasso.pred.gt","lasso.pred.gtAt1","lasso.pred.gtAt0", "ptHt", "ptSt")
  # lasso predictions for ptHt and ptSt only provide probability of At = 1
  # now add ptHt and ptSt for observed At
  data_withpred$ptHtobs = data_withpred[,At]*data_withpred$ptHt + (1-data_withpred[,At])*(1-data_withpred$ptHt)
  data_withpred$ptStobs = data_withpred[,At]*data_withpred$ptSt + (1-data_withpred[,At])*(1-data_withpred$ptSt)
  
  
  return(list(data_withpred, varselect_gt, varselect_ptHt, varselect_ptSt))
  
}

pesudo_outcomecal = function(data_withpred) {
  # Input:
  # the outcome of function simple_lasso
  
  # Output:
  # produce dataset that contains a new column "yDR". This is the pseudo outcome we hope to get
  
  At = data_withpred$action
  ptSt = data_withpred$ptSt
  y = data_withpred$outcome
  gt = data_withpred$lasso.pred.gt
  gtAt1 = data_withpred$lasso.pred.gtAt1
  gtAt0 = data_withpred$lasso.pred.gtAt0
  
  # wt
  # wt = data_withpred$ptStobs/data_withpred$ptHtobs
  wt = data_withpred$ptStobs/(data_withpred$prob*At + (1-data_withpred$prob)* (1-At))
  data_withpred$yDR = wt*(At - ptSt)*(y - gt)/(ptSt * (1-ptSt)) + (gtAt1 - gtAt0)
  
  return(data_withpred)
}



