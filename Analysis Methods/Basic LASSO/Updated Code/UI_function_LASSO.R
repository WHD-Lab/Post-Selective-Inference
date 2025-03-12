# function for basic LASSO DR-WCLS

DR_WCLS_LASSO = function(data, fold, ID, time, Ht, St, At, outcome, method_pesu, 
                                 lam = NULL, noise_scale = NULL, splitrat = 0.8, virtualenv_path, 
                                 beta = NULL, level = 0.9, core_num = NULL){
  # data: raw data without pesudo-outcome, ptSt.
  # fold: # of folds hope to split when generating pesudo outcome
  # ID: the name of column where participants' ID are stored
  # time: the name of column where time in study are stored
  # Ht: a vector that contains column names of Ht variables
  # St: a vector that contains column names of St variables; St should be a subset of Ht
  # At: column names of treatment (At)
  # outcome: column names of outcome variable
  # method_pesu: the machines learning method used when generate estimates of the nuisance parameters, and those values will be used to calculate the 
  #             pesudo outcome
  # lam: the value of lambda used to compute beta. If it is not provided, the default value will be used
  # noise_scale: Scale of Gaussian noise added to objective. Default is sqrt((1 - splitrat)/splitrat*NT)*sd(y).
  #             where omega is drawn from IID normals with standard deviation noise_scale
  # splitrat: the corresponding to the data splitting rate. Details can read "Exact Selective Inference with Randomization" page 15 equation (10).
  #           This value will be used only when user doesn't provide the noise_scale.
  # virtualenv_path: Python virtual environment path
  # beta: the true coefficient value (if simulation is conducted)
  # level: the CI significant level
  # core_num: the number of cores will be used for parallel calculation
  
  if(method_pesu == "CVLASSO") {
    ps = pesudo_outcome_generator_CVlasso(fold, ID, data, Ht, St, At, outcome, core_num)
  }
  
  if(method_pesu == "RandomForest") {
    ps = pesudo_outcome_generator_rf_v2(fold, ID, data, Ht, St, At, outcome, core_num)
  }
  
  my_formula = as.formula(paste("yDR ~ ", paste(St, collapse = " + ")))
  
  if(is.null(lam) & is.null(noise_scale)) {
    select = variable_selection_PY_penal_int(ps, ID, my_formula, splitrat=splitrat, virtualenv_path= virtualenv_path, beta = beta)
  }
  
  if(!is.null(lam) & is.null(noise_scale)) {
    select = variable_selection_PY_penal_int(ps, ID, my_formula, lam = lam, splitrat = splitrat, virtualenv_path= virtualenv_path, beta = beta)
  }
  
  if(is.null(lam) & !is.null(noise_scale)) {
    select = variable_selection_PY_penal_int(ps, ID, my_formula, noise_scale = noise_scale, splitrat = splitrat, virtualenv_path= virtualenv_path, beta = beta)
  }
  
  if(!is.null(lam) & !is.null(noise_scale)) {
    select = variable_selection_PY_penal_int(ps, ID, my_formula, lam = lam, noise_scale = noise_scale,
                                             splitrat = splitrat, virtualenv_path= virtualenv_path, beta = beta)
  }
   
  AsyNormbeta_shared = joint_dist_Penal_Int_shared(E = select$E, NE = select$NE, pes_outcome = "yDR", data = ps, id = ID, time = time)
  PQR_shared = PQR_Pint_shared(AsyNormbeta_shared, select)
  
  CI_per_select_var = function(ej) {
    AsyNormbeta_ej = joint_dist_Penal_Int_ej(AsyNormbeta_shared, ej)
    PQR_ej = PQR_Pint_ej(PQR_shared, AsyNormbeta_ej, AsyNormbeta_shared)
    condition = conditional_dist(PQR_shared, PQR_ej, AsyNormbeta_shared, AsyNormbeta_ej, select)
    CI = pivot_split_update(PQR_shared, PQR_ej, condition, AsyNormbeta_shared, AsyNormbeta_ej, select, level = level, pes_outcome = "yDR", data = ps,
                       id = ID, time = time)
    
    return(CI)
  }
  
  require(parallel)
  
  if(is.null(core_num)) {cl = makeCluster(core_num)} else {cl = makeCluster(detectCores())}
  
  # the list that stores the ej vectors for selected variables
  ejs = list()
  for(i in 1:length(select$E)) {
    vec_ej = rep(0,length(select$E))
    vec_ej[i] = 1
    
    ejs[[i]] = vec_ej
  }
  
  # List all objects in the global environment
  all_functions <- ls(envir = globalenv())
  
  # Filter the list to keep only the function names (not variables)
  function_names <- all_functions[sapply(all_functions, function(x) is.function(get(x)))]
  
  clusterExport(cl, varlist = function_names, envir = environment())
  # do parallel calculation
  results = parLapply(cl, ejs, CI_per_select_var)
  
  stopCluster(cl)
  
  # now makes the returned outcomes look like a table
  final_results = data.frame(
    E = character(),
    GEE_est = numeric(),
    post_beta = numeric(),
    pvalue = numeric(),
    lowCI = numeric(),
    upperCI = numeric(),
    prop_low = numeric(),
    prop_up = numeric()
  )
  
  for(i in 1:length(results)) {final_results[i,] = unlist(results[[i]])}
  
  return(final_results)
}




