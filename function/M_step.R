#function to fit transition-specific cox models
fit_cox_mstate <- function(data_mstate,
                           covar_tran,
                           tran_state,
                           tran_state_cure){
  
  data_split <- split(data_mstate, data_mstate[[ "fromto"]])
  
  cox_models <- lapply(seq_along(tran_state), function(i){
    
    tr <- tran_state[i]
    
    rhs <- covar_tran[[i]]
    
    if(tr %in% tran_state_cure){
      cure_var<-paste0("cure",tr)
      rhs <- c(rhs, cure_var)
    }
    
    formula_cox <- as.formula(
      paste("Surv(Tstart,Tstop,status) ~", paste(rhs, collapse = "+"))
    )
    
    coxph(
      formula_cox,
      weights = weight,
      data = data_split[[tr]],
      method = "breslow"
    )
    
  })
  
  names(cox_models) <- tran_state
  
  return(cox_models)
}

#function for M step in the EM algorithm
M_step <- function(data_longlong, data_logit, covar_tran,covar_cure,tran_state,tran_state_cure){
  
  data_mstate <- data_longlong[data_longlong$weight != 0,]
  fit_cox<- fit_cox_mstate(data_mstate,
                           covar_tran,
                           tran_state,
                           tran_state_cure)
  
  
  fit_alpha <- suppressWarnings(glm(
    reformulate(covar_cure, response = "cure"),
    weights = weight,
    data = data_logit,
    family = binomial(link="logit")
  ))
  
  list(
    fit_cox = fit_cox,
    fit_alpha = fit_alpha
  )
}
