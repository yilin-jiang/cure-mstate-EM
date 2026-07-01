EM_mstate_cure <- function(  data_long0,
                             id,
                             covar_cure,
                             covar_tran,
                             tran_state,
                             tran_noncure,
                             tran_cure,
                             pi_init, # can be NULL when start is given
                             model0, # can be NULL when start is given
                             start=NULL,
                           epsilon_tol=1e-4,
                           max_iter=200,
                           t_cure=NA){
  
  start_time <- Sys.time()
  tran_cured<-paste0(tran_cure,".2")
  tran_state_cure<-tran_state[which(tran_noncure %in% tran_cure)]
  tran_correction<-tran_noncure[which(!tran_noncure %in% tran_cure)]
  
  if (!is.null(start)) {
    fit_cox   <- start$fit_cox
    fit_alpha <- start$fit_alpha
    data_long0$fromto<-paste0(data_long0$from,data_long0$to)
    data_longlong <- create_datalonglong(
      data_long = data_long0,
      transition = tran_state_cure,
      id = id
    )
    data_longlong<-update_pi(fit_alpha= fit_alpha,cov= covar_cure,data_longlong = data_longlong)
    data_longlong <- calculate_hazard_dt(
      model = fit_cox,
      tran = tran_state,
      data_long = data_longlong
    )
    data_longlong <-add_likelihood(data_longlong)
    data_longlong<-calculate_weight_mstate(data_longlong=data_longlong,
                                           id,
                                           tran_noncure,
                                           tran_cure=tran_cured,tran_correction,t_cure)
    
    data_logit <- create_datalogit(
      data_longlong,
      id = id,
      cov = covar_cure
    )
    } else {
  
  init <- initialize_em(
    model0,
    data_long0,
    pi_init,
    tran_state,
    tran_state_cure,
    tran_noncure,
    tran_cured,
    tran_correction,
    covar_cure,
    id,
    t_cure
  )
  
  data_longlong <- init$data_longlong
  data_logit <- init$data_logit
    }
  stepM <- M_step(data_longlong, data_logit, covar_tran,covar_cure,tran_state,tran_state_cure)
  
  fit_cox <- stepM$fit_cox
  fit_alpha <- stepM$fit_alpha
  
  
  loglik <- calculate_loglik_ob(data_longlong, id)
  
  iter <- 0
  epsilon <- 1
  loglik_trace <- loglik
 
  
  while(epsilon > epsilon_tol & iter < max_iter){
    
    stepE <- E_step(
      data_longlong,
      fit_alpha,
      fit_cox, 
      tran_state,
      tran_noncure,
      tran_cured,
      tran_correction,
      id,
      t_cure
    )
    data_longlong<-stepE$data_longlong
    data_logit<-stepE$data_logit
    
    stepM <- M_step(
      data_longlong,
      data_logit,
      covar_tran,
      covar_cure,
      tran_state,
      tran_state_cure
    )
    
    fit_cox <- stepM$fit_cox
    fit_alpha <- stepM$fit_alpha
    
    loglik_new <- calculate_loglik_ob(
      data_longlong,
      id
    )
    
    epsilon <- loglik_new - loglik
    loglik <- loglik_new
    
    iter <- iter + 1
    loglik_trace <- c(loglik_trace, loglik)
    
    cat("Iter:",iter,"LogLik:",loglik,"\n")
  }
  
  runtime <- Sys.time() - start_time
  
  list(
    cox_coef = lapply(fit_cox, coef),
    logit_coef = coef(fit_alpha),
    fit_cox=fit_cox,
    fit_logit=fit_alpha,
    cox_elf = data_longlong,
    logit_elf = data_logit,
    iterations = iter,
    loglik_trace = loglik_trace,
    runtime = runtime
  )
}
