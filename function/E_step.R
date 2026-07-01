E_step <- function(data_longlong, fit_alpha, fit_cox, tran_state,tran_noncure,tran_cured, tran_correction,id,t_cure){
  
  covar_cure<-all.vars(fit_alpha$formula)[-1]
 
   data_longlong <- update_pi(
    fit_alpha = fit_alpha,
    cov = covar_cure,
    data_longlong = data_longlong
  )
  data_longlong<-calculate_hazard_dt(model=fit_cox,tran=tran_state,data_long=data_longlong)
  
  data_longlong <- add_likelihood(data_longlong)
  
  data_longlong <- calculate_weight_mstate(
    data_longlong=data_longlong,
    id,
    tran_noncure =tran_noncure,
    tran_cure = tran_cured,
    tran_correction = tran_correction,t_cure
  )
  data_logit=create_datalogit(data_longlong=data_longlong,id,cov=covar_cure)
  
  list(
    data_longlong = data_longlong,
    data_logit = data_logit
  )
}
