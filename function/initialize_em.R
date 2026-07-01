initialize_em <- function(model0, data_long0, pi_init, tran_state,tran_state_cure,tran_noncure,tran_cured,tran_correction,covar_cure,id,t_cure){
  
  data_long <- calculate_hazard_dt(
    model = model0,
    tran = tran_state,
    data_long = data_long0
  )
  
  data_longlong <- create_datalonglong(
    data_long = data_long,
    transition = tran_state_cure,
    id = id
  )
  
  data_longlong <- add_likelihood(data_longlong)
  
  data_longlong$pi <- pi_init
  
  data_longlong<-calculate_weight_mstate(data_longlong=data_longlong,
                                         id,
                                        tran_noncure,
                                        tran_cure=tran_cured,tran_correction,t_cure)
  
  data_logit <- create_datalogit(
    data_longlong,
    id = id,
    cov = covar_cure
  )
  
  list(
    data_longlong = data_longlong,
    data_logit = data_logit
  )
}
