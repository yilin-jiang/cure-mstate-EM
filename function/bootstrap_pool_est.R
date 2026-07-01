pool_est=function(alpha_est,beta_est,tran_state,no.iter){
  names(alpha_est)=paste0(names(alpha_est),"_alpha")
  betas_est=c()
  for (tran_index in 1:length(tran_state)){
    beta=beta_est[[tran_index]]
    names(beta)=paste0(names(beta),"_tran",tran_state[tran_index])
    betas_est=c(betas_est,beta)
  }
  return(c(alpha_est,betas_est,"no.iter"=no.iter))
}
