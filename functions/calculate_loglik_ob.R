calculate_loglik_ob=function(data_longlong,id){
  #data_longlong should contain id, pi, Prob_D_g1,Prob_D_g0
  
  data_likelihood=data_longlong[,c(id,"pi","Prob_D_g1","Prob_D_g0")]
  data_likelihood=data_likelihood[!duplicated(data_likelihood),]
  rownames(data_likelihood)=1:nrow(data_likelihood)
  data_likelihood$lik_ob=data_likelihood$pi*data_likelihood$Prob_D_g1+(1-data_likelihood$pi)*data_likelihood$Prob_D_g0
  data_likelihood$loglik_ob=log(data_likelihood$lik_ob)
  Loglikelihood_ob=sum(data_likelihood$loglik_ob)
  
  return(Loglikelihood_ob)
}