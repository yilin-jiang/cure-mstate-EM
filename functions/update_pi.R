update_pi=function(fit_alpha,cov,data_longlong){

  data_longlong$pi = predict(fit_alpha, newdata=data_longlong[,cov],
                          type="response")
  return(data_longlong)

}