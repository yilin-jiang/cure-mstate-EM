add_likelihood=function(data){
  #data should include column hazard, status, cumHaz
  data$lik=(data$hazard^data$status)*exp(-data$cumHaz)
  return(data)
}