gendataboot=function(data,seed){
  set.seed(seed)
  N=nrow(data)
  id_sample=sample(1:N,N,replace=TRUE)
  data_boot=data[id_sample,]
  rownames(data_boot)=1:N
  return(data_boot)
}