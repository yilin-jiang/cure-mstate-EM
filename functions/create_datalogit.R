create_datalogit=function(data_longlong,id,cov){
  #data_longlong:data with column id, Expect_G_D,cov
  #id: colname of subject id
  #cov: cov in data_longlong,used for weighted logistic model
 
  data_logit=data_longlong[,c(id,"Expect_G_D",cov)]
  data_logit=data_logit[!duplicated(data_logit),]
  colnames(data_logit)[2]="weight"
  data_logit$cure=1 #noncured
  
  data_logit_cure0=data_logit
  data_logit_cure0$cure=0
  data_logit_cure0$weight=1-data_logit_cure0$weight
  data_logit=rbind(data_logit,data_logit_cure0)
  order_index=order(data_logit[,id],data_logit$cure)
  data_logit=data_logit[with(data_logit,order_index),]
  rownames(data_logit)=1:nrow(data_logit)
  
  return(data_logit)
}