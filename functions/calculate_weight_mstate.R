calculate_weight_mstate=function(data_longlong,id,tran_noncure,tran_cure,tran_correction){
  #data: data in long-long format with column id, trans, lik,pi
  #id: column names that represents subject id
  #tran_noncure:trans with which subjects could possibly be noncured
  ##tran_cure:trans with which subjects could possibly be cured
  #tran_correction: trans that if an subject had an event, he must be noncured 
 

  patientid=unique(data_longlong[,id])
  N=length(patientid)
  
  data_longlong$Prob_D_g1=NA
  data_longlong$Prob_D_g0=NA
  data_longlong$Expect_G_D=NA
  data_longlong$weight=NA
  
  for (pt in 1:N){
    row_index1=which(data_longlong[,id]==patientid[pt])
    row_index_g1=which(data_longlong[,id]==patientid[pt] & data_longlong$trans %in% tran_noncure)
    row_index_g0=which(data_longlong[,id]==patientid[pt] & data_longlong$trans %in% tran_cure)
    data_longlong$Prob_D_g1[row_index1]=prod(data_longlong$lik[row_index_g1])
    data_longlong$Prob_D_g0[row_index1]=prod(data_longlong$lik[row_index_g0])
    if (sum(data_longlong[which(data_longlong[,id]==patientid[pt]&data_longlong$trans %in% tran_correction),]$status)>=1) {data_longlong$Prob_D_g0[row_index1]=0} #check if status=1 for tran_correction
    data_longlong$Expect_G_D[row_index1]=data_longlong$Prob_D_g1[row_index1]*data_longlong$pi[row_index1]/(data_longlong$Prob_D_g1[row_index1]*data_longlong$pi[row_index1]+data_longlong$Prob_D_g0[row_index1]*(1-data_longlong$pi[row_index1]))
    
    data_longlong$weight[row_index_g1]=data_longlong$Expect_G_D[row_index_g1]
    data_longlong$weight[row_index_g0]=1-data_longlong$Expect_G_D[row_index_g0]
    
  }
  return(data_longlong)
}