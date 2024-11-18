calculate_hazard=function(model,tran,data_long,data_mstate,id,...){
  
  #model: a list of coxph models on data (separate form of mstate model)
  #tran: transition vector that includes all possible from-to pairs, in the same order as in the model list
  #data_long: data in long or extended long format with column id, from, to, time, status and covariates
  #data_mstate: the data that is referred to in model fit
  #id: column names that represents subject id
  
  data_long$cumHaz=NA
  data_long$hazard=NA
  
  patientid=unique(data_long[,id])
  N=length(patientid)
  
  for (pt in 1:N){
    row_index=which(data_long[,id]==patientid[pt])
    no.row=length(row_index)
    for (i in 1 :no.row){
      fromto=paste0(data_long$from[row_index[i]],data_long$to[row_index[i]])
      model_number<<-which(tran==fromto)
      
      model_fit=survfit(model[[model_number]],newdata=data_long[row_index[i],],se.fit=FALSE)
      Haz=data.frame(time=model_fit$time,cumhaz=model_fit$cumhaz)
      Haz=rbind(c(0,0),Haz)
      Haz=Haz[!duplicated(Haz$cumhaz),]
      Haz$hazard=diff(c(0,Haz$cumhaz))
      
      #find the location of the timepoint
      Haz$timedif_Tstop=data_long[row_index[i],"Tstop"]-Haz$time
      t_tstop=max(which(Haz$timedif_Tstop>=0))
      Haz$timedif_Tstart=data_long[row_index[i],"Tstart"]-Haz$time
      t_tstart=max(which(Haz$timedif_Tstart>=0))
      
      data_long[row_index[i],"cumHaz"]=Haz[t_tstop,"cumhaz"]-Haz[t_tstart,"cumhaz"]
      data_long[row_index[i],"hazard"]=ifelse(data_long$status[row_index[i]]==0,1,Haz[t_tstop,"hazard"])
    }
  }
  return(data_long)
}