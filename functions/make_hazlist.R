make_hazlist=function(tran.index,data_cov,time_list,tran.nc2,curecol){
  #tran.nc2 is the transition list requiring to change cure column
  #curecol is the corresponding cure column names based on tran.nc2 list
  data=data_cov
  if(tran.index %in% tran.nc2){
    curecolname=curecol[which(tran.nc2==tran.index)]
    data[,curecolname]=1}
  
  model=cox_full[[round(tran.index)]]
  model_fit=survfit(model,newdata=data,se.fit=FALSE)
  
  Haz_list=data.frame(time=model_fit$time,Haz=model_fit$cumhaz)
  Haz_list=rbind(c(0,0),Haz_list)
  Haz_list$trans=round(tran.index)
  
  insert_index=match(Haz_list$time,time_list)
  
  time_frame_trans=data.frame(time=time_list,Haz=NA,trans=NA)
  time_frame_trans[insert_index,]=Haz_list
  time_frame_trans$Haz=mstate:::NAfix(time_frame_trans$Haz)
  time_frame_trans$trans=round(tran.index)
  
  return(time_frame_trans)
}