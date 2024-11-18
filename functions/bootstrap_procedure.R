boot_pro=function(boot_id,data,data_coxfit,time,status,tmat,var,tran_cure,alpha_fit,beta_fit,subjectid,max.iter=80){
  data_boot=gendataboot(data=data,seed=boot_id)
  data_long_boot= msprep(data = data_boot, trans = tmat, time = c(NA, time), 
                          status = c(NA,status), 
                          keep = var)
  data_long_boot$fromto=paste0(data_long_boot$from,data_long_boot$to)
  tran<<-as.character(sort(as.numeric(unique(data_long_boot$fromto)))) #Models in beta_fit must be in the same order as tran
  data_longlong_boot=create_datalonglong(data_long=data_long_boot,transition=tran_cure,id=subjectid)
  result_boot=em(data_longlong_boot=data_longlong_boot,data_coxfit=data_coxfit,subjectid=subjectid,
     alpha_fit=alpha_fit,beta_fit=beta_fit,var=var,tran=tran,tran_cure=tran_cure,eps=0.001,max.iter=max.iter)
  return(result_boot=data.frame(t(result_boot),"seed"=boot_id))
  }

gendataboot=function(data,seed){
  set.seed(seed)
  N=nrow(data)
  id_sample=sample(1:N,N,replace=TRUE)
  data_boot=data[id_sample,]
  rownames(data_boot)=1:N
  return(data_boot)
}

em=function(data_longlong_boot,data_coxfit,subjectid,alpha_fit,beta_fit,var,tran,tran_cure,eps,max.iter){
  beta_fit_boot=beta_fit
  alpha_fit_boot=alpha_fit
  var_cure=paste0("cure",tran_cure)
  data_mstate<<-data_coxfit[data_coxfit$weight!=0,]
  
  tran_set=unique(data_longlong_boot$trans)
  tran_noncure=tran_set[which(as.numeric(tran_set) %% 1 ==0)]
  tran_cure=tran_set[!tran_set %in% tran_noncure]
  tran_correction=tran_noncure[which(!as.numeric(tran_noncure) %in% as.integer(tran_cure))]


  no.iter=0
  epsilon=1

  while (epsilon >eps && no.iter<max.iter){
  #E STEP#
  data_longlong_boot=update_pi(fit_alpha=alpha_fit_boot,cov=var,data_longlong = data_longlong_boot)
  data_longlong_boot=calculate_hazard(model=beta_fit_boot,tran=tran,data_long=data_longlong_boot,data_mstate=data_mstate,id=subjectid)
  data_longlong_boot=add_likelihood(data=data_longlong_boot)
  data_longlong_boot=calculate_weight_mstate(data_longlong=data_longlong_boot,id=subjectid,
                                        tran_noncure=tran_noncure,
                                        tran_cure=tran_cure,tran_correction=tran_correction)
  data_logit=create_datalogit(data_longlong=data_longlong_boot,id=subjectid,cov=var)
  
  #M STEP#
  data_mstate<<-data_longlong_boot[data_longlong_boot$weight!=0,]
  
  for(model_number in 1:length(tran)){
    beta_fit_boot[[model_number]]=coxph(beta_fit[[model_number]]$formula,weights=weight,data=data_mstate[data_mstate$fromto == tran[model_number],],method="efron")
  }
  
  alpha_fit_boot=glm(alpha_fit$formula,weights=weight,data=data_logit,family=binomial(link="logit"))  
  
  #calculate observed log-likelihood#
  Loglik_ob_new=calculate_loglik_ob(data_longlong=data_longlong_boot,id=subjectid)
  if (no.iter ==0){
    Loglik_ob=Loglik_ob_new-1}
  epsilon=Loglik_ob_new-Loglik_ob
  Loglik_ob=Loglik_ob_new
  no.iter=no.iter+1

  
  print(Loglik_ob)
  print(no.iter)
  print(epsilon)
  }

  result_boot=extract_est(alpha_fit_boot=alpha_fit_boot,beta_fit_boot=beta_fit_boot,tran=tran,no.iter=no.iter,Loglik_ob=Loglik_ob,epsilon=epsilon)
  return(result_boot)
  }

extract_est=function(alpha_fit_boot,beta_fit_boot,tran,no.iter,Loglik_ob,epsilon){
  alpha_est=coef(alpha_fit_boot)
  names(alpha_est)=paste0(names(alpha_est),"_alpha")
  beta_est=c()
  for (tran_index in 1:length(tran)){
    beta=coef(beta_fit_boot[[tran_index]])
    names(beta)=paste0(names(beta),"_tran",tran[tran_index])
    beta_est=c(beta_est,beta)
  }
  return(c(alpha_est,beta_est,"no.iter"=no.iter,"Loglik_ob"=Loglik_ob,"epsilon"=epsilon))
}
