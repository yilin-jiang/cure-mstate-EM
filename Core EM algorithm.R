library(mstate)
data("ebmt4")
ebmt=ebmt4

#load functions
source("calculate_hazard.R")
source("create_datalonglong.R")
source("add_likelihood.R")
source("calculate_weight_mstate.R")
source("create_datalogit.R")
source("calculate_loglik_ob.R")
source("update_pi.R")

########## initial values #########
###### pi: probability of being noncured ######
relap_rawprop=proportions(table(ebmt$rel.s))[2]

correction=0.15

pi=relap_rawprop+correction

#### prepare data in long format ####
#transition 5 to 6 does not need to be involved in EM algorithm
tmat0=transMat(x = list(c(2, 3, 5, 6), c(4, 5, 6), c(4, 5, 6), c(5, 6),
                       c(), c()), names = c("Tx", "Rec", "AE", "Rec+AE", "Rel", "Death"))
tmat0

data_long0= msprep(data = ebmt, trans = tmat0, time = c(NA, "rec", "ae","recae", "rel", "srv"), 
                   status = c(NA, "rec.s", "ae.s", "recae.s","rel.s", "srv.s"), keep = c("match", "proph", "year", "agecl"))
events(data_long0)

##### betas #####
#### standard mstate model ####
tran=c("12","13","15","16","24","25","26","34","35","36","45","46") 

cox_tran12=coxph(Surv(Tstart,Tstop,status)~year+agecl+proph+match,data=data_long0[data_long0$trans=="1",],method="breslow")
cox_tran13=coxph(Surv(Tstart,Tstop,status)~year+agecl+proph+match,data=data_long0[data_long0$trans=="2",],method="breslow")
cox_tran15=coxph(Surv(Tstart,Tstop,status)~year+agecl+proph+match,data=data_long0[data_long0$trans=="3",],method="breslow")
cox_tran16=coxph(Surv(Tstart,Tstop,status)~year+agecl+proph+match,data=data_long0[data_long0$trans=="4",],method="breslow")
cox_tran24=coxph(Surv(Tstart,Tstop,status)~year+agecl+proph+match,data=data_long0[data_long0$trans=="5",],method="breslow")
cox_tran25=coxph(Surv(Tstart,Tstop,status)~year+agecl+proph+match,data=data_long0[data_long0$trans=="6",],method="breslow")
cox_tran26=coxph(Surv(Tstart,Tstop,status)~year+agecl+proph+match,data=data_long0[data_long0$trans=="7",],method="breslow")
cox_tran34=coxph(Surv(Tstart,Tstop,status)~year+agecl+proph+match,data=data_long0[data_long0$trans=="8",],method="breslow")
cox_tran35=coxph(Surv(Tstart,Tstop,status)~year+agecl+proph+match,data=data_long0[data_long0$trans=="9",],method="breslow")
cox_tran36=coxph(Surv(Tstart,Tstop,status)~year+agecl+proph+match,data=data_long0[data_long0$trans=="10",],method="breslow")
cox_tran45=coxph(Surv(Tstart,Tstop,status)~year+agecl+proph+match,data=data_long0[data_long0$trans=="11",],method="breslow")
cox_tran46=coxph(Surv(Tstart,Tstop,status)~year+agecl+proph+match,data=data_long0[data_long0$trans=="12",],method="breslow")

model0=list(cox_tran12,cox_tran13,cox_tran15,cox_tran16,cox_tran24,cox_tran25,cox_tran26,cox_tran34,cox_tran35,cox_tran36,cox_tran45,cox_tran46)

########## EM algorithm #########
###### E step: calculate the expected complete data likelihood based on the observed data and recent parameter estimates#####
##### calculate hazard #######
data_long=calculate_hazard(model=model0,tran=tran,data_long=data_long0,data_mstate=data_long0,id="id")

##### create data in extended long format#####
data_longlong=create_datalonglong(data_long=data_long,transition=c("12","13","16","24","26","34","36","46"),id="id")

##### calculate data likelihood contribution per row #####
data_longlong=add_likelihood(data=data_longlong)

##### calculate conditional probabilities and assign weights#####
##### G = 1 means noncured; 0 means cured########################
##### P(Di|Gi=1) ：product of likelihood from trans =1 to 12#####
##### P(Di|Gi=0) ：0 when status=1 for trans in c(3,6,9,11), otherwise product of likelihood from trans =1.2,2.2,4.2,5.2,7.2,8.2,10.2,12.2#####
##### E(Gi|Di)=P(Gi=1|Di)=P(Di|Gi=1)*pi/(P(Di|Gi=1)*pi+P(Di|Gi=0)*(1-pi))#####

data_longlong$pi=pi
data_longlong=calculate_weight_mstate(data_longlong=data_longlong,id="id",
                                      tran_noncure=c("1","2","3","4","5","6","7","8","9","10","11","12"),
                                      tran_cure=c("1.2","2.2","4.2","5.2","7.2","8.2","10.2","12.2"),tran_correction=c("3","6","9","11"))

##### create data_logit for logistic regression ######
data_logit=create_datalogit(data_longlong=data_longlong,id="id",cov=c("year","agecl","proph","match"))

###### M step:maximize expected data likelihood with respect to parameters, alphas & betas #########

#### fit weighted mstate model ####
##### assume same coefficient estimates and proportional transition hazard for transition 1 vs 12, 2 vs 22, 4 vs 42, 5 vs 52, 7 vs 72, 8 vs 82, 10 vs 102, 122 vs 122 #####
data_mstate=data_longlong[data_longlong$weight!=0,]

cox_tran12=coxph(Surv(Tstart,Tstop,status)~year+agecl+proph+match+cure12,weights=weight,data=data_mstate[data_mstate$fromto=="12",],method="breslow")
cox_tran13=coxph(Surv(Tstart,Tstop,status)~year+agecl+proph+match+cure13,weights=weight,data=data_mstate[data_mstate$fromto=="13",],method="breslow")
cox_tran15=coxph(Surv(Tstart,Tstop,status)~year+agecl+proph+match,weights=weight,data=data_mstate[data_mstate$fromto=="15",],method="breslow")
cox_tran16=coxph(Surv(Tstart,Tstop,status)~year+agecl+proph+match+cure16,weights=weight,data=data_mstate[data_mstate$fromto=="16",],method="breslow")
cox_tran24=coxph(Surv(Tstart,Tstop,status)~year+agecl+proph+match+cure24,weights=weight,data=data_mstate[data_mstate$fromto=="24",],method="breslow")
cox_tran25=coxph(Surv(Tstart,Tstop,status)~year+agecl+proph+match,weights=weight,data=data_mstate[data_mstate$fromto=="25",],method="breslow")
cox_tran26=coxph(Surv(Tstart,Tstop,status)~year+agecl+proph+match,weights=weight,data=data_mstate[data_mstate$fromto=="26",],method="breslow")
cox_tran34=coxph(Surv(Tstart,Tstop,status)~year+agecl+proph+match+cure34,weights=weight,data=data_mstate[data_mstate$fromto=="34",],method="breslow")
cox_tran35=coxph(Surv(Tstart,Tstop,status)~year+agecl+proph+match,weights=weight,data=data_mstate[data_mstate$fromto=="35",],method="breslow")
cox_tran36=coxph(Surv(Tstart,Tstop,status)~year+agecl+proph+match+cure36,weights=weight,data=data_mstate[data_mstate$fromto=="36",],method="breslow")
cox_tran45=coxph(Surv(Tstart,Tstop,status)~year+agecl+proph+match,weights=weight,data=data_mstate[data_mstate$fromto=="45",],method="breslow")
cox_tran46=coxph(Surv(Tstart,Tstop,status)~year+agecl+proph+match+cure46,weights=weight,data=data_mstate[data_mstate$fromto=="46",],method="breslow")

cox_mstate=list(cox_tran12,cox_tran13,cox_tran15,cox_tran16,cox_tran24,cox_tran25,cox_tran26,cox_tran34,cox_tran35,cox_tran36,cox_tran45,cox_tran46)

#### fit weighted logistic model ####
fit_alpha=glm(cure~year+agecl+proph+match,weights=weight,data=data_logit,family=binomial(link="logit"))

########calculate observed log-likelihood#########
Loglik_ob=calculate_loglik_ob(data_longlong=data_longlong,"id")
Loglik_ob 

#### EM loop ####
loglik_result=c(Loglik_ob)
no.iter=0
epsilon=1

while(epsilon >0.0001){
  
  #E step
  data_longlong=update_pi(fit_alpha=fit_alpha,cov=c("year","agecl","proph","match"),data_longlong = data_longlong)
  
  data_longlong=calculate_hazard(model=cox_mstate,tran=tran,data_mstate=data_mstate,data_long=data_longlong,id="id")
  
  data_longlong=add_likelihood(data=data_longlong)
  
  
  data_longlong=calculate_weight_mstate(data_longlong=data_longlong,id="id",
                                        tran_noncure=c("1","2","3","4","5","6","7","8","9","10","11","12"),
                                        tran_cure=c("1.2","2.2","4.2","5.2","7.2","8.2","10.2","12.2"),tran_correction=c("3","6","9","11"))
  
  data_logit=create_datalogit(data_longlong=data_longlong,id="id",cov=c("year","agecl","proph","match"))
  
  #M step
  data_mstate=data_longlong[data_longlong$weight!=0,]
  cox_tran12=coxph(Surv(Tstart,Tstop,status)~year+agecl+proph+match+cure12,weights=weight,data=data_mstate[data_mstate$fromto=="12",],method="breslow")
  cox_tran13=coxph(Surv(Tstart,Tstop,status)~year+agecl+proph+match+cure13,weights=weight,data=data_mstate[data_mstate$fromto=="13",],method="breslow")
  cox_tran15=coxph(Surv(Tstart,Tstop,status)~year+agecl+proph+match,weights=weight,data=data_mstate[data_mstate$fromto=="15",],method="breslow")
  cox_tran16=coxph(Surv(Tstart,Tstop,status)~year+agecl+proph+match+cure16,weights=weight,data=data_mstate[data_mstate$fromto=="16",],method="breslow")
  cox_tran24=coxph(Surv(Tstart,Tstop,status)~year+agecl+proph+match+cure24,weights=weight,data=data_mstate[data_mstate$fromto=="24",],method="breslow")
  cox_tran25=coxph(Surv(Tstart,Tstop,status)~year+agecl+proph+match,weights=weight,data=data_mstate[data_mstate$fromto=="25",],method="breslow")
  cox_tran26=coxph(Surv(Tstart,Tstop,status)~year+agecl+proph+match,weights=weight,data=data_mstate[data_mstate$fromto=="26",],method="breslow")
  cox_tran34=coxph(Surv(Tstart,Tstop,status)~year+agecl+proph+match+cure34,weights=weight,data=data_mstate[data_mstate$fromto=="34",],method="breslow")
  cox_tran35=coxph(Surv(Tstart,Tstop,status)~year+agecl+proph+match,weights=weight,data=data_mstate[data_mstate$fromto=="35",],method="breslow")
  cox_tran36=coxph(Surv(Tstart,Tstop,status)~year+agecl+proph+match+cure36,weights=weight,data=data_mstate[data_mstate$fromto=="36",],method="breslow")
  cox_tran45=coxph(Surv(Tstart,Tstop,status)~year+agecl+proph+match,weights=weight,data=data_mstate[data_mstate$fromto=="45",],method="breslow")
  cox_tran46=coxph(Surv(Tstart,Tstop,status)~year+agecl+proph+match+cure46,weights=weight,data=data_mstate[data_mstate$fromto=="46",],method="breslow")
  
  cox_mstate=list(cox_tran12,cox_tran13,cox_tran15,cox_tran16,cox_tran24,cox_tran25,cox_tran26,cox_tran34,cox_tran35,cox_tran36,cox_tran45,cox_tran46)
  
  
  fit_alpha=glm(cure~year+agecl+proph+match,weights=weight,data=data_logit,family=binomial(link="logit"))
  
  ########calculate observed log-likelihood#########
  Loglik_ob_new=calculate_loglik_ob(data_longlong=data_longlong,"id")
  
  epsilon=Loglik_ob_new-Loglik_ob
  Loglik_ob=Loglik_ob_new
  
  no.iter=no.iter+1
  loglik_result=append(loglik_result,Loglik_ob_new)
  print(Loglik_ob)
  print(no.iter)
}

#save results from the last iteration
saveRDS(cox_mstate,file="cox_mstate.rda")
saveRDS(fit_alpha,file="fit_alpha.rda")
write.csv(data_longlong,"data_longlong.csv",row.names=FALSE)