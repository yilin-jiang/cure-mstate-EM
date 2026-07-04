
library(GenKern)
library(MultiCure)
library(mstate)
#load functions
files <- list.files("../function", full.names = TRUE)
invisible(lapply(files, source))

n_sim_start=1
n_sim_stop=500

###############
n_sim <- n_sim_stop-n_sim_start+1

scenario_results <- vector("list", n_sim)

set.seed(123)  # for reproducibility
seeds <- sample(1e6, 500)

for (i in n_sim_start:n_sim_stop) {
  set.seed(seeds[i])
  # --- simulate data ---
  NONE <- SimulateMultiCure(type = "NoMissingness")
  t_cure<-50
  ## ---- ELF method ----
  ########## initialize pi #########
  relap_rawprop=proportions(table(NONE$delta_R))[2]
  
  correction=0.1
  
  pi_init=1-relap_rawprop-correction #probability of being cured 
  
  ##### initialize betas #####
  #### standard mstate model ####
  tmat0<-trans.illdeath()
  data_long0<-msprep(data=NONE,trans=tmat0,time=c(NA,"Y_R","Y_D"),
                     status=c(NA,"delta_R","delta_D"),keep=c("X1","X2"))
  
  cox_tran12=coxph(Surv(Tstart,Tstop,status)~X1+X2,data=data_long0[data_long0$trans=="1",],method="breslow")
  cox_tran13=coxph(Surv(Tstart,Tstop,status)~X1+X2,data=data_long0[data_long0$trans=="2",],method="breslow")
  cox_tran23=coxph(Surv(Tstart,Tstop,status)~X1+X2,data=data_long0[data_long0$trans=="3",],method="breslow")
  
  model0=list(cox_tran12,cox_tran13,cox_tran23)
  
  # Run EM algorithm
  result_elf <- EM_mstate_cure(
    data_long0=data_long0 ,
    id="id",
    covar_cure=c("X1","X2"),
    covar_tran=list(c("X1","X2"),c("X1","X2"),c("X1","X2")),
    tran_state=c("12","13","23"),
    tran_noncure = c("1","2","3"),
    tran_cure = c("2"),
    pi_init = pi_init,
    model0 = model0,
    start=NULL,
    epsilon_tol =-1,
    max_iter = 100,
    t_cure=t_cure
  )

  ## ---Multicure method ---#
  Cov <- data.frame(X1 = NONE$X1, X2 = NONE$X2)
  VARS <- names(Cov)
  
  TransCov <- list(
    Trans13 = VARS,
    Trans24 = VARS,
    Trans14 = VARS,
    Trans34 = VARS,
    PNonCure = VARS
  )
  
  datWIDE <- data.frame(
    Y_R = NONE$Y_R,
    Y_D = NONE$Y_D,
    delta_R = NONE$delta_R,
    delta_D = NONE$delta_D,
    G = NONE$G
  )
  
  # --- model fit ---
  start_time <- Sys.time()
  fit_multicure <- MultiCure(
    iternum = 100,
    datWIDE = datWIDE,
    Cov = Cov,
    ASSUME = "ProportionalHazard",
    TransCov = TransCov,
    BASELINE = "cox"
  )
  
  
  beta_multicure<-fit_multicure$beta[c(1:4,7:9)]
  beta_multicure[5]<--beta_multicure[5]#transform beta0 to gamma_13
  names(beta_multicure)<-c("X1_12","X2_12","X1_13","X2_13","gamma_13","X1_23","X2_23")
  
  alpha_multicure<--fit_multicure$alpha #transform since G=1 indicates noncure in Beesley2019
  names(alpha_multicure)<-c("intercept","X1_cure","X2_cure")
  
  runtime_multicure <- Sys.time() - start_time
  result_multicure<-list(cox_coef=beta_multicure,alpha_coef=alpha_multicure,
                         runtime=runtime_multicure)
  
  scenario_results[[i]] <- list(
    ELF = result_elf,
    Multicure = result_multicure,
    seed=seeds[i]
  )
  
  cat("Finished simulation:", i, "\n")
}


## save scenario result
saveRDS(
  scenario_results,
  paste0("results_zerotail/results_zerotail.rds")
)