run_scenario <- function(scen,n_sim,base_seed,full_factorial){
  
  seeds <-  base_seed + scen*100000 + 1:n_sim
  
  scenario_results <- vector("list", n_sim)
  
  for(i in 1:n_sim){
    
    cat("Scenario", scen, "Dataset", i, "\n")
    
    set.seed(seeds[i])
    
    ## ---- Generate data ----
    data<-generate_data(n=2000,alpha0=full_factorial$alpha0[scen],scale_cens = full_factorial$cen_scale[scen])
    
    ## ---- ELF method ----
    ########## initialize pi #########
    relap_rawprop=proportions(table(data$wide_data_obs$status2))[2]
    
    correction=0.1
    
    pi_init=1-relap_rawprop-correction #probability of being cured 
    
    ##### initialize betas #####
    #### standard mstate model ####
    data_long0<-data$long_data_obs
    
    cox_tran12=coxph(Surv(Tstart,Tstop,status)~X+V,data=data_long0[data_long0$trans=="1",],method="breslow")
    cox_tran13=coxph(Surv(Tstart,Tstop,status)~X+V,data=data_long0[data_long0$trans=="2",],method="breslow")
    cox_tran23=coxph(Surv(Tstart,Tstop,status)~X+V,data=data_long0[data_long0$trans=="3",],method="breslow")
    
    model0=list(cox_tran12,cox_tran13,cox_tran23)
    
    
    # Run EM algorithm
    result_elf <- EM_mstate_cure(
      data_long0=data_long0 ,
      id="id",
      covar_cure=c("X","V"),
      covar_tran=list(c("X","V"),c("X","V"),c("X","V")),
      tran_state=c("12","13","23"),
      tran_noncure = c("1","2","3"),
      tran_cure = c("2"),
      pi_init = pi_init,
      model0 = model0,
      epsilon_tol = 1e-4,
      max_iter = 200,
      t_cure=NA
    )
    
    ## ---- Multicure method ----
    data_wide<-data$wide_data_obs
    Cov <- data_wide[,c("X","V")]
    VARS <- names(Cov)
    TransCov <-list(Trans13 = VARS, Trans24 = VARS, Trans14 = VARS, Trans34 = VARS, PNonCure = VARS)
    datWIDE <- data_wide[,c("T2","T3","status2","status3")]
    colnames(datWIDE)<-c("Y_R", "Y_D", "delta_R" , "delta_D")
    datWIDE$G<-ifelse(datWIDE$delta_R==1,1,NA) # This takes value 1 for known non-cured, 0 for "known" cured and NA for unknown cur'e status
    iternum=result_elf$iterations
    
    start_time <- Sys.time()
    fit_multicure <- MultiCure(iternum = iternum, datWIDE, Cov, ASSUME = "ProportionalHazard", TransCov = TransCov, BASELINE = "cox")
    
    beta_multicure<-fit_multicure$beta[c(1:4,7:9)]
    beta_multicure[5]<--beta_multicure[5]#transform beta0 to gamma_13
    names(beta_multicure)<-c("X_12","V_12","X_13","V_13","gamma_13","X_23","V_23")
    
    alpha_multicure<--fit_multicure$alpha #transform since G=1 indicates noncure in Beesley2019
    names(alpha_multicure)<-c("intercept","X_cure","V_cure")
    
    runtime_multicure <- Sys.time() - start_time
    result_multicure<-list(cox_coef=beta_multicure,alpha_coef=alpha_multicure,
                           runtime=runtime_multicure)
    
    scenario_results[[i]] <- list(
      ELF = result_elf,
      Multicure = result_multicure,
      seed=seeds[i]
    )
    
  }
  
  ## save scenario result
  saveRDS(
    scenario_results,
    paste0("simulation 2 results/results_scenario_", scen, ".rds")
  )
  
}