
load("scenarios_v1.RData")

#load packages
library(dplyr)
library(tidyr)

extract_elf_coef <- function(rep) {
  c12 <- rep[["ELF"]][["cox_coef"]][["12"]]
  c13 <- rep[["ELF"]][["cox_coef"]][["13"]]
  c23 <- rep[["ELF"]][["cox_coef"]][["23"]]
  c_cure<-rep[["ELF"]][["logit_coef"]]
  
  data.frame(
    X_12 = unname(c12["X"]),
    V_12 = unname(c12["V"]),
    X_13 = unname(c13["X"]),
    V_13 = unname(c13["V"]),
    gamma_13 = unname(c13["cure13"]),
    X_23 = unname(c23["X"]),
    V_23 = unname(c23["V"]),
    intercept_cure=unname(c_cure["(Intercept)"]),
    X_cure=unname(c_cure["X"]),
    V_cure=unname(c_cure["V"]),
    method="ELF"
  )
}

extract_multicure_coef <- function(rep) {
  cox_coef <- rep[["Multicure"]][["cox_coef"]]
  alpha_coef<-rep[["Multicure"]][["logit_coef"]]

  c(cox_coef,
    alpha_coef,
    method="Multicure"
  )
}

summarize_sim <- function(est_vec, true) {
  nsim <- length(est_vec)
  diffs <- est_vec - true
  
  # Performance measures
  bias   <- mean(diffs)
  rmse   <- sqrt(mean(diffs^2))
  emp_se <- sd(est_vec)
  mean_est <- mean(est_vec)
  
  # Monte Carlo errors
  mce_bias   <- sd(diffs) / sqrt(nsim)
  mce_emp_se <- emp_se / sqrt(2 * (nsim - 1))
  mce_rmse   <- sd(diffs^2) / (2 * rmse * sqrt(nsim))  # delta method
  
  data.frame(true_value = true,mean_est = mean_est,bias, rmse, emp_se, mce_bias, mce_emp_se, mce_rmse)
}

#load simulation results per scenario
for (i in 1:9) {
  assign(
    paste0("results_scenario_", i),
    readRDS(paste0("simulation 2 results/results_scenario_", i, ".rds"))
  )
}
for (scen in 1:9){

elf_coef<-do.call(rbind, lapply(get(paste0("results_scenario_",scen)), extract_elf_coef))
multicure_coef<-do.call(rbind, lapply(get(paste0("results_scenario_",scen)), extract_multicure_coef))
colnames(multicure_coef)[8]<-"intercept_cure"

all_results<-rbind(elf_coef, multicure_coef)
all_results<- all_results %>%
  mutate(across(-method, ~ as.numeric(.)))%>%
  mutate(method = as.factor(method))


true_coef<-data.frame("X_12"=0.5, "V_12"=-0.3,"X_13"=0.3,"V_13"=0.2,"gamma_13"=-0.4,      
                      "X_23"=0.4, "V_23"=0.1, "intercept_cure"=full_factorial$alpha0[which(full_factorial$scen_num==scen)],
                      "X_cure"=0.5,"V_cure"=-0.2)



all_summaries <- list()
methods<-unique(all_results$method)
params <- c("X_12", "V_12", "X_13", "V_13", "gamma_13", 
            "X_23", "V_23", "intercept_cure", "X_cure", "V_cure")

for (m in methods) {
  df_method <- all_results %>% filter(method == m)
  
  for (p in params) {
    est_vec <- df_method[[p]]  # vector of 500 repetitions
    true_val <- true_coef[[p]]
    res <- summarize_sim(est_vec, true_coef[[p]])
    res$method <- m
    res$parameter <- p
    all_summaries[[length(all_summaries) + 1]] <- res
  }
}

tidy_summary <- do.call(rbind, all_summaries)
tidy_summary <- tidy_summary[, c("method", "parameter", "true_value", "mean_est","bias", "rmse", "emp_se",
                                 "mce_bias", "mce_emp_se", "mce_rmse")]

write.table(tidy_summary,file=paste0("simulation 2 results/performance_",scen,".csv"),sep = ",",row.names=FALSE)

}

##### computation time comparison####
scen <- 1:9

runtime_list <- lapply(scen, function(s) {
  scen_obj <- get(paste0("results_scenario_", s))  # get actual scenario object
  
  do.call(rbind, lapply(scen_obj, function(repetition) {
    data.frame(
      ELF = repetition[["ELF"]]$runtime,
      Multicure = repetition[["Multicure"]]$runtime
    )
  }))
})

all_runtime <- do.call(rbind, runtime_list)

# convert to long format
all_runtime_long <- all_runtime %>%
  pivot_longer(
    cols = c(ELF, Multicure),
    names_to = "method",
    values_to = "runtime"
  )

runtime_summary <- all_runtime_long %>%
  group_by(method) %>%
  summarise(
    mean_runtime = mean(runtime),
    sd_runtime = sd(runtime),
    median_runtime = median(runtime),
    IQR_runtime = IQR(runtime),
    n = n()
  )

write.table(runtime_summary,file="simulation 2 results/runtime_summary.csv",sep = ",",row.names=FALSE)
