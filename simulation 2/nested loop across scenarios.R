setwd("../simulation 2")
source("../function/nested loop plot with mce.R")
library(data.table)
library(ggplot2)
library(latex2exp)

#load scenarios for simulation
load("scenarios_v1.RData")

###Combine results across scenarios####
all_tables <- list()

for (scen_num in 1:9) {
  
  summary_scenario_table <- read.table(
    paste0("simulation 2 results/performance_", scen_num, ".csv"),
    sep = ",",
    header = TRUE
  )
  
  summary_scenario_table$method <- factor(summary_scenario_table$method)
  summary_scenario_table$parameter <- factor(summary_scenario_table$parameter)
  summary_scenario_table <- as.data.table(summary_scenario_table)
  
  # add scenario column
  summary_scenario_table[, scen_num := scen_num]
  
  # store in list
  all_tables[[scen_num]] <- summary_scenario_table
}

# combine all tables
all_results <- rbindlist(all_tables)

all_results<-merge(all_results,as.data.table(full_factorial),by="scen_num")

####NLP across scenarios#####
estim_list <- c("bias", "rmse", "emp_se")
mce_list   <- c("mce_bias", "mce_rmse", "mce_emp_se")
ylab_list<-c("Bias","RMSE","Empirical SE")
parameters <- c("X_12","V_12","X_13","V_13","gamma_13","intercept_cure","X_cure","V_cure")
param_ylab<-c(
  "$X_{12}$", 
  "$V_{12}$", 
  "$X_{13}$", 
  "$V_{13}$", 
  "$\\gamma_{13}$", 
  "$intercept_{cure}$", 
  "$X_{cure}$", 
  "$V_{cure}$"
)

plot_list <- list()

for (i in seq_along(estim_list)) {
  estim <- estim_list[i]
  mce   <- mce_list[i]
  
  for (j in seq_along(parameters)) {
    param<-parameters[j]
    dat_sub <- all_results[parameter == param]
    dat_sub[, `:=`(
      "cure proportion" = factor(cure_prop),
      "censoring rate" = factor(censor_rate),
      method = factor(method),
      parameter = factor(parameter)
    )]
    ylab<-paste0(ylab_list[i]," for ", param_ylab[j])
    p <- ggplot_nlp_mce(
      dat = dat_sub,
      estim = estim,
      mce = mce,
      method_var = "method",
      true = 0,
      step_factors = c("cure proportion","censoring rate")   
    ) +
      labs(
        y = TeX(ylab)
      )+
      theme(axis.text.x = element_text(color = "black",size = 13),
            axis.text.y=element_text(color = "black",size = 13),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12),
            panel.grid.major.x = ggplot2::element_blank(),
            panel.grid.minor.x = ggplot2::element_blank())
    
    point_layer_idx <- which(sapply(p$layers, function(l) inherits(l$geom, "GeomPoint")))
    # Update the size of the points 
      p$layers[[point_layer_idx]]$aes_params$size <- 4  

    errorbar_layer_idx <- which(sapply(p$layers, function(l) inherits(l$geom, "GeomErrorbar")))
    p$layers[[errorbar_layer_idx]]$aes_params$linewidth <- 1 
  
      plot_list[[paste(param, estim, sep = "_")]] <- p
  }
}

for (name in names(plot_list)) {
  ggsave(
    filename = paste0("simulation 2 results/plot_", name, ".pdf"),
    plot = plot_list[[name]],
    width = 8,
    height = 8
  )
}

