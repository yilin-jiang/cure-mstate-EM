#load R packages and functions needed
library(dplyr)
library(tidyr)
library(ggplot2)
library(latex2exp)

source("../function/nested loop plot with mce.R")

#functions needed
process_simulation_results <- function(file_path, 
                                       file_name) {
  # Load data
  results <- readRDS(paste0(file_path, file_name))
  
  # Define extraction functions
  extract_elf_coef_s1 <- function(rep) {
    c12 <- rep[["ELF"]][["cox_coef"]][["12"]]
    c13 <- rep[["ELF"]][["cox_coef"]][["13"]]
    c23 <- rep[["ELF"]][["cox_coef"]][["23"]]
    c_cure <- rep[["ELF"]][["logit_coef"]]
    
    data.frame(
      X1_12 = unname(c12["X1"]),
      X2_12 = unname(c12["X2"]),
      X1_13 = unname(c13["X1"]),
      X2_13 = unname(c13["X2"]),
      gamma_13 = unname(c13["cure13"]),
      X1_23 = unname(c23["X1"]),
      X2_23 = unname(c23["X2"]),
      intercept_cure = unname(c_cure["(Intercept)"]),
      X1_cure = unname(c_cure["X1"]),
      X2_cure = unname(c_cure["X2"]),
      method = "ELF"
    )
  }
  
  extract_multicure_coef <- function(rep) {
    cox_coef <- rep[["Multicure"]][["cox_coef"]]
    alpha_coef <- rep[["Multicure"]][["logit_coef"]]
    
    # Combine into a named numeric vector
    c(cox_coef, alpha_coef, method = "Multicure")
  }
  
  # Apply extraction
  elf_coef <- do.call(rbind, lapply(results, extract_elf_coef_s1))
  multicure_coef <- do.call(rbind, lapply(results, extract_multicure_coef))
  
  # Format data
  colnames(multicure_coef)[8] <- "intercept_cure"
  
  all_results <- rbind(elf_coef, multicure_coef)
  all_results <- all_results %>%
    mutate(across(-method, ~ as.numeric(.))) %>%
    mutate(method = as.factor(method))
  
  # Set true coefficients
  params <- c("X1_12", "X2_12", "X1_13", "X2_13", "gamma_13", 
              "X1_23", "X2_23", "intercept_cure", "X1_cure", "X2_cure")
  beta_true<-c(0.5,0.5,0.5,0.5,0,0.5,0.5)
  alpha_true<-c(-0.5,-0.5,-0.5) #since G=1 represents cured in our study but uncured in Beesley 2019, the coefficient signs are flipped
  true_coef <- c(beta_true, alpha_true)
  names(true_coef) <- params
  
  # Summary function
  summarize_sim <- function(est_vec, true) {
    nsim <- length(est_vec)
    diffs <- est_vec - true
    
    bias   <- mean(diffs)
    rmse   <- sqrt(mean(diffs^2))
    emp_se <- sd(est_vec)
    mean_est <- mean(est_vec)
    
    mce_bias   <- sd(diffs) / sqrt(nsim)
    mce_emp_se <- emp_se / sqrt(2 * (nsim - 1))
    mce_rmse   <- sd(diffs^2) / (2 * rmse * sqrt(nsim))
    
    data.frame(true_value = true, mean_est = mean_est,
               bias, rmse, emp_se, mce_bias, mce_emp_se, mce_rmse)
  }
  
  # Compute summaries
  methods <- unique(all_results$method)
  all_summaries <- list()
  
  for (m in methods) {
    df_method <- all_results %>% filter(method == m)
    
    for (p in params) {
      est_vec <- df_method[[p]]
      true_val <- true_coef[[p]]
      res <- summarize_sim(est_vec, true_val)
      res$method <- m
      res$parameter <- p
      all_summaries[[length(all_summaries) + 1]] <- res
    }
  }
  
  tidy_summary <- do.call(rbind, all_summaries)
  tidy_summary <- tidy_summary[, c("method", "parameter", "true_value", "mean_est", 
                                   "bias", "rmse", "emp_se",
                                   "mce_bias", "mce_emp_se", "mce_rmse")]
  
  # Save CSV
  output_csv <- paste0(file_path, "performance_", 
                       gsub("\\.rds$", "", file_name), ".csv")
  write.table(tidy_summary, file = output_csv, sep = ",", row.names = FALSE)
  return(tidy_summary)
}
create_nested_plots <- function(summary_table, path,scen) {
  level = c("X1_12", "X2_12", "X1_13", "X2_13", "gamma_13", 
            "intercept_cure", "X1_cure", "X2_cure")
  summary_table$method <- factor(summary_table$method)
  summary_table$parameter <- factor(summary_table$parameter, levels = level)
  summary_table <- data.table::as.data.table(summary_table)
  summary_table <- summary_table[!is.na(parameter)]
  estimates <- list(
    list(estim = "bias", mce = "mce_bias", ylab = "Bias"),
    list(estim = "rmse", mce = "mce_rmse", ylab = "RMSE"),
    list(estim = "emp_se", mce = "mce_emp_se", ylab = "Empirical SE")
  )
  
  for (e in estimates) {
    p_mce <- ggplot_nlp_mce(
      dat = summary_table,
      estim = e$estim,
      mce = e$mce,
      method_var = "method",
      true = 0,  # Bias/RMSE/SE have true = 0 in this context
      step_factors = "parameter"
    ) + labs(y = e$ylab, x = "Parameter") +
      scale_x_continuous(name = "", labels = NULL)+theme(
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor.x = ggplot2::element_blank()
      )

    keep_layers <- sapply(p_mce$layers, function(l) {
      is_text <- inherits(l$geom, "GeomText")
      is_step_with_group <- inherits(l$geom, "GeomStep") && 
        !is.null(l$mapping$group) && 
        rlang::as_label(l$mapping$group) == "step_ID"
      
      !(is_text || is_step_with_group)
    })

    p_mce$layers <- p_mce$layers[keep_layers]
    
    point_layer_idx <- which(sapply(p_mce$layers, function(l) inherits(l$geom, "GeomPoint")))
    # Update the size of the points 
      p_mce$layers[[point_layer_idx]]$aes_params$size <- 4  
    
    errorbar_layer_idx <- which(sapply(p_mce$layers, function(l) inherits(l$geom, "GeomErrorbar")))
    p_mce$layers[[errorbar_layer_idx]]$aes_params$linewidth <- 1 
    
    estim_col <- rlang::as_name(p_mce$mapping$y) 
    mce_col <- paste0("mce_",e$estim)
    
    max_val <- max(p_mce$data[[estim_col]] + 1.96 * p_mce$data[[mce_col]])
    upper_bound <- max_val * 1.1  # 10% margin above the data
    min_val <- min(p_mce$data[[estim_col]] - 1.96 * p_mce$data[[mce_col]])
    lower_bound <- ifelse(e$estim=="bias",min_val * 1.1,0)
    
    p_mce <- p_mce + ggplot2::coord_cartesian(ylim = c(lower_bound, upper_bound))
    
    p_mce <- p_mce +
      scale_x_continuous(
        breaks = seq(1.5, 8.5, by = 1),  
        labels = TeX(c(
          "$X_{1,12}$", 
          "$X_{2,12}$", 
          "$X_{1,13}$", 
          "$X_{2,13}$", 
          "$\\gamma_{13}$", 
          "$intercept_{cure}$", 
          "$X_{1,cure}$", 
          "$X_{2,cure}$"
        ))
      ) +
      theme(axis.text.x = element_text(color = "black",size = 13),
            axis.text.y=element_text(color = "black",size = 13),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12))

    ggsave(
      filename = paste0(path, "nested_", e$estim, "_",scen,".pdf"),
      plot = p_mce,
      width = 10,
      height = 8
    )
  }
}

get_runtime_summary <- function(file_path, file_name) {
  results <- readRDS(paste0(file_path, file_name))
  wide_df <- do.call(rbind, lapply(results, function(repetition) {
    data.frame(
      ELF       = repetition[["ELF"]]$runtime,
      Multicure = repetition[["Multicure"]]$runtime
    )
  }))
  long_df <- pivot_longer(wide_df, cols = c(ELF, Multicure),
                          names_to = "method", values_to = "runtime")
  summary_df <- long_df %>%
    group_by(method) %>%
    summarise(mean_runtime = mean(runtime),
              sd_runtime = sd(runtime),
              median_runtime = median(runtime),
              IQR_runtime = IQR(runtime),
              n = n(),
              .groups = "drop")
  return(summary = summary_df)
}

###zero tail results ###
#performance measures
summary_zero<-process_simulation_results(
  file_path = "results_zerotail/",
  file_name = "results_zerotail.rds"
)

create_nested_plots(summary_table = summary_zero, 
                        path = "results_zerotail/",scen="zerotail")

#run time comparison
runtime_zerotail<-get_runtime_summary( file_path = "results_zerotail/",
                               file_name = "results_zerotail.rds")
write.table(runtime_zerotail,file="results_zerotail/runtime_summary_zerotail.csv",sep = ",",row.names=FALSE)

### nozero tail results ###
#performance measures
summary_nozero<-process_simulation_results(
  file_path = "results_nozerotail/",
  file_name = "results_nozerotail.rds"
)

create_nested_plots(summary_table = summary_nozero, 
                    path = "results_nozerotail/",
                    scen="nozerotail")

#run time comparison
runtime_nozerotail<-get_runtime_summary( file_path = "results_nozerotail/",
                                       file_name = "results_nozerotail.rds")
write.table(runtime_nozerotail,file="results_nozerotail/runtime_summary_nozerotail.csv",sep = ",",row.names=FALSE)
